//! cas 多项式层（P0-M2）：有序稀疏分布式多项式，泛型于 [`Ring`] 系数。
//!
//! 设计见 `cas/DESIGN.md` D3：
//! - **有序**：项按 `PolyRing` 携带的单项式序（缺省 DegRevLex）严格升序，
//!   [`Poly::terms`] 的次序即序列化/打印/代码生成次序（D6 四方同序）；
//! - **稀疏**：零系数项不存在，指数向量平坦存储（`exps[i*nvars..]`）；
//! - **确定性**：构造即规范化（合并同类项、去零、按项序排序）；
//!   乘法经哈希合并后**必经排序出口**，内部哈希迭代序不进入输出。
//!
//! M2 范围：加/减/乘/幂/偏导/求值。gcd、除法、因式分解属 M3/M4/M5；
//! ℚ 上"本原 ℤ 表示 + 容度"的去分母优化同属 M3。

use cas_domain::{Rational, Ring};
use std::cmp::Ordering;

mod factor;
mod factor_mv;
mod gcd;
use std::collections::HashMap;
use std::fmt;
use std::sync::Arc;

/// 单项式序。变元位置即 `PolyRing::vars` 的下标。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MonOrder {
    /// 字典序：左起第一个非零差决定。
    Lex,
    /// 分次反字典序（Gröbner 友好，缺省）：先比总次数；同次时右起
    /// 第一个非零差为负者更小。
    DegRevLex,
}

pub(crate) fn cmp_monomials(order: MonOrder, a: &[u32], b: &[u32]) -> Ordering {
    debug_assert_eq!(a.len(), b.len());
    match order {
        MonOrder::Lex => {
            for (&x, &y) in a.iter().zip(b.iter()) {
                if x != y {
                    return x.cmp(&y);
                }
            }
            Ordering::Equal
        }
        MonOrder::DegRevLex => {
            let da: u64 = a.iter().map(|&e| e as u64).sum();
            let db: u64 = b.iter().map(|&e| e as u64).sum();
            da.cmp(&db).then_with(|| {
                for i in (0..a.len()).rev() {
                    let d = a[i] as i64 - b[i] as i64;
                    if d != 0 {
                        return d.cmp(&0);
                    }
                }
                Ordering::Equal
            })
        }
    }
}

/// 多项式环身份 =（变元表, 项序）。相同身份的 Poly 才可相互运算。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PolyRing {
    pub vars: Vec<Box<str>>,
    pub order: MonOrder,
}

impl PolyRing {
    pub fn new<I: IntoIterator<Item: Into<Box<str>>>>(vars: I, order: MonOrder) -> Arc<Self> {
        let vars: Vec<Box<str>> = vars.into_iter().map(Into::into).collect();
        assert!(!vars.is_empty(), "PolyRing 至少一个变元");
        Arc::new(PolyRing { vars, order })
    }

    pub fn nvars(&self) -> usize {
        self.vars.len()
    }
}

/// 稀疏有序多项式。`exps` 长度恒为 `nvars × nterms`，第 i 项指数为
/// `exps[i*nvars..(i+1)*nvars]`；项按环的项序严格升序，系数零的项不存在。
#[derive(Clone, PartialEq, Eq)]
pub struct Poly<C: Ring> {
    ring: Arc<PolyRing>,
    exps: Vec<u32>,
    coefs: Vec<C>,
}

impl<C: Ring> Poly<C> {
    /// 由（指数向量, 系数）序列构造，构造即规范化：合并同类项、去零、排序。
    pub fn from_terms<I: IntoIterator<Item = (Vec<u32>, C)>>(
        ring: Arc<PolyRing>,
        terms: I,
    ) -> Self {
        let nvars = ring.nvars();
        let mut acc: HashMap<Vec<u32>, C> = HashMap::new();
        for (e, c) in terms {
            assert_eq!(e.len(), nvars, "指数向量长度须等于变元数 {nvars}");
            if c.is_zero() {
                continue;
            }
            let slot = acc.entry(e).or_insert_with(Ring::zero);
            *slot = slot.add(&c);
        }
        let mut items: Vec<(Vec<u32>, C)> = acc.into_iter().filter(|(_, c)| !c.is_zero()).collect();
        let order = ring.order;
        items.sort_by(|(a, _), (b, _)| cmp_monomials(order, a, b));
        let exps: Vec<u32> = items.iter().flat_map(|(e, _)| e.iter().copied()).collect();
        let coefs: Vec<C> = items.into_iter().map(|(_, c)| c).collect();
        Poly { ring, exps, coefs }
    }

    pub fn zero(ring: Arc<PolyRing>) -> Self {
        Poly {
            ring,
            exps: Vec::new(),
            coefs: Vec::new(),
        }
    }

    pub fn constant(ring: Arc<PolyRing>, c: C) -> Self {
        if c.is_zero() {
            Self::zero(ring)
        } else {
            let n = ring.nvars();
            Poly {
                ring,
                exps: vec![0; n],
                coefs: vec![c],
            }
        }
    }

    pub fn ring(&self) -> &Arc<PolyRing> {
        &self.ring
    }

    pub fn nterms(&self) -> usize {
        self.coefs.len()
    }

    pub fn is_zero(&self) -> bool {
        self.coefs.is_empty()
    }

    pub fn is_constant(&self) -> bool {
        self.coefs.len() <= 1 && self.exps.iter().all(|&e| e == 0)
    }

    /// 总次数（零多项式为 0）。
    pub fn degree(&self) -> u64 {
        self.nvars_terms()
            .map(|(e, _)| e.iter().map(|&x| x as u64).sum())
            .max()
            .unwrap_or(0)
    }

    /// 稳定项迭代器：次序 = 项序 = 序列化/打印/代码生成次序（D6）。
    pub fn terms(&self) -> impl Iterator<Item = (&[u32], &C)> {
        let n = self.ring.nvars();
        self.coefs
            .iter()
            .enumerate()
            .map(move |(i, c)| (&self.exps[i * n..(i + 1) * n], c))
    }

    fn nvars_terms(&self) -> impl Iterator<Item = (&[u32], &C)> {
        self.terms()
    }

    fn exps_of(&self, i: usize) -> &[u32] {
        let n = self.ring.nvars();
        &self.exps[i * n..(i + 1) * n]
    }

    pub fn add(&self, other: &Self) -> Self {
        self.assert_same_ring(other);
        // 两序列均已按项序升序 → 双指针归并
        let (mut i, mut j, mut items) = (0, 0, Vec::with_capacity(self.nterms() + other.nterms()));
        while i < self.nterms() && j < other.nterms() {
            match cmp_monomials(self.ring.order, self.exps_of(i), other.exps_of(j)) {
                Ordering::Less => {
                    items.push((self.exps_of(i).to_vec(), self.coefs[i].clone()));
                    i += 1;
                }
                Ordering::Greater => {
                    items.push((other.exps_of(j).to_vec(), other.coefs[j].clone()));
                    j += 1;
                }
                Ordering::Equal => {
                    let c = self.coefs[i].add(&other.coefs[j]);
                    if !c.is_zero() {
                        items.push((self.exps_of(i).to_vec(), c));
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        items.extend(self.terms_from(i).map(|(e, c)| (e.to_vec(), c.clone())));
        items.extend(other.terms_from(j).map(|(e, c)| (e.to_vec(), c.clone())));
        Self::from_sorted(self.ring.clone(), items)
    }

    fn terms_from(&self, start: usize) -> impl Iterator<Item = (&[u32], &C)> {
        let n = self.ring.nvars();
        (start..self.nterms()).map(move |i| (&self.exps[i * n..(i + 1) * n], &self.coefs[i]))
    }

    /// 由**已按项序升序**的项构造（归并路径专用，跳过排序）。
    fn from_sorted(ring: Arc<PolyRing>, items: Vec<(Vec<u32>, C)>) -> Self {
        let exps: Vec<u32> = items.iter().flat_map(|(e, _)| e.iter().copied()).collect();
        let coefs: Vec<C> = items.into_iter().map(|(_, c)| c).collect();
        Poly { ring, exps, coefs }
    }

    pub fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    pub fn neg(&self) -> Self {
        Poly {
            ring: self.ring.clone(),
            exps: self.exps.clone(),
            coefs: self.coefs.iter().map(Ring::neg).collect(),
        }
    }

    pub fn mul(&self, other: &Self) -> Self {
        self.assert_same_ring(other);
        if self.is_zero() || other.is_zero() {
            return Self::zero(self.ring.clone());
        }
        let mut acc: HashMap<Vec<u32>, C> = HashMap::with_capacity(self.nterms() * other.nterms());
        for i in 0..self.nterms() {
            let (ea, ca) = (self.exps_of(i), &self.coefs[i]);
            for j in 0..other.nterms() {
                let key: Vec<u32> = ea
                    .iter()
                    .zip(other.exps_of(j))
                    .map(|(&x, &y)| x + y)
                    .collect();
                let slot = acc.entry(key).or_insert_with(Ring::zero);
                *slot = slot.add(&ca.mul(&other.coefs[j]));
            }
        }
        let mut items: Vec<(Vec<u32>, C)> = acc.into_iter().filter(|(_, c)| !c.is_zero()).collect();
        let order = self.ring.order;
        // 哈希迭代不进入输出（D6 纪律）：出口必经排序
        items.sort_by(|(a, _), (b, _)| cmp_monomials(order, a, b));
        Self::from_sorted(self.ring.clone(), items)
    }

    pub fn pow(&self, e: u32) -> Self {
        match e {
            0 => Self::constant(self.ring.clone(), C::one()),
            1 => self.clone(),
            _ => {
                // 平方求幂
                let mut acc = Self::constant(self.ring.clone(), C::one());
                let mut base = self.clone();
                let mut k = e;
                while k > 0 {
                    if k & 1 == 1 {
                        acc = acc.mul(&base);
                    }
                    k >>= 1;
                    if k > 0 {
                        base = base.mul(&base);
                    }
                }
                acc
            }
        }
    }

    /// 对第 `idx` 个变元求偏导。
    pub fn deriv(&self, idx: usize) -> Self {
        assert!(idx < self.ring.nvars(), "变元下标越界");
        let items: Vec<(Vec<u32>, C)> = self
            .terms()
            .filter(|(e, _)| e[idx] > 0)
            .map(|(e, c)| {
                let mut e2 = e.to_vec();
                e2[idx] -= 1;
                let c2 = c.mul(&C::from_i64_coeff(e[idx] as i64));
                (e2, c2)
            })
            .collect();
        Self::from_terms(self.ring.clone(), items)
    }

    /// 在给定点求值（点分量按变元次序）。长度不符属编程错误。
    pub fn eval(&self, vals: &[C]) -> C {
        assert_eq!(vals.len(), self.ring.nvars(), "求值点分量数须等于变元数");
        let mut sum = C::zero();
        for (e, c) in self.terms() {
            let mut m = c.clone();
            for (j, &power) in e.iter().enumerate() {
                m = m.mul(&vals[j].pow_u32(power));
            }
            sum = sum.add(&m);
        }
        sum
    }

    fn assert_same_ring(&self, other: &Self) {
        assert_eq!(
            &*self.ring, &*other.ring,
            "多项式须属于同一 PolyRing（变元表+项序）"
        );
    }
}

// 系数从 u64 构造的小助手已并入 Ring::from_i64_coeff。

impl Poly<Rational> {
    /// f64 求值（系数精确转 f64 后按项求和）——数值交叉验证用。
    pub fn eval_f64(&self, vals: &[f64]) -> f64 {
        assert_eq!(vals.len(), self.ring.nvars());
        let mut sum = 0.0;
        for (e, c) in self.terms() {
            let mut m = c.to_f64();
            for (j, &power) in e.iter().enumerate() {
                m *= vals[j].powi(power as i32);
            }
            sum += m;
        }
        sum
    }
}

impl<C: Ring + fmt::Debug> fmt::Debug for Poly<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let mut first = true;
        for (e, c) in self.terms() {
            if !first {
                write!(f, " + ")?;
            }
            first = false;
            write!(f, "({c:?})")?;
            for (j, &power) in e.iter().enumerate() {
                if power > 0 {
                    write!(f, "*{}^{}", self.ring.vars[j], power)?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cas_domain::Integer;

    fn ring6() -> Arc<PolyRing> {
        PolyRing::new(["q1", "q2", "q3", "p1", "p2", "p3"], MonOrder::DegRevLex)
    }

    fn rat(n: i64, d: i64) -> Rational {
        Rational::from_ints(&Integer::from_i64(n), &Integer::from_i64(d)).unwrap()
    }

    fn term(_ring: &Arc<PolyRing>, e: [u32; 6], c: Rational) -> (Vec<u32>, Rational) {
        (e.to_vec(), c)
    }

    #[test]
    fn 构造即规范化() {
        let ring = ring6();
        // 乱序 + 重复同类项 + 零系数 → 合并排序去零
        let p = Poly::from_terms(
            ring.clone(),
            [
                term(&ring, [0, 0, 0, 0, 0, 0], rat(1, 2)),
                term(&ring, [2, 0, 0, 0, 0, 0], rat(1, 3)),
                term(&ring, [2, 0, 0, 0, 0, 0], rat(1, 6)),
                term(&ring, [1, 0, 0, 0, 0, 0], rat(0, 1)),
            ],
        );
        assert_eq!(p.nterms(), 2);
        // grevlex：常数（次数 0）在最前
        let (e0, c0) = p.terms().next().unwrap();
        assert!(e0.iter().all(|&x| x == 0) && c0 == &rat(1, 2));
        // 同类项已合并
        let (_, c1) = p.terms().nth(1).unwrap();
        assert_eq!(*c1, rat(1, 2));
    }

    #[test]
    fn 加乘分配律与交换律() {
        let ring = ring6();
        let mk = |seed: &mut u64| -> Poly<Rational> {
            // 确定性小随机多项式
            let mut items = Vec::new();
            for _ in 0..6 {
                let r = xorshift(seed);
                let e: Vec<u32> = (0..6).map(|k| ((r >> (k * 5)) % 3) as u32).collect();
                let c = rat((r % 19) as i64 - 9, ((r >> 8) % 7 + 1) as i64);
                items.push((e, c));
            }
            Poly::from_terms(ring.clone(), items)
        };
        let mut s = 12345u64;
        for _ in 0..20 {
            let (a, b, c) = (mk(&mut s), mk(&mut s), mk(&mut s));
            assert!(a.add(&b) == b.add(&a), "加法交换律");
            assert!(a.mul(&b) == b.mul(&a), "乘法交换律");
            assert!(a.mul(&b).mul(&c) == a.mul(&b.mul(&c)), "结合律");
            assert!(a.add(&b).mul(&c) == a.mul(&c).add(&b.mul(&c)), "分配律");
        }
    }

    #[test]
    fn 幂偏导与求值() {
        let ring = ring6();
        // (1/2 q1)^4 = 1/16 q1^4
        let p = Poly::from_terms(ring.clone(), [term(&ring, [1, 0, 0, 0, 0, 0], rat(1, 2))]).pow(4);
        assert_eq!(p.nterms(), 1);
        assert_eq!(p.terms().next().unwrap().1, &rat(1, 16));

        // d/dq1 (3 q1^2 p2) = 6 q1 p2
        let f = Poly::from_terms(ring.clone(), [term(&ring, [2, 0, 0, 0, 1, 0], rat(3, 1))]);
        let df = f.deriv(0);
        assert_eq!(
            df.terms().next().unwrap(),
            (&[1, 0, 0, 0, 1, 0][..], &rat(6, 1))
        );

        // 求值一致性：(f*g)(x) == f(x)*g(x)
        let g = Poly::from_terms(ring.clone(), [term(&ring, [0, 1, 0, 0, 0, 0], rat(-1, 3))]);
        let pt: Vec<Rational> = [2, -3, 5, 7, 11, 13].map(|v| rat(v, 1)).to_vec();
        let lhs = f.mul(&g).eval(&pt);
        let rhs = f.eval(&pt).mul(&g.eval(&pt));
        assert_eq!(lhs, rhs);

        // f64 求值与精确求值一致
        let ptf: Vec<f64> = [2.0, -3.0, 5.0, 7.0, 11.0, 13.0].to_vec();
        let exact = f.mul(&g).eval(&pt);
        let flt = f.mul(&g).eval_f64(&ptf);
        assert!((exact.to_f64() - flt).abs() < 1e-12);
    }

    #[test]
    fn 项序_grevlex约定() {
        let a = [2u32, 0, 0, 0, 0, 0]; // q1^2，总次 2
        let b = [0u32, 0, 0, 0, 0, 1]; // p3，总次 1
        let c = [1u32, 0, 0, 0, 0, 1]; // q1 p3，总次 2
        // grevlex：先比总次数（p3 < q1^2）
        assert_eq!(cmp_monomials(MonOrder::DegRevLex, &b, &a), Ordering::Less);
        // 同总次：右起第一个非零差为负者更小 → q1^2 < q1*p3
        assert_eq!(cmp_monomials(MonOrder::DegRevLex, &a, &c), Ordering::Less);
        // 字典序：左起第一个非零差为正者更大 → q1^2 > q1*p3
        assert_eq!(cmp_monomials(MonOrder::Lex, &a, &c), Ordering::Greater);
    }

    fn xorshift(s: &mut u64) -> u64 {
        let mut x = *s | 1;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        *s = x;
        x
    }
}
