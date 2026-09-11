//! 多变元除法与 gcd（M3，设计 D3）。
//!
//! - [`Poly::div_rem`]：域系数（ℚ）下的主单项式除法（多除子），
//!   Gröbner 归约的原语；余式的主单项式不被任何除式主单项式整除。
//! - [`Poly::gcd`]：递归**本原 PRS**——把 ℚ[x₁..xₙ] 视作 R[x₁]
//!   （R = ℚ[x₂..xₙ]），伪除法（只乘不除系数，环上精确）+ 每步取本原
//!   部分控制系数膨胀；系数环上的内容与公因子递归到一元（Euclid）。
//! - 规范化：本原（ℚ 内容归一）且主项系数为正——gcd 在相差 ℚ 常数
//!   意义下唯一，规范化后完全确定。

use crate::{MonOrder, Poly, PolyRing, cmp_monomials};
use cas_domain::{Integer, Rational};
use std::collections::HashMap;
use std::sync::Arc;

impl PolyRing {
    /// 去掉第 idx 个变元后的系数环（同项序），供递归 gcd。
    fn tail(&self, idx: usize) -> Arc<PolyRing> {
        let vars: Vec<Box<str>> = self
            .vars
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != idx)
            .map(|(_, v)| v.clone())
            .collect();
        PolyRing::new(vars, self.order)
    }
}

impl Poly<Rational> {
    /// 主单项式除法。确定性：每步消去被除式当前**最大**主单项式，
    /// 多除子按给定次序匹配。
    pub fn div_rem(&self, divs: &[&Self]) -> (Vec<Self>, Self) {
        for d in divs {
            assert!(!d.is_zero(), "除式为零");
            assert_eq!(&*self.ring, &*d.ring, "多项式须属于同一 PolyRing");
        }
        let n = self.ring.nvars();
        let order = self.ring.order;
        // 各除式的主单项式（项序最大项 = 升序存储的末项；注意 terms() 是升序）
        let divs_lt: Vec<(Vec<u32>, Rational)> = divs
            .iter()
            .map(|d| {
                let (e, c) = d.terms().last().expect("除式非零");
                (e.to_vec(), c.clone())
            })
            .collect();
        let mut qs: Vec<HashMap<Vec<u32>, Rational>> = vec![HashMap::new(); divs.len()];
        let mut rem: HashMap<Vec<u32>, Rational> = HashMap::new();
        let mut p: HashMap<Vec<u32>, Rational> =
            self.terms().map(|(e, c)| (e.to_vec(), c.clone())).collect();

        while let Some((lt_e, lt_c)) = leading(&p, order) {
            let mut matched = None;
            for (i, (le, lc)) in divs_lt.iter().enumerate() {
                if (0..n).all(|j| lt_e[j] >= le[j]) {
                    let qe: Vec<u32> = lt_e.iter().zip(le).map(|(&a, &b)| a - b).collect();
                    let qc = lt_c.div(lc).expect("主系数非零");
                    matched = Some((i, qe, qc));
                    break;
                }
            }
            match matched {
                Some((i, qe, qc)) => {
                    let qslot = qs[i].entry(qe.clone()).or_insert_with(Rational::zero);
                    *qslot = qslot.add(&qc);
                    for (de, dc) in divs[i].terms() {
                        let key: Vec<u32> = de.iter().zip(&qe).map(|(&a, &b)| a + b).collect();
                        let slot = p.entry(key).or_insert_with(Rational::zero);
                        *slot = slot.sub(&qc.mul(dc));
                    }
                }
                None => {
                    let slot = rem.entry(lt_e.clone()).or_insert_with(Rational::zero);
                    *slot = slot.add(&lt_c);
                    p.remove(&lt_e);
                    continue;
                }
            }
            p.retain(|_, c| !c.is_zero());
        }

        let qs: Vec<Poly<Rational>> = qs
            .into_iter()
            .map(|q| Poly::from_terms(self.ring.clone(), q))
            .collect();
        let r = Poly::from_terms(self.ring.clone(), rem);
        (qs, r)
    }

    /// 精确除法：g 整除 self 时返回商，否则 None。
    pub fn exact_div(&self, g: &Self) -> Option<Self> {
        let (mut qs, r) = self.div_rem(&[g]);
        if r.is_zero() && qs.len() == 1 {
            qs.pop()
        } else {
            None
        }
    }

    /// 最大公因子（本原规范化：ℚ 内容归一、主项系数为正）。
    /// gcd(0, g) = normalize(g)；gcd(c₁, c₂) = 1（ℚ 上常数互素）。
    pub fn gcd(&self, other: &Self) -> Self {
        assert_eq!(&*self.ring, &*other.ring, "多项式须属于同一 PolyRing");
        self.gcd_raw(other).normalize_primitive()
    }

    fn gcd_raw(&self, other: &Self) -> Self {
        if self.is_zero() {
            return other.clone();
        }
        if other.is_zero() {
            return self.clone();
        }
        if self.ring.nvars() == 1 {
            // ℚ[x]：Euclid
            let (mut a, mut b) = (self.clone(), other.clone());
            while !b.is_zero() {
                let (_, r) = a.div_rem(&[&b]);
                a = b;
                b = r;
            }
            return a;
        }
        // 多变元：视作 R[x₁]，R = ℚ[x₂..xₙ]，递归本原 PRS
        let ring0 = self.ring.tail(0);
        let ua = Uni::from_poly(self, &ring0);
        let ub = Uni::from_poly(other, &ring0);
        let ca = ua.content();
        let cb = ub.content();
        let c = ca.gcd_raw(&cb);
        let mut a = ua.exact_div_content(&ca);
        let mut b = ub.exact_div_content(&cb);
        while !b.is_zero() {
            let r = a.pseudo_rem(&b);
            if r.is_zero() {
                a = b;
                break;
            }
            let cr = r.content();
            let r = r.exact_div_content(&cr);
            a = b;
            b = r;
        }
        // a 为主变元视角的本原 gcd；乘回内容公因子并回嵌全部变元
        let g = a.to_poly(self.ring.clone());
        let c_full = embed_tail(&c, &self.ring.clone());
        g.mul(&c_full)
    }

    /// 本原规范化：除以 ℚ 内容使分子整系、分子 gcd 为 1，主项（项序最大
    /// 项）系数取正。零多项式原样。
    pub fn normalize_primitive(&self) -> Self {
        if self.is_zero() {
            return self.clone();
        }
        let mut d_lcm = Integer::one();
        for (_, c) in self.terms() {
            d_lcm = lcm(&d_lcm, &c.den());
        }
        let mut m = Integer::zero();
        for (_, c) in self.terms() {
            let scaled = c.num().mul(&d_lcm.div_exact(&c.den()));
            m = if m.is_zero() { scaled } else { m.gcd(&scaled) };
        }
        if m.is_zero() {
            m = Integer::one();
        }
        // p' = p · d_lcm / m；符号按主项（项序最大，即迭代末项）
        let mut terms: Vec<(Vec<u32>, Rational)> = self
            .terms()
            .map(|(e, c)| {
                let scaled = c.mul(&Rational::from_integer(&d_lcm));
                let v = scaled.div(&Rational::from_integer(&m)).expect("内容整除");
                (e.to_vec(), v)
            })
            .collect();
        if let Some((_, lc)) = terms.last() {
            if lc.is_negative() {
                for (_, c) in &mut terms {
                    *c = c.neg();
                }
            }
        }
        Poly::from_terms(self.ring.clone(), terms)
    }
}

/// map 中按项序最大的非零项。
fn leading(p: &HashMap<Vec<u32>, Rational>, order: MonOrder) -> Option<(Vec<u32>, Rational)> {
    p.iter()
        .filter(|(_, c)| !c.is_zero())
        .max_by(|(a, _), (b, _)| cmp_monomials(order, a, b))
        .map(|(e, c)| (e.clone(), c.clone()))
}

fn lcm(a: &Integer, b: &Integer) -> Integer {
    if a.is_zero() || b.is_zero() {
        return Integer::zero();
    }
    a.div_exact(&a.gcd(b)).mul(b)
}

/// 系数环（少一个变元）上的多项式回嵌到全环：指数向量前插 0。
fn embed_tail(p: &Poly<Rational>, ring: &Arc<PolyRing>) -> Poly<Rational> {
    let items: Vec<(Vec<u32>, Rational)> = p
        .terms()
        .map(|(e, c)| {
            let mut full = Vec::with_capacity(e.len() + 1);
            full.push(0);
            full.extend_from_slice(e);
            (full, c.clone())
        })
        .collect();
    Poly::from_terms(ring.clone(), items)
}

/// 主变元视角的稀疏一元表示：系数为剩余变元上的 Poly（ring0）。
#[derive(Clone)]
pub(crate) struct Uni {
    /// 按主变元次数降序；系数非零。
    terms: Vec<(u32, Poly<Rational>)>,
    ring0: Arc<PolyRing>,
}

impl Uni {
    pub(crate) fn from_poly(p: &Poly<Rational>, ring0: &Arc<PolyRing>) -> Self {
        let mut by_deg: HashMap<u32, Vec<(Vec<u32>, Rational)>> = HashMap::new();
        for (e, c) in p.terms() {
            by_deg
                .entry(e[0])
                .or_default()
                .push((e[1..].to_vec(), c.clone()));
        }
        let mut terms: Vec<(u32, Poly<Rational>)> = by_deg
            .into_iter()
            .map(|(d, items)| (d, Poly::from_terms(ring0.clone(), items)))
            .filter(|(_, p)| !p.is_zero())
            .collect();
        terms.sort_by_key(|(d, _)| std::cmp::Reverse(*d));
        Uni {
            terms,
            ring0: ring0.clone(),
        }
    }

    fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    fn deg(&self) -> u32 {
        self.terms[0].0
    }

    fn lc(&self) -> &Poly<Rational> {
        &self.terms[0].1
    }

    /// 系数（剩余变元多项式）的公因子——递归 gcd。
    pub(crate) fn content(&self) -> Poly<Rational> {
        let mut acc = Poly::zero(self.ring0.clone());
        for (_, c) in &self.terms {
            acc = acc.gcd_raw(c);
        }
        acc
    }

    /// 各系数被内容整除（构造保证整除）。
    pub(crate) fn exact_div_content(&self, c: &Poly<Rational>) -> Self {
        let terms: Vec<(u32, Poly<Rational>)> = self
            .terms
            .iter()
            .map(|(d, p)| (*d, p.exact_div(c).expect("内容必整除")))
            .collect();
        Uni {
            terms,
            ring0: self.ring0.clone(),
        }
    }

    /// 伪余式：r ← lc(b)·r − lc(r)·x^Δ·b，循环至 deg(r) < deg(b)。
    /// 只用环乘法（无系数除法），在 R[x₁] 上精确。
    fn pseudo_rem(&self, b: &Uni) -> Self {
        let db = b.deg();
        let mut r = self.clone();
        while !r.is_zero() && r.deg() >= db {
            let shift = r.deg() - db;
            let lcb = b.lc().clone();
            let lcr = r.lc().clone();
            let mut acc: HashMap<u32, Poly<Rational>> = HashMap::new();
            for (d, c) in &r.terms {
                let slot = acc
                    .entry(*d)
                    .or_insert_with(|| Poly::zero(self.ring0.clone()));
                *slot = slot.add(&lcb.mul(c));
            }
            for (d, c) in &b.terms {
                let slot = acc
                    .entry(d + shift)
                    .or_insert_with(|| Poly::zero(self.ring0.clone()));
                *slot = slot.sub(&lcr.mul(c));
            }
            let mut terms: Vec<(u32, Poly<Rational>)> =
                acc.into_iter().filter(|(_, p)| !p.is_zero()).collect();
            terms.sort_by_key(|(d, _)| std::cmp::Reverse(*d));
            // 注：朴素伪除 lc(b)^k 的系数膨胀在 2 变元高次（deg ≳ 6）时
            // 不可算。内层本原化可治但把递归 gcd 扇出放大到小输入超时，
            // 正解为子结果式 PRS（挂 M5b，见 DESIGN §18）。
            r = Uni {
                terms,
                ring0: self.ring0.clone(),
            };
        }
        r
    }

    /// 回嵌全环：指数向量 = [主变元次数, 剩余指数…]。
    pub(crate) fn to_poly(&self, ring: Arc<PolyRing>) -> Poly<Rational> {
        let mut items: Vec<(Vec<u32>, Rational)> = Vec::new();
        for (d, coef) in &self.terms {
            for (e, c) in coef.terms() {
                let mut full = Vec::with_capacity(ring.nvars());
                full.push(*d);
                full.extend_from_slice(e);
                items.push((full, c.clone()));
            }
        }
        Poly::from_terms(ring, items)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::MonOrder;

    fn ringn(names: &[&str]) -> Arc<PolyRing> {
        PolyRing::new(names.iter().copied(), MonOrder::DegRevLex)
    }

    fn rat(n: i64, d: i64) -> Rational {
        Rational::from_ints(&Integer::from_i64(n), &Integer::from_i64(d)).unwrap()
    }

    fn t(e: Vec<u32>, c: (i64, i64)) -> (Vec<u32>, Rational) {
        (e, rat(c.0, c.1))
    }

    #[test]
    fn 一元_gcd() {
        let ring = ringn(&["x"]);
        // (x^2 - 1) 与 (x - 1) 的 gcd = x - 1
        let f = Poly::from_terms(ring.clone(), [t(vec![2], (1, 1)), t(vec![0], (-1, 1))]);
        let g = Poly::from_terms(ring.clone(), [t(vec![1], (1, 1)), t(vec![0], (-1, 1))]);
        let d = f.gcd(&g);
        assert_eq!(d.nterms(), 2);
        assert_eq!(d.terms().next().unwrap(), (&[0][..], &rat(-1, 1)));
        assert_eq!(d.terms().nth(1).unwrap(), (&[1][..], &rat(1, 1)));

        // 常数与任何式互素 → 1
        let c = Poly::from_terms(ring.clone(), [t(vec![0], (4, 6))]);
        assert!(g.gcd(&c).is_constant());

        // gcd(f·h, g·h) = h·gcd(f,g)（h 本原）
        let h = Poly::from_terms(ring.clone(), [t(vec![3], (2, 1)), t(vec![1], (1, 1))]);
        let d1 = f.mul(&h).gcd(&g.mul(&h));
        let d2 = f.gcd(&g).mul(&h).normalize_primitive();
        assert_eq!(d1, d2);
    }

    #[test]
    fn 多元_gcd_六变元() {
        let ring = ringn(&["q1", "q2", "q3", "p1", "p2", "p3"]);
        let mk = |e: [u32; 6], c: (i64, i64)| t(e.to_vec(), c);
        // h = p1*q2 + q3（本原），f、g 与 h 互素
        let h = Poly::from_terms(
            ring.clone(),
            [
                mk([0, 1, 0, 1, 0, 0], (1, 1)),
                mk([0, 0, 1, 0, 0, 0], (1, 1)),
            ],
        );
        let f = Poly::from_terms(
            ring.clone(),
            [
                mk([1, 0, 0, 0, 0, 0], (1, 1)),
                mk([0, 0, 0, 0, 1, 0], (-2, 3)),
            ],
        );
        let g = Poly::from_terms(
            ring.clone(),
            [
                mk([0, 0, 0, 0, 0, 1], (5, 1)),
                mk([2, 0, 0, 1, 0, 0], (1, 2)),
            ],
        );
        let d = f.mul(&h).gcd(&g.mul(&h));
        assert_eq!(d, h);
        // 互质对 → 1
        assert!(f.gcd(&g).is_constant());
    }

    #[test]
    fn 精确除法与除法不变量() {
        let ring = ringn(&["x", "y"]);
        let mk = |e: [u32; 2], c: (i64, i64)| t(e.to_vec(), c);
        let f = Poly::from_terms(
            ring.clone(),
            [mk([2, 0], (1, 1)), mk([1, 1], (3, 2)), mk([0, 2], (-1, 5))],
        );
        let g = Poly::from_terms(ring.clone(), [mk([1, 0], (2, 1)), mk([0, 1], (1, 1))]);
        // f·g 精确除以 g == f
        let q = f.mul(&g).exact_div(&g);
        assert_eq!(q.as_ref(), Some(&f));
        // 不整除 → None（x 除以 y）
        let x = Poly::from_terms(ring.clone(), [mk([1, 0], (1, 1))]);
        let y = Poly::from_terms(ring.clone(), [mk([0, 1], (1, 1))]);
        assert_eq!(x.exact_div(&y), None);
        // 除法不变量：p = q·d + r
        let p = f.mul(&g);
        let (qs, r) = p.div_rem(&[&f]);
        let lhs = qs[0].mul(&f).add(&r);
        assert_eq!(lhs, p);
        assert!(r.is_zero());
    }

    #[test]
    fn 随机性质() {
        // gcd(f·h, g·h) = gcd(f,g)·h（h 本原）× 随机组
        let ring = ringn(&["q1", "q2", "q3", "p1", "p2", "p3"]);
        let mut xs = 4242u64;
        let mut nxt = move || {
            xs ^= xs << 13;
            xs ^= xs >> 7;
            xs ^= xs << 17;
            xs
        };
        fn rand_poly(
            nxt: &mut dyn FnMut() -> u64,
            deg: u64,
            ring: &Arc<PolyRing>,
        ) -> Poly<Rational> {
            let items: Vec<(Vec<u32>, Rational)> = (0..4)
                .map(|_| {
                    let r = nxt();
                    let e: Vec<u32> = (0..6)
                        .map(|k| ((r >> (k * 4)) % (deg + 1) % 3) as u32)
                        .collect();
                    (e, rat((r % 13) as i64 - 6, 1))
                })
                .collect();
            Poly::from_terms(ring.clone(), items)
        }
        for _ in 0..10 {
            let f = rand_poly(&mut nxt, 2, &ring);
            let g = rand_poly(&mut nxt, 2, &ring);
            // h 取两项整数系数（本原）
            let h = Poly::from_terms(
                ring.clone(),
                [
                    t(vec![1, 0, 0, 0, 0, 0], (1, 1)),
                    t(vec![0, 0, 0, 1, 1, 0], ((nxt() % 7) as i64, 1)),
                ],
            );
            if f.is_zero() || g.is_zero() || h.is_zero() {
                continue;
            }
            let d1 = f.mul(&h).gcd(&g.mul(&h));
            let d2 = f.gcd(&g).mul(&h).normalize_primitive();
            assert_eq!(d1, d2, "f={f:?} g={g:?} h={h:?}");
            // gcd 确为公因子
            assert!(f.mul(&h).exact_div(&d1).is_some());
            assert!(g.mul(&h).exact_div(&d1).is_some());
        }
    }

    #[test]
    fn 本原规范化() {
        let ring = ringn(&["x", "y"]);
        // (2/3)x + 4 → 本原化 = x + 6（内容 2/3 除尽）
        let p = Poly::from_terms(ring.clone(), [t(vec![1, 0], (2, 3)), t(vec![0, 0], (4, 1))]);
        let n = p.normalize_primitive();
        let items: Vec<_> = n.terms().map(|(e, c)| (e.to_vec(), c.clone())).collect();
        assert_eq!(
            items,
            vec![(vec![0, 0], rat(6, 1)), (vec![1, 0], rat(1, 1))]
        );
        // 主项为负 → 整体取正
        let p = Poly::from_terms(ring.clone(), [t(vec![1, 0], (-1, 1))]);
        let n = p.normalize_primitive();
        assert_eq!(n.terms().next().unwrap().1, &rat(1, 1));
    }
}
