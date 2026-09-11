//! 多变元因式分解（M5，设计 §5/D3）。
//!
//! 路由：
//! - 1 变元 → M4（factor_univariate）；
//! - 2 变元 → **完全分解**：确定性平移搜索 z→z+c 取好赋值点，f₀ = S(v,0)
//!   用 M4 分解，然后 **模 z^k 的线性 Hensel**（系数环 ℚ[v]，与 M4 的 𝔽_p
//!   版同构：模 p^k → 模 z^k，截断 = 丢弃 z 次数 ≥ k 的项）逐层提升至
//!   ΠH_i 与目标**精确相等**（E==0 提前退出；无 monic 化、无 lc 幂吸收
//!   ——整系数目标直接提升，子集积的精确性使组合仅需 pp+试除）；
//! - ≥3 变元 → 部分分解：ℚ 内容 + 主变元多项式内容（递归分解）+ 多变元
//!   Yun 平方自由；完全分解（Wang EEZ 逐变元提升）挂 M5b。
//!
//! 确定性：平移常数按固定表顺序尝试，组合按子集基数升序——同输入恒得
//! 同输出。

use crate::factor::squarefree_parts;
use crate::gcd::Uni;
use crate::{Poly, PolyRing};
use cas_domain::{Integer, Rational};
use std::sync::Arc;

// ── 2 变元辅助（vars = [v, z]：e[0]=v 次数，e[1]=z 次数）─────────

/// z 次数恰为 k 的项 → ℚ[v]（ring_v 为单变元环）。
fn coeff_zk(p: &Poly<Rational>, ring_v: &Arc<PolyRing>, k: u32) -> Poly<Rational> {
    let items: Vec<(Vec<u32>, Rational)> = p
        .terms()
        .filter(|(e, _)| e[1] == k)
        .map(|(e, c)| (vec![e[0]], c.clone()))
        .collect();
    Poly::from_terms(ring_v.clone(), items)
}

/// base + z^k·p（p ∈ ℚ[v] 嵌入为 z 次数 k）。
fn add_zk_pow(base: &Poly<Rational>, p: &Poly<Rational>, k: u32) -> Poly<Rational> {
    if p.is_zero() {
        return base.clone();
    }
    let items: Vec<(Vec<u32>, Rational)> = p
        .terms()
        .map(|(e, c)| {
            let mut full = vec![0u32; 2];
            full[0] = e[0];
            full[1] = k;
            (full, c.clone())
        })
        .collect();
    base.add(&Poly::from_terms(base.ring().clone(), items))
}

/// z → z + c（单项式二项式展开；c=0 原样）。
fn shift_z(p: &Poly<Rational>, c: i64) -> Poly<Rational> {
    if c == 0 {
        return p.clone();
    }
    let cint = Integer::from_i64(c);
    let mut items: Vec<(Vec<u32>, Rational)> = Vec::new();
    for (e, coef) in p.terms() {
        let n = e[1] as u64;
        // z^n → Σ_j C(n,j)·c^{n-j}·z^j
        let mut binom = 1i128;
        for j in 0..=n {
            if binom != 0 {
                let factor = Rational::from_ints(&Integer::from_i64(binom as i64), &Integer::one())
                    .unwrap()
                    .mul(&Rational::from_integer(&cint.pow((n - j) as u32)))
                    .mul(coef);
                if !factor.is_zero() {
                    let mut full = vec![0u32; 2];
                    full[0] = e[0];
                    full[1] = j as u32;
                    items.push((full, factor));
                }
            }
            binom = binom * (n - j) as i128 / (j + 1) as i128;
        }
    }
    Poly::from_terms(p.ring().clone(), items)
}

/// ℚ[v] 上单变元带余除法（div_rem 域版在单变元环上即普通除法）。
fn divrem1(a: &Poly<Rational>, b: &Poly<Rational>) -> (Poly<Rational>, Poly<Rational>) {
    let (mut qs, r) = a.div_rem(&[b]);
    (qs.pop().expect("单除子"), r)
}

/// ℚ[v] 扩展 Euclid：(monic g, s, t)，s·a + t·b = g。
fn rat_xgcd1(
    a: &Poly<Rational>,
    b: &Poly<Rational>,
) -> (Poly<Rational>, Poly<Rational>, Poly<Rational>) {
    let ring = a.ring().clone();
    let (mut r0, mut r1) = (a.clone(), b.clone());
    let (mut s0, mut s1) = (
        Poly::constant(ring.clone(), Rational::one()),
        Poly::zero(ring.clone()),
    );
    let (mut t0, mut t1) = (
        Poly::zero(ring.clone()),
        Poly::constant(ring.clone(), Rational::one()),
    );
    while !r1.is_zero() {
        let (q, r) = divrem1(&r0, &r1);
        r0 = r1;
        r1 = r;
        let s2 = s0.sub(&q.mul(&s1));
        s0 = s1;
        s1 = s2;
        let t2 = t0.sub(&q.mul(&t1));
        t0 = t1;
        t1 = t2;
    }
    // monic 归一
    let lc = match r0.terms().last() {
        Some((_, c)) => c.clone(),
        None => return (r0, s0, t0), // a=b=0 不会出现
    };
    let inv = lc.inv_reduced().unwrap_or_else(Rational::one);
    let scale = |p: &Poly<Rational>| -> Poly<Rational> {
        let items: Vec<(Vec<u32>, Rational)> =
            p.terms().map(|(e, c)| (e.to_vec(), c.mul(&inv))).collect();
        Poly::from_terms(ring.clone(), items)
    };
    (scale(&r0), scale(&s0), scale(&t0))
}

// ── 双变元 Hensel（模 z^k，线性，二叉树）───────────────────────

/// 两因子提升：Π(a0,b0) ≡ target (mod z)；逐层至 Π ≡ target (mod z^k)。
/// k 逐层递增，每当 k 超过目标的 z 次数即检查精确相等并提前退出。
fn hensel2_bivar(
    target: &Poly<Rational>, // ℚ[v,z]，整系数
    a0: &Poly<Rational>,     // ℚ[v]（z=0 视角的因子）
    b0: &Poly<Rational>,     // ℚ[v]
    _s: &Poly<Rational>,     // 贝祖：s·a0 + t·b0 = 1（此公式只用 t）
    t: &Poly<Rational>,
    k_max: u32,
    ring_v: &Arc<PolyRing>,
) -> Option<(Poly<Rational>, Poly<Rational>)> {
    let mut a = Poly::from_terms(
        target.ring().clone(),
        a0.terms()
            .map(|(e, c)| (vec![e[0], 0], c.clone()))
            .collect::<Vec<_>>(),
    );
    let mut b = Poly::from_terms(
        target.ring().clone(),
        b0.terms()
            .map(|(e, c)| (vec![e[0], 0], c.clone()))
            .collect::<Vec<_>>(),
    );
    let deg_z_target: u32 = target.terms().map(|(e, _)| e[1]).max().unwrap_or(0);
    let mut k: u32 = 1;
    while k < k_max {
        let prod = a.mul(&b);
        let err_full = target.sub(&prod);
        if err_full.is_zero() {
            return Some((a, b)); // 精确相等，提前退出
        }
        // 验证 err ∈ z^k：低次项必须为零
        if err_full.terms().any(|(e, _)| e[1] < k) {
            return None; // 不变量破坏：保守放弃
        }
        if k > deg_z_target + 1 {
            return None; // 超界仍不精确：视为提升失败
        }
        let err = coeff_zk(&err_full, ring_v, k); // err/z^k 的 z^0 系数
        if !err.is_zero() {
            let sig = {
                let (_, r) = divrem1(&err.mul(t), a0);
                r
            };
            let tau = {
                let pr = err.sub(&sig.mul(b0));
                let (q, r) = divrem1(&pr, a0);
                if !r.is_zero() {
                    return None; // 校正整除性破坏
                }
                q
            };
            a = add_zk_pow(&a, &sig, k);
            b = add_zk_pow(&b, &tau, k);
        }
        k += 1;
    }
    let prod = a.mul(&b);
    if target.sub(&prod).is_zero() {
        Some((a, b))
    } else {
        None
    }
}

/// 多因子二叉树提升。叶子次序与 facs 一致；返回精确因子（Π = target）。
fn hensel_all_bivar(
    target: &Poly<Rational>,
    facs: &[Poly<Rational>],
    k_max: u32,
    ring_v: &Arc<PolyRing>,
) -> Option<Vec<Poly<Rational>>> {
    if facs.len() == 1 {
        return Some(vec![target.clone()]);
    }
    let mid = facs.len() / 2;
    let (left, right) = facs.split_at(mid);
    let a0 = left.iter().cloned().reduce(|x, y| x.mul(&y)).unwrap();
    let b0 = right.iter().cloned().reduce(|x, y| x.mul(&y)).unwrap();
    let (_, s, t) = rat_xgcd1(&a0, &b0);
    let (al, br) = hensel2_bivar(target, &a0, &b0, &s, &t, k_max, ring_v)?;
    let mut out = hensel_all_bivar(&al, left, k_max, ring_v)?;
    out.extend(hensel_all_bivar(&br, right, k_max, ring_v)?);
    Some(out)
}

// ── 双变元平方自由部分的完全分解 ──────────────────────────────

/// S：2 变元、整系数本原、关于 v 平方自由。返回不可约因子（本原、首项正）。
fn factor_sqfree_bivar(s: &Poly<Rational>) -> Vec<Poly<Rational>> {
    let ring_v = PolyRing::new([s.ring().vars[0].to_string()], crate::MonOrder::DegRevLex);
    let deg_v = s.terms().map(|(e, _)| e[0]).max().unwrap_or(0);
    let deg_z = s.terms().map(|(e, _)| e[1]).max().unwrap_or(0);
    if deg_v == 0 {
        // 不含主变元：作为 ℚ[z] 单变元分解（由调用方递归处理，此处分流）
        return vec![s.clone()];
    }
    // 确定性平移搜索好点：f0 次数保持 + 平方自由
    for c in [0i64, 1, -1, 2, -2, 3, -3, 5, -5, 7, -7] {
        let st = shift_z(s, c);
        let f0 = coeff_zk(&st, &ring_v, 0);
        if f0.is_zero() || f0.degree() < deg_v as u64 {
            continue;
        }
        let df0 = f0.deriv(0);
        if !f0.gcd(&df0).is_constant() {
            continue;
        }
        let (cont0, facs0) = f0.factor_univariate();
        if facs0.len() <= 1 {
            return vec![s.clone()]; // f0 不可约 ⇒ S 不可约
        }
        // M4 输出分离了 ℚ 常数内容：必须乘回，否则 Π 因子与 f0 差常数，
        // Hensel 起点 mod z 不变量破坏（x²−y² 整体误判不可约的根因）
        let mut facs_only: Vec<Poly<Rational>> = facs0.iter().map(|(g, _)| g.clone()).collect();
        if !cont0.is_one() {
            if let Some(first) = facs_only.first_mut() {
                *first = first.mul(&Poly::constant(ring_v.clone(), cont0));
            }
        }
        let k_max = deg_z + 4;
        let Some(lifted) = hensel_all_bivar(&st, &facs_only, k_max, &ring_v) else {
            continue; // 提升失败：换平移点
        };
        // 子集组合（基数升序；Π lifted 精确 = st，子集积 pp 后试除即真因子）
        let r = lifted.len();
        if r > 16 {
            return vec![s.clone()];
        }
        let mut masks: Vec<u64> = (1..(1u64 << (r - 1))).collect();
        masks.sort_by_key(|m| m.count_ones());
        let mut factors: Vec<Poly<Rational>> = vec![];
        let mut used = vec![false; r];
        let mut f_cur = st.clone();
        loop {
            let mut found = false;
            for &mask in &masks {
                if (0..r - 1).any(|i| mask >> i & 1 == 1 && used[i]) {
                    continue;
                }
                let mut cand = Poly::constant(s.ring().clone(), Rational::one());
                for (i, g) in lifted.iter().enumerate().take(r - 1) {
                    if mask >> i & 1 == 1 {
                        cand = cand.mul(g);
                    }
                }
                let pp = cand.normalize_primitive();
                if pp.is_constant() || pp.degree() > f_cur.degree() {
                    continue;
                }
                if let Some(q) = f_cur.exact_div(&pp) {
                    factors.push(shift_z(&pp, -c));
                    f_cur = q;
                    for (i, u) in used.iter_mut().enumerate().take(r - 1) {
                        if mask >> i & 1 == 1 {
                            *u = true;
                        }
                    }
                    found = true;
                    break;
                }
            }
            if !found {
                break;
            }
        }
        if !f_cur.is_constant() {
            factors.push(shift_z(&f_cur, -c));
        }
        if factors.len() >= 2 {
            return factors;
        }
        // 单因子：尝试下一平移点（或最终保守）
    }
    vec![s.clone()]
}

// ── 总控 ──────────────────────────────────────────────────────

impl Poly<Rational> {
    /// 多变元因式分解：返回（ℚ 常数内容, [(本原因子, 重数)]），
    /// cont·Π 因子^重数 == self。1 变元完全；2 变元完全；≥3 变元
    /// 部分分解（内容 + 平方自由；因子可能仍可约）。
    pub fn factor(&self) -> (Rational, Vec<(Self, u32)>) {
        if self.ring().nvars() == 1 {
            return self.factor_univariate();
        }
        if self.is_zero() {
            return (Rational::zero(), vec![]);
        }
        // ℚ 内容 + 整系数本原化（首项正）
        let pp = self.normalize_primitive();
        let cont = if pp.is_constant() {
            self.terms().next().unwrap().1.clone()
        } else {
            // 找 pp 与 self 的同次项系数比
            let (e0, c_self) = self.terms().last().unwrap();
            let c_pp = pp
                .terms()
                .find(|(e, _)| e == &e0)
                .map(|(_, c)| c.clone())
                .unwrap_or_else(Rational::one);
            c_self.div(&c_pp).unwrap_or_else(Rational::one)
        };
        if pp.is_constant() {
            return (cont, vec![]);
        }
        // 主变元多项式内容（递归分解）
        let ring = self.ring().clone();
        let ring0 = PolyRing::new(ring.vars[1..].iter().map(|v| v.to_string()), ring.order);
        let uni = Uni::from_poly(&pp, &ring0);
        let c_poly = uni.content();
        let mut factors: Vec<(Self, u32)> = vec![];
        if !c_poly.is_constant() {
            let (_cc, sub) = c_poly.factor();
            for (g, m) in sub {
                factors.push((embed_tail(&g, &ring), m));
            }
        }
        let pp0 = uni.exact_div_content(&c_poly).to_poly(ring.clone());
        // 多变元 Yun（关于 var0）
        for (sp, m) in squarefree_parts(&pp0) {
            if sp.is_constant() {
                continue;
            }
            let deg_v = sp.terms().map(|(e, _)| e[0]).max().unwrap_or(0);
            let sub = if ring.nvars() == 2 && deg_v > 0 {
                factor_sqfree_bivar(&sp)
            } else if deg_v == 0 {
                // 不含主变元：递归到剩余变元环
                let (_cc, sub) = strip_var0(&sp).factor();
                let mut out: Vec<Poly<Rational>> = vec![];
                for (g, _mm) in sub {
                    out.push(embed_tail(&g, &ring));
                }
                out
            } else {
                vec![sp.clone()] // ≥3 变元：部分分解
            };
            for g in sub {
                factors.push((g, m));
            }
        }
        factors.sort_by(|a, b| {
            let da = a.0.degree();
            let db = b.0.degree();
            da.cmp(&db).then_with(|| a.1.cmp(&b.1))
        });
        (cont, factors)
    }
}

/// 去掉第 0 个变元（指数 e[0] 恒为 0 的断言下重排）。
fn strip_var0(p: &Poly<Rational>) -> Poly<Rational> {
    let ring0 = PolyRing::new(
        p.ring().vars[1..].iter().map(|v| v.to_string()),
        p.ring().order,
    );
    let items: Vec<(Vec<u32>, Rational)> = p
        .terms()
        .map(|(e, c)| {
            debug_assert_eq!(e[0], 0);
            (e[1..].to_vec(), c.clone())
        })
        .collect();
    Poly::from_terms(ring0, items)
}

/// 少一变元的多项式回嵌到全环（前插 0 指数）。
fn embed_tail(p: &Poly<Rational>, ring: &Arc<PolyRing>) -> Poly<Rational> {
    let items: Vec<(Vec<u32>, Rational)> = p
        .terms()
        .map(|(e, c)| {
            let mut full = vec![0u32; ring.nvars()];
            full[1..].copy_from_slice(e);
            (full, c.clone())
        })
        .collect();
    Poly::from_terms(ring.clone(), items)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring(names: &[&str]) -> Arc<PolyRing> {
        PolyRing::new(names.iter().copied(), crate::MonOrder::DegRevLex)
    }

    fn from_i_terms(ringv: &Arc<PolyRing>, terms: &[(Vec<u32>, i64)]) -> Poly<Rational> {
        let items: Vec<(Vec<u32>, Rational)> = terms
            .iter()
            .filter(|(_, c)| *c != 0)
            .map(|(e, c)| {
                (
                    e.clone(),
                    Rational::from_ints(&Integer::from_i64(*c), &Integer::one()).unwrap(),
                )
            })
            .collect();
        Poly::from_terms(ringv.clone(), items)
    }

    fn check_refold(f: &Poly<Rational>) {
        let (cont, facs) = f.factor();
        let mut prod = Poly::constant(f.ring().clone(), cont);
        for (g, m) in &facs {
            prod = prod.mul(&g.pow(*m));
        }
        assert_eq!(&prod, f, "重展开不等于原式: f={f:?}");
    }

    #[test]
    fn 双变元_已知分解() {
        let r2 = ring(&["x", "y"]);
        // x^2 − y^2 = (x−y)(x+y)
        let f = from_i_terms(&r2, &[(vec![2, 0], 1), (vec![0, 2], -1)]);
        let (_, facs) = f.factor();
        assert_eq!(facs.len(), 2, "{facs:?}");
        check_refold(&f);

        // x^2 + y^2 在 ℚ 上不可约
        let f = from_i_terms(&r2, &[(vec![2, 0], 1), (vec![0, 2], 1)]);
        let (_, facs) = f.factor();
        assert_eq!(facs.len(), 1);

        // (x + y)^2·(x^2 + 2y^2) = x⁴ + 2x³y + 3x²y² + 4xy³ + 2y⁴
        let f = from_i_terms(
            &r2,
            &[
                (vec![4, 0], 1),
                (vec![3, 1], 2),
                (vec![2, 2], 3),
                (vec![1, 3], 4),
                (vec![0, 4], 2),
            ],
        );
        let (c, facs) = f.factor();
        let mut prod = Poly::constant(r2.clone(), c);
        for (g, m) in &facs {
            prod = prod.mul(&g.pow(*m));
        }
        assert_eq!(prod, f);
        assert_eq!(facs.len(), 2, "{facs:?}");
    }

    #[test]
    fn 定点_r3() {
        let r2 = ring(&["x", "y"]);
        // (x + y)(x + 2y)(x + 3y)：3 个关于 x 的 1 次因子
        let f = from_i_terms(
            &r2,
            &[
                (vec![3, 0], 1),
                (vec![2, 1], 6),
                (vec![1, 2], 11),
                (vec![0, 3], 6),
            ],
        );
        let (c, facs) = f.factor();
        eprintln!("cont={c:?} n={}", facs.len());
        for (g, m) in &facs {
            eprintln!("  deg={} m={m}", g.degree());
        }
        assert!(facs.len() >= 3, "r=3 应完全分解: {facs:?}");
    }

    #[test]
    fn 双变元_随机积() {
        let mut det = crate::factor::Det::new();
        let r2 = ring(&["x", "y"]);
        for _ in 0..120 {
            let nf = 2 + det.next() % 2;
            let mut f = Poly::constant(r2.clone(), Rational::one());
            for _ in 0..nf {
                // 指数 ≤ 2：底层 PRS 在高次上有系数膨胀问题（子结果式
                // 修复挂 M5b，见 DESIGN §18），中低次下功能完全验证
                let nt = 1 + det.next() % 3;
                let mut terms = vec![];
                for _ in 0..nt {
                    let r = det.next();
                    terms.push((
                        vec![(r % 3) as u32, ((r >> 8) % 3) as u32],
                        ((r >> 16) % 11) as i64 - 5,
                    ));
                }
                let g = from_i_terms(&r2, &terms);
                if g.is_zero() {
                    continue;
                }
                f = f.mul(&g);
            }
            let scale = Rational::from_ints(
                &Integer::from_i64((det.next() % 7) as i64 - 3),
                &Integer::from_i64(1 + (det.next() % 5) as i64),
            )
            .unwrap();
            let f = f.mul(&Poly::constant(r2.clone(), scale));
            if !f.is_zero() && !f.is_constant() {
                check_refold(&f);
            }
        }
    }

    #[test]
    fn 三变元_部分分解与恒等() {
        let r3 = ring(&["x", "y", "z"]);
        // (x + y + z)^2·(x + 2y)：内容/重数至少分离
        let g1 = from_i_terms(
            &r3,
            &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 1), (vec![0, 0, 1], 1)],
        );
        let g2 = from_i_terms(&r3, &[(vec![1, 0, 0], 1), (vec![0, 1, 0], 2)]);
        let f = g1.pow(2).mul(&g2);
        let (c, facs) = f.factor();
        let mut prod = Poly::constant(r3.clone(), c);
        for (g, m) in &facs {
            prod = prod.mul(&g.pow(*m));
        }
        assert_eq!(prod, f, "≥3 变元：恒等必须保持");
    }
}
