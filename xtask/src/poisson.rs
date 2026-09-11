//! poisson — 6 变元泊松括号（qiao 正规化流水线的核心负载，CONTEXT.md §三）。
//!
//! 变元布局与 MATLAB poly_poisson 一致：q1 q2 q3 p1 p2 p3，
//! {f,g} = Σ_k (∂f/∂q_k · ∂g/∂p_k − ∂f/∂p_k · ∂g/∂q_k)。
//!
//! 检验：反对称/双线性/Jacobi 恒等式（精确，规范形相等）+
//! 有限差分数值交叉验证（独立于 deriv 代码路径）+ 计时。

use crate::util::Lcg;
use cas_domain::{Integer, Rational};
use cas_poly::{MonOrder, Poly, PolyRing};
use std::process::ExitCode;
use std::sync::Arc;
use std::time::Instant;

fn ring6() -> Arc<PolyRing> {
    PolyRing::new(["q1", "q2", "q3", "p1", "p2", "p3"], MonOrder::DegRevLex)
}

fn rat(n: i64, d: i64) -> Rational {
    Rational::from_ints(&Integer::from_i64(n), &Integer::from_i64(d)).unwrap()
}

/// 泊松括号：k = 0..3 为 q 下标，k+3 为对应 p。
fn poisson(f: &Poly<Rational>, g: &Poly<Rational>) -> Poly<Rational> {
    let mut acc = Poly::zero(f.ring().clone());
    for k in 0..3usize {
        let t1 = f.deriv(k).mul(&g.deriv(k + 3));
        let t2 = f.deriv(k + 3).mul(&g.deriv(k));
        acc = acc.add(&t1).sub(&t2);
    }
    acc
}

/// 随机稀疏多项式：总次 ≤ deg，nterms 项，小有理系数。
fn rand_poly(rng: &mut Lcg, ring: &Arc<PolyRing>, nterms: usize, deg: u32) -> Poly<Rational> {
    let mut items = Vec::with_capacity(nterms);
    for _ in 0..nterms {
        let e: Vec<u32> = (0..6)
            .map(|_| (rng.below((deg + 1) as u64) % 3) as u32)
            .collect();
        let total: u32 = e.iter().sum();
        let e: Vec<u32> = if total > deg {
            e.iter().map(|&x| x.min(1)).collect()
        } else {
            e
        };
        let c = rat(rng.below(19) as i64 - 9, rng.below(7) as i64 + 1);
        items.push((e, c));
    }
    Poly::from_terms(ring.clone(), items)
}

fn rand_point(rng: &mut Lcg) -> [f64; 6] {
    let mut v = [0.0f64; 6];
    for x in &mut v {
        let mag = 0.2 + rng.uniform(0.0, 1.0);
        *x = if rng.below(2) == 0 { mag } else { -mag };
    }
    v
}

pub(crate) fn run(args: &[String]) -> ExitCode {
    let seed: u64 = crate::arg(args, "--seed", "20260911")
        .parse()
        .unwrap_or(20260911);
    let terms: usize = crate::arg(args, "--terms", "120").parse().unwrap_or(120);
    let ring = ring6();
    let mut rng = Lcg::new(seed);
    let mut ok = true;

    println!("泊松括号 {{f,g}}：6 变元（q1..q3, p1..p3），grevlex，精确有理系数");

    // ── 精确恒等式 ──
    for t in 0..10 {
        let f = rand_poly(&mut rng, &ring, 6 + t, 3);
        let g = rand_poly(&mut rng, &ring, 6 + t, 3);
        let h = rand_poly(&mut rng, &ring, 5 + t, 3);

        // 反对称：{f,g} == -{g,f}
        let fg = poisson(&f, &g);
        let gf = poisson(&g, &f).neg();
        if fg != gf {
            eprintln!("反对称失败 @t={t}");
            ok = false;
        }
        // 双线性：{2f+h, g} == 2{f,g} + {h,g}
        let two = Poly::constant(ring.clone(), rat(2, 1));
        let lhs = poisson(&two.mul(&f).add(&h), &g);
        let rhs = two.mul(&poisson(&f, &g)).add(&poisson(&h, &g));
        if lhs != rhs {
            eprintln!("双线性失败 @t={t}");
            ok = false;
        }
        // Jacobi：{f,{g,h}} + {g,{h,f}} + {h,{f,g}} == 0
        let j = poisson(&f, &poisson(&g, &h))
            .add(&poisson(&g, &poisson(&h, &f)))
            .add(&poisson(&h, &poisson(&f, &g)));
        if !j.is_zero() {
            eprintln!("Jacobi 失败 @t={t}");
            ok = false;
        }
    }
    println!(
        "恒等式    : 反对称/双线性/Jacobi × 10 组随机多项式 {}",
        if ok { "全部通过" } else { "存在失败" }
    );

    // ── 数值交叉验证：{f,g}(x) vs 中心差分（独立于 deriv 路径）──
    let f = rand_poly(&mut rng, &ring, terms, 6);
    let g = rand_poly(&mut rng, &ring, terms.max(2) * 2 / 3, 5);
    let mut max_rel = 0.0f64;
    for _ in 0..5 {
        let x = rand_point(&mut rng);
        let exact = poisson(&f, &g).eval_f64(&x);
        let h = 1e-6;
        let mut fd = 0.0;
        for k in 0..3 {
            let mut xp = x;
            let mut xm = x;
            xp[k] += h;
            xm[k] -= h;
            let mut yp = x;
            let mut ym = x;
            yp[k + 3] += h;
            ym[k + 3] -= h;
            fd += (f.eval_f64(&xp) - f.eval_f64(&xm)) / (2.0 * h)
                * (g.eval_f64(&yp) - g.eval_f64(&ym))
                / (2.0 * h);
            let mut xp2 = x;
            let mut xm2 = x;
            xp2[k + 3] += h;
            xm2[k + 3] -= h;
            let mut yp2 = x;
            let mut ym2 = x;
            yp2[k] += h;
            ym2[k] -= h;
            fd -= (f.eval_f64(&xp2) - f.eval_f64(&xm2)) / (2.0 * h)
                * (g.eval_f64(&yp2) - g.eval_f64(&ym2))
                / (2.0 * h);
        }
        let rel = (exact - fd).abs() / (exact.abs() + fd.abs() + 1.0);
        max_rel = max_rel.max(rel);
    }
    let fd_ok = max_rel < 1e-5;
    println!(
        "数值交叉  : {{f,g}}(x) vs 中心差分（h=1e-6），5 点最大相对差 {max_rel:.2e} {}",
        if fd_ok { "通过" } else { "失败" }
    );
    ok &= fd_ok;

    // ── 计时（基准负载）──
    let t = Instant::now();
    let b = poisson(&f, &g);
    let d = t.elapsed();
    println!(
        "计时      : f(6 次 {} 项) × g(5 次 {} 项) → {{f,g}} {} 项，耗时 {d:?}",
        f.nterms(),
        g.nterms(),
        b.nterms()
    );
    if ok {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn 泊松括号_单对单() {
        // {q1, p1} = 1；{q1^2, p1} = 2 q1
        let ring = ring6();
        let q1 = Poly::from_terms(ring.clone(), [(vec![1, 0, 0, 0, 0, 0], rat(1, 1))]);
        let p1 = Poly::from_terms(ring.clone(), [(vec![0, 0, 0, 1, 0, 0], rat(1, 1))]);
        let b = poisson(&q1, &p1);
        assert!(b.is_constant() && b.terms().next().unwrap().1 == &rat(1, 1));

        let q1_2 = q1.pow(2);
        let b = poisson(&q1_2, &p1);
        assert_eq!(
            b.terms().next().unwrap(),
            (&[1, 0, 0, 0, 0, 0][..], &rat(2, 1))
        );

        // {p1, q1} = -1
        let b = poisson(&p1, &q1);
        assert_eq!(b.terms().next().unwrap().1, &rat(-1, 1));
    }
}
