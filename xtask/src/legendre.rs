//! legendre — 复刻 MATLAB `Code02_Legendre_expr.m`：把三体距离倒数 `1/|r0−q|`
//! 展开成 `(q1,q2,q3)` 的次数 ≤ N 的多项式系数表。
//!
//! 与 MATLAB 逐行对应：
//!
//! ```matlab
//! P_list(1,1) = 1;  P_list(2,1) = dot(qv,r0v)/q/r0;      % P_0, P_1(t)
//! for i = 2:N
//!     n = i-1;
//!     P_list(i+1) = (2n+1)/(n+1)·t·P_list(i) − n/(n+1)·P_list(i-1);
//! end
//! P_list(i) = P_list(i)·(q/r0)^(i-1);
//! P = Σ P_list(i)/r0;  simplify; expand;
//! P = subs(P, q, sqrt(q1²+q2²+q3²));  expand;
//! Le_poly = expr2poly(P, [q1 q2 q3 p1 p2 p3]);
//! ```
//!
//! 输出：按幂次排序的 `pow1..pow6 | 系数` 文本表。
//! `--eval rx,ry,rz` 时把系数在给定几何上求值（`r0 = |r|`、`q` 已代入），
//! 便于与 MATLAB 侧 `double(subs(coeff,{rx,ry,rz},…))` 逐项对拍。
//!
//! 用法：
//!   xtask legendre --n 15 --out /tmp/Leg15.txt [--eval 1,0,0]

use cas_expr::{Context, Expr};
use std::process::ExitCode;

fn opt<'a>(args: &'a [String], name: &str, dflt: &'a str) -> &'a str {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .map(String::as_str)
        .unwrap_or(dflt)
}

/// 构造 `1/|r0−q|` 的展开并返回按幂次排序的系数表。
///
/// 变量顺序与 MATLAB `expr2poly(P, [q1 q2 q3 p1 p2 p3])` 一致（后三个不出现）。
pub fn legendre_table(ctx: &Context, n_max: u32) -> Vec<(Vec<u32>, Expr)> {
    let rx = ctx.sym("rx");
    let ry = ctx.sym("ry");
    let rz = ctx.sym("rz");
    let q1 = ctx.sym("q1");
    let q2 = ctx.sym("q2");
    let q3 = ctx.sym("q3");
    let r0 = ctx.sym("r0");
    // p1..p3 只是占位：与 MATLAB `expr2poly(P, [q1 q2 q3 p1 p2 p3])` 的变量表一致，
    // 保证输出幂次是 6 维（它们在 Legendre 段不出现，幂次恒为 0）。
    let _p1 = ctx.sym("p1");
    let _p2 = ctx.sym("p2");
    let _p3 = ctx.sym("p3");

    // 直接算齐次多项式 A_n := |q|^n·P_n(t)（t = (q·r0)/(|q|·r0)），避免 sqrt：
    //   A_0 = 1
    //   A_1 = (q·r0)/r0
    //   A_{n+1} = (2n+1)/(n+1)·((q·r0)/r0)·A_n − n/(n+1)·|q|²·A_{n−1}
    // 与 MATLAB 的 `P_list(i)·(q/r0)^{i-1}` + `subs(q,sqrt(...))` 数学等价，
    // 但全程只有 q 的正幂与 r0 的负幂（后者属于系数），不会产生 sqrt 的非多项式项。
    let dot = ctx.add(&[
        ctx.mul(&[q1.clone(), rx.clone()]),
        ctx.mul(&[q2.clone(), ry.clone()]),
        ctx.mul(&[q3.clone(), rz.clone()]),
    ]);
    let r0_inv = ctx.pow(&r0, &ctx.int(-1));
    let qt = ctx.mul(&[dot, r0_inv.clone()]); // (q·r0)/r0
    let q_sq = ctx.add(&[
        ctx.pow(&q1, &ctx.int(2)),
        ctx.pow(&q2, &ctx.int(2)),
        ctx.pow(&q3, &ctx.int(2)),
    ]);

    let mut a_list: Vec<Expr> = Vec::with_capacity(n_max as usize + 1);
    a_list.push(ctx.int(1));
    if n_max >= 1 {
        a_list.push(qt.clone());
    }
    for i in 2..=n_max {
        let n = (i - 1) as i64;
        let c1 = ctx.rational(2 * n + 1, n + 1).expect("有理数");
        let c2 = ctx.rational(-n, n + 1).expect("有理数");
        let term1 = ctx.mul(&[c1, qt.clone(), a_list[i as usize - 1].clone()]);
        let term2 = ctx.mul(&[c2, q_sq.clone(), a_list[i as usize - 2].clone()]);
        a_list.push(ctx.add(&[term1, term2]));
    }

    // P = (1/r0)·Σ_n (1/r0)^n·A_n
    let mut sum = ctx.int(0);
    for (n, a) in a_list.iter().enumerate() {
        let scaled = if n == 0 {
            a.clone()
        } else {
            ctx.mul(&[a.clone(), ctx.pow(&r0, &ctx.int(-(n as i64)))])
        };
        sum = ctx.add(&[sum, scaled]);
    }
    let p = ctx.mul(&[sum, r0_inv]);
    let p = ctx.expand(&p);

    let mut table = ctx
        .monomial_coeffs(&p, &["q1", "q2", "q3", "p1", "p2", "p3"])
        .unwrap_or_else(|e| panic!("monomial_coeffs 失败：{e}"));
    table.sort_by(|a, b| b.0.cmp(&a.0));
    table
}

pub fn run(args: &[String]) -> ExitCode {
    let n: u32 = opt(args, "--n", "15").parse().unwrap_or(15);
    let out = opt(args, "--out", "");
    let eval_at = opt(args, "--eval", "");

    let ctx = Context::new();
    let table = legendre_table(&ctx, n);

    let mut lines: Vec<String> = Vec::with_capacity(table.len());
    for (pow, coef) in &table {
        let head = pow
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let body = if eval_at.is_empty() {
            ctx.inspect(|i| i.dump(coef.raw_id()))
        } else {
            let nums: Vec<f64> = eval_at
                .split(',')
                .map(|s| s.trim().parse::<f64>().unwrap_or(0.0))
                .collect();
            let (rx, ry, rz) = (
                nums.first().copied().unwrap_or(1.0),
                nums.get(1).copied().unwrap_or(0.0),
                nums.get(2).copied().unwrap_or(0.0),
            );
            let r0 = (rx * rx + ry * ry + rz * rz).sqrt();
            let vals: Vec<(&str, f64)> = vec![("rx", rx), ("ry", ry), ("rz", rz), ("r0", r0)];
            match ctx.eval_float(coef, &vals) {
                Some(v) => format!("{v:.17e}"),
                None => "NaN".to_string(),
            }
        };
        lines.push(format!("{head} | {body}"));
    }
    let text = lines.join("\n") + "\n";
    if out.is_empty() {
        print!("{text}");
    } else {
        if let Err(e) = std::fs::write(out, &text) {
            eprintln!("写 {out} 失败：{e}");
            return ExitCode::FAILURE;
        }
        eprintln!("已写 {out}（{} 项，N={n}）", table.len());
    }
    ExitCode::SUCCESS
}
