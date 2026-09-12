//! hamilton — 复刻 MATLAB `Code03_Hamilton_expr.m`：把二次型与三体势能拼成 `H_poly`。
//!
//! 与 MATLAB 逐行对应：
//!
//! ```matlab
//! H = (fv')*qv + 1/2*(pv')*pv + (pv')*Cpq*qv + 1/2*(qv')*Cqq*qv;
//! H = expand(simplify(H));  H_poly = expr2poly(H);
//! % 三体：Le_poly 分别把 (rx,ry,rz,r0) 换成 (re*/rm*/rs*)，各乘 -mu_e/-mu_m/-mu_s
//! H_poly = poly_simplify([H_poly; Pe_poly; Pm_poly; Ps_poly]);
//! ```
//!
//! 输出 828 行（N=15）= 12（二次型）+ 816（三体按幂次合并后）；行序为**幂次向量升序字典序**
//! （实测 oracle `Hamilton_poly15.mat` 的首行为 [0,0,0,0,0,0]、次行 [0,0,0,0,0,2]）。
//!
//! 用法：
//!   xtask hamilton --n 15 --out /tmp/Ham15.txt [--eval "rex,rey,rez,...,mu_s"]

use crate::legendre::legendre_table;
use cas_expr::{Context, Expr};
use std::process::ExitCode;

fn opt<'a>(args: &'a [String], name: &str, dflt: &'a str) -> &'a str {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .map(String::as_str)
        .unwrap_or(dflt)
}

/// 二次型与线性项：`f·q + ½ p·p + p·Cpq·q + ½ q·Cqq·q`。
fn quadratic_terms(ctx: &Context) -> Vec<(Vec<u32>, Expr)> {
    let q: Vec<Expr> = (1..=3).map(|i| ctx.sym(&format!("q{i}"))).collect();
    let p: Vec<Expr> = (1..=3).map(|i| ctx.sym(&format!("p{i}"))).collect();
    let f: Vec<Expr> = (1..=3).map(|i| ctx.sym(&format!("f{i}"))).collect();
    let cpq: Vec<Expr> = (1..=9).map(|i| ctx.sym(&format!("Cpq{i}"))).collect();
    let cqq: Vec<Expr> = (1..=9).map(|i| ctx.sym(&format!("Cqq{i}"))).collect();

    let mut h = ctx.int(0);
    // f·q
    for (fi, qi) in f.iter().zip(q.iter()) {
        h = ctx.add(&[h, ctx.mul(&[fi.clone(), qi.clone()])]);
    }
    // ½ p·p
    let mut pp = ctx.int(0);
    for pi in p.iter() {
        pp = ctx.add(&[pp, ctx.mul(&[pi.clone(), pi.clone()])]);
    }
    h = ctx.add(&[h, ctx.mul(&[ctx.rational(1, 2).unwrap(), pp])]);
    // p·Cpq·q（MATLAB 列主序：Cpq(i,j) = cpq[3*(j-1) + (i-1)]，1-based 名字 Cpq1..9）
    let mut pcq = ctx.int(0);
    for i in 0..3 {
        for j in 0..3 {
            pcq = ctx.add(&[
                pcq,
                ctx.mul(&[p[i].clone(), cpq[3 * j + i].clone(), q[j].clone()]),
            ]);
        }
    }
    h = ctx.add(&[h, pcq]);
    // ½ q·Cqq·q
    let mut qcq = ctx.int(0);
    for i in 0..3 {
        for j in 0..3 {
            let w = if i == j {
                ctx.rational(1, 2).unwrap()
            } else {
                ctx.int(1)
            };
            qcq = ctx.add(&[
                qcq,
                ctx.mul(&[w, q[i].clone(), cqq[3 * j + i].clone(), q[j].clone()]),
            ]);
        }
    }
    h = ctx.add(&[h, qcq]);

    let h = ctx.expand(&h);
    ctx.monomial_coeffs(&h, &["q1", "q2", "q3", "p1", "p2", "p3"])
        .expect("二次型应能提取")
}

/// 三体势能：把 Legendre 表逐体替换 `(rx,ry,rz,r0)→(r*_*, …, r*0)` 并乘 `-mu_*`。
fn three_body_terms(ctx: &Context, n_max: u32) -> Vec<(Vec<u32>, Expr)> {
    let le = legendre_table(ctx, n_max);
    let mut out = Vec::with_capacity(le.len() * 3);
    for (mu, names) in [
        ("mu_e", ["rex", "rey", "rez", "re0"]),
        ("mu_m", ["rmx", "rmy", "rmz", "rm0"]),
        ("mu_s", ["rsx", "rsy", "rsz", "rs0"]),
    ] {
        let rx = ctx.sym(names[0]);
        let ry = ctx.sym(names[1]);
        let rz = ctx.sym(names[2]);
        let r0 = ctx.sym(names[3]);
        let m = ctx.sym(mu);
        let neg_m = ctx.mul(&[ctx.int(-1), m]);
        for (pow, coef) in &le {
            let c = ctx.subst(
                coef,
                &[
                    ("rx", rx.clone()),
                    ("ry", ry.clone()),
                    ("rz", rz.clone()),
                    ("r0", r0.clone()),
                ],
            );
            let c = ctx.mul(&[neg_m.clone(), c]);
            out.push((pow.clone(), c));
        }
    }
    out
}

/// 全表：二次型 + 三体，按幂次合并（系数相加）并按幂次升序排序。
pub fn hamilton_table(ctx: &Context, n_max: u32) -> Vec<(Vec<u32>, Expr)> {
    let mut all = quadratic_terms(ctx);
    all.extend(three_body_terms(ctx, n_max));

    // 按幂次分组求和（poly_simplify 的对应物）
    let mut groups: std::collections::BTreeMap<Vec<u32>, Vec<Expr>> =
        std::collections::BTreeMap::new();
    for (pow, c) in all {
        groups.entry(pow).or_default().push(c);
    }
    let mut out = Vec::with_capacity(groups.len());
    for (pow, terms) in groups {
        let coef = ctx.add(&terms);
        out.push((pow, coef));
    }
    out
}

pub fn run(args: &[String]) -> ExitCode {
    let n: u32 = opt(args, "--n", "15").parse().unwrap_or(15);
    let out = opt(args, "--out", "");
    let eval_at = opt(args, "--eval", "");

    let ctx = Context::new();
    let table = hamilton_table(&ctx, n);

    let mut lines = Vec::with_capacity(table.len());
    for (pow, coef) in &table {
        let head = pow
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let body = if eval_at.is_empty() {
            ctx.inspect(|i| i.dump(coef.raw_id()))
        } else {
            // 12 个数：rex,rey,rez, rmx,rmy,rmz, rsx,rsy,rsz, mu_e,mu_m,mu_s
            let v: Vec<f64> = eval_at
                .split(',')
                .map(|s| s.trim().parse::<f64>().unwrap_or(0.0))
                .collect();
            let g = |i: usize| v.get(i).copied().unwrap_or(0.0);
            let norm = |a: (f64, f64, f64)| (a.0 * a.0 + a.1 * a.1 + a.2 * a.2).sqrt();
            let re = (g(0), g(1), g(2));
            let rm = (g(3), g(4), g(5));
            let rs = (g(6), g(7), g(8));
            let vals: Vec<(&str, f64)> = vec![
                ("rex", re.0),
                ("rey", re.1),
                ("rez", re.2),
                ("re0", norm(re)),
                ("rmx", rm.0),
                ("rmy", rm.1),
                ("rmz", rm.2),
                ("rm0", norm(rm)),
                ("rsx", rs.0),
                ("rsy", rs.1),
                ("rsz", rs.2),
                ("rs0", norm(rs)),
                ("mu_e", g(9)),
                ("mu_m", g(10)),
                ("mu_s", g(11)),
            ];
            match ctx.eval_float(coef, &vals) {
                Some(x) => format!("{x:.17e}"),
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
