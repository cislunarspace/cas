//! oracle — sympy 差分 harness（M2 硬门槛：随机 expand 语义一致）。
//!
//! 协议（制表符分隔，与 cas/oracle/sympy_expand.py 对应；plain 输出无制表符，
//! 按列切分安全）：
//!   请求行: E\t<i>\t<plain 表达式>\t<x y z w 的有理点，形如 3/7 -2/5 ...>
//!   响应行: V\t<i>\t<sympy expand 后在该点的值 num/den>；非有限值输出 ERR
//!
//! 响应行序与请求行序一致（管道顺序），按序归位到 (case, point)。
//! 判据：随机 ℚ 点交叉精确求值（设计 §8）——our expand 与 sympy expand
//! 在同点值相等 ⇔ 语义一致；字符串形态差异不作失败。

use crate::util::Lcg;
use cas::prelude::*;
use cas_domain::{Integer, Rational};
use std::io::{BufRead, BufReader, Write};
use std::process::{Command, ExitCode, Stdio};
use std::time::Instant;

const SCRIPT: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../oracle/sympy_expand.py");
const VARS: [&str; 4] = ["x", "y", "z", "w"];
const PTS_PER_CASE: usize = 3;

fn rand_rat(rng: &mut Lcg) -> Rational {
    let n = rng.below(19) as i64 - 9;
    let d = rng.below(9) as i64 + 1;
    Rational::from_ints(&Integer::from_i64(n), &Integer::from_i64(d)).unwrap()
}

/// 精确表达式生成（无浮点/函数——oracle 语料限定 int/rat/sym）。
fn gen_expr(ctx: &Context, rng: &mut Lcg, depth: u32) -> Expr {
    if depth == 0 || rng.below(3) == 0 {
        match rng.below(3) {
            0 => ctx.int(rng.below(19) as i64 - 9),
            1 => ctx
                .rational(rng.below(15) as i64 - 7, rng.below(8) as i64 + 2)
                .unwrap(),
            _ => ctx.sym(VARS[rng.below(VARS.len() as u64) as usize]),
        }
    } else {
        match rng.below(5) {
            0 | 1 => gen_expr(ctx, rng, depth - 1) + gen_expr(ctx, rng, depth - 1),
            2 => gen_expr(ctx, rng, depth - 1) * gen_expr(ctx, rng, depth - 1),
            3 => gen_expr(ctx, rng, depth - 1) / gen_expr(ctx, rng, depth - 1),
            _ => gen_expr(ctx, rng, depth - 1).pow(rng.below(7) as i64 - 3),
        }
    }
}

fn parse_nd(s: &str) -> Option<Rational> {
    let (n, d) = s.split_once('/')?;
    let n = Integer::parse(n)?;
    let d = Integer::parse(d)?;
    Rational::from_ints(&n, &d)
}

pub(crate) fn run(args: &[String]) -> ExitCode {
    let cases: usize = crate::arg(args, "--cases", "1000").parse().unwrap_or(1000);
    let seed: u64 = crate::arg(args, "--seed", "42").parse().unwrap_or(42);
    let ctx = Context::new();
    let mut rng = Lcg::new(seed);

    // 生成语料与点；点须使原式可精确求值（避开 0 的负幂）
    let mut corpus: Vec<(Expr, Vec<[Rational; 4]>)> = Vec::with_capacity(cases);
    let mut texts: Vec<String> = Vec::with_capacity(cases);
    while corpus.len() < cases {
        let e = gen_expr(&ctx, &mut rng, 4);
        let mut pts = Vec::new();
        'pts: for _ in 0..PTS_PER_CASE {
            for _ in 0..8 {
                let p: Vec<Rational> = VARS.iter().map(|_| rand_rat(&mut rng)).collect();
                let vals: Vec<(&str, Rational)> =
                    VARS.iter().zip(&p).map(|(n, v)| (*n, v.clone())).collect();
                if ctx.eval_rational(&e, &vals).is_some() {
                    pts.push([p[0].clone(), p[1].clone(), p[2].clone(), p[3].clone()]);
                    continue 'pts;
                }
            }
            break 'pts; // 该式取不到可求值点，弃式重生成
        }
        if pts.len() == PTS_PER_CASE {
            texts.push(cas::plain(&ctx, &e));
            corpus.push((e, pts));
        }
    }

    // 我方：expand 后逐点精确求值
    let t = Instant::now();
    let mut our_vals: Vec<Vec<Rational>> = Vec::with_capacity(corpus.len());
    for (e, pts) in &corpus {
        let g = ctx.expand(e);
        let mut row = Vec::with_capacity(pts.len());
        for p in pts {
            let vals: Vec<(&str, Rational)> = VARS
                .iter()
                .zip(p.iter())
                .map(|(n, v)| (*n, v.clone()))
                .collect();
            row.push(ctx.eval_rational(&g, &vals).expect("点已预检可求值"));
        }
        our_vals.push(row);
    }
    let d_ours = t.elapsed();

    // sympy oracle
    let mut child = Command::new("python3")
        .arg(SCRIPT)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap_or_else(|e| {
            eprintln!("启动 python3 {SCRIPT} 失败: {e}（需要 sympy）");
            std::process::exit(2);
        });
    let mut stdin = child.stdin.take().expect("stdin");
    let stdout = child.stdout.take().expect("stdout");

    let mut req = String::new();
    for (i, (text, pts)) in texts.iter().zip(corpus.iter()).enumerate() {
        for p in &pts.1 {
            let pt_str = p
                .iter()
                .map(|r| format!("{}/{}", r.num(), r.den()))
                .collect::<Vec<_>>()
                .join(" ");
            req.push_str("E\t");
            req.push_str(&i.to_string());
            req.push('\t');
            req.push_str(text);
            req.push('\t');
            req.push_str(&pt_str);
            req.push('\n');
        }
    }
    let t = Instant::now();
    // 请求量可远超管道缓冲，写入放独立线程，避免与响应读取互相阻塞
    let writer = std::thread::spawn(move || {
        stdin.write_all(req.as_bytes()).expect("写 oracle 请求");
    });

    let expected = cases * PTS_PER_CASE;
    let mut sympy_vals: Vec<Vec<Option<Rational>>> = vec![Vec::new(); cases];
    let mut line_no = 0usize;
    for line in BufReader::new(stdout).lines() {
        let line = line.expect("读 oracle 响应");
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() != 3 || cols[0] != "V" {
            continue;
        }
        let case = line_no / PTS_PER_CASE;
        let val = if cols[2] == "ERR" {
            None
        } else {
            parse_nd(cols[2])
        };
        if let Some(v) = sympy_vals.get_mut(case) {
            v.push(val);
        }
        line_no += 1;
        if line_no == expected {
            break;
        }
    }
    let d_sympy = t.elapsed();
    writer.join().expect("写线程");
    let _ = child.wait();

    // 语义比较
    let (mut agree, mut mismatch, mut skip) = (0usize, 0usize, 0usize);
    for (i, ours) in our_vals.iter().enumerate() {
        for (k, ov) in ours.iter().enumerate() {
            match sympy_vals
                .get(i)
                .and_then(|v| v.get(k))
                .and_then(|o| o.as_ref())
            {
                Some(sv) => {
                    if sv == ov {
                        agree += 1;
                    } else {
                        mismatch += 1;
                        if mismatch <= 3 {
                            eprintln!(
                                "语义不一致 case#{i} 点{k}: ours={ov:?} sympy={sv:?}\n  expr: {}",
                                texts[i]
                            );
                        }
                    }
                }
                None => skip += 1,
            }
        }
    }

    println!("\n── sympy expand 语义对拍 ──");
    println!("语料      : {cases} 条 × {PTS_PER_CASE} 点（seed {seed}）");
    println!("一致      : {agree}；不一致: {mismatch}；跳过: {skip}");
    println!(
        "耗时      : 我方 expand+精确求值 {d_ours:?}（{:.0} 条/s）；sympy {d_sympy:?}（{:.0} 条/s）",
        cases as f64 / d_ours.as_secs_f64().max(1e-9),
        cases as f64 / d_sympy.as_secs_f64().max(1e-9),
    );
    if mismatch == 0 && skip == 0 {
        println!("判定      : 通过（语义一致率 100%）");
        ExitCode::SUCCESS
    } else {
        println!("判定      : 失败（mismatch={mismatch}, skip={skip}）");
        ExitCode::FAILURE
    }
}
