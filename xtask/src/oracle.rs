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

/// 单变元随机积：2–3 个 [-6,6] 整系数因子（次数 1..=6）× 有理内容。
fn gen_univar_product(ctx: &Context, rng: &mut Lcg) -> Expr {
    let x = ctx.sym("x");
    let nf = 2 + rng.below(2);
    let mut acc = ctx.int(1);
    for _ in 0..nf {
        let deg = 1 + rng.below(6);
        let mut f = ctx.int(0);
        for k in 0..=deg {
            let c = ctx.int(rng.below(13) as i64 - 6);
            f = f + c * x.clone().pow(k as i64);
        }
        acc = acc * f;
    }
    let den = 1 + rng.below(5) as i64;
    let num = rng.below(9) as i64 - 4;
    acc * ctx.rational(num, den).unwrap()
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
    let op = crate::arg(args, "--op", "expand");
    let op_cancel = op == "cancel";
    let op_factor = op == "factor";
    let tag = if op_cancel {
        "C"
    } else if op_factor {
        "F"
    } else {
        "E"
    };
    let ctx = Context::new();
    let mut rng = Lcg::new(seed);

    // 生成语料与点；点须使原式可精确求值（避开 0 的负幂）
    let mut corpus: Vec<(Expr, Vec<[Rational; 4]>)> = Vec::with_capacity(cases);
    let mut texts: Vec<String> = Vec::with_capacity(cases);
    while corpus.len() < cases {
        let e = if op_factor {
            gen_univar_product(&ctx, &mut rng)
        } else {
            gen_expr(&ctx, &mut rng, 4)
        };
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
        let g = if op_cancel {
            ctx.cancel(e)
        } else if op_factor {
            ctx.factor(e)
        } else {
            ctx.expand(e)
        };
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
        if op_factor {
            // factor 对拍只需表达式；点列以 "-" 占位
            req.push_str(&format!("F\t{i}\t{text}\t-\n"));
            continue;
        }
        for p in &pts.1 {
            let pt_str = p
                .iter()
                .map(|r| format!("{}/{}", r.num(), r.den()))
                .collect::<Vec<_>>()
                .join(" ");
            req.push_str(tag);
            req.push('\t');
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

    let per_case: usize = if op_factor { 1 } else { PTS_PER_CASE };
    let expected = cases * per_case;
    let mut sympy_factor: Vec<String> = vec![String::new(); cases];
    let mut sympy_vals: Vec<Vec<Option<Rational>>> = vec![Vec::new(); cases];
    let mut line_no = 0usize;
    for line in BufReader::new(stdout).lines() {
        let line = line.expect("读 oracle 响应");
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() != 3 || cols[0] != "V" {
            continue;
        }
        let case = line_no / per_case;
        if op_factor {
            if let Some(v) = sympy_factor.get_mut(case) {
                *v = cols[2].to_string();
            }
            line_no += 1;
            if line_no == expected {
                break;
            }
            continue;
        }
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

    // 比较：expand/cancel 走点值；factor 走 重数:次数 多重集 + 重展开点值
    let (mut agree, mut mismatch, mut skip) = (0usize, 0usize, 0usize);
    if op_factor {
        // 我方从 factor 结果读 (重数, 次数)：Mul 参数逐个归类
        let multiset_of = |g: &Expr| -> Option<String> {
            ctx.inspect(|i| {
                let mut ms: Vec<String> = Vec::new();
                collect_ms_pub(i, g.raw_id(), &mut ms).map(|()| {
                    ms.sort();
                    ms.join(",")
                })
            })
        };
        for (i, ((e, _), text)) in corpus.iter().zip(&texts).enumerate() {
            let g = ctx.factor(e);
            let ours = multiset_of(&g);
            match (ours, sympy_factor.get(i).map(String::as_str)) {
                (Some(o), Some(sy)) if !sy.is_empty() => {
                    if o == *sy {
                        agree += 1;
                    } else {
                        mismatch += 1;
                        if mismatch <= 3 {
                            eprintln!("多重集不一致 case#{i}: ours={o} sympy={sy}\n  expr: {text}");
                        }
                    }
                }
                _ => skip += 1,
            }
        }
    } else {
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
    }

    println!("\n── sympy {op} 语义对拍 ──");
    println!(
        "语料      : {cases} 条 × {}（seed {seed}）",
        if op_factor { "多重集" } else { "3 点" }
    );
    println!("一致      : {agree}；不一致: {mismatch}；跳过: {skip}");
    println!(
        "耗时      : 我方 {op} {d_ours:?}（{:.0} 条/s）；sympy {d_sympy:?}（{:.0} 条/s）",
        cases as f64 / d_ours.as_secs_f64().max(1e-9),
        cases as f64 / d_sympy.as_secs_f64().max(1e-9),
    );
    if mismatch == 0 && skip == 0 {
        println!("判定      : 通过（一致率 100%）");
        ExitCode::SUCCESS
    } else {
        println!("判定      : 失败（mismatch={mismatch}, skip={skip}）");
        ExitCode::FAILURE
    }
}

/// 从 factor 结果收集 (重数:次数) 多重集（oracle factor 模式用）。
fn collect_ms_pub(i: &cas_expr::Inspector<'_>, id: u32, ms: &mut Vec<String>) -> Option<()> {
    use cas_expr::Kind;
    fn deg_of(i: &cas_expr::Inspector<'_>, id: u32) -> Option<u64> {
        match i.kind(id) {
            Kind::Int(_) | Kind::Rat(_) | Kind::Float(_) => Some(0),
            Kind::Sym(_) => Some(1),
            Kind::Pow { base, exp } => {
                let b = deg_of(i, base)?;
                let e = match i.kind(exp) {
                    Kind::Int(v) => v.to_i64().map(|v| v as u64),
                    _ => None,
                }?;
                Some(b.saturating_mul(e))
            }
            Kind::Mul(args) => args.iter().map(|&a| deg_of(i, a)).sum::<Option<u64>>(),
            Kind::Add(args) => {
                let mut m: Option<u64> = None;
                for &a in args {
                    let d = deg_of(i, a)?;
                    m = Some(m.map_or(d, |v| v.max(d)));
                }
                m
            }
            Kind::Fn { .. } => None,
        }
    }
    match i.kind(id) {
        Kind::Mul(args) => {
            for &a in args {
                collect_ms_pub(i, a, ms)?;
            }
            Some(())
        }
        Kind::Pow { base, exp } => {
            let d = deg_of(i, base)?;
            let m = match i.kind(exp) {
                Kind::Int(v) => v.to_i64(),
                _ => None,
            }?;
            ms.push(format!("{m}:{d}"));
            Some(())
        }
        Kind::Sym(_) => {
            ms.push("1:1".to_string());
            Some(())
        }
        Kind::Int(_) | Kind::Rat(_) | Kind::Float(_) => Some(()),
        Kind::Add(_) => {
            ms.push(format!("1:{}", deg_of(i, id)?));
            Some(())
        }
        Kind::Fn { .. } => None,
    }
}
