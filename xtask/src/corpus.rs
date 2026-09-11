//! corpus — 真实语料分析（hpoly_generated.rs，只读）。
//!
//! 基础档：转写 → 解析 → 往返/跨上下文/确定性/吞吐。
//! --eval：数值求值交叉验证——对原文跑一个**独立实现**的 f64 求值器，
//! 与 cas 的 parse→eval、parse→print→parse→eval 两条通道对拍。

use crate::util::{Lcg, fnv64};
use cas::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::process::ExitCode;
use std::time::Instant;

pub(crate) fn run(args: &[String]) -> ExitCode {
    let file = crate::arg(
        args,
        "--file",
        "../dynamics/normalization/rust/src/hpoly_generated.rs",
    );
    let do_eval = crate::flag(args, "--eval");
    let src = fs::read_to_string(&file).unwrap_or_else(|e| {
        eprintln!("读取语料失败: {file}: {e}\n（在 cas/ 目录下运行，或用 --file 指定绝对路径）");
        std::process::exit(2);
    });
    let (params, raw) = extract(&src);
    let corpus: Vec<String> = raw.iter().map(|s| transcribe(s)).collect();
    if corpus.is_empty() {
        eprintln!("未提取到任何 h[i] 表达式：{file}");
        return ExitCode::from(2);
    }
    println!("语料: {file}");
    println!("参数绑定: {} 个；表达式: {} 条", params.len(), corpus.len());

    let ctx = Context::new();
    let t = Instant::now();
    let mut fails: Vec<(usize, String, &str)> = Vec::new();
    let mut exprs = Vec::with_capacity(corpus.len());
    for (i, e) in corpus.iter().enumerate() {
        match cas::parse(&ctx, e) {
            Ok(x) => exprs.push((i, x)),
            Err(err) => fails.push((i, err.to_string(), e.as_str())),
        }
    }
    let d_parse = t.elapsed();
    let parsed = exprs.len();

    let t = Instant::now();
    let texts: Vec<(usize, String)> = exprs
        .iter()
        .map(|(i, e)| (*i, cas::plain(&ctx, e)))
        .collect();
    let d_print = t.elapsed();

    let t = Instant::now();
    let mut rt_fail: Vec<usize> = Vec::new();
    for ((i, e), (_, txt)) in exprs.iter().zip(&texts) {
        match cas::parse(&ctx, txt) {
            Ok(back) => {
                if !(back == *e) {
                    rt_fail.push(*i);
                }
            }
            Err(err) => {
                eprintln!("往返解析失败 #{i}: {err} | {txt}");
                rt_fail.push(*i);
            }
        }
    }
    let d_rt = t.elapsed();

    let mut xc_fail: Vec<usize> = Vec::new();
    let ctx2 = Context::new();
    for (i, txt) in &texts {
        match cas::parse(&ctx2, txt) {
            Ok(e2) => {
                if cas::plain(&ctx2, &e2) != *txt {
                    xc_fail.push(*i);
                }
            }
            Err(_) => xc_fail.push(*i),
        }
    }

    let orig_chars: usize = corpus.iter().map(|s| s.len()).sum();
    let canon_chars: usize = texts.iter().map(|(_, t)| t.len()).sum();
    let digest = fnv64(&texts.iter().map(|(_, t)| t.as_str()).collect::<String>());

    println!("\n── M1 档：解析/往返/确定性 ──");
    println!(
        "解析      : 成功 {parsed}/{}, 失败 {}",
        corpus.len(),
        fails.len()
    );
    for (i, msg, e) in fails.iter().take(5) {
        println!("  失败 #{i}: {msg} | {e}");
    }
    println!(
        "往返      : parse∘print 同节点 {} 条通过, {} 条失败",
        parsed - rt_fail.len(),
        rt_fail.len()
    );
    println!(
        "跨上下文  : plain 字节一致 {} 条通过, {} 条失败",
        parsed - xc_fail.len(),
        xc_fail.len()
    );
    println!(
        "字符量    : 原文 {orig_chars} → 规范形 {canon_chars} 字节（{:.1}%）",
        100.0 * canon_chars as f64 / orig_chars.max(1) as f64
    );
    println!(
        "arena     : 唯一节点 {} 个（hash-consing 跨表达式共享）",
        ctx.node_count()
    );
    println!(
        "耗时      : 解析 {d_parse:?}（{:.0} 条/s）、打印 {d_print:?}（{:.0} 条/s）、往返 {d_rt:?}",
        parsed as f64 / d_parse.as_secs_f64().max(1e-9),
        parsed as f64 / d_print.as_secs_f64().max(1e-9),
    );
    println!("确定性摘要 FNV64(规范形全文): {digest:016x}");

    let mut ok = fails.is_empty() && rt_fail.is_empty() && xc_fail.is_empty();

    if do_eval {
        ok &= eval_cross_check(&ctx, &corpus, &exprs, &texts, &params);
    }

    if ok {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    }
}

// ── 数值交叉验证 ──────────────────────────────────────────────

fn eval_cross_check(
    ctx: &Context,
    corpus: &[String],
    exprs: &[(usize, Expr)],
    texts: &[(usize, String)],
    params: &[String],
) -> bool {
    println!("\n── --eval 档：数值求值交叉验证（独立求值器对拍）──");
    let mut rng = Lcg::new(20260911);
    let rounds = 5usize;
    let (mut max_rel, mut worst) = (0.0f64, String::new());
    let mut a_fail = 0usize; // 直接通道 (parse→eval) 违例
    let mut c_fail = 0usize; // 往返通道 (parse→print→parse→eval) 违例
    let tol: f64 = 1e-6;

    // 往返节点（canonical 文本再解析），供通道 C
    let roundtrip: Vec<(usize, Expr)> = texts
        .iter()
        .map(|(i, t)| (*i, cas::parse(ctx, t).expect("M1 档已保证可解析")))
        .collect();

    for _ in 0..rounds {
        let pt: Vec<(&str, f64)> = params
            .iter()
            .map(|p| (p.as_str(), param_value(p, &mut rng)))
            .collect();
        let mut pv: HashMap<&str, f64> = HashMap::new();
        for (k, v) in &pt {
            pv.insert(k, *v);
        }
        for (j, text) in corpus.iter().enumerate() {
            let b = match eval_indep(text, &pv) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("独立求值失败 #{j}: {e}");
                    a_fail += 1;
                    continue;
                }
            };
            let vals: Vec<(&str, f64)> = pt.clone();
            // 通道 A：cas 直接解析后求值
            if let Some(a) = ctx.eval_float(&exprs[j].1, &vals) {
                let rel = rel_diff(a, b);
                if rel > tol {
                    a_fail += 1;
                }
                if rel > max_rel {
                    max_rel = rel;
                    worst = format!("#{j} A {a:.6e} vs indep {b:.6e}");
                }
            } else {
                a_fail += 1;
            }
            // 通道 C：规范化往返后求值
            if let Some(c) = ctx.eval_float(&roundtrip[j].1, &vals) {
                let rel = rel_diff(c, b);
                if rel > tol {
                    c_fail += 1;
                }
                if rel > max_rel {
                    max_rel = rel;
                    worst = format!("#{j} C {c:.6e} vs indep {b:.6e}");
                }
            } else {
                c_fail += 1;
            }
        }
    }

    let total = rounds * corpus.len();
    println!(
        "点数      : {rounds} 组 × {} 条；容差 rel {tol:e}",
        corpus.len()
    );
    println!(
        "直接通道  : {} 违例；往返通道: {} 违例（共 {total} 次比较/通道）",
        a_fail, c_fail
    );
    println!("最大相对差: {max_rel:e}（{worst}）");
    a_fail == 0 && c_fail == 0
}

fn rel_diff(a: f64, b: f64) -> f64 {
    (a - b).abs() / (a.abs() + b.abs() + 1.0)
}

/// 现实量级的参数采样（按 hpoly 参数命名前缀）。
fn param_value(name: &str, rng: &mut Lcg) -> f64 {
    if name.starts_with("mu_") {
        rng.uniform(0.001, 1.0)
    } else if name.ends_with('0') {
        rng.uniform(0.05, 1.0) // re0/rm0/rs0
    } else if name.starts_with('r') {
        rng.uniform(-0.3, 0.3)
    } else if name.starts_with('f') {
        rng.uniform(0.001, 0.01)
    } else {
        rng.uniform(-0.5, 0.5) // Cpq*/Cqq*
    }
}

/// 独立 f64 求值器：对转写后的 plain 文本直接递归下降求值。
/// 与 cas-parse/cas-eval 完全独立的代码路径——交叉验证的意义所在。
/// 文法子集（语料实测）：字面量、36 标识符、+ - * /、^正整数、括号、一元负。
fn eval_indep(src: &str, vals: &HashMap<&str, f64>) -> Result<f64, String> {
    let cs: Vec<char> = src.chars().collect();
    let mut p = 0usize;
    let v = parse_expr(&cs, &mut p, vals)?;
    skip_ws(&cs, &mut p);
    if p != cs.len() {
        return Err(format!("多余字符@{p}: {}", cs[p]));
    }
    Ok(v)
}

fn skip_ws(cs: &[char], p: &mut usize) {
    while *p < cs.len() && cs[*p].is_whitespace() {
        *p += 1;
    }
}

fn parse_expr(cs: &[char], p: &mut usize, vals: &HashMap<&str, f64>) -> Result<f64, String> {
    let mut l = parse_term(cs, p, vals)?;
    loop {
        skip_ws(cs, p);
        match cs.get(*p) {
            Some('+') => {
                *p += 1;
                l += parse_term(cs, p, vals)?;
            }
            Some('-') => {
                *p += 1;
                l -= parse_term(cs, p, vals)?;
            }
            _ => return Ok(l),
        }
    }
}

fn parse_term(cs: &[char], p: &mut usize, vals: &HashMap<&str, f64>) -> Result<f64, String> {
    let mut l = parse_unary(cs, p, vals)?;
    loop {
        skip_ws(cs, p);
        match cs.get(*p) {
            Some('*') => {
                *p += 1;
                l *= parse_unary(cs, p, vals)?;
            }
            Some('/') => {
                *p += 1;
                let d = parse_unary(cs, p, vals)?;
                if d == 0.0 {
                    return Err("除零".into());
                }
                l /= d;
            }
            _ => return Ok(l),
        }
    }
}

fn parse_unary(cs: &[char], p: &mut usize, vals: &HashMap<&str, f64>) -> Result<f64, String> {
    skip_ws(cs, p);
    match cs.get(*p) {
        Some('-') => {
            *p += 1;
            Ok(-parse_unary(cs, p, vals)?)
        }
        Some('+') => {
            *p += 1;
            parse_unary(cs, p, vals)
        }
        _ => parse_postfix(cs, p, vals),
    }
}

fn parse_postfix(cs: &[char], p: &mut usize, vals: &HashMap<&str, f64>) -> Result<f64, String> {
    let a = parse_atom(cs, p, vals)?;
    skip_ws(cs, p);
    if cs.get(*p) == Some(&'^') {
        *p += 1;
        let e = parse_unary(cs, p, vals)?;
        return Ok(a.powf(e));
    }
    Ok(a)
}

fn parse_atom(cs: &[char], p: &mut usize, vals: &HashMap<&str, f64>) -> Result<f64, String> {
    skip_ws(cs, p);
    match cs.get(*p) {
        Some('(') => {
            *p += 1;
            let v = parse_expr(cs, p, vals)?;
            skip_ws(cs, p);
            if cs.get(*p) != Some(&')') {
                return Err(format!("缺 ) @{p}"));
            }
            *p += 1;
            Ok(v)
        }
        Some(c) if c.is_ascii_digit() => {
            let start = *p;
            while *p < cs.len() && (cs[*p].is_ascii_digit() || cs[*p] == '.') {
                *p += 1;
            }
            if matches!(cs.get(*p), Some('e') | Some('E')) {
                *p += 1;
                if matches!(cs.get(*p), Some('+') | Some('-')) {
                    *p += 1;
                }
                while *p < cs.len() && cs[*p].is_ascii_digit() {
                    *p += 1;
                }
            }
            let s: String = cs[start..*p].iter().collect();
            s.parse::<f64>().map_err(|e| e.to_string())
        }
        Some(c) if c.is_ascii_alphabetic() || *c == '_' => {
            let start = *p;
            while *p < cs.len() && (cs[*p].is_ascii_alphanumeric() || cs[*p] == '_') {
                *p += 1;
            }
            let name: String = cs[start..*p].iter().collect();
            vals.get(name.as_str())
                .copied()
                .ok_or_else(|| format!("未绑定参数 {name}"))
        }
        other => Err(format!("意外字符 @{p}: {other:?}")),
    }
}

// ── 语料提取与转写 ────────────────────────────────────────────

/// 从 ExportHpoly 生成的 Rust 源里提取 `let X = p[k];` 参数名与
/// `h[i] = <表达式>;`（支持跨行续接）。
fn extract(src: &str) -> (Vec<String>, Vec<String>) {
    let mut params: Vec<String> = Vec::new();
    let mut exprs: Vec<String> = Vec::new();
    let mut cur: Option<String> = None;
    for line in src.lines() {
        let l = line.trim();
        if l.starts_with("let ") && l.contains("p[") {
            if let Some(name) = l.split_whitespace().nth(1) {
                params.push(name.to_string());
            }
            continue;
        }
        let after_h = l
            .strip_prefix("h[")
            .and_then(|r| r.split_once(']'))
            .map(|(_, rest)| rest);
        if let Some(rest) = after_h {
            if let Some(prev) = cur.take() {
                exprs.push(prev);
            }
            let e = rest.trim_start().trim_start_matches('=').trim_start();
            cur = Some(e.to_string());
        } else if cur.is_some() {
            if let Some(c) = cur.as_mut() {
                c.push(' ');
                c.push_str(l);
            }
        }
        if cur.as_ref().is_some_and(|c| c.ends_with(';')) {
            let done = cur.take().unwrap();
            exprs.push(done.trim_end_matches(';').trim().to_string());
        }
    }
    if let Some(prev) = cur {
        exprs.push(prev);
    }
    (params, exprs)
}

/// Rust 形态 → cas plain 文法：`.powi(N)` → `^N`，其余透传。
pub(crate) fn transcribe(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if s.is_char_boundary(i) && s[i..].starts_with(".powi(") {
            let open = i + ".powi(".len();
            if let Some(rel) = s[open..].find(')') {
                let close = open + rel;
                out.push('^');
                out.push_str(&s[open..close]);
                i = close + 1;
                continue;
            }
        }
        let ch_len = s[i..].chars().next().map_or(1, char::len_utf8);
        out.push_str(&s[i..i + ch_len]);
        i += ch_len;
    }
    out
}
