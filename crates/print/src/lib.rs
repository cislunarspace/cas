//! 输出层：plain 与 LaTeX。
//!
//! plain 是往返格式（D6 承诺）：对任意规范形表达式，`parse(print(e)) == e`
//! 逐节点成立，输出字节完全由表达式全序决定。LaTeX 面向阅读，无往返承诺。

use cas_domain::{Integer, Rational};
use cas_expr::{Context, Expr, Inspector, Kind};

pub fn plain(ctx: &Context, e: &Expr) -> String {
    let mut out = String::new();
    ctx.inspect(|ins| fmt_plain(ins, e.raw_id(), &mut out, 0));
    out
}

pub fn latex(ctx: &Context, e: &Expr) -> String {
    let mut out = String::new();
    ctx.inspect(|ins| fmt_latex(ins, e.raw_id(), &mut out, 0));
    out
}

const MAX_DEPTH: u32 = 10_000;

fn guard(depth: u32) {
    assert!(depth <= MAX_DEPTH, "表达式嵌套超过 {MAX_DEPTH} 层");
}

// ── 数值格式化 ──────────────────────────────────────────────────

enum Num<'a> {
    I(&'a Integer),
    R(&'a Rational),
    F(f64),
}

fn num_of<'a>(k: &Kind<'a>) -> Option<Num<'a>> {
    match k {
        Kind::Int(v) => Some(Num::I(v)),
        Kind::Rat(r) => Some(Num::R(r)),
        Kind::Float(v) => Some(Num::F(*v)),
        _ => None,
    }
}

/// 浮点文本：保证含小数点或指数，可被解析器原样读回。
/// 大/小量级切换 `{:e}` 科学计数（Rust 的 Display 永不用指数，1e300 会打
/// 出 301 位数字）；Display 与 LowerExp 都是最短往返表示，切换无精度损失。
fn float_str(v: f64) -> String {
    let mag = v.abs();
    let mut s = if v != 0.0 && (mag >= 1e16 || mag < 1e-4) {
        format!("{v:e}")
    } else {
        format!("{v}")
    };
    if let Some(epos) = s.find('e') {
        if !s[..epos].contains('.') {
            s.insert_str(epos, ".0");
        }
    } else if !s.contains('.') {
        s.push_str(".0");
    }
    s
}

fn num_plain(n: &Num<'_>) -> String {
    match n {
        Num::I(v) => v.to_string(),
        Num::R(r) => r.to_string(),
        Num::F(v) => float_str(*v),
    }
}

fn num_is_neg(n: &Num<'_>) -> bool {
    match n {
        Num::I(v) => v.is_negative(),
        Num::R(r) => r.is_negative(),
        Num::F(v) => *v < 0.0, // arena 中 Float 恒非负，防御性保留
    }
}

fn num_abs_plain(n: &Num<'_>) -> String {
    let mut s = num_plain(n);
    if s.starts_with('-') {
        s.remove(0);
    }
    s
}

fn num_latex(n: &Num<'_>) -> String {
    match n {
        Num::I(v) => v.to_string(),
        Num::R(r) if r.is_negative() => {
            format!("-\\frac{{{}}}{{{}}}", r.num().neg(), r.den())
        }
        Num::R(r) => format!("\\frac{{{}}}{{{}}}", r.num(), r.den()),
        Num::F(v) => float_latex(*v),
    }
}

fn float_latex(v: f64) -> String {
    let s = float_str(v);
    match s.find('e') {
        Some(p) => {
            let (mant, exp) = s.split_at(p);
            format!("{mant}\\times10^{{{}}}", &exp[1..])
        }
        None => s,
    }
}

// ── plain ──────────────────────────────────────────────────────

fn fmt_plain(ins: &Inspector<'_>, id: u32, out: &mut String, depth: u32) {
    guard(depth);
    match ins.kind(id) {
        Kind::Int(_) | Kind::Rat(_) | Kind::Float(_) | Kind::Sym(_) => atom_plain(ins, id, out),
        Kind::Fn { head, args } => {
            out.push_str(head);
            out.push('(');
            for (i, &a) in args.iter().enumerate() {
                if i > 0 {
                    out.push_str(", ");
                }
                fmt_plain(ins, a, out, depth + 1);
            }
            out.push(')');
        }
        Kind::Pow { base, exp } => {
            let base_paren = match ins.kind(base) {
                Kind::Rat(_) | Kind::Mul(_) | Kind::Add(_) | Kind::Pow { .. } => true,
                k => num_of(&k).is_some_and(|n| num_is_neg(&n)),
            };
            if base_paren {
                out.push('(');
                fmt_plain(ins, base, out, depth + 1);
                out.push(')');
            } else {
                fmt_plain(ins, base, out, depth + 1);
            }
            out.push('^');
            // 指数裸写仅限单值类（整数/非负浮点/符号/函数）；其余括起
            let exp_paren = !matches!(
                ins.kind(exp),
                Kind::Int(_) | Kind::Float(_) | Kind::Sym(_) | Kind::Fn { .. }
            );
            if exp_paren {
                out.push('(');
                fmt_plain(ins, exp, out, depth + 1);
                out.push(')');
            } else {
                fmt_plain(ins, exp, out, depth + 1);
            }
        }
        Kind::Mul(args) => mul_plain(ins, args, out, depth),
        Kind::Add(args) => add_plain(ins, args, out, depth),
    }
}

fn atom_plain(ins: &Inspector<'_>, id: u32, out: &mut String) {
    match ins.kind(id) {
        Kind::Int(v) => out.push_str(&v.to_string()),
        Kind::Rat(r) => out.push_str(&r.to_string()),
        Kind::Float(v) => out.push_str(&float_str(v)),
        Kind::Sym(n) => out.push_str(n),
        _ => unreachable!("atom"),
    }
}

/// 拆出 Mul 的精确数值系数（浮点因子不算系数，按普通因子打印）。
fn split_coeff<'a>(ins: &'a Inspector<'_>, args: &'a [u32]) -> (Option<u32>, &'a [u32]) {
    match ins.kind(args[0]) {
        Kind::Int(_) | Kind::Rat(_) => (Some(args[0]), &args[1..]),
        _ => (None, args),
    }
}

fn coeff_is_neg_one(ins: &Inspector<'_>, c: u32) -> bool {
    matches!(ins.kind(c), Kind::Int(v) if v.to_i64() == Some(-1))
}

fn factors_plain(ins: &Inspector<'_>, fs: &[u32], out: &mut String, depth: u32) {
    for (i, &f) in fs.iter().enumerate() {
        if i > 0 {
            out.push('*');
        }
        if matches!(ins.kind(f), Kind::Add(_)) {
            out.push('(');
            fmt_plain(ins, f, out, depth + 1);
            out.push(')');
        } else {
            fmt_plain(ins, f, out, depth + 1);
        }
    }
}

fn mul_plain(ins: &Inspector<'_>, args: &[u32], out: &mut String, depth: u32) {
    let (coeff, rest) = split_coeff(ins, args);
    match coeff {
        None => factors_plain(ins, rest, out, depth),
        Some(c) => {
            let n = num_of(&ins.kind(c)).expect("已判 Int/Rat");
            if rest.is_empty() {
                out.push_str(&num_plain(&n));
            } else if coeff_is_neg_one(ins, c) {
                out.push('-');
                factors_plain(ins, rest, out, depth);
            } else {
                out.push_str(&num_plain(&n));
                out.push('*');
                factors_plain(ins, rest, out, depth);
            }
        }
    }
}

/// 项是否负系数领衔（仅精确数值参与符号合并；浮点恒非负）。
fn neg_led(ins: &Inspector<'_>, a: u32) -> bool {
    let k = ins.kind(a);
    if let Some(n) = num_of(&k) {
        return num_is_neg(&n);
    }
    if let Kind::Mul(args) = k {
        if let Some(c) = split_coeff(ins, args).0 {
            if let Some(n) = num_of(&ins.kind(c)) {
                return num_is_neg(&n);
            }
        }
    }
    false
}

/// 负系数项按幅值打印（符号由外层 " - " 承担）；正项整项打印。
fn term_body_plain(ins: &Inspector<'_>, a: u32, out: &mut String, depth: u32) {
    match ins.kind(a) {
        k if num_of(&k).is_some() => {
            let n = num_of(&k).expect("数值");
            if num_is_neg(&n) {
                out.push_str(&num_abs_plain(&n));
            } else {
                out.push_str(&num_plain(&n));
            }
        }
        Kind::Mul(args) => {
            let (coeff, rest) = split_coeff(ins, args);
            match coeff {
                Some(c) if num_of(&ins.kind(c)).is_some_and(|n| num_is_neg(&n)) => {
                    let n = num_of(&ins.kind(c)).expect("数值");
                    if coeff_is_neg_one(ins, c) {
                        factors_plain(ins, rest, out, depth);
                    } else {
                        out.push_str(&num_abs_plain(&n));
                        out.push('*');
                        factors_plain(ins, rest, out, depth);
                    }
                }
                _ => fmt_plain(ins, a, out, depth),
            }
        }
        _ => fmt_plain(ins, a, out, depth),
    }
}

fn add_plain(ins: &Inspector<'_>, args: &[u32], out: &mut String, depth: u32) {
    for (i, &a) in args.iter().enumerate() {
        let neg = neg_led(ins, a);
        if i == 0 {
            if neg {
                out.push('-');
            }
        } else {
            out.push_str(if neg { " - " } else { " + " });
        }
        term_body_plain(ins, a, out, depth);
    }
}

// ── latex ──────────────────────────────────────────────────────

fn fmt_latex(ins: &Inspector<'_>, id: u32, out: &mut String, depth: u32) {
    guard(depth);
    match ins.kind(id) {
        Kind::Int(v) => out.push_str(&v.to_string()),
        Kind::Rat(r) => out.push_str(&num_latex(&Num::R(r))),
        Kind::Float(v) => out.push_str(&float_latex(v)),
        Kind::Sym(n) => out.push_str(n),
        Kind::Fn { head, args } => match head {
            "sqrt" if args.len() == 1 => {
                out.push_str("\\sqrt{");
                fmt_latex(ins, args[0], out, depth + 1);
                out.push('}');
            }
            "abs" if args.len() == 1 => {
                out.push_str("\\left|");
                fmt_latex(ins, args[0], out, depth + 1);
                out.push_str("\\right|");
            }
            "sin" | "cos" | "tan" | "exp" | "log" => {
                out.push_str(&format!("\\{head}\\left("));
                for (i, &a) in args.iter().enumerate() {
                    if i > 0 {
                        out.push_str(", ");
                    }
                    fmt_latex(ins, a, out, depth + 1);
                }
                out.push_str("\\right)");
            }
            _ => {
                out.push_str(&format!("\\operatorname{{{head}}}\\left("));
                for (i, &a) in args.iter().enumerate() {
                    if i > 0 {
                        out.push_str(", ");
                    }
                    fmt_latex(ins, a, out, depth + 1);
                }
                out.push_str("\\right)");
            }
        },
        Kind::Pow { base, exp } => {
            let base_paren = match ins.kind(base) {
                Kind::Rat(_) | Kind::Mul(_) | Kind::Add(_) | Kind::Pow { .. } => true,
                k => num_of(&k).is_some_and(|n| num_is_neg(&n)),
            };
            if base_paren {
                out.push_str("\\left(");
                fmt_latex(ins, base, out, depth + 1);
                out.push_str("\\right)");
            } else {
                fmt_latex(ins, base, out, depth + 1);
            }
            out.push_str("^{");
            fmt_latex(ins, exp, out, depth + 1);
            out.push('}');
        }
        Kind::Mul(args) => mul_latex(ins, args, out, depth),
        Kind::Add(args) => add_latex(ins, args, out, depth),
    }
}

fn factors_latex(ins: &Inspector<'_>, fs: &[u32], out: &mut String, depth: u32) {
    for (i, &f) in fs.iter().enumerate() {
        if i > 0 {
            out.push_str(" \\cdot ");
        }
        if matches!(ins.kind(f), Kind::Add(_)) {
            out.push_str("\\left(");
            fmt_latex(ins, f, out, depth + 1);
            out.push_str("\\right)");
        } else {
            fmt_latex(ins, f, out, depth + 1);
        }
    }
}

fn mul_latex(ins: &Inspector<'_>, args: &[u32], out: &mut String, depth: u32) {
    let (coeff, rest) = split_coeff(ins, args);
    match coeff {
        None => factors_latex(ins, rest, out, depth),
        Some(c) => {
            let n = num_of(&ins.kind(c)).expect("已判 Int/Rat");
            if rest.is_empty() {
                out.push_str(&num_latex(&n));
            } else if coeff_is_neg_one(ins, c) {
                out.push('-');
                factors_latex(ins, rest, out, depth);
            } else {
                out.push_str(&num_latex(&n));
                out.push_str(" \\cdot ");
                factors_latex(ins, rest, out, depth);
            }
        }
    }
}

fn term_body_latex(ins: &Inspector<'_>, a: u32, out: &mut String, depth: u32) {
    match ins.kind(a) {
        k if num_of(&k).is_some() => {
            let n = num_of(&k).expect("数值");
            if num_is_neg(&n) {
                let abs = match &n {
                    Num::I(v) => v.abs().to_string(),
                    Num::R(r) => r.abs().to_string(),
                    Num::F(v) => float_latex(*v),
                };
                out.push_str(&abs);
            } else {
                out.push_str(&num_latex(&n));
            }
        }
        Kind::Mul(args) => {
            let (coeff, rest) = split_coeff(ins, args);
            match coeff {
                Some(c) if num_of(&ins.kind(c)).is_some_and(|n| num_is_neg(&n)) => {
                    let n = num_of(&ins.kind(c)).expect("数值");
                    if coeff_is_neg_one(ins, c) {
                        factors_latex(ins, rest, out, depth);
                    } else {
                        let abs = match n {
                            Num::I(v) => v.abs().to_string(),
                            Num::R(r) => r.abs().to_string(),
                            Num::F(v) => float_latex(v),
                        };
                        out.push_str(&abs);
                        out.push_str(" \\cdot ");
                        factors_latex(ins, rest, out, depth);
                    }
                }
                _ => fmt_latex(ins, a, out, depth),
            }
        }
        _ => fmt_latex(ins, a, out, depth),
    }
}

fn add_latex(ins: &Inspector<'_>, args: &[u32], out: &mut String, depth: u32) {
    for (i, &a) in args.iter().enumerate() {
        let neg = neg_led(ins, a);
        if i == 0 {
            if neg {
                out.push('-');
            }
        } else {
            out.push_str(if neg { " - " } else { " + " });
        }
        term_body_latex(ins, a, out, depth);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ctx_xy() -> (Context, Expr, Expr) {
        let ctx = Context::new();
        let x = ctx.sym("x");
        let y = ctx.sym("y");
        (ctx, x, y)
    }

    #[test]
    fn plain_基础规范形() {
        let (ctx, x, y) = ctx_xy();
        assert_eq!(plain(&ctx, &(x.clone() + x.clone())), "2*x");
        assert_eq!(plain(&ctx, &(x.clone() * x.clone())), "x^2");
        assert_eq!(plain(&ctx, &(x.clone().pow(2) * x.clone().pow(3))), "x^5");
        assert_eq!(plain(&ctx, &(x.clone() + y.clone())), "x + y");
        assert_eq!(plain(&ctx, &(y.clone() - x.clone())), "y - x");
        assert_eq!(plain(&ctx, &(x.clone() - x.clone())), "0");
        assert_eq!(plain(&ctx, &(x.clone() * y.clone()).pow(2)), "x^2*y^2");
        assert_eq!(plain(&ctx, &(x.clone() + y.clone()).pow(2)), "(x + y)^2");
    }

    #[test]
    fn plain_数值() {
        let (ctx, x, _) = ctx_xy();
        assert_eq!(plain(&ctx, &ctx.rational(1, 3).unwrap()), "1/3");
        assert_eq!(plain(&ctx, &(x.clone() / ctx.int(2))), "1/2*x");
        assert_eq!(
            plain(&ctx, &(ctx.rational(3, 2).unwrap() * x.clone())),
            "3/2*x"
        );
        assert_eq!(
            plain(&ctx, &ctx.int(2).pow(100)),
            "1267650600228229401496703205376"
        );
        assert_eq!(plain(&ctx, &ctx.int(2).pow(-2)), "1/4");
        assert_eq!(plain(&ctx, &ctx.float(0.0)), "0.0");
        assert_eq!(plain(&ctx, &ctx.float(1e300)), "1.0e300");
        assert_eq!(plain(&ctx, &(ctx.float(1.5) * x.clone())), "1.5*x");
        assert_eq!(plain(&ctx, &(-ctx.float(2.5))), "-2.5");
    }

    #[test]
    fn plain_函数与幂() {
        let (ctx, x, y) = ctx_xy();
        let a = ctx.sym("a");
        let b = ctx.sym("b");
        assert_eq!(
            plain(&ctx, &ctx.call("sin", std::slice::from_ref(&x))),
            "sin(x)"
        );
        assert_eq!(
            plain(
                &ctx,
                &(ctx.call("sin", std::slice::from_ref(&x))
                    * ctx.call("sin", std::slice::from_ref(&x)))
            ),
            "sin(x)^2"
        );
        assert_eq!(
            plain(&ctx, &x.clone().pow(&(a.clone() + b.clone()))),
            "x^(a + b)"
        );
        assert_eq!(plain(&ctx, &x.clone().pow(&a).pow(2)), "x^(2*a)");
        assert_eq!(
            plain(&ctx, &ctx.pow(&ctx.rational(1, 2).unwrap(), &x.clone())),
            "(1/2)^x"
        );
        assert_eq!(plain(&ctx, &ctx.pow(&ctx.int(-2), &x.clone())), "(-2)^x");
        assert_eq!(
            plain(&ctx, &ctx.call("atan2", &[y.clone(), x.clone()])),
            "atan2(y, x)"
        );
    }

    #[test]
    fn latex_基础() {
        let (ctx, x, y) = ctx_xy();
        assert_eq!(latex(&ctx, &(x.clone() + y.clone())), "x + y");
        assert_eq!(
            latex(&ctx, &(x.clone() / ctx.int(2))),
            "\\frac{1}{2} \\cdot x"
        );
        assert_eq!(
            latex(&ctx, &(x.clone() + y.clone()).pow(2)),
            "\\left(x + y\\right)^{2}"
        );
        assert_eq!(
            latex(&ctx, &ctx.call("sqrt", std::slice::from_ref(&x))),
            "\\sqrt{x}"
        );
        assert_eq!(latex(&ctx, &(-x.clone())), "-x");
        assert_eq!(
            latex(&ctx, &(y.clone() - ctx.int(3) * x.clone())),
            "y - 3 \\cdot x"
        );
        assert_eq!(latex(&ctx, &ctx.float(1e300)), "1.0\\times10^{300}");
    }
}
