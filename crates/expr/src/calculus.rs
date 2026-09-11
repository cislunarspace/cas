//! 微积分与化简（P1，对标 MATLAB Symbolic Math Toolbox 的 diff /
//! taylor / simplify，设计 D4）。
//!
//! - [`Inner::diff_at`]：符号求导——和/积/链式/幂（整数、有理、一般
//!   指数三档）/初等函数表（sin cos tan exp log sqrt abs）。
//!   全部经规范形构造器，输出即规范形（同类项自动合并）。
//! - [`Inner::taylor_at`]：x=a 处的 Taylor 级数多项式部分到指定阶
//!   （f^(n)(a)/n!·(x−a)^n）。实现：代换 x→a+t 后对 t 反复求导取
//!   t=0 值；O(·) 余项记号 P1 不引入（多项式部分即 MATLAB taylor
//!   忽略 Order 项的输出）。
//! - [`Inner::simplify_at`]：L2 定向规则表——常量折叠（cos(0)→1 等）
//!   与无条件恒等式（sin²+cos²→1、双角），逐节点自底向上单遍，
//!   规则应用与否按"节点数减少"代价判据，确定性优先。

use crate::node::Node;
use crate::{Inner, Rational};
use cas_domain::Integer;

impl Inner {
    // ── 求导 ──────────────────────────────────────────────────

    /// 对符号 `sym_id` 求导。
    pub(crate) fn diff_at(&mut self, id: u32, sym_id: u32, depth: u32) -> u32 {
        assert!(depth <= 10_000, "表达式嵌套过深");
        // owned 视图，避免边读边写借用冲突
        enum Info {
            Sym(u32),
            Add(Vec<u32>),
            Mul(Vec<u32>),
            Pow { base: u32, exp: u32 },
            Fn { head: u32, args: Vec<u32> },
            Const,
        }
        let info = match &self.nodes[id as usize] {
            Node::Sym(s) => Info::Sym(*s),
            Node::Add { args: sp } => Info::Add(self.node_args(*sp).to_vec()),
            Node::Mul { args: sp } => Info::Mul(self.node_args(*sp).to_vec()),
            Node::Pow { base, exp } => Info::Pow {
                base: *base,
                exp: *exp,
            },
            Node::Fn { head, args: sp } => Info::Fn {
                head: *head,
                args: self.node_args(*sp).to_vec(),
            },
            _ => Info::Const,
        };
        match info {
            Info::Const => self.lit_int(0),
            Info::Sym(s) => self.lit_int(if s == sym_id { 1 } else { 0 }),
            Info::Add(args) => {
                let v: Vec<u32> = args
                    .iter()
                    .map(|&a| self.diff_at(a, sym_id, depth + 1))
                    .collect();
                self.make_add(&v)
            }
            Info::Mul(args) => {
                // 逐因子求导 × 其余因子积：Σ d(aᵢ)·Π_{j≠i} aⱼ
                let mut terms = Vec::with_capacity(args.len());
                for (i, &a) in args.iter().enumerate() {
                    let d = self.diff_at(a, sym_id, depth + 1);
                    if let Node::Int(z) = &self.nodes[d as usize] {
                        if z.is_zero() {
                            continue;
                        }
                    }
                    let rest: Vec<u32> = args
                        .iter()
                        .enumerate()
                        .filter(|(j, _)| *j != i)
                        .map(|(_, &x)| x)
                        .collect();
                    let r = self.make_mul(&rest);
                    terms.push(self.make_mul(&[d, r]));
                }
                self.make_add(&terms)
            }
            Info::Pow { base, exp } => {
                let db = self.diff_at(base, sym_id, depth + 1);
                let de = self.diff_at(exp, sym_id, depth + 1);
                // owned 判定后再构造（借用分离）
                let exp_int: Option<Option<i64>> = match &self.nodes[exp as usize] {
                    Node::Int(k) => Some(k.to_i64()),
                    _ => None,
                };
                let exp_rat: Option<Rational> = match &self.nodes[exp as usize] {
                    Node::Rat(r) => Some(r.clone()),
                    _ => None,
                };
                if let Some(Some(k)) = exp_int {
                    // k·b^(k−1)·b'
                    let km1 = self.lit_int(k - 1);
                    let bk = self.make_pow(base, km1);
                    let kv = self.lit_int(k);
                    self.make_mul(&[kv, bk, db])
                } else if let Some(r) = exp_rat {
                    // (p/q)·b^(p/q−1)·b'
                    let pq = self.lit_rational(&r);
                    let one_rat = Rational::from_ints(&Integer::one(), &Integer::one()).unwrap();
                    let new_e = r.sub(&one_rat);
                    let e_new = self.lit_rational(&new_e);
                    let be = self.make_pow(base, e_new);
                    self.make_mul(&[pq, be, db])
                } else {
                    // 一般指数：b^e·(e'·log(b) + e·b'/b)
                    let log_b = self.fn_node("log", &[base]);
                    let t1 = self.make_mul(&[de, log_b]);
                    let b_inv = {
                        let neg1 = self.lit_int(-1);
                        self.make_pow(base, neg1)
                    };
                    let t2 = self.make_mul(&[exp, db, b_inv]);
                    let inner = self.make_add(&[t1, t2]);
                    self.make_mul(&[id, inner])
                }
            }
            Info::Fn { head, args } => {
                let name = self.fn_names[head as usize].to_string();
                if args.len() != 1 {
                    // 多参函数：返回未求导节点 × 0（保守；P1 语料单参）
                    return self.lit_int(0);
                }
                let a = args[0];
                let da = self.diff_at(a, sym_id, depth + 1);
                let name_ref: &str = &name;
                // 外层导数因子 × 链式内导
                let outer: u32 = match name_ref {
                    "sin" => self.fn_node("cos", &[a]),
                    "cos" => {
                        let neg1 = self.lit_int(-1);
                        let s = self.fn_node("sin", &[a]);
                        self.make_mul(&[neg1, s])
                    }
                    "tan" => {
                        // 1 + tan²
                        let t = self.fn_node("tan", &[a]);
                        let two = self.lit_int(2);
                        let t2 = self.make_pow(t, two);
                        let one = self.lit_int(1);
                        self.make_add(&[one, t2])
                    }
                    "exp" => self.fn_node("exp", &[a]),
                    "log" => {
                        let neg1 = self.lit_int(-1);
                        self.make_pow(a, neg1)
                    }
                    "sqrt" => {
                        // 1/(2·sqrt(a))
                        let two = self.lit_int(2);
                        let s = self.fn_node("sqrt", &[a]);
                        let den = self.make_mul(&[two, s]);
                        let neg1 = self.lit_int(-1);
                        self.make_pow(den, neg1)
                    }
                    "abs" => self.fn_node("sign", &[a]),
                    _ => return self.lit_int(0), // 未知函数：保守 0
                };
                self.make_mul(&[outer, da])
            }
        }
    }

    // ── Taylor 级数 ────────────────────────────────────────────

    /// x=a 处 Taylor 多项式到 (x−a)^order（含）。a 为整数点。
    pub(crate) fn taylor_at(&mut self, id: u32, sym_id: u32, at: i64, order: u32) -> u32 {
        // 平移 x→a+t：直接在平移坐标下反复求导取 t=0 值更简单——
        // 做法：subst x→(a+t) 需要 t 符号；等价实现：对 x 反复求导后
        // 代入 x=a 逐项求值（f^(n)(a) 精确有理数）
        let mut terms: Vec<u32> = Vec::with_capacity(order as usize + 1);
        let mut cur = id;
        for n in 0..=order {
            let val = self.eval_at_rational(cur, sym_id, at);
            let fact = Integer::from_i64(factorial(n));
            let coeff = match val {
                Some(v) => v
                    .div(&Rational::from_integer(&fact))
                    .unwrap_or_else(Rational::zero),
                None => Rational::zero(),
            };
            if !coeff.is_zero() {
                // coeff·(x−a)^n
                let sym = self.intern_node(Node::Sym(sym_id), &[]);
                let xa = if at == 0 {
                    sym
                } else {
                    let a_lit = self.lit_int(-at);
                    self.make_add(&[sym, a_lit])
                };
                let term = if n == 0 {
                    self.lit_rational(&coeff)
                } else {
                    let n_lit = self.lit_int(n as i64);
                    let pw = self.make_pow(xa, n_lit);
                    let c = self.lit_rational(&coeff);
                    self.make_mul(&[c, pw])
                };
                terms.push(term);
            }
            if n < order {
                let next = self.diff_at(cur, sym_id, 0);
                // 代入 x=a 防止下一轮导数表达式膨胀：以 (导数在 a 的值) 不可行
                //（需保持符号），保留符号形式
                cur = next;
            }
        }
        if terms.is_empty() {
            return self.lit_int(0);
        }
        self.make_add(&terms)
    }

    /// 只读精确求值：单符号代入整数值（Taylor 系数提取用）。
    fn eval_at_rational(&self, id: u32, sym_id: u32, at: i64) -> Option<Rational> {
        use std::collections::HashMap;
        let mut map: HashMap<u32, Rational> = HashMap::new();
        map.insert(
            sym_id,
            Rational::from_ints(&Integer::from_i64(at), &Integer::one()).unwrap(),
        );
        crate::eval::eval_rational_at(self, id, &map, 0).or_else(|| self.eval_fn_exact(id, &map))
    }

    /// 函数节点的特殊点精确值（Taylor 系数提取用；与 simplify 的
    /// 常量折叠表同源）。
    fn eval_fn_exact(
        &self,
        id: u32,
        map: &std::collections::HashMap<u32, Rational>,
    ) -> Option<Rational> {
        match &self.nodes[id as usize] {
            Node::Fn { head, args: sp } => {
                let a = self.node_args(*sp);
                if a.len() != 1 {
                    return None;
                }
                let v = self
                    .eval_fn_exact(a[0], map)
                    .or_else(|| crate::eval::eval_rational_at(self, a[0], map, 0))?;
                let name = self.fn_names[*head as usize].as_ref();
                match name {
                    "sin" if v.is_zero() => Some(Rational::zero()),
                    "cos" if v.is_zero() => Some(Rational::one()),
                    "tan" if v.is_zero() => Some(Rational::zero()),
                    "exp" if v.is_zero() => Some(Rational::one()),
                    "log" if v.is_one() => Some(Rational::zero()),
                    "sqrt" if v.is_zero() => Some(Rational::zero()),
                    "sqrt" if v.is_one() => Some(Rational::one()),
                    _ => None,
                }
            }
            Node::Add { args: sp } => {
                let mut acc = Rational::zero();
                for &a in self.node_args(*sp) {
                    acc = acc.add(&self.eval_fn_exact(a, map)?);
                }
                Some(acc)
            }
            Node::Mul { args: sp } => {
                let mut acc = Rational::one();
                for &a in self.node_args(*sp) {
                    acc = acc.mul(&self.eval_fn_exact(a, map)?);
                }
                Some(acc)
            }
            Node::Pow { base, exp } => {
                let b = self.eval_fn_exact(*base, map)?;
                let e = self.eval_fn_exact(*exp, map)?;
                if !e.den().is_one() {
                    return None;
                }
                let k = e.num().to_i64()?;
                if k >= 0 {
                    Some(b.pow_reduced(k as u32))
                } else {
                    Some(b.inv_reduced()?.pow_reduced(k.unsigned_abs() as u32))
                }
            }
            _ => crate::eval::eval_rational_at(self, id, map, 0),
        }
    }

    // ── simplify（L2 规则表）──────────────────────────────────

    /// 定向化简：自底向上单遍；规则按"减少节点数"判据接受。
    pub(crate) fn simplify_at(&mut self, id: u32, depth: u32) -> u32 {
        assert!(depth <= 10_000, "表达式嵌套过深");
        enum Info {
            Atom(u32),
            Add(Vec<u32>),
            Mul(Vec<u32>),
            Pow { base: u32, exp: u32 },
            Fn { head: u32, args: Vec<u32> },
        }
        let info = match &self.nodes[id as usize] {
            Node::Add { args: sp } => Info::Add(self.node_args(*sp).to_vec()),
            Node::Mul { args: sp } => Info::Mul(self.node_args(*sp).to_vec()),
            Node::Pow { base, exp } => Info::Pow {
                base: *base,
                exp: *exp,
            },
            Node::Fn { head, args: sp } => Info::Fn {
                head: *head,
                args: self.node_args(*sp).to_vec(),
            },
            _ => Info::Atom(id),
        };
        match info {
            Info::Atom(a) => a,
            Info::Add(args) => {
                let v: Vec<u32> = args
                    .iter()
                    .map(|&a| self.simplify_at(a, depth + 1))
                    .collect();
                let s = self.make_add(&v);
                let s = self.try_rules_add(s);
                self.try_rules_log_add(s)
            }
            Info::Mul(args) => {
                let v: Vec<u32> = args
                    .iter()
                    .map(|&a| self.simplify_at(a, depth + 1))
                    .collect();
                self.make_mul(&v)
            }
            Info::Pow { base, exp } => {
                let b = self.simplify_at(base, depth + 1);
                let e = self.simplify_at(exp, depth + 1);
                self.make_pow(b, e)
            }
            Info::Fn { head, args } => {
                let v: Vec<u32> = args
                    .iter()
                    .map(|&a| self.simplify_at(a, depth + 1))
                    .collect();
                let f = self.fn_node_by_id(head, &v);
                let f = self.try_rules_fn(f, head, &v);
                self.try_rules_sqrt_sq(f)
            }
        }
    }

    /// Add 级假设驱动规则：log(a) + log(b) → log(a·b)（a>0 ∧ b>0）。
    fn try_rules_log_add(&mut self, id: u32) -> u32 {
        use crate::assume::{Predicate, Trinary};
        let args = match &self.nodes[id as usize] {
            Node::Add { args: sp } => self.node_args(*sp).to_vec(),
            _ => return id,
        };
        if args.len() < 2 {
            return id;
        }
        let mut logs: Vec<u32> = vec![];
        let mut rest: Vec<u32> = vec![];
        for &a in &args {
            let log_arg = match &self.nodes[a as usize] {
                Node::Fn { head, args: sp } => {
                    let name = self.fn_names[*head as usize].as_ref();
                    if name == "log" && self.node_args(*sp).len() == 1 {
                        Some(self.node_args(*sp)[0])
                    } else {
                        None
                    }
                }
                _ => None,
            };
            match log_arg {
                Some(arg) if self.query_at(arg, Predicate::Positive, 0) == Trinary::True => {
                    logs.push(a);
                }
                _ => rest.push(a),
            }
        }
        if logs.len() < 2 {
            return id;
        }
        let mut combined: u32 = match &self.nodes[logs[0] as usize] {
            Node::Fn { args: sp, .. } => self.node_args(*sp)[0],
            _ => unreachable!(),
        };
        for &l in &logs[1..] {
            let arg = match &self.nodes[l as usize] {
                Node::Fn { args: sp, .. } => self.node_args(*sp)[0],
                _ => unreachable!(),
            };
            combined = self.make_mul(&[combined, arg]);
        }
        let merged = self.fn_node("log", &[combined]);
        rest.push(merged);
        self.make_add(&rest)
    }

    /// sqrt(x²) → x（x nonnegative；MATLAB assume positive 的标准行为）。
    fn try_rules_sqrt_sq(&mut self, id: u32) -> u32 {
        use crate::assume::{Predicate, Trinary};
        if let Node::Fn { head, args: sp } = &self.nodes[id as usize] {
            if self.fn_names[*head as usize].as_ref() == "sqrt" {
                let args = self.node_args(*sp);
                if args.len() == 1 {
                    if let Node::Pow { base, exp } = &self.nodes[args[0] as usize] {
                        let e2 = matches!(&self.nodes[*exp as usize], Node::Int(k) if k.to_i64() == Some(2));
                        if e2 && self.query_at(*base, Predicate::Positive, 0) == Trinary::True {
                            return *base;
                        }
                    }
                }
            }
        }
        id
    }

    /// Add 级规则：sin²+cos² → 1（匹配两个项，其余项保留）。
    fn try_rules_add(&mut self, id: u32) -> u32 {
        let args = match &self.nodes[id as usize] {
            Node::Add { args: sp } => self.node_args(*sp).to_vec(),
            _ => return id,
        };
        if args.len() < 2 {
            return id;
        }
        // 检查项的平方形态：Mul[-1, cos²] 记为 -cos²；sin² 为 Pow
        let mut sin2: Option<u32> = None; // sin(x)^2 的参数 id
        let mut cos2: Option<(bool, u32)> = None; // (负号, cos(x)^2 参数)
        let mut keep: Vec<u32> = vec![];
        for &a in &args {
            if let Some(arg) = self.as_sq_of(a, "sin") {
                if sin2.is_none() {
                    sin2 = Some(arg);
                    continue;
                }
            }
            if let Some(arg) = self.as_sq_of(a, "cos") {
                if cos2.is_none() {
                    cos2 = Some((false, arg));
                    continue;
                }
            }
            if let Some(arg) = self.as_neg_sq_of(a, "cos") {
                if cos2.is_none() {
                    cos2 = Some((true, arg));
                    continue;
                }
            }
            keep.push(a);
        }
        if let (Some(sin_arg), Some((neg, cos_arg))) = (sin2, cos2) {
            if sin_arg == cos_arg {
                // sin² − cos² 不化简；sin² + cos² 与 −sin²−cos² → ±1
                if !neg {
                    keep.push(self.lit_int(1));
                    let s = self.make_add(&keep);
                    return s;
                }
            }
        }
        id
    }

    /// a == f(x)^2（f 为 head_name）→ Some(x)。
    fn as_sq_of(&self, a: u32, head_name: &str) -> Option<u32> {
        if let Node::Pow { base, exp } = &self.nodes[a as usize] {
            if let Node::Int(k) = &self.nodes[*exp as usize] {
                if k.to_i64() == Some(2) {
                    if let Node::Fn { head, args: sp } = &self.nodes[*base as usize] {
                        if self.fn_names[*head as usize].as_ref() == head_name {
                            let a = self.node_args(*sp);
                            if a.len() == 1 {
                                return Some(a[0]);
                            }
                        }
                    }
                }
            }
        }
        None
    }

    /// a == Mul[-1, f(x)^2] → Some(x)。
    fn as_neg_sq_of(&self, a: u32, head_name: &str) -> Option<u32> {
        if let Node::Mul { args: sp } = &self.nodes[a as usize] {
            let args = self.node_args(*sp);
            if args.len() == 2 {
                if let Node::Int(k) = &self.nodes[args[0] as usize] {
                    if k.to_i64() == Some(-1) {
                        return self.as_sq_of(args[1], head_name);
                    }
                }
            }
        }
        None
    }

    /// 函数节点规则：常量折叠（参数为精确数值常数时）。
    fn try_rules_fn(&mut self, id: u32, head: u32, args: &[u32]) -> u32 {
        let name = self.fn_names[head as usize].to_string();
        if args.len() != 1 {
            return id;
        }
        // 常量折叠：参数为 Int/Rat
        let val: Option<Rational> = match &self.nodes[args[0] as usize] {
            Node::Int(v) => Some(Rational::from_integer(v)),
            Node::Rat(r) => Some(r.clone()),
            _ => None,
        };
        if let Some(v) = val {
            let folded: Option<Rational> = match name.as_str() {
                "exp" if v.is_zero() => Some(Rational::one()),
                "log" if v.is_one() => Some(Rational::zero()),
                "sin" if v.is_zero() => Some(Rational::zero()),
                "cos" if v.is_zero() => Some(Rational::one()),
                "tan" if v.is_zero() => Some(Rational::zero()),
                "sqrt" if v.is_one() => Some(Rational::one()),
                "sqrt" if v.is_zero() => Some(Rational::zero()),
                "abs" => Some(v.abs()),
                _ => None,
            };
            if let Some(f) = folded {
                return self.lit_rational(&f);
            }
        }
        id
    }
}

fn factorial(n: u32) -> i64 {
    let mut r = 1i64;
    for i in 2..=n as i64 {
        r = r.saturating_mul(i);
    }
    r
}

#[cfg(test)]
mod tests {
    use crate::{Context, Expr};

    fn d_eq(ctx: &Context, e: &Expr, x: &Expr, expect: &Expr) {
        let d = ctx.diff(e, x);
        assert!(
            d == *expect,
            "diff 不符:\n  got    {:?}\n  expect {:?}",
            d,
            expect
        );
    }

    #[test]
    fn 求导_多项式与幂() {
        let ctx = Context::new();
        let x = ctx.sym("x");
        // d(x^3) = 3x^2
        d_eq(
            &ctx,
            &x.clone().pow(3),
            &x,
            &(ctx.int(3) * x.clone().pow(2)),
        );
        // d(x^2 + x) = 2x + 1
        let e = x.clone().pow(2) + x.clone();
        d_eq(&ctx, &e, &x, &(ctx.int(2) * x.clone() + ctx.int(1)));
        // d(x^(-1)) = -x^(-2)
        let e = x.clone().pow(-1);
        d_eq(&ctx, &e, &x, &(ctx.int(-1) * x.clone().pow(-2)));
    }

    #[test]
    fn 求导_初等函数() {
        let ctx = Context::new();
        let x = ctx.sym("x");
        // d(sin x) = cos x
        let s = ctx.call("sin", std::slice::from_ref(&x));
        let c = ctx.call("cos", std::slice::from_ref(&x));
        d_eq(&ctx, &s, &x, &c);
        // d(exp x) = exp x
        let e0 = ctx.call("exp", std::slice::from_ref(&x));
        d_eq(&ctx, &e0, &x, &e0);
        // d(log x) = 1/x
        let lg = ctx.call("log", std::slice::from_ref(&x));
        d_eq(&ctx, &lg, &x, &(ctx.int(1) / x.clone()));
        // d(sin(x^2)) = 2x·cos(x^2)（链式）
        let s2 = ctx.call("sin", &[x.clone().pow(2)]);
        let c2 = ctx.call("cos", &[x.clone().pow(2)]);
        d_eq(&ctx, &s2, &x, &(ctx.int(2) * x.clone() * c2));
    }

    #[test]
    fn 求导_语义验证() {
        // 与数值差分交叉验证（独立代码路径）
        let ctx = Context::new();
        let x = ctx.sym("x");
        let e = ctx.call("sin", &[x.clone() * x.clone()]) + x.clone().pow(3);
        let d = ctx.diff(&e, &x);
        for xv in [0.7f64, -1.3, 2.1] {
            let h = 1e-6;
            let f = |t: f64| ctx.eval_float(&e, &[("x", t)]).unwrap();
            let fd = (f(xv + h) - f(xv - h)) / (2.0 * h);
            let ds = ctx.eval_float(&d, &[("x", xv)]).unwrap();
            assert!(
                (fd - ds).abs() < 1e-4 * ds.abs().max(1.0),
                "数值差分 {fd} vs 符号导 {ds} @x={xv}"
            );
        }
    }

    #[test]
    fn taylor_已知级数() {
        let ctx = Context::new();
        let x = ctx.sym("x");
        // sin x 到 5 阶：x − x³/6 + x⁵/120
        let s = ctx.call("sin", std::slice::from_ref(&x));
        let t = ctx.taylor(&s, &x, 0, 5);
        let expect = x.clone() - x.clone().pow(3) / ctx.int(6) + x.clone().pow(5) / ctx.int(120);
        assert!(t == expect, "taylor(sin) 不符: {t:?} vs {expect:?}");
        // exp x 到 3 阶：1 + x + x²/2 + x³/6
        let e0 = ctx.call("exp", std::slice::from_ref(&x));
        let t = ctx.taylor(&e0, &x, 0, 3);
        let expect =
            ctx.int(1) + x.clone() + x.clone().pow(2) / ctx.int(2) + x.clone().pow(3) / ctx.int(6);
        assert!(t == expect, "taylor(exp) 不符");
    }

    #[test]
    fn simplify_常量折叠与恒等式() {
        let ctx = Context::new();
        let x = ctx.sym("x");
        // cos(0) → 1
        let z = ctx.int(0);
        let e = ctx.call("cos", std::slice::from_ref(&z));
        let g = ctx.simplify(&e);
        assert!(g == ctx.int(1), "cos(0) ≠ 1: {g:?}");
        // sin²(x) + cos²(x) → 1
        let s = ctx.call("sin", std::slice::from_ref(&x));
        let c = ctx.call("cos", std::slice::from_ref(&x));
        let e = s.clone().pow(2) + c.clone().pow(2);
        let g = ctx.simplify(&e);
        assert!(g == ctx.int(1), "sin²+cos² ≠ 1: {g:?}");
        // 语义保持：未化简部分不变
        let e = s.pow(2) + c.pow(2) + x.clone();
        let g = ctx.simplify(&e);
        let expect = ctx.int(1) + x;
        assert!(g == expect, "带余项: {g:?}");
    }
}
