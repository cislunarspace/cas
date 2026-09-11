//! 规范形构造器（L0"构造即规范"的全部逻辑，见 DESIGN.md D1/D4）。
//!
//! 不变量由 [`Inner::make_add`] / [`Inner::make_mul`] / [`Inner::make_pow`]
//! 建立：嵌套扁平化、参数按全序排序、精确数值折叠、同类项/同底幂合并。
//! 浮点字面量不参与任何算术折叠（D2 的域边界）；负浮点在 `Context::float`
//! 入口已拆成 `-1 * |v|`，因此 arena 里的 Float 节点值恒非负。

use crate::Inner;
use crate::hash::hash_node;
use crate::node::{Node, Span};
use cas_domain::{Integer, Rational};

/// 数值幂折叠的规模守卫：预计结果位数超过此值则保留未折叠形态。
/// 上限 2^16 位（约 8 KB 常数）：既挡住 `2^999999999` 这类字面量炸弹，
/// 也截断"折叠结果的再折叠"链（如 ((p/q)^-253)^-256 会滚出几十万比特的
/// 常数，其打印/解析/归约都是平方级负担——perf_probe 实测单例 14 s）。
/// 更大的精确常数仍可通过解析器字面量进入（解析路径不做幂折叠）。
const MAX_FOLD_BITS: u64 = 1 << 16;

impl Inner {
    // ── interning 基础 ────────────────────────────────────────────

    pub(crate) fn intern_node(&mut self, node: Node, args: &[u32]) -> u32 {
        let h = hash_node(&node, args, &self.hashes);
        if let Some(cands) = self.intern.get(&h) {
            for &c in cands {
                if self.nodes[c as usize] == node && self.stored_args(c) == args {
                    return c;
                }
            }
        }
        let off = self.args_slab.len() as u32;
        let len = u16::try_from(args.len()).expect("参数超过 65535 个");
        self.args_slab.extend_from_slice(args);
        let id = self.nodes.len() as u32;
        let node = match node {
            Node::Fn { head, .. } => Node::Fn {
                head,
                args: Span { off, len },
            },
            Node::Mul { .. } => Node::Mul {
                args: Span { off, len },
            },
            Node::Add { .. } => Node::Add {
                args: Span { off, len },
            },
            other => other,
        };
        self.nodes.push(node);
        self.hashes.push(h);
        self.intern.entry(h).or_default().push(id);
        id
    }

    /// 读节点自身携带的参数区段。
    fn stored_args(&self, id: u32) -> &[u32] {
        match &self.nodes[id as usize] {
            Node::Fn { args: sp, .. } | Node::Mul { args: sp } | Node::Add { args: sp } => {
                self.span_slice(*sp)
            }
            _ => &[],
        }
    }

    pub(crate) fn node_args(&self, sp: Span) -> &[u32] {
        self.span_slice(sp)
    }

    fn span_slice(&self, sp: Span) -> &[u32] {
        let lo = sp.off as usize;
        &self.args_slab[lo..lo + sp.len as usize]
    }

    pub(crate) fn lit_int(&mut self, v: i64) -> u32 {
        self.intern_node(Node::Int(Integer::from_i64(v)), &[])
    }

    pub(crate) fn lit_integer(&mut self, v: &Integer) -> u32 {
        self.intern_node(Node::Int(v.clone()), &[])
    }

    /// 有理数字面量：分母为 1 时归一为 Int 节点（Float/Rat 节点规范的一部分）。
    pub(crate) fn lit_rational(&mut self, r: &Rational) -> u32 {
        if r.den().is_one() {
            let n = r.num();
            return self.lit_integer(&n);
        }
        self.intern_node(Node::Rat(r.clone()), &[])
    }

    /// 浮点字面量：调用方保证 `v` 非负且非 NaN；`-0.0` 在此归一为 `0.0`。
    pub(crate) fn lit_float(&mut self, v: f64) -> u32 {
        debug_assert!(!v.is_nan() && v >= 0.0);
        let v = if v == 0.0 { 0.0 } else { v };
        self.intern_node(Node::Float { bits: v.to_bits() }, &[])
    }

    pub(crate) fn sym_node(&mut self, name: &str) -> u32 {
        let sid = match self.syms.get(name) {
            Some(&s) => s,
            None => {
                let sid = self.sym_names.len() as u32;
                self.sym_names.push(name.into());
                self.syms.insert(name.into(), sid);
                sid
            }
        };
        self.intern_node(Node::Sym(sid), &[])
    }

    pub(crate) fn fn_node(&mut self, head: &str, args: &[u32]) -> u32 {
        let hid = match self.fn_heads.get(head) {
            Some(&h) => h,
            None => {
                let hid = self.fn_names.len() as u32;
                self.fn_names.push(head.into());
                self.fn_heads.insert(head.into(), hid);
                hid
            }
        };
        self.fn_node_by_id(hid, args)
    }

    /// 按已有 head 索引重建函数节点（subst/expand 等变换路径用）。
    pub(crate) fn fn_node_by_id(&mut self, head: u32, args: &[u32]) -> u32 {
        self.intern_node(
            Node::Fn {
                head,
                args: Span { off: 0, len: 0 },
            },
            args,
        )
    }

    // ── 规范形构造 ────────────────────────────────────────────────

    pub(crate) fn make_add(&mut self, args: &[u32]) -> u32 {
        // 1. 扁平化嵌套 Add
        let mut flat: Vec<u32> = Vec::with_capacity(args.len());
        for &a in args {
            match &self.nodes[a as usize] {
                Node::Add { args: sp } => flat.extend_from_slice(self.node_args(*sp)),
                _ => flat.push(a),
            }
        }

        // 2. 归类。浮点（含浮点系数的 Mul 整体）不作算术合并，原样保留。
        enum Arg {
            Num(Rational),
            Atom, // Sym/Fn/Pow/Add/无精确系数的 Mul：整体作 key
            MulExact { coeff: Rational, rest: Vec<u32> },
        }
        let mut acc = Rational::zero();
        let mut terms: Vec<(u32, Rational)> = Vec::new();
        let mut floats: Vec<u32> = Vec::new();
        for a in flat {
            let arg = match &self.nodes[a as usize] {
                Node::Int(v) => Arg::Num(Rational::from_integer(v)),
                Node::Rat(r) => Arg::Num(r.clone()),
                Node::Float { .. } => {
                    floats.push(a);
                    continue;
                }
                Node::Sym(_) | Node::Fn { .. } | Node::Pow { .. } | Node::Add { .. } => Arg::Atom,
                Node::Mul { args: sp } => {
                    let margs = self.node_args(*sp);
                    match margs.split_first() {
                        Some((first, rest)) => match &self.nodes[*first as usize] {
                            Node::Int(v) if rest.is_empty() => Arg::Num(Rational::from_integer(v)),
                            Node::Int(v) => Arg::MulExact {
                                coeff: Rational::from_integer(v),
                                rest: rest.to_vec(),
                            },
                            Node::Rat(r) if rest.is_empty() => Arg::Num(r.clone()),
                            Node::Rat(r) => Arg::MulExact {
                                coeff: r.clone(),
                                rest: rest.to_vec(),
                            },
                            _ => Arg::Atom, // 浮点居首的 Mul：浮点系数，不合并
                        },
                        None => Arg::Atom, // 不可达：Mul 至少两参
                    }
                }
            };
            match arg {
                Arg::Num(v) => acc = acc.add(&v),
                Arg::Atom => merge_term(&mut terms, a, Rational::one()),
                Arg::MulExact { coeff, rest } => {
                    let key = self.make_mul(&rest);
                    merge_term(&mut terms, key, coeff);
                }
            }
        }

        // 3. 重建：零系数项丢弃；系数 1 裸 key；否则 系数*key
        let mut out: Vec<u32> = floats;
        for (key, coeff) in terms {
            if coeff.is_zero() {
                continue;
            }
            let term = if coeff.is_one() {
                key
            } else {
                let c = self.lit_rational(&coeff);
                self.make_mul(&[c, key])
            };
            out.push(term);
        }
        if !acc.is_zero() {
            let c = self.lit_rational(&acc);
            out.push(c);
        }

        // 4. 全序排序后建节点
        let me = &*self;
        out.sort_by(|&x, &y| me.cmp_ids(x, y));
        match out.len() {
            0 => self.lit_int(0),
            1 => out[0],
            _ => self.intern_node(
                Node::Add {
                    args: Span { off: 0, len: 0 },
                },
                &out,
            ),
        }
    }

    pub(crate) fn make_mul(&mut self, args: &[u32]) -> u32 {
        // 1. 扁平化嵌套 Mul
        let mut flat: Vec<u32> = Vec::with_capacity(args.len());
        for &a in args {
            match &self.nodes[a as usize] {
                Node::Mul { args: sp } => flat.extend_from_slice(self.node_args(*sp)),
                _ => flat.push(a),
            }
        }

        // 2. 精确数值进累乘系数；其余按底归类，指数收集待求和。
        //    约定与 sympy 一致：x^a · x^b → x^(a+b) 无条件成立（x=0 的退化点除外）。
        enum Fac {
            Exact(Rational),
            Power { base: u32, exp: u32 },
        }
        let one = self.lit_int(1);
        let mut coeff = Rational::one();
        let mut bases: Vec<(u32, Vec<u32>)> = Vec::new();
        for a in flat {
            let fac = match &self.nodes[a as usize] {
                Node::Int(v) => Fac::Exact(Rational::from_integer(v)),
                Node::Rat(r) => Fac::Exact(r.clone()),
                Node::Pow { base, exp } => Fac::Power {
                    base: *base,
                    exp: *exp,
                },
                _ => Fac::Power { base: a, exp: one }, // Sym/Fn/Add/Float
            };
            match fac {
                Fac::Exact(v) => {
                    coeff = coeff.mul(&v);
                    if coeff.is_zero() {
                        return self.lit_int(0);
                    }
                }
                Fac::Power { base, exp } => match bases.iter_mut().find(|(b, _)| *b == base) {
                    Some((_, es)) => es.push(exp),
                    None => bases.push((base, vec![exp])),
                },
            }
        }

        // 3. 同底合并：指数求和；和为零的因子消失；折叠回数值的并进系数
        let mut out: Vec<u32> = Vec::new();
        for (base, exps) in bases {
            let e = self.make_add(&exps);
            if matches!(&self.nodes[e as usize], Node::Int(v) if v.is_zero()) {
                continue;
            }
            let f = self.make_pow(base, e);
            let folded = match &self.nodes[f as usize] {
                Node::Int(v) => Some(Rational::from_integer(v)),
                Node::Rat(r) => Some(r.clone()),
                _ => None,
            };
            match folded {
                Some(v) => coeff = coeff.mul(&v),
                None => out.push(f),
            }
        }
        if coeff.is_zero() {
            return self.lit_int(0);
        }
        if !coeff.is_one() {
            let c = self.lit_rational(&coeff);
            out.push(c);
        }

        let me = &*self;
        out.sort_by(|&x, &y| me.cmp_ids(x, y));
        match out.len() {
            0 => self.lit_int(1),
            1 => out[0],
            _ => self.intern_node(
                Node::Mul {
                    args: Span { off: 0, len: 0 },
                },
                &out,
            ),
        }
    }

    pub(crate) fn make_pow(&mut self, base: u32, exp: u32) -> u32 {
        enum E {
            NotInt,
            Zero,
            One,
            Int(i64),
        }
        let e = match &self.nodes[exp as usize] {
            Node::Int(v) if v.is_zero() => E::Zero,
            Node::Int(v) if v.is_one() => E::One,
            Node::Int(v) => match v.to_i64() {
                Some(i) => E::Int(i),
                None => E::NotInt,
            },
            _ => E::NotInt,
        };
        match e {
            E::Zero => return self.lit_int(1), // 约定 0^0 = 1（与 sympy 一致）
            E::One => return base,
            _ => {}
        }

        enum B {
            Int(Integer),
            Rat(Rational),
            One,
            Zero,
            Pow { base: u32, exp: u32 },
            Mul(Vec<u32>),
            Other,
        }
        let b = match &self.nodes[base as usize] {
            Node::Int(v) if v.is_one() => B::One,
            Node::Int(v) if v.is_zero() => B::Zero,
            Node::Int(v) => B::Int(v.clone()),
            Node::Rat(r) if r.is_one() => B::One,
            Node::Rat(r) => B::Rat(r.clone()),
            Node::Pow { base, exp } => B::Pow {
                base: *base,
                exp: *exp,
            },
            Node::Mul { args: sp } => B::Mul(self.node_args(*sp).to_vec()),
            _ => B::Other,
        };
        // 底为精确 1：1^任何 = 1
        if matches!(b, B::One) {
            return self.lit_int(1);
        }

        if let E::Int(ei) = e {
            match b {
                B::Zero => {
                    if ei > 0 {
                        return self.lit_int(0);
                    }
                    // 0 的负幂无定义：保留节点
                }
                B::Int(v) => {
                    if let Some(id) = self.fold_int_pow(&v, ei) {
                        return id;
                    }
                }
                B::Rat(r) => {
                    if let Some(id) = self.fold_rat_pow(&r, ei) {
                        return id;
                    }
                }
                B::Pow { base: b2, exp: e2 } => {
                    // (b^e)^k = b^(e·k)：整数外幂单值、主值一致，安全折叠
                    let prod = self.make_mul(&[exp, e2]);
                    return self.make_pow(b2, prod);
                }
                B::Mul(fs) if ei > 0 => {
                    // 非负整数幂分配到因子：(a·b)^k = a^k·b^k
                    let pw: Vec<u32> = fs.iter().map(|&f| self.make_pow(f, exp)).collect();
                    return self.make_mul(&pw);
                }
                _ => {}
            }
        }
        self.intern_node(Node::Pow { base, exp }, &[])
    }

    /// 整数底的整数幂折叠；返回 None 表示保留未折叠形态。
    fn fold_int_pow(&mut self, v: &Integer, ei: i64) -> Option<u32> {
        // ±1 特判：不增长，直接按奇偶折叠
        if v == &Integer::from_i64(-1) {
            return Some(self.lit_int(if ei.rem_euclid(2) == 0 { 1 } else { -1 }));
        }
        let mag = ei.unsigned_abs();
        if v.bit_len().max(1).saturating_mul(mag) > MAX_FOLD_BITS {
            return None;
        }
        let m = u32::try_from(mag).ok()?;
        if ei >= 0 {
            Some(self.lit_integer(&v.pow(m)))
        } else {
            // v 非零（零底已分流）
            let r = Rational::from_ints(&Integer::one(), &v.pow(m)).expect("非零底");
            Some(self.lit_rational(&r))
        }
    }

    fn fold_rat_pow(&mut self, r: &Rational, ei: i64) -> Option<u32> {
        let mag = ei.unsigned_abs();
        let bits = r.num().bit_len().max(r.den().bit_len()).max(1);
        if bits.saturating_mul(mag) > MAX_FOLD_BITS {
            return None;
        }
        let m = u32::try_from(mag).ok()?;
        // r 既约 ⇒ 幂与逆元仍既约，走跳过 gcd 的快速路径
        if ei >= 0 {
            Some(self.lit_rational(&r.pow_reduced(m)))
        } else {
            Some(self.lit_rational(&r.inv_reduced()?.pow_reduced(m)))
        }
    }
}

fn merge_term(terms: &mut Vec<(u32, Rational)>, key: u32, coeff: Rational) {
    match terms.iter_mut().find(|(k, _)| *k == key) {
        Some((_, c)) => *c = c.add(&coeff),
        None => terms.push((key, coeff)),
    }
}
