//! expr ↔ poly 桥（M3：expand 的多项式快路径，设计 D3）。
//!
//! 快路径只接受纯多项式子树（Int/Rat/Sym/非负整幂/加乘）：提取为
//! `Poly<Rational>`（变元序 = 名字节序，与创建序无关），在多项式域完成
//! 分配/求幂后重建表达式。重建经规范形构造器，故快路径与朴素 DAG 分配
//! 必然产出**同一节点**（同构性由测试钉死）。
//!
//! 含浮点、函数、负幂的子树返回 None，由 DAG 路径接管。

use crate::Inner;
use crate::node::Node;
use cas_domain::Rational;
use cas_poly::{MonOrder, Poly, PolyRing};
use std::collections::HashMap;
use std::sync::Arc;

/// 快路径结果项数上限（与 DAG 路径的 MAX_EXPAND_TERMS 对齐）。
const MAX_FAST_TERMS: u128 = 1_000_000;

/// 收集子树全部符号（符号表 id，去重）。
fn collect_syms(inner: &Inner, id: u32, out: &mut Vec<u32>, depth: u32) {
    assert!(depth <= 10_000, "表达式嵌套过深");
    match &inner.nodes[id as usize] {
        Node::Sym(s) => {
            if !out.contains(s) {
                out.push(*s);
            }
        }
        Node::Pow { base, exp } => {
            collect_syms(inner, *base, out, depth + 1);
            collect_syms(inner, *exp, out, depth + 1);
        }
        Node::Fn { args: sp, .. } | Node::Mul { args: sp } | Node::Add { args: sp } => {
            for &a in inner.node_args(*sp) {
                collect_syms(inner, a, out, depth + 1);
            }
        }
        _ => {}
    }
}

/// 提取为多项式。返回（ring，变元 sym-id 表，多项式）。
pub(crate) fn to_poly(inner: &Inner, id: u32) -> Option<(Arc<PolyRing>, Vec<u32>, Poly<Rational>)> {
    let mut syms: Vec<u32> = Vec::new();
    collect_syms(inner, id, &mut syms, 0);
    // 变元序 = 名字节序升序（D6；跨构造路径确定性）
    let mut named: Vec<(u32, String)> = syms
        .iter()
        .map(|&s| (s, inner.sym_names[s as usize].to_string()))
        .collect();
    named.sort_by(|a, b| a.1.cmp(&b.1));
    let ring = PolyRing::new(named.iter().map(|(_, n)| n.clone()), MonOrder::DegRevLex);
    let var_ids: Vec<u32> = named.iter().map(|&(s, _)| s).collect();
    let vi: HashMap<u32, usize> = named
        .iter()
        .enumerate()
        .map(|(i, &(s, _))| (s, i))
        .collect();
    let p = poly_at(inner, id, &ring, &vi, 0)?;
    Some((ring, var_ids, p))
}

fn poly_at(
    inner: &Inner,
    id: u32,
    ring: &Arc<PolyRing>,
    vi: &HashMap<u32, usize>,
    depth: u32,
) -> Option<Poly<Rational>> {
    assert!(depth <= 10_000, "表达式嵌套过深");
    let nvars = ring.nvars();
    match &inner.nodes[id as usize] {
        Node::Int(v) => Some(Poly::constant(ring.clone(), Rational::from_integer(v))),
        Node::Rat(r) => Some(Poly::constant(ring.clone(), r.clone())),
        Node::Sym(s) => {
            let j = vi.get(s)?;
            let mut e = vec![0u32; nvars];
            e[*j] = 1;
            Some(Poly::from_terms(ring.clone(), [(e, Rational::one())]))
        }
        Node::Pow { base, exp } => {
            let k = match &inner.nodes[*exp as usize] {
                Node::Int(v) => v.to_i64(),
                _ => None,
            }?;
            if k < 0 || k > u32::MAX as i64 {
                return None;
            }
            let b = poly_at(inner, *base, ring, vi, depth + 1)?;
            // 项数上界 C(t+k-1, t-1)（t 个单项带重复取 k 个）——防 (x+y+z+w)^256 类
            if est_terms(b.nterms() as u128, k as u128) > MAX_FAST_TERMS {
                return None;
            }
            Some(b.pow(k as u32))
        }
        Node::Mul { args: sp } => {
            let mut acc = Poly::constant(ring.clone(), Rational::one());
            for &a in inner.node_args(*sp) {
                acc = acc.mul(&poly_at(inner, a, ring, vi, depth + 1)?);
            }
            Some(acc)
        }
        Node::Add { args: sp } => {
            let mut acc = Poly::zero(ring.clone());
            for &a in inner.node_args(*sp) {
                acc = acc.add(&poly_at(inner, a, ring, vi, depth + 1)?);
            }
            Some(acc)
        }
        Node::Float { .. } | Node::Fn { .. } => None,
    }
}

/// C(t+k-1, t-1) 的饱和估计（t = 基底项数，k = 幂次）。
fn est_terms(t: u128, k: u128) -> u128 {
    if t <= 1 || k == 0 {
        return 1;
    }
    // C(t+k-1, m)，m = min(k, t-1)
    let m = k.min(t - 1);
    let n = t + k - 1;
    let mut r: u128 = 1;
    for i in 0..m {
        r = r.saturating_mul(n - i) / (i + 1);
        if r > MAX_FAST_TERMS {
            return r;
        }
    }
    r
}

/// 多项式重建为规范形表达式（项 → Mul[系数, Pow(sym, e)…] → 整批 make_add）。
pub(crate) fn from_poly(
    inner: &mut Inner,
    ring: &Arc<PolyRing>,
    var_ids: &[u32],
    poly: &Poly<Rational>,
) -> u32 {
    let mut terms = Vec::with_capacity(poly.nterms());
    for (e, c) in poly.terms() {
        let coeff = inner.lit_rational(c);
        let mut factors: Vec<u32> = Vec::with_capacity(ring.nvars() + 1);
        factors.push(coeff);
        for (j, &power) in e.iter().enumerate() {
            if power > 0 {
                let sym = inner.intern_node(Node::Sym(var_ids[j]), &[]);
                let lit = inner.lit_int(power as i64);
                factors.push(inner.make_pow(sym, lit));
            }
        }
        let t = if factors.len() == 1 {
            factors[0]
        } else {
            inner.make_mul(&factors)
        };
        terms.push(t);
    }
    if terms.is_empty() {
        return inner.lit_int(0);
    }
    inner.make_add(&terms)
}
