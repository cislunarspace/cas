//! 求值（只读遍历）。
//!
//! - `eval_rational_at`：ℚ 上精确求值。Float 字面量与函数节点不参与
//!   （D2 域边界：浮点不进精确算术；初等函数的精确求值属 P1）。
//! - `eval_float_at`：f64 数值求值，支持初等函数头（sin/cos/tan/exp/log/
//!   sqrt/abs），用于数值交叉验证（差分测试的语义判据通道）。

use crate::node::Node;
use crate::{Inner, Rational};
use std::collections::HashMap;

const MAX_DEPTH: u32 = 10_000;

pub(crate) fn eval_rational_at(
    inner: &Inner,
    id: u32,
    vals: &HashMap<u32, Rational>,
    depth: u32,
) -> Option<Rational> {
    assert!(depth <= MAX_DEPTH, "表达式嵌套过深");
    match &inner.nodes[id as usize] {
        Node::Int(v) => Some(Rational::from_integer(v)),
        Node::Rat(r) => Some(r.clone()),
        Node::Float { .. } => None,
        Node::Sym(s) => vals.get(s).cloned(),
        Node::Pow { base, exp } => {
            let k = match &inner.nodes[*exp as usize] {
                Node::Int(v) => v.to_i64(),
                _ => None,
            }?;
            if k.unsigned_abs() > u32::MAX as u64 {
                return None; // 指数天文数字：不 folding 也不求值
            }
            let b = eval_rational_at(inner, *base, vals, depth + 1)?;
            if k >= 0 {
                Some(b.pow_reduced(k as u32))
            } else {
                let inv = b.inv_reduced()?;
                Some(inv.pow_reduced(k.unsigned_abs() as u32))
            }
        }
        Node::Mul { args: sp } => {
            let mut acc = Rational::one();
            for &a in inner.node_args(*sp) {
                acc = acc.mul(&eval_rational_at(inner, a, vals, depth + 1)?);
            }
            Some(acc)
        }
        Node::Add { args: sp } => {
            let mut acc = Rational::zero();
            for &a in inner.node_args(*sp) {
                acc = acc.add(&eval_rational_at(inner, a, vals, depth + 1)?);
            }
            Some(acc)
        }
        Node::Fn { .. } => None,
    }
}

pub(crate) fn eval_float_at(
    inner: &Inner,
    id: u32,
    vals: &HashMap<u32, f64>,
    depth: u32,
) -> Option<f64> {
    assert!(depth <= MAX_DEPTH, "表达式嵌套过深");
    match &inner.nodes[id as usize] {
        Node::Int(v) => v
            .to_i64()
            .map(|x| x as f64)
            .or_else(|| v.to_bigint().to_string().parse::<f64>().ok()),
        Node::Rat(r) => Some(r.to_f64()),
        Node::Float { bits } => Some(f64::from_bits(*bits)),
        Node::Sym(s) => vals.get(s).copied(),
        Node::Pow { base, exp } => {
            let b = eval_float_at(inner, *base, vals, depth + 1)?;
            let e = eval_float_at(inner, *exp, vals, depth + 1)?;
            let r = b.powf(e);
            if r.is_finite() { Some(r) } else { None }
        }
        Node::Mul { args: sp } => {
            let mut acc = 1.0;
            for &a in inner.node_args(*sp) {
                acc *= eval_float_at(inner, a, vals, depth + 1)?;
            }
            Some(acc)
        }
        Node::Add { args: sp } => {
            let mut acc = 0.0;
            for &a in inner.node_args(*sp) {
                acc += eval_float_at(inner, a, vals, depth + 1)?;
            }
            Some(acc)
        }
        Node::Fn { head, args: sp } => {
            let args = inner.node_args(*sp);
            if args.len() != 1 {
                return None;
            }
            let v = eval_float_at(inner, args[0], vals, depth + 1)?;
            let name = &inner.fn_names[*head as usize];
            match name.as_ref() {
                "sin" => Some(v.sin()),
                "cos" => Some(v.cos()),
                "tan" => Some(v.tan()),
                "exp" => Some(v.exp()),
                "log" if v > 0.0 => Some(v.ln()),
                "sqrt" if v >= 0.0 => Some(v.sqrt()),
                "abs" => Some(v.abs()),
                _ => None,
            }
        }
    }
}
