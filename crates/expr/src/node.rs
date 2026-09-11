//! arena 节点定义。
//!
//! 规范形不变量（构造器保证，见 canonical.rs）：
//! - `Int`/`Rat` 数值已折叠、既约；`Rat` 分母 > 1（分母为 1 归一为 `Int`）；
//! - `Float` 仅作字面量，`-0.0` 已归一为 `0.0`，禁止 NaN；不参与算术折叠；
//! - `Add`/`Mul` 为 n 元：已扁平化（无同名嵌套）、参数按全序升序、
//!   同类项已合并、数值系数居首、至少两个参数；
//! - `Pow` 的指数不为 0/1，数值底的整数幂已折叠（受规模守卫限制）。

use cas_domain::{Integer, Rational};

/// 节点在 arena 中的参数区段：`args_slab[off..off + len]`。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Span {
    pub off: u32,
    pub len: u16,
}

/// 相等性只比内容：带参数的节点不比 `Span`（同一内容可在 slab 的不同位置，
/// 参数内容由 intern 侧的 `stored_args` 比较覆盖）。
#[derive(Clone, Debug, Eq)]
pub(crate) enum Node {
    Int(Integer),
    Rat(Rational),
    /// f64 位模式（负值与 -0.0 已在入口归一；多精度 prec 字段随 P1 引入）。
    Float {
        bits: u64,
    },
    /// 符号表索引（`sym_names`）。
    Sym(u32),
    /// 函数应用：head 为函数名表索引，参数保持语义次序（不排序）。
    Fn {
        head: u32,
        args: Span,
    },
    Pow {
        base: u32,
        exp: u32,
    },
    Mul {
        args: Span,
    },
    Add {
        args: Span,
    },
}

impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Node::Int(a), Node::Int(b)) => a == b,
            (Node::Rat(a), Node::Rat(b)) => a == b,
            (Node::Float { bits: a, .. }, Node::Float { bits: b, .. }) => a == b,
            (Node::Sym(a), Node::Sym(b)) => a == b,
            (Node::Pow { base: a, exp: b }, Node::Pow { base: c, exp: d }) => a == c && b == d,
            // 参数内容在 intern_node 里经 stored_args 比对
            (Node::Fn { head: a, .. }, Node::Fn { head: b, .. }) => a == b,
            (Node::Mul { .. }, Node::Mul { .. }) => true,
            (Node::Add { .. }, Node::Add { .. }) => true,
            _ => false,
        }
    }
}
