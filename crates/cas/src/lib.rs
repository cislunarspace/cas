//! # symcas - Fast, deterministic computer algebra system in Rust
//!
//! Core features: canonical expression representation (arena + hash-consing),
//! deterministic total order, plain / LaTeX formatting, and roundtrip parsing.
//! Polynomial kernel, differentiation, series expansion, and simplification.
//!
//! # Example 1: Construction and Canonical Form
//!
//! ```
//! use symcas::prelude::*;
//!
//! let ctx = Context::new();
//! let x = ctx.sym("x");
//! let y = ctx.sym("y");
//!
//! // Operator overloading constructs canonical forms: flattening, ordering,
//! // exact arithmetic folding, and term combination happen immediately.
//! let e = (x.clone() + y.clone()).pow(3) - x.clone().pow(3) - y.clone().pow(3);
//! assert_eq!(symcas::plain(&ctx, &e), "(x + y)^3 - x^3 - y^3");
//!
//! // 同一数学内容的两条构造路径落在同一节点（hash-consing）
//! let a = x.clone() + x.clone();
//! let b = ctx.int(2) * x.clone();
//! assert!(a == b);
//! ```
//!
//! # 例 2：解析与往返
//!
//! ```
//! use symcas::prelude::*;
//!
//! let ctx = Context::new();
//! let x = ctx.sym("x");
//!
//! let e = symcas::parse(&ctx, "x + x + 1/2*x").unwrap();
//! assert_eq!(symcas::plain(&ctx, &e), "5/2*x");
//!
//! // Plain output parses back to the exact same expression node (roundtrip guarantee)
//! let back = symcas::parse(&ctx, "5/2*x").unwrap();
//! assert!(back == e);
//!
//! // sym! 宏批量声明符号
//! let (u, v) = sym!(&ctx, u, v);
//! let _ = u + v;
//! ```

pub use cas_expr::{Context, Expr, Integer, Rational};

pub use cas_parse::ParseError;

/// 解析 plain 文法（见 `cas-parse` 文档）。
pub fn parse(ctx: &Context, src: &str) -> Result<Expr, ParseError> {
    cas_parse::parse(ctx, src)
}

/// plain 输出：往返格式，字节由确定性全序决定。
pub fn plain(ctx: &Context, e: &Expr) -> String {
    cas_print::plain(ctx, e)
}

/// LaTeX 输出：面向阅读，无往返承诺。
pub fn latex(ctx: &Context, e: &Expr) -> String {
    cas_print::latex(ctx, e)
}

pub mod prelude {
    pub use crate::{latex, parse, plain};
    pub use cas_expr::{Context, Expr, sym};
}
