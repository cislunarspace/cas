//! # cas——Rust 通用符号代数引擎
//!
//! P0-M1 骨架：表达式规范形（arena + hash-consing）、确定性全序、
//! plain/LaTeX 输出、可往返解析。多项式内核（M2）、求导与级数（P1）、
//! Gröbner（P2）、代码生成（P3）按 `cas/DESIGN.md` 分期落地。
//!
//! # 例 1：构造与规范形
//!
//! ```
//! use cas::prelude::*;
//!
//! let ctx = Context::new();
//! let x = ctx.sym("x");
//! let y = ctx.sym("y");
//!
//! // 运算符重载构造规范形：扁平化、排序、数值折叠、同类合并即刻发生
//! let e = (x.clone() + y.clone()).pow(3) - x.clone().pow(3) - y.clone().pow(3);
//! assert_eq!(cas::plain(&ctx, &e), "(x + y)^3 - x^3 - y^3");
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
//! use cas::prelude::*;
//!
//! let ctx = Context::new();
//! let x = ctx.sym("x");
//!
//! let e = cas::parse(&ctx, "x + x + 1/2*x").unwrap();
//! assert_eq!(cas::plain(&ctx, &e), "5/2*x");
//!
//! // plain 输出可解析回同一表达式（D6 往返承诺）
//! let back = cas::parse(&ctx, "5/2*x").unwrap();
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
