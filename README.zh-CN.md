# symcas

快速、确定性的 Rust 计算机代数系统。

`symcas` 提供带哈希共享的规范表达式 arena、基于任意精度有理数的精确算术、有序稀疏多元多项式、符号微分、级数展开、定向化简、可往返的 plain 解析，以及 LaTeX 渲染。

[English](README.md) | **简体中文**

[![CI](https://github.com/cislunarspace/cas/actions/workflows/ci.yml/badge.svg)](https://github.com/cislunarspace/cas/actions/workflows/ci.yml)
[![Crates.io](https://img.shields.io/crates/v/symcas)](https://crates.io/crates/symcas)
[![Docs.rs](https://docs.rs/symcas/badge.svg)](https://docs.rs/symcas)
[![License](https://img.shields.io/crates/l/symcas)](#许可证)

## 目录

- [核心特性](#核心特性)
- [架构](#架构)
- [快速上手](#快速上手)
- [延伸阅读](#延伸阅读)
- [许可证](#许可证)

## 核心特性

- **Canonical Construction & O(1) Equality**：表达式在构造时自动规范化：展平、项排序、精确算术折叠、同类项合并。同一 `Context` 内等价的表达式共享相同的 arena 句柄（`Expr::raw_id`）。
- **Exact Numeric Domain**：基于 `num-bigint` / `num-integer` 的任意精度整数与有理数，小整数使用 `i64` 内联优化。
- **Multivariate Polynomials**：系数环泛型的有序稀疏分布多项式，支持加法、乘法、精确除法、多项式余式序列（PRS）、多元 GCD 与无平方分解。
- **Symbolic Calculus**：初等函数（`sin`、`cos`、`tan`、`exp`、`log`、`sqrt`、`abs`）的符号微分、泰勒级数展开、带规模上界的代数展开与消去。
- **Assumptions & Predicates**：三值逻辑（`True`、`False`、`Unknown`）与去重谓词（`Positive`、`NonZero`、`Integer` 等）。
- **Deterministic Roundtrip**：机器可读的 `plain` 格式解析后得到完全相同的表达式节点（`parse(plain(e)) == e`），并提供面向发表的 `latex` 输出。

## 架构

系统由多个模块化 crate 组成：

| Crate | Role |
|---|---|
| `symcas` | 高层门面 crate：prelude、`parse`、`plain` 与 `latex` API |
| `cas-domain` | 数值核心：`Integer`、`Rational`、`Ring` 与 `Field` trait |
| `cas-poly` | 多元多项式代数、单项式序（`Lex`、`DegRevLex`）与 GCD |
| `cas-expr` | 表达式 arena、DAG 规范化、微积分与假设 |
| `cas-print` | 确定性的 `plain` 与 `latex` 打印器 |
| `cas-parse` | `plain` 数学语法的递归下降解析器 |

## 快速上手

将 `symcas` 添加到你的 `Cargo.toml`：

```toml
[dependencies]
symcas = "0.1"
```

需要 Rust 1.85+（由 Cargo.toml 中的 `rust-version` 强制约束）。

### 基本流程

每次使用 symcas 都遵循同一套标准流程：

1. **创建 `Context`**（持有全部表达式的 arena）：`Context::new()`。
2. **创建原子**：符号 `ctx.sym("x")`、整数 `ctx.int(2)`、精确有理数 `ctx.int(1) / ctx.int(2)`；`sym!(&ctx, x, y)` 可批量声明符号。
3. **构造或解析表达式**：用 `+`、`-`、`*`、`/`、`.pow(n)` 与 `ctx.call("sin", &[x.clone()])` 组合原子，或用 `symcas::parse(&ctx, "sin(x)^2")` 读取文本；构造即完成规范化。
4. **做变换**：`ctx.expand(&e)`、`ctx.simplify(&e)`、`ctx.diff(&e, &x)`、`ctx.taylor(&e, &x, 0, 4)`、`ctx.subst(&e, &[("x", v)])`、`ctx.cancel(&e)`、`ctx.factor(&e)`。
5. **输出结果**：`symcas::plain(&ctx, &e)`（可往返的机器格式）或 `symcas::latex(&ctx, &e)`（面向发表的排版格式）。

下面的每个示例都是这套流程在某一场景下的应用。

### 1. 表达式与规范形

*适用于判断写法不同的两个公式是否为同一数学对象：*

$$ (x + y)^{2} = x^{2} + 2xy + y^{2} \qquad x + x = 2x $$

```rust
use symcas::prelude::*;

let ctx = Context::new();
let x = ctx.sym("x");
let y = ctx.sym("y");

// Arithmetic operations construct canonical forms on the fly
let e = (x.clone() + y.clone()).pow(3) - x.clone().pow(3) - y.clone().pow(3);
assert_eq!(symcas::plain(&ctx, &e), "(x + y)^3 - x^3 - y^3");

// Hash-consing: expressions with identical canonical structure share identical node IDs
let a = x.clone() + x.clone();
let b = ctx.int(2) * x.clone();
assert!(a == b);
```

### 2. 解析与往返保证

*适用于把随手写下的输入归一为唯一确定的规范形式：*

$$ x + x + \frac{1}{2}x \;\to\; \frac{5}{2}x $$

```rust
use symcas::prelude::*;

let ctx = Context::new();
let e = symcas::parse(&ctx, "x + x + 1/2*x").unwrap();
assert_eq!(symcas::plain(&ctx, &e), "5/2*x");

// Parsing plain output yields the exact same canonical node
let back = symcas::parse(&ctx, "5/2*x").unwrap();
assert!(back == e);
```

### 3. 展开、化简与替换

*适用于验证手工推导的恒等式、精简冗余的表达：*

$$ (x + y)^{2} = x^{2} + y^{2} + 2xy \qquad \sin^{2}x + \cos^{2}x = 1 $$

```rust
use symcas::prelude::*;

let ctx = Context::new();
let (x, y) = sym!(&ctx, x, y);

// Expand products and powers
let e = (x.clone() + y.clone()).pow(2);
let expanded = ctx.expand(&e);
assert_eq!(symcas::plain(&ctx, &expanded), "x^2 + y^2 + 2*x*y");

// Algebraic substitution
let sub = ctx.subst(&expanded, &[("x", ctx.int(1))]);
assert_eq!(symcas::plain(&ctx, &sub), "1 + y^2 + 2*y");

// Trigonometric simplification
let trig = ctx.call("sin", &[x.clone()]).pow(2) + ctx.call("cos", &[x.clone()]).pow(2);
let sim = ctx.simplify(&trig);
assert_eq!(symcas::plain(&ctx, &sim), "1");
```

### 4. 微分与泰勒展开

*适用于计算不便手算的导数、用级数做局部近似：*

$$ \frac{\mathrm{d}}{\mathrm{d}x}\left(\sin(x)\,x^{2}\right) = 2x\sin(x) + x^{2}\cos(x) $$

$$ e^{x} = 1 + x + \frac{x^{2}}{2} + \frac{x^{3}}{6} + \frac{x^{4}}{24} + \cdots $$

```rust
use symcas::prelude::*;

let ctx = Context::new();
let x = ctx.sym("x");

// Symbolic derivative d/dx (sin(x) * x^2)
let f = ctx.call("sin", &[x.clone()]) * x.clone().pow(2);
let df = ctx.diff(&f, &x);
assert_eq!(symcas::plain(&ctx, &df), "2*x*sin(x) + x^2*cos(x)");

// Taylor series around x = 0 to order 4 for exp(x)
let ex = ctx.call("exp", &[x.clone()]);
let s = ctx.taylor(&ex, &x, 0, 4);
assert_eq!(symcas::plain(&ctx, &s), "1 + x + 1/24*x^4 + 1/6*x^3 + 1/2*x^2");
```

### 5. 有理式约简与多项式操作

*适用于约去有理式公因式、因式分解暴露多项式零点：*

$$ \frac{x^{2} - 1}{x - 1} = x + 1 \qquad x^{2} - 4 = (x - 2)(x + 2) $$

```rust
use symcas::prelude::*;

let ctx = Context::new();
let x = ctx.sym("x");

// Rational cancellation via polynomial GCD
let frac = (x.clone().pow(2) - ctx.int(1)) / (x.clone() - ctx.int(1));
let reduced = ctx.cancel(&frac);
assert_eq!(symcas::plain(&ctx, &reduced), "1 + x");

// Factoring polynomials
let poly = x.clone().pow(2) - ctx.int(4);
let factored = ctx.factor(&poly);
assert_eq!(symcas::plain(&ctx, &factored), "(-2 + x)*(2 + x)");
```

## 延伸阅读

- [DESIGN.md](DESIGN.md)：设计文档（中文），涵盖目标、决策 D1 至 D6 与验收标准
- [MATLAB-PARITY.md](MATLAB-PARITY.md)：与 MATLAB Symbolic Math Toolbox 的功能对齐状态
- [docs.rs](https://docs.rs/symcas) 上的 API 参考文档

## 许可证

根据以下任一许可证授权：

- Apache License, Version 2.0（[LICENSE-APACHE](LICENSE-APACHE) 或 <http://www.apache.org/licenses/LICENSE-2.0>）
- MIT License（[LICENSE-MIT](LICENSE-MIT) 或 <http://opensource.org/licenses/MIT>）

由你选择。
