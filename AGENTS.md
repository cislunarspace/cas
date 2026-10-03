# AGENTS.md

symcas 的 agent 协作约定:Rust workspace 实现的确定性计算机代数系统(CAS)。

## Project Overview

提供规范化表达式 arena(hash-consing)、任意精度精确有理数算术、有序稀疏多元多项式、符号微分与级数展开、确定性 plain 往返解析与 LaTeX 输出。核心能力:构造即规范化(O(1) 相等判断)、三值逻辑假设谓词、`parse(plain(e)) == e` 的往返一致性。

## Architecture & Data Flow

七个 crate 组成 workspace,依赖方向自底向上:`cas-domain`(数值核心)← `cas-poly`(多项式代数)← `cas-expr`(arena、规范化、微积分、假设)← `cas-print`/`cas-parse` ← `cas`(门面:prelude、parse、plain、latex)。`xtask` 提供开发辅助任务。所有表达式存活于 `Context` arena,等价表达式共享同一 handle。

## Key Directories

| 路径 | 用途 |
|---|---|
| `crates/domain` | `Integer`/`Rational`/`Ring`/`Field` 数值核心 |
| `crates/poly` | 多项式代数:序(Lex/DegRevLex)、PRS、GCD、无平方分解 |
| `crates/expr` | 表达式 arena、DAG 规范化、微积分、假设 |
| `crates/print` | 确定性 plain 与 latex 打印器 |
| `crates/parse` | plain 语法的递归下降解析器 |
| `crates/cas` | 门面 crate:对外 API |
| `xtask` | 开发辅助任务 |
| `MATLAB-PARITY.md` | MATLAB 对齐说明 |

## Development Commands

构建:`cargo build --workspace`;测试:`cargo test --workspace --locked`;lint:`cargo fmt --all -- --check`、`cargo clippy --workspace --all-targets --locked -- -D warnings`。与 CI(.github/workflows/ci.yml)一致。

## 约定

改动需有测试(含 proptest 属性测试);提交信息用中文;遵循既有命名与 crate 边界,不跨层反向依赖。
