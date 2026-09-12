//! cas 开发任务入口（设计见 cas/DESIGN.md §7/§9）。
//!
//! 子命令：
//!   corpus           — 转写并分析真实语料（hpoly_generated.rs，只读）；
//!                      --eval 追加数值求值交叉验证（独立 f64 求值器对拍）
//!   poisson          — 6 变元泊松括号（qiao 正规化的核心负载）：
//!                      反对称/双线性/Jacobi 精确检验 + 有限差分数值交叉验证
//!   oracle           — sympy 差分 harness：随机精确表达式 expand 语义对拍
//!   bench            — 基准：expand (1+x+y+z+w)^20 等（M2 跟踪指标）
//!
//! 用法（在 cas/ 下，建议 --release）：
//!   cargo run --release -p xtask -- <子命令> [选项]

mod bench;
mod corpus;
mod hamilton;
mod legendre;
mod oracle;
mod poisson;
mod util;

use std::env;
use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("corpus") => corpus::run(&args[1..]),
        Some("poisson") => poisson::run(&args[1..]),
        Some("legendre") => legendre::run(&args[1..]),
        Some("hamilton") => hamilton::run(&args[1..]),
        Some("oracle") => oracle::run(&args[1..]),
        Some("bench") => bench::run(&args[1..]),
        _ => {
            eprintln!(
                "用法: cargo run -p xtask -- corpus|poisson|oracle|bench [选项]\n\
                 \x20 corpus  [--file <路径>] [--eval]\n\
                 \x20 poisson [--seed N] [--terms N]\n\
                 \x20 oracle  [--cases N] [--seed N]\n\
                 \x20 bench"
            );
            ExitCode::from(2)
        }
    }
}

pub(crate) fn arg(args: &[String], name: &str, default: &str) -> String {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1))
        .cloned()
        .unwrap_or_else(|| default.to_string())
}

pub(crate) fn flag(args: &[String], name: &str) -> bool {
    args.iter().any(|a| a == name)
}
