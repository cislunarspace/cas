//! bench — M2 跟踪指标的简版基准（criterion 套件随后随 benches/ 落地）。
//!
//! expand (1+x+y+z+w)^20（设计 §7 的 expand 基准家族）与 sympy 单次对照。

use cas_expr::Kind;
use std::process::ExitCode;
use std::time::Instant;
use symcas::prelude::*;

pub(crate) fn run(_args: &[String]) -> ExitCode {
    let ctx = Context::new();
    let (x, y, z, w) = sym!(&ctx, x, y, z, w);

    let base = ctx.int(1) + x + y + z + w;
    for k in [10u32, 16, 20] {
        let e = base.clone().pow(k);
        let t = Instant::now();
        let g = ctx.expand(&e);
        let d = t.elapsed();
        let n = ctx.inspect(|i| match i.kind(g.raw_id()) {
            Kind::Add(args) => args.len(),
            _ => 1,
        });
        let expect = binom(k as u64 + 4, 4);
        println!(
            "expand (1+x+y+z+w)^{k:<2} : {n:6} 项（二项式校验 {}）{d:?}",
            if n == expect { "✓" } else { "✗" }
        );
    }

    // sympy 单次对照
    let out = std::process::Command::new("python3")
        .arg("-c")
        .arg(
            "import time, sympy; x,y,z,w=sympy.symbols('x y z w');\
             t=time.perf_counter(); e=sympy.expand((1+x+y+z+w)**20);\
             print(f'{time.perf_counter()-t:.3f}s {len(e.args)} 项')",
        )
        .output();
    match out {
        Ok(o) if o.status.success() => {
            println!(
                "sympy  (1+x+y+z+w)^20 : {}",
                String::from_utf8_lossy(&o.stdout).trim()
            );
        }
        _ => println!("sympy 对照不可用（跳过）"),
    }
    ExitCode::SUCCESS
}

fn binom(n: u64, k: u64) -> usize {
    let mut r = 1u64;
    for i in 0..k {
        r = r * (n - i) / (i + 1);
    }
    r as usize
}
