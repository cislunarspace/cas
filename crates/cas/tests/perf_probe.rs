//! 全管线性能探针（#[ignore]，手动运行）：build / plain / parse / alt 分相计时。
//! 运行：cargo test -p cas --test perf_probe -- --ignored --nocapture

use cas::prelude::*;
use cas_expr::test_gen::{self, Lcg, Shape};
use std::time::{Duration, Instant};

fn rand_shape(rng: &mut Lcg, depth: u32) -> Shape {
    if depth == 0 || rng.next_u64() % 3 == 0 {
        match rng.next_u64() % 4 {
            0 => Shape::Int((rng.next_u64() % 2001) as i64 - 1000),
            1 => Shape::Rat((rng.next_u64() % 101) as i64 - 50, rng.next_u64() % 20 + 1),
            2 => Shape::Float(match rng.next_u64() % 4 {
                0 => -0.0,
                1 => 1e300,
                2 => 1.5,
                _ => (rng.next_u64() % 2_000_001) as f64 / 1000.0 - 1000.0,
            }),
            _ => Shape::Sym((rng.next_u64() % test_gen::SYMBOL_NAMES.len() as u64) as usize),
        }
    } else {
        match rng.next_u64() % 4 {
            0 => {
                let exp = match rng.next_u64() % 3 {
                    0 => Shape::Int((rng.next_u64() % 65) as i64 - 32),
                    1 => Shape::Rat(rng.next_u64() as i64 % 20 + 1, rng.next_u64() % 20 + 1),
                    _ => rand_shape(rng, depth - 1),
                };
                Shape::Pow(Box::new(rand_shape(rng, depth - 1)), Box::new(exp))
            }
            1 => Shape::Fn {
                head: (rng.next_u64() % test_gen::HEAD_NAMES.len() as u64) as usize,
                args: vec![rand_shape(rng, depth - 1)],
            },
            k => {
                let n = (rng.next_u64() % 3 + 2) as usize;
                let items: Vec<Shape> = (0..n).map(|_| rand_shape(rng, depth - 1)).collect();
                if k == 2 {
                    Shape::Mul(items)
                } else {
                    Shape::Add(items)
                }
            }
        }
    }
}

#[test]
#[ignore]
fn 探针_全管线计时() {
    let mut rng = Lcg::new(20260911);
    let (mut t_build, mut t_plain, mut t_parse, mut t_alt) = (
        Duration::ZERO,
        Duration::ZERO,
        Duration::ZERO,
        Duration::ZERO,
    );
    let mut worst: Vec<(u128, String)> = Vec::new();
    for case in 0..3000u32 {
        let s = rand_shape(&mut rng, 6);
        let ctx = Context::new();

        let t = Instant::now();
        let e = test_gen::build(&ctx, &s);
        let d_build = t.elapsed();
        t_build += d_build;

        let t = Instant::now();
        let txt = cas::plain(&ctx, &e);
        let d_plain = t.elapsed();
        t_plain += d_plain;

        let t = Instant::now();
        let back = cas::parse(&ctx, &txt).expect("往返解析");
        let d_parse = t.elapsed();
        t_parse += d_parse;
        assert!(back == e);

        let t = Instant::now();
        let mut rng2 = Lcg::new(0x9e37_79b9 ^ case as u64);
        let alt = test_gen::build_alt(&ctx, &s, &mut rng2);
        let d_alt = t.elapsed();
        t_alt += d_alt;
        assert!(alt == e);

        // 单例耗时（此前误用累计值，已修）
        let us = (d_build + d_plain + d_parse + d_alt).as_micros();
        if us > 3000 {
            worst.push((
                us,
                format!(
                    "case {case}: 文本长度={} 节点={} shape={s:?}",
                    txt.len(),
                    ctx.node_count()
                ),
            ));
        }
    }
    worst.sort_by_key(|(us, _)| std::cmp::Reverse(*us));
    for (us, msg) in worst.iter().take(3) {
        println!("[{us}µs] {msg}");
    }
    println!("总计: build={t_build:?} plain={t_plain:?} parse={t_parse:?} alt={t_alt:?}");
}
