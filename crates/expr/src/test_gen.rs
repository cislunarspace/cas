//! 测试语料生成（feature = "test-gen"，仅测试依赖）。
//!
//! 提供：随机表达式形状 [`Shape`] 与 proptest 策略 [`shape_strategy`]；
//! 规范顺序构造 [`build`]；等价构造路径 [`build_alt`]（洗牌 + 恒等式包裹
//! + 重新分组）。后两者必须产出同一规范形——这是规范形不变量的直接测试。

use crate::{Context, Expr};
use proptest::prelude::*;

pub const SYMBOL_NAMES: [&str; 5] = ["x", "y", "z", "a", "b_2"];
pub const HEAD_NAMES: [&str; 4] = ["sin", "cos", "exp", "f"];

#[derive(Clone, Debug)]
pub enum Shape {
    Int(i64),
    Rat(i64, u64),
    Float(f64),
    Sym(usize),
    Fn { head: usize, args: Vec<Shape> },
    Pow(Box<Shape>, Box<Shape>),
    Mul(Vec<Shape>),
    Add(Vec<Shape>),
}

pub fn shape_strategy(max_depth: u32) -> impl Strategy<Value = Shape> {
    let leaf = prop_oneof![
        (-1000i64..=1000).prop_map(Shape::Int),
        (-50i64..=50, 1u64..=20).prop_map(|(n, d)| Shape::Rat(n, d)),
        prop_oneof![
            Just(-0.0f64),
            Just(0.0),
            Just(1.5),
            Just(-2.25),
            Just(1e300),
            (-1e6f64..1e6f64),
        ]
        .prop_map(Shape::Float),
        (0..SYMBOL_NAMES.len()).prop_map(Shape::Sym),
    ];
    // 字面指数限定小量级：随机语料若带 ±1000 的字面指数会折叠出几千比特的
    // 常数——那是 M2 数值基准的领地，不是规范形属性测试要覆盖的负载。
    let exp = prop_oneof![
        (-32i64..=32).prop_map(Shape::Int),
        (1u64..=20, 1u64..=20).prop_map(|(n, d)| Shape::Rat(n as i64, d)),
        leaf.clone(),
    ];
    leaf.clone().prop_recursive(max_depth, 64, 6, move |inner| {
        prop_oneof![
            (inner.clone(), exp.clone()).prop_map(|(b, e)| Shape::Pow(Box::new(b), Box::new(e))),
            (0..HEAD_NAMES.len(), inner.clone()).prop_map(|(head, a)| Shape::Fn {
                head,
                args: vec![a]
            }),
            proptest::collection::vec(inner.clone(), 2..=4).prop_map(Shape::Mul),
            proptest::collection::vec(inner, 2..=4).prop_map(Shape::Add),
        ]
    })
}

/// 确定性小随机源：不引入 rand 依赖，同种子同序列。
pub struct Lcg(u64);

impl Lcg {
    pub fn new(seed: u64) -> Self {
        Lcg(seed | 1)
    }

    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    pub fn shuffle<T>(&mut self, v: &mut [T]) {
        for i in (1..v.len()).rev() {
            let j = (self.next_u64() % (i as u64 + 1)) as usize;
            v.swap(i, j);
        }
    }
}

/// 按形状的规范顺序构造。
pub fn build(ctx: &Context, s: &Shape) -> Expr {
    match s {
        Shape::Int(v) => ctx.int(*v),
        Shape::Rat(n, d) => ctx.rational(*n, *d as i64).expect("分母非零"),
        Shape::Float(v) => ctx.float(*v),
        Shape::Sym(i) => ctx.sym(SYMBOL_NAMES[*i]),
        Shape::Fn { head, args } => {
            let built: Vec<Expr> = args.iter().map(|a| build(ctx, a)).collect();
            ctx.call(HEAD_NAMES[*head], &built)
        }
        Shape::Pow(b, e) => ctx.pow(&build(ctx, b), &build(ctx, e)),
        Shape::Mul(fs) => {
            let built: Vec<Expr> = fs.iter().map(|f| build(ctx, f)).collect();
            ctx.mul(&built)
        }
        Shape::Add(ts) => {
            let built: Vec<Expr> = ts.iter().map(|t| build(ctx, t)).collect();
            ctx.add(&built)
        }
    }
}

/// 等价构造路径：Add/Mul 子项洗牌 + 两两折叠改变结合顺序；子表达式按
/// 概率作恒等包裹（`*1`、`+0`、`^1`，规范形中全部折叠回去）。
pub fn build_alt(ctx: &Context, s: &Shape, rng: &mut Lcg) -> Expr {
    match s {
        Shape::Add(ts) => {
            let mut built: Vec<Expr> = ts.iter().map(|t| build_alt(ctx, t, rng)).collect();
            rng.shuffle(&mut built);
            let mut acc = built.pop().expect("Add 至少两项");
            while let Some(e) = built.pop() {
                acc = ctx.add(&[acc, e]);
            }
            maybe_wrap(ctx, acc, rng)
        }
        Shape::Mul(fs) => {
            let mut built: Vec<Expr> = fs.iter().map(|f| build_alt(ctx, f, rng)).collect();
            rng.shuffle(&mut built);
            let mut acc = built.pop().expect("Mul 至少两项");
            while let Some(e) = built.pop() {
                acc = ctx.mul(&[acc, e]);
            }
            maybe_wrap(ctx, acc, rng)
        }
        Shape::Fn { head, args } => {
            let built: Vec<Expr> = args.iter().map(|a| build_alt(ctx, a, rng)).collect();
            maybe_wrap(ctx, ctx.call(HEAD_NAMES[*head], &built), rng)
        }
        Shape::Pow(b, e) => maybe_wrap(
            ctx,
            ctx.pow(&build_alt(ctx, b, rng), &build_alt(ctx, e, rng)),
            rng,
        ),
        leaf => maybe_wrap(ctx, build(ctx, leaf), rng),
    }
}

fn maybe_wrap(ctx: &Context, e: Expr, rng: &mut Lcg) -> Expr {
    let one = ctx.int(1);
    match rng.next_u64() % 4 {
        0 => ctx.mul(&[e, one]),
        1 => ctx.add(&[e, ctx.int(0)]),
        2 => ctx.pow(&e, &one),
        _ => e,
    }
}
