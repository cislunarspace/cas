//! 规范形不变量的属性测试（M1 硬门槛之一：构造路径无关，500 例）。

use cas_expr::Context;
use cas_expr::test_gen::{self, Lcg, shape_strategy};
use proptest::prelude::*;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(500))]

    /// 同一形状的两条构造路径（规范顺序 vs 洗牌/重分组/恒等包裹）
    /// 必须落在同一 arena 节点——hash-consing 规范形的直接检验。
    #[test]
    fn 构造路径无关(s in shape_strategy(5)) {
        let ctx = Context::new();
        let a = test_gen::build(&ctx, &s);
        let mut rng = Lcg::new(0x9e37_79b9 ^ (a.raw_id() as u64));
        let b = test_gen::build_alt(&ctx, &s, &mut rng);
        prop_assert!(a == b, "路径差异: a={a:?} b={b:?}");
        prop_assert_eq!(a.raw_id(), b.raw_id());
    }

    /// 加法与乘法的交换律在节点层面成立。
    #[test]
    fn 交换律(s1 in shape_strategy(3), s2 in shape_strategy(3)) {
        let ctx = Context::new();
        let a = test_gen::build(&ctx, &s1);
        let b = test_gen::build(&ctx, &s2);
        prop_assert!(ctx.add(&[a.clone(), b.clone()]) == ctx.add(&[b.clone(), a.clone()]));
        prop_assert!(ctx.mul(&[a.clone(), b.clone()]) == ctx.mul(&[b.clone(), a.clone()]));
    }

    /// 相同形状在两个独立 Context 中结构相等（跨上下文确定性）。
    #[test]
    fn 跨上下文结构相等(s in shape_strategy(4)) {
        let c1 = Context::new();
        let c2 = Context::new();
        let a = test_gen::build(&c1, &s);
        let b = test_gen::build(&c2, &s);
        prop_assert!(a == b, "跨上下文差异: a={a:?} b={b:?}");
    }

    /// 恒等式包裹在任意嵌套深度下都折叠回去。
    #[test]
    fn 恒等式包裹折叠(s in shape_strategy(4)) {
        let ctx = Context::new();
        let a = test_gen::build(&ctx, &s);
        let one = ctx.int(1);
        let wrapped = ctx.mul(&[ctx.add(&[ctx.pow(&a.clone(), &one), ctx.int(0)]), one]);
        prop_assert!(wrapped == a);
    }
}
