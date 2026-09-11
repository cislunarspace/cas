//! M1 硬门槛：print∘parse 往返 1000 例 + 跨上下文字节一致性。

use cas::prelude::*;
use cas_expr::test_gen;
use proptest::prelude::*;

proptest! {
    #![proptest_config(ProptestConfig::with_cases(1000))]

    /// 对随机规范形 e：parse(print(e)) 与 e 是同一节点；
    /// 且相同形状在另一个 Context 重建，plain 输出逐字节一致（D6）。
    #[test]
    fn print_parse_往返(s in test_gen::shape_strategy(5)) {
        let ctx = Context::new();
        let e = test_gen::build(&ctx, &s);
        let txt = cas::plain(&ctx, &e);
        let back = cas::parse(&ctx, &txt)
            .unwrap_or_else(|err| panic!("往返解析失败: {err}\n原文: {txt}\n形状: {s:?}"));
        prop_assert!(back == e, "往返差异:\n  原: {e:?}\n  回: {back:?}\n  文本: {txt}");

        let ctx2 = Context::new();
        let e2 = test_gen::build(&ctx2, &s);
        let txt2 = cas::plain(&ctx2, &e2);
        prop_assert_eq!(txt2, txt, "跨上下文输出漂移");
    }

    /// 连续两次解析产出同一节点（解析器自身的确定性）。
    #[test]
    fn 解析确定性(s in test_gen::shape_strategy(4)) {
        let ctx = Context::new();
        let e = test_gen::build(&ctx, &s);
        let txt = cas::plain(&ctx, &e);
        let p1 = cas::parse(&ctx, &txt).unwrap();
        let p2 = cas::parse(&ctx, &txt).unwrap();
        prop_assert_eq!(p1.raw_id(), p2.raw_id());
    }
}

#[test]
fn 黄金快照_规范形输出() {
    let ctx = Context::new();
    let (x, y, z) = sym!(&ctx, x, y, z);

    // 手工钉死一小组输出字节：任何全序/打印变更都会在此暴露（D6 快照防线）
    assert_eq!(cas::plain(&ctx, &(x.clone() + x.clone())), "2*x");
    assert_eq!(
        cas::plain(&ctx, &(x.clone() + y.clone() + z.clone())),
        "x + y + z"
    );
    assert_eq!(
        cas::plain(
            &ctx,
            &(ctx.int(3) * x.clone().pow(2) * y.clone()
                + ctx.int(3) * x.clone() * y.clone().pow(2))
        ),
        "3*x*y^2 + 3*x^2*y"
    );
    assert_eq!(
        cas::plain(
            &ctx,
            &((x.clone() + y.clone()).pow(3) - x.clone().pow(3) - y.clone().pow(3))
        ),
        "(x + y)^3 - x^3 - y^3"
    );
    assert_eq!(cas::plain(&ctx, &(y.clone() - x.clone())), "y - x");
    assert_eq!(
        cas::plain(&ctx, &(x.clone() * (x.clone() + y.clone()))),
        "x*(x + y)"
    );
    assert_eq!(
        cas::plain(
            &ctx,
            &(ctx.call("sin", std::slice::from_ref(&x))
                + ctx.call("sin", std::slice::from_ref(&x)))
        ),
        "2*sin(x)"
    );
    assert_eq!(
        cas::latex(&ctx, &(ctx.int(3) * x.clone() * y.clone())),
        "3 \\cdot x \\cdot y"
    );
}
