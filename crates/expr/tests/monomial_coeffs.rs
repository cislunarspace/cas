//! `Context::monomial_coeffs` 的行为测试：对拍规格来自
//! `math-tools/matlab/poly-operator/poly_operator/expr2poly.m`（嵌套 `coeffs(...,'All')`）。

use cas_expr::{CasErrorKind, Context};
use std::collections::BTreeMap;

fn map_of(pairs: Vec<(Vec<u32>, cas_expr::Expr)>, ctx: &Context) -> BTreeMap<Vec<u32>, String> {
    pairs
        .into_iter()
        .map(|(k, v)| (k, ctx.inspect(|i| i.dump(v.raw_id()))))
        .collect()
}

#[test]
fn expands_trinomial_cube() {
    let ctx = Context::new();
    let x = ctx.sym("x");
    let y = ctx.sym("y");
    let one = ctx.int(1);
    let e = ctx.expand(&ctx.pow(&ctx.add(&[one, x.clone(), y.clone()]), &ctx.int(3)));

    let pairs = ctx.monomial_coeffs(&e, &["x", "y"]).expect("应能提取");
    let m = map_of(pairs, &ctx);

    // (1+x+y)^3 = Σ 多项系数；共 10 个单项式
    assert_eq!(
        m.len(),
        10,
        "项数应为 10，实际 {:?}",
        m.keys().collect::<Vec<_>>()
    );

    // 逐幂次核对系数（dump 形如 int(6)）
    let expect = |exps: [u32; 2], coef: i64| {
        let k = vec![exps[0], exps[1]];
        let got = m.get(&k).unwrap_or_else(|| panic!("缺幂次 {k:?}"));
        assert_eq!(got, &format!("int({coef})"), "幂次 {k:?} 系数不符");
    };
    expect([3, 0], 1);
    expect([2, 1], 3);
    expect([2, 0], 3);
    expect([1, 2], 3);
    expect([1, 1], 6);
    expect([1, 0], 3);
    expect([0, 3], 1);
    expect([0, 2], 3);
    expect([0, 1], 3);
    expect([0, 0], 1);
}

#[test]
fn merges_like_terms_and_keeps_free_symbols_as_coefficients() {
    let ctx = Context::new();
    let x = ctx.sym("x");
    let a = ctx.sym("a");
    // a*x + 2*x + a^2  →  x 的系数为 a+2，另有一项常数 a^2
    let e = ctx.expand(&ctx.add(&[
        ctx.mul(&[a.clone(), x.clone()]),
        ctx.mul(&[ctx.int(2), x.clone()]),
        ctx.pow(&a, &ctx.int(2)),
    ]));
    let m = map_of(ctx.monomial_coeffs(&e, &["x"]).expect("应能提取"), &ctx);
    assert_eq!(m.len(), 2);
    assert!(
        m.get(&vec![1]).unwrap().contains("sym(a)"),
        "x 的系数应含自由符号 a，实际 {}",
        m[&vec![1]]
    );
    assert!(m.contains_key(&vec![0]), "应有常数项");
}

#[test]
fn rejects_non_polynomial_use_of_variable() {
    let ctx = Context::new();
    let x = ctx.sym("x");
    let e = ctx.call("sqrt", &[ctx.add(&[ctx.int(1), x.clone()])]);
    let err = ctx.monomial_coeffs(&e, &["x"]).expect_err("sqrt(x) 应报错");
    assert_eq!(err.kind, CasErrorKind::NonPolynomial);
}

#[test]
fn rejects_negative_exponent() {
    let ctx = Context::new();
    let x = ctx.sym("x");
    let e = ctx.pow(&x, &ctx.int(-2));
    let err = ctx.monomial_coeffs(&e, &["x"]).expect_err("x^-2 应报错");
    assert_eq!(err.kind, CasErrorKind::NegativeExponent);
}

#[test]
fn rejects_unknown_symbol() {
    let ctx = Context::new();
    let x = ctx.sym("x");
    let err = ctx
        .monomial_coeffs(&x, &["nope"])
        .expect_err("未知变量应报错");
    assert_eq!(err.kind, CasErrorKind::UnknownSymbol);
}

#[test]
fn matches_legendre_recurrence_shape() {
    // 勒让德多项式 P_3(x) = (5x^3 - 3x)/2：验证有理系数与幂次提取
    let ctx = Context::new();
    let x = ctx.sym("x");
    let p3 = ctx.mul(&[
        ctx.rational(1, 2).unwrap(),
        ctx.add(&[
            ctx.mul(&[ctx.int(5), ctx.pow(&x, &ctx.int(3))]),
            ctx.mul(&[ctx.int(-3), x.clone()]),
        ]),
    ]);
    let e = ctx.expand(&p3);
    let m = map_of(ctx.monomial_coeffs(&e, &["x"]).expect("应能提取"), &ctx);
    assert_eq!(m.len(), 2);
    // 规范形会把数值系数折叠（1/2·5 → 5/2），这正是要保留的行为
    assert_eq!(m[&vec![3]], "rat(5/2)");
    assert_eq!(m[&vec![1]], "rat(-3/2)");
}
