//! cas 表达式层（P0-M1）：arena + hash-consing、规范形构造、确定性全序。
//!
//! # 核心不变量
//!
//! - **规范形即身份**：构造经 intern 表去重，同一 Context 内结构相等 ⇔
//!   索引相等（[`Expr::raw_id`] 相同），相等性 O(1)。
//! - **构造即规范**：`Add`/`Mul` 在构造器内完成扁平化、全序排序、精确数值
//!   折叠、同类项/同底幂合并；任何时刻取出的 `Expr` 都是规范形。
//! - **确定性**：全序见 [`order`] 模块文档（`ORDER_VERSION = 1`）；全库哈希
//!   用固定键（禁 `RandomState`），哈希迭代不进入输出。
//! - **浮点边界**：Float 仅作字面量（值恒非负，负号由系数 -1 携带），
//!   不参与任何算术折叠；`-0.0` 归一为 `0.0`；NaN 拒绝。
//!
//! arena 只增不减（P0 不回收，P4 评估压缩）；`Context` 按问题作用域创建。
//! 详见 `cas/DESIGN.md` D1/D6。

use std::cell::{Ref, RefCell};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fmt;
use std::rc::Rc;

mod canonical;
mod eval;
mod factor;
mod hash;
mod node;
mod order;
mod poly_bridge;
#[cfg(feature = "test-gen")]
pub mod test_gen;
mod transform;

use hash::FxBuild;
pub(crate) use node::Node;

pub use cas_domain::{Integer, Rational};

/// arena 与全部查表的宿主。类型名公开是为 `IntoExpr` 等签名可提及；
/// 字段全部 crate 私有，外部无法构造或窥视。
pub struct Inner {
    pub(crate) nodes: Vec<Node>,
    pub(crate) hashes: Vec<u64>,
    pub(crate) intern: HashMap<u64, Vec<u32>, FxBuild>,
    pub(crate) args_slab: Vec<u32>,
    pub(crate) syms: HashMap<Box<str>, u32, FxBuild>,
    pub(crate) sym_names: Vec<Box<str>>,
    pub(crate) fn_heads: HashMap<Box<str>, u32, FxBuild>,
    pub(crate) fn_names: Vec<Box<str>>,
}

impl Inner {
    fn new() -> Self {
        Inner {
            nodes: vec![Node::Int(Integer::zero())], // id 0 占位，永不引用
            hashes: vec![0],
            intern: HashMap::with_hasher(FxBuild::default()),
            args_slab: Vec::new(),
            syms: HashMap::with_hasher(FxBuild::default()),
            sym_names: Vec::new(),
            fn_heads: HashMap::with_hasher(FxBuild::default()),
            fn_names: Vec::new(),
        }
    }
}

/// 表达式构造与求值的宿主。按问题作用域创建；`Clone` 得到共享同一 arena
/// 的另一个入口（不复制数据）。
#[derive(Clone)]
pub struct Context {
    inner: Rc<RefCell<Inner>>,
}

impl Context {
    pub fn new() -> Self {
        Context {
            inner: Rc::new(RefCell::new(Inner::new())),
        }
    }

    fn wrap(&self, id: u32) -> Expr {
        Expr {
            ctx: Rc::clone(&self.inner),
            id,
        }
    }

    fn with(&self, f: impl FnOnce(&mut Inner) -> u32) -> Expr {
        let id = f(&mut self.inner.borrow_mut());
        self.wrap(id)
    }

    /// 声明符号。名字须匹配 `[A-Za-z_][A-Za-z0-9_]*`（违反属编程错误，直接 panic）。
    pub fn sym(&self, name: &str) -> Expr {
        validate_name(name, "符号");
        self.with(|c| c.sym_node(name))
    }

    pub fn int(&self, v: i64) -> Expr {
        self.with(|c| c.lit_int(v))
    }

    pub fn integer(&self, v: &Integer) -> Expr {
        self.with(|c| c.lit_integer(v))
    }

    /// 浮点字面量。负值返回 `-1 * |v|` 的规范形（Float 节点恒非负）；
    /// NaN 属编程错误。
    pub fn float(&self, v: f64) -> Expr {
        assert!(!v.is_nan(), "Float 字面量禁止 NaN");
        if v < 0.0 {
            self.with(|c| {
                let f = c.lit_float(-v);
                let neg = c.lit_int(-1);
                c.make_mul(&[neg, f])
            })
        } else {
            self.with(|c| c.lit_float(v))
        }
    }

    /// 有理数字面量；分母为零返回 `None`。
    pub fn rational(&self, n: i64, d: i64) -> Option<Expr> {
        let r = Rational::from_ints(&Integer::from_i64(n), &Integer::from_i64(d))?;
        Some(self.with(|c| c.lit_rational(&r)))
    }

    /// n 元加法（运算符重载的批量化入口）。
    pub fn add(&self, args: &[Expr]) -> Expr {
        let ids: Vec<u32> = args.iter().map(|e| e.id).collect();
        self.with(|c| c.make_add(&ids))
    }

    /// n 元乘法。
    pub fn mul(&self, args: &[Expr]) -> Expr {
        let ids: Vec<u32> = args.iter().map(|e| e.id).collect();
        self.with(|c| c.make_mul(&ids))
    }

    pub fn pow(&self, base: &Expr, exp: &Expr) -> Expr {
        self.with(|c| c.make_pow(base.id, exp.id))
    }

    /// 函数应用。参数保持语义次序（不排序）；head 命名规则同符号。
    pub fn call(&self, head: &str, args: &[Expr]) -> Expr {
        validate_name(head, "函数");
        let ids: Vec<u32> = args.iter().map(|e| e.id).collect();
        self.with(|c| c.fn_node(head, &ids))
    }

    /// 表达式全序（同一 Context 内；跨 Context 用 `Expr::eq` 的结构比较）。
    pub fn cmp_expr(&self, a: &Expr, b: &Expr) -> Ordering {
        let inner = self.inner.borrow();
        inner.cmp_ids(a.id, b.id)
    }

    pub fn eq_expr(&self, a: &Expr, b: &Expr) -> bool {
        self.cmp_expr(a, b) == Ordering::Equal
    }

    /// arena 规模（节点数；诊断与后续基准用）。
    pub fn node_count(&self) -> usize {
        self.inner.borrow().nodes.len()
    }

    /// 只读检视。持有 Inspector 期间不可构造新表达式（构造方独占借用 arena）。
    pub fn inspect<R>(&self, f: impl FnOnce(&Inspector<'_>) -> R) -> R {
        f(&Inspector {
            inner: self.inner.borrow(),
        })
    }

    // ── L1 显式变换与求值 ──────────────────────────────────────

    /// 展开（L1）：分配乘积 over 和式、展开非负整数幂（带规模守卫）、
    /// 函数参数内部展开。每步经规范形构造器，输出为规范形。
    pub fn expand(&self, e: &Expr) -> Expr {
        self.with(|c| c.expand_at(e.id, 0))
    }

    /// 有理函数约化（L1）：分子分母的多项式公因子（gcd）约去。
    /// 非多项式因子原样保留；无可约化时返回与输入同一节点。
    pub fn cancel(&self, e: &Expr) -> Expr {
        self.with(|c| c.cancel_at(e.id))
    }

    /// 多项式因式分解（M4/M5）：`cont · Π 因子^重数`，重建为规范形。
    /// 1/2 变元完全分解，≥3 变元部分（内容+平方自由）；非多项式输入返回原节点。
    pub fn factor(&self, e: &Expr) -> Expr {
        self.with(|c| c.factor_at(e.id))
    }

    /// 代换：按符号名替换子表达式（替换值须属同一 Context）。
    /// 重建经规范形构造器，`subst(e, x→x)` 与 `e` 同节点。
    pub fn subst(&self, e: &Expr, map: &[(&str, Expr)]) -> Expr {
        let m = {
            let inner = self.inner.borrow();
            map.iter()
                .filter_map(|(name, val)| inner.syms.get(*name).map(|&s| (s, val.id)))
                .collect::<std::collections::HashMap<u32, u32>>()
        };
        self.with(|c| c.subst_at(e.id, &m, 0))
    }

    /// ℚ 上精确求值。Float 字面量与函数节点返回 `None`（D2 域边界）。
    pub fn eval_rational(&self, e: &Expr, vals: &[(&str, Rational)]) -> Option<Rational> {
        let inner = self.inner.borrow();
        let m = vals
            .iter()
            .filter_map(|(n, v)| inner.syms.get(*n).map(|&s| (s, v.clone())))
            .collect::<std::collections::HashMap<u32, Rational>>();
        eval::eval_rational_at(&inner, e.id, &m, 0)
    }

    /// f64 数值求值（差分测试的语义判据通道）。支持初等函数头。
    pub fn eval_float(&self, e: &Expr, vals: &[(&str, f64)]) -> Option<f64> {
        let inner = self.inner.borrow();
        let m = vals
            .iter()
            .filter_map(|(n, v)| inner.syms.get(*n).map(|&s| (s, *v)))
            .collect::<std::collections::HashMap<u32, f64>>();
        eval::eval_float_at(&inner, e.id, &m, 0)
    }
}

impl Default for Context {
    fn default() -> Self {
        Self::new()
    }
}

fn validate_name(name: &str, what: &str) {
    let mut chars = name.chars();
    let valid = matches!(chars.next(), Some(c) if c.is_ascii_alphabetic() || c == '_')
        && chars.all(|c| c.is_ascii_alphanumeric() || c == '_');
    assert!(
        valid && !name.is_empty(),
        "{what}名非法: {name:?}，须匹配 [A-Za-z_][A-Za-z0-9_]*"
    );
}

/// 表达式句柄：arena 索引 + 共享上下文引用。
///
/// `Clone` 廉价（引用计数 + 4 字节索引）。同一 Context 内结构相等 ⇔ 索引相等；
/// 跨 Context 的 `==` 走结构比较。不实现 `Hash`/`Ord`：句柄的序只在所属
/// Context 内有意义，用 [`Context::cmp_expr`]。
#[derive(Clone)]
pub struct Expr {
    ctx: Rc<RefCell<Inner>>,
    id: u32,
}

impl Expr {
    /// arena 索引。仅在同一 Context 内有意义（相等判据/诊断用）。
    pub fn raw_id(&self) -> u32 {
        self.id
    }

    fn bin(self, other: Expr, f: impl FnOnce(&mut Inner, u32, u32) -> u32) -> Expr {
        let ctx = Rc::clone(&self.ctx);
        let id = f(&mut ctx.borrow_mut(), self.id, other.id);
        Expr { ctx, id }
    }

    /// 幂。指数可传 `Expr`/`&Expr`/整数字面量/浮点字面量。
    pub fn pow<E: IntoExpr>(self, exp: E) -> Expr {
        let ctx = Rc::clone(&self.ctx);
        let mut inner = ctx.borrow_mut();
        let eid = exp.into_id(&mut inner);
        let id = inner.make_pow(self.id, eid);
        drop(inner);
        Expr { ctx, id }
    }
}

/// 四种（自有/引用 × 自有/引用）组合的二元运算符转发。
macro_rules! forward_binop {
    ($tr:ident, $method:ident, $inner:expr) => {
        impl std::ops::$tr<Expr> for Expr {
            type Output = Expr;
            fn $method(self, rhs: Expr) -> Expr {
                self.bin(rhs, $inner)
            }
        }
        impl std::ops::$tr<&Expr> for Expr {
            type Output = Expr;
            fn $method(self, rhs: &Expr) -> Expr {
                self.bin(rhs.clone(), $inner)
            }
        }
        impl std::ops::$tr<Expr> for &Expr {
            type Output = Expr;
            fn $method(self, rhs: Expr) -> Expr {
                self.clone().bin(rhs, $inner)
            }
        }
        impl std::ops::$tr<&Expr> for &Expr {
            type Output = Expr;
            fn $method(self, rhs: &Expr) -> Expr {
                self.clone().bin(rhs.clone(), $inner)
            }
        }
    };
}

forward_binop!(Add, add, |c, a, b| c.make_add(&[a, b]));
forward_binop!(Sub, sub, |c, a, b| {
    let neg1 = c.lit_int(-1);
    let neg = c.make_mul(&[neg1, b]);
    c.make_add(&[a, neg])
});
forward_binop!(Mul, mul, |c, a, b| c.make_mul(&[a, b]));
forward_binop!(Div, div, |c, a, b| {
    let neg1 = c.lit_int(-1);
    let inv = c.make_pow(b, neg1);
    c.make_mul(&[a, inv])
});

impl std::ops::Neg for Expr {
    type Output = Expr;
    fn neg(self) -> Expr {
        let ctx = Rc::clone(&self.ctx);
        let id = {
            let mut c = ctx.borrow_mut();
            let neg1 = c.lit_int(-1);
            c.make_mul(&[neg1, self.id])
        };
        Expr { ctx, id }
    }
}

impl std::ops::Neg for &Expr {
    type Output = Expr;
    fn neg(self) -> Expr {
        (*self).clone().neg()
    }
}

impl PartialEq for Expr {
    fn eq(&self, other: &Self) -> bool {
        if Rc::ptr_eq(&self.ctx, &other.ctx) {
            return self.id == other.id;
        }
        let x = self.ctx.borrow();
        let y = other.ctx.borrow();
        order::deep_eq(&x, self.id, &y, other.id)
    }
}

impl Eq for Expr {}

impl fmt::Debug for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(
            &Inspector {
                inner: self.ctx.borrow(),
            }
            .dump(self.id),
        )
    }
}

/// [`Expr::pow`] 可接受的指数类型。
pub trait IntoExpr {
    fn into_id(self, c: &mut Inner) -> u32;
}

impl IntoExpr for Expr {
    fn into_id(self, _c: &mut Inner) -> u32 {
        self.id
    }
}

impl IntoExpr for &Expr {
    fn into_id(self, _c: &mut Inner) -> u32 {
        self.id
    }
}

impl IntoExpr for i64 {
    fn into_id(self, c: &mut Inner) -> u32 {
        c.lit_int(self)
    }
}

impl IntoExpr for i32 {
    fn into_id(self, c: &mut Inner) -> u32 {
        c.lit_int(self as i64)
    }
}

impl IntoExpr for u32 {
    fn into_id(self, c: &mut Inner) -> u32 {
        c.lit_int(self as i64)
    }
}

impl IntoExpr for f64 {
    fn into_id(self, c: &mut Inner) -> u32 {
        assert!(self >= 0.0, "pow 的浮点指数须非负（负浮点请用 Expr 形态）");
        c.lit_float(self)
    }
}

/// 只读检视器：在 [`Context::inspect`] 闭包内遍历表达式结构。
pub struct Inspector<'a> {
    inner: Ref<'a, Inner>,
}

/// 节点的只读视图。复合节点给子节点索引，调用方继续经 Inspector 递归。
#[derive(Debug)]
pub enum Kind<'a> {
    Int(&'a Integer),
    Rat(&'a Rational),
    Float(f64),
    Sym(&'a str),
    Fn { head: &'a str, args: &'a [u32] },
    Pow { base: u32, exp: u32 },
    Mul(&'a [u32]),
    Add(&'a [u32]),
}

impl Inspector<'_> {
    pub fn kind(&self, id: u32) -> Kind<'_> {
        match &self.inner.nodes[id as usize] {
            Node::Int(v) => Kind::Int(v),
            Node::Rat(v) => Kind::Rat(v),
            Node::Float { bits, .. } => Kind::Float(f64::from_bits(*bits)),
            Node::Sym(s) => Kind::Sym(&self.inner.sym_names[*s as usize]),
            Node::Fn { head, args } => Kind::Fn {
                head: &self.inner.fn_names[*head as usize],
                args: self.inner.node_args(*args),
            },
            Node::Pow { base, exp } => Kind::Pow {
                base: *base,
                exp: *exp,
            },
            Node::Mul { args } => Kind::Mul(self.inner.node_args(*args)),
            Node::Add { args } => Kind::Add(self.inner.node_args(*args)),
        }
    }

    pub fn cmp(&self, a: u32, b: u32) -> Ordering {
        self.inner.cmp_ids(a, b)
    }

    /// 诊断转储：`add[mul[int(3), pow[sym(x), int(2)]], ...]` 风格。
    pub fn dump(&self, id: u32) -> String {
        let mut s = String::new();
        self.dump_at(id, &mut s, 0);
        s
    }

    fn dump_at(&self, id: u32, out: &mut String, depth: u32) {
        assert!(depth <= 10_000, "表达式嵌套过深");
        let sub = |c: &Self, cid: u32, o: &mut String, d: u32| c.dump_at(cid, o, d);
        match self.kind(id) {
            Kind::Int(v) => out.push_str(&format!("int({v})")),
            Kind::Rat(r) => out.push_str(&format!("rat({r})")),
            Kind::Float(v) => out.push_str(&format!("flt({v})")),
            Kind::Sym(n) => out.push_str(&format!("sym({n})")),
            Kind::Fn { head, args } => {
                out.push_str(&format!("fn({head}, ["));
                for (i, &a) in args.iter().enumerate() {
                    if i > 0 {
                        out.push_str(", ");
                    }
                    sub(self, a, out, depth + 1);
                }
                out.push_str("])");
            }
            Kind::Pow { base, exp } => {
                out.push_str("pow(");
                sub(self, base, out, depth + 1);
                out.push_str(", ");
                sub(self, exp, out, depth + 1);
                out.push(')');
            }
            Kind::Mul(args) => {
                out.push_str("mul[");
                for (i, &a) in args.iter().enumerate() {
                    if i > 0 {
                        out.push_str(", ");
                    }
                    sub(self, a, out, depth + 1);
                }
                out.push(']');
            }
            Kind::Add(args) => {
                out.push_str("add[");
                for (i, &a) in args.iter().enumerate() {
                    if i > 0 {
                        out.push_str(", ");
                    }
                    sub(self, a, out, depth + 1);
                }
                out.push(']');
            }
        }
    }
}

/// 批量声明符号：`let (x, y, z) = sym!(&ctx, x, y, z);`（单个名字直接返回 Expr）。
#[macro_export]
macro_rules! sym {
    ($ctx:expr, $name:ident) => {
        $ctx.sym(stringify!($name))
    };
    ($ctx:expr, $($name:ident),+ $(,)?) => {
        ($($ctx.sym(stringify!($name))),+)
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn 规范形_同类合并与数值折叠() {
        let ctx = Context::new();
        let x = ctx.sym("x");
        let y = ctx.sym("y");

        // x + x 与 2*x 是同一节点
        let a = x.clone() + x.clone();
        let b = ctx.int(2) * x.clone();
        assert!(a == b && a.raw_id() == b.raw_id());

        // (a+b)+c 扁平化，c+(b+a) 同一规范形
        let c = ctx.sym("c");
        let p = (x.clone() + y.clone()) + c.clone();
        let q = c + (y + x);
        assert_eq!(p.raw_id(), q.raw_id());

        // 数值折叠
        let n = ctx.int(2) + ctx.int(3);
        assert!(
            ctx.inspect(
                |i| matches!(i.kind(n.raw_id()), Kind::Int(v) if *v == Integer::from_i64(5))
            )
        );
    }

    #[test]
    fn 规范形_同底幂合并() {
        let ctx = Context::new();
        let x = ctx.sym("x");
        let a = ctx.sym("a");
        let b = ctx.sym("b");

        // x * x == x^2
        let p = x.clone() * x.clone();
        let q = x.clone().pow(2);
        assert_eq!(p.raw_id(), q.raw_id());

        // x^2 * x^3 == x^5
        let p = x.clone().pow(2) * x.clone().pow(3);
        let q = x.clone().pow(5);
        assert_eq!(p.raw_id(), q.raw_id());

        // x^a * x^b == x^(a+b)
        let p = x.clone().pow(&a) * x.clone().pow(&b);
        let q = x.clone().pow(&a + &b);
        assert_eq!(p.raw_id(), q.raw_id());

        // (x^a)^2 == x^(2*a)；(x^a)^b 保持嵌套（非整数外幂）
        let p = x.clone().pow(&a).pow(2);
        let q = x.clone().pow(&ctx.int(2) * &a);
        assert_eq!(p.raw_id(), q.raw_id());
        let nested = x.clone().pow(&a).pow(&b);
        assert!(ctx.inspect(|i| matches!(i.kind(nested.raw_id()), Kind::Pow { .. })));
    }

    #[test]
    fn 规范形_数值幂折叠与守卫() {
        let ctx = Context::new();
        let p = ctx.int(2).pow(100);
        assert!(ctx.inspect(|i| matches!(i.kind(p.raw_id()), Kind::Int(v)
            if *v == Integer::parse("1267650600228229401496703205376").unwrap())));

        // 负整数幂 → 有理数
        let p = ctx.int(2).pow(-2);
        assert!(ctx.inspect(|i| matches!(i.kind(p.raw_id()), Kind::Rat(_))));

        // 0^0 = 1，0^3 = 0，0^-1 无定义：保留 Pow 节点
        assert_eq!(ctx.int(0).pow(0).raw_id(), ctx.int(1).raw_id());
        assert_eq!(ctx.int(0).pow(3).raw_id(), ctx.int(0).raw_id());
        let zero_neg_pow = ctx.int(0).pow(-1);
        assert!(ctx.inspect(|i| matches!(i.kind(zero_neg_pow.raw_id()), Kind::Pow { .. })));

        // 规模守卫：超大指数不折叠
        let huge = ctx.int(2).pow(1 << 30);
        assert!(ctx.inspect(|i| matches!(i.kind(huge.raw_id()), Kind::Pow { .. })));
    }

    #[test]
    fn 规范形_除法与负浮点() {
        let ctx = Context::new();
        let x = ctx.sym("x");

        // x / 2 == 1/2 * x
        let p = x.clone() / ctx.int(2);
        let q = ctx.mul(&[ctx.rational(1, 2).unwrap(), x.clone()]);
        assert_eq!(p.raw_id(), q.raw_id());

        // 负浮点 = -1 * |v|
        let p = ctx.float(-2.5);
        let q = -ctx.float(2.5);
        assert_eq!(p.raw_id(), q.raw_id());

        // x/x == 1（sympy 同约定的退化点）
        assert_eq!((x.clone() / x.clone()).raw_id(), ctx.int(1).raw_id());
    }

    #[test]
    fn 全序_数值与符号() {
        let ctx = Context::new();
        let x = ctx.sym("x");
        let y = ctx.sym("y");
        let two = ctx.int(2);
        let half = ctx.rational(1, 2).unwrap();
        let three = ctx.int(3);
        let f1 = ctx.float(1.5);

        // 精确数按值：1/2 < 2 < 3；Float 排在精确数之后
        assert_eq!(ctx.cmp_expr(&half, &two), Ordering::Less);
        assert_eq!(ctx.cmp_expr(&two, &three), Ordering::Less);
        assert_eq!(ctx.cmp_expr(&three, &f1), Ordering::Less);

        // 符号按名字节序，与创建先后无关
        let big = ctx.sym("zz");
        assert_eq!(ctx.cmp_expr(&x, &y), Ordering::Less);
        assert_eq!(ctx.cmp_expr(&y, &big), Ordering::Less);

        // 数值在 Mul 中居首
        let m = ctx.int(3) * x.clone() * y.clone();
        let first = ctx.inspect(|i| match i.kind(m.raw_id()) {
            Kind::Mul(args) => args[0],
            k => panic!("期望 Mul: {k:?}"),
        });
        assert!(ctx.cmp_expr(&ctx.wrap(first), &x) == Ordering::Less);
    }

    #[test]
    fn 全序_复合节点字典序() {
        let ctx = Context::new();
        let x = ctx.sym("x");
        let y = ctx.sym("y");

        let x2 = x.clone().pow(2);
        let xy = x.clone() * y.clone();
        let sum = x.clone() + y.clone();

        // Sym < Pow < Mul < Add
        assert_eq!(ctx.cmp_expr(&x, &x2), Ordering::Less);
        assert_eq!(ctx.cmp_expr(&x2, &xy), Ordering::Less);
        assert_eq!(ctx.cmp_expr(&xy, &sum), Ordering::Less);
    }

    // ── L1 变换与求值 ─────────────────────────────────────────

    #[test]
    fn 展开_黄金快照() {
        // 断言用与手工构造表达式的节点恒等（比字符串更强）
        let ctx = Context::new();
        let x = ctx.sym("x");
        let y = ctx.sym("y");

        let e = (x.clone() + y.clone()) * (x.clone() - y.clone());
        let g = ctx.expand(&e);
        let expect = x.clone().pow(2) - y.clone().pow(2);
        assert_eq!(g.raw_id(), expect.raw_id());

        let e = (x.clone() + y.clone()).pow(2);
        let g = ctx.expand(&e);
        let expect = x.clone().pow(2) + y.clone().pow(2) + ctx.int(2) * x.clone() * y.clone();
        assert_eq!(g.raw_id(), expect.raw_id());

        // 立方：交叉项系数正确（含 3 与 1/3 的往返核对）
        let e = (x.clone() + y.clone()).pow(3);
        let g = ctx.expand(&e);
        let expect = x.clone().pow(3)
            + y.clone().pow(3)
            + ctx.int(3) * x.clone() * y.clone().pow(2)
            + ctx.int(3) * x.clone().pow(2) * y.clone();
        assert_eq!(g.raw_id(), expect.raw_id());
    }

    #[test]
    fn 展开_语义一致与幂等() {
        let ctx = Context::new();
        let x = ctx.sym("x");
        let y = ctx.sym("y");
        let z = ctx.sym("z");

        let e = ((x.clone() + y.clone()).pow(3) - z.clone() * (x.clone() + y.clone()))
            * (x.clone() - z.clone());
        let g = ctx.expand(&e);
        // 幂等
        let g2 = ctx.expand(&g);
        assert_eq!(g.raw_id(), g2.raw_id());
        // 与原式在随机点精确同值
        for (xv, yv, zv) in [(2, -3, 5), (-1, 4, 7), (0, 9, -11)] {
            let pts = [
                ("x", Rational::from_integer(&Integer::from_i64(xv))),
                ("y", Rational::from_integer(&Integer::from_i64(yv))),
                ("z", Rational::from_integer(&Integer::from_i64(zv))),
            ];
            let a = ctx.eval_rational(&e, &pts);
            let b = ctx.eval_rational(&g, &pts);
            assert_eq!(a, b, "展开改变语义");
        }

        // （规模守卫在展开_四元二十次幂与 huge 用例覆盖）
        let huge = (x.clone() + y.clone()).pow(20_000);
        assert!(ctx.inspect(|i| matches!(i.kind(huge.raw_id()), Kind::Pow { .. })));
    }

    #[test]
    fn 展开_四元二十次幂() {
        // 设计基准负载：项数 C(23,3) = 1771
        let ctx = Context::new();
        let (x, y, z, w) = crate::sym!(&ctx, x, y, z, w);
        let base = ctx.int(1) + x + y + z + w;
        let g = ctx.expand(&base.pow(20));
        let n = ctx.inspect(|i| match i.kind(g.raw_id()) {
            Kind::Add(args) => args.len(),
            _ => 0,
        });
        assert_eq!(n, 10_626); // C(24,4)
    }

    #[test]
    fn 代换() {
        let ctx = Context::new();
        let x = ctx.sym("x");
        let y = ctx.sym("y");

        let e = x.clone() + y.clone();
        let g = ctx.subst(&e, &[("x", x.clone().pow(2))]);
        let expect = x.clone().pow(2) + y.clone();
        assert_eq!(g.raw_id(), expect.raw_id());

        // 恒等代换回原节点（intern 去重）
        let same = ctx.subst(&e, &[("x", x.clone())]);
        assert_eq!(same.raw_id(), e.raw_id());

        // 对换
        let p = x.clone() * y.clone();
        let q = ctx.subst(&p, &[("x", y.clone()), ("y", x.clone())]);
        assert_eq!(q.raw_id(), p.raw_id()); // x*y 交换后仍同节点
    }

    #[test]
    fn 展开快慢路径同构() {
        // M3 快路径（poly 域）与朴素 DAG 分配必须产出同一节点
        let ctx = Context::new();
        let (x, y, z, w) = crate::sym!(&ctx, x, y, z, w);
        let syms = [x.clone(), y.clone(), z.clone(), w.clone()];

        let mut xs = 777u64;
        let mut nxt = move || {
            xs ^= xs << 13;
            xs ^= xs >> 7;
            xs ^= xs << 17;
            xs
        };
        fn gen_poly(
            ctx: &Context,
            syms: &[Expr],
            nxt: &mut impl FnMut() -> u64,
            depth: u32,
        ) -> Expr {
            if depth == 0 || nxt() % 3 == 0 {
                match nxt() % 3 {
                    0 => ctx.int((nxt() % 19) as i64 - 9),
                    1 => ctx
                        .rational((nxt() % 15) as i64 - 7, (nxt() % 8) as i64 + 2)
                        .unwrap(),
                    _ => syms[(nxt() % syms.len() as u64) as usize].clone(),
                }
            } else {
                match nxt() % 3 {
                    0 => gen_poly(ctx, syms, nxt, depth - 1) + gen_poly(ctx, syms, nxt, depth - 1),
                    1 => gen_poly(ctx, syms, nxt, depth - 1) * gen_poly(ctx, syms, nxt, depth - 1),
                    _ => gen_poly(ctx, syms, nxt, depth - 1).pow((nxt() % 5) as i64),
                }
            }
        }
        let slow_of = |ctx: &Context, e: &Expr| -> u32 {
            let mut inner = ctx.inner.borrow_mut();
            inner.expand_impl(e.raw_id(), 0, false)
        };

        for _ in 0..60 {
            let e = gen_poly(&ctx, &syms, &mut nxt, 4);
            let fast = ctx.expand(&e);
            let slow = slow_of(&ctx, &e);
            assert_eq!(fast.raw_id(), slow, "快慢路径不同构: {e:?}");
        }

        // 混合表达式（函数/浮点使快路径部分让位）仍须同构
        let e = ctx.call("sin", std::slice::from_ref(&x)) * (y.clone() + z.clone()).pow(5)
            + ctx.float(1.5) * (x.clone() + w.clone()).pow(3);
        assert_eq!(ctx.expand(&e).raw_id(), slow_of(&ctx, &e));
    }

    #[test]
    fn 因式分解() {
        let ctx = Context::new();
        let x = ctx.sym("x");

        // x^2 − 1 → (x − 1)(x + 1)（节点恒等）
        let e = x.clone().pow(2) - ctx.int(1);
        let g = ctx.factor(&e);
        let expect = (x.clone() - ctx.int(1)) * (x.clone() + ctx.int(1));
        assert_eq!(g.raw_id(), expect.raw_id());

        // (x^2 − 1)^2 → (x − 1)^2 (x + 1)^2
        let e = x.clone().pow(2) - ctx.int(1);
        let e = e.pow(2);
        let g = ctx.factor(&e);
        let expect = (x.clone() - ctx.int(1)).pow(2) * (x.clone() + ctx.int(1)).pow(2);
        assert_eq!(g.raw_id(), expect.raw_id());

        // x^4 + 4 在 ℚ 上可约（Sophie Germain）
        let e = x.clone().pow(4) + ctx.int(4);
        let g = ctx.factor(&e);
        let expect = (x.clone().pow(2) - ctx.int(2) * x.clone() + ctx.int(2))
            * (x.clone().pow(2) + ctx.int(2) * x.clone() + ctx.int(2));
        assert_eq!(g.raw_id(), expect.raw_id());

        // x^4 + 1 不可约：重建后应与原式同节点
        let e = x.clone().pow(4) + ctx.int(1);
        let g = ctx.factor(&e);
        let expect = x.clone().pow(4) + ctx.int(1);
        assert_eq!(g.raw_id(), expect.raw_id());

        // 双变元：x^2 − y^2 → 两个一次因子（节点形态随 pp 归一约定，
        // 断言用 dump 形态 + 随机点语义恒等）
        let y = ctx.sym("y");
        let e = x.clone().pow(2) - y.clone().pow(2);
        let g = ctx.factor(&e);
        let s = ctx.inspect(|i| i.dump(g.raw_id()));
        assert_eq!(
            s,
            "mul[int(-1), add[sym(x), sym(y)], add[sym(y), mul[int(-1), sym(x)]]]"
        );
        for (xv, yv) in [(3, 2), (-5, 7)] {
            let pts = [
                ("x", Rational::from_integer(&Integer::from_i64(xv))),
                ("y", Rational::from_integer(&Integer::from_i64(yv))),
            ];
            assert_eq!(ctx.eval_rational(&e, &pts), ctx.eval_rational(&g, &pts));
        }
    }

    #[test]
    fn 精确与数值求值() {
        let ctx = Context::new();
        let x = ctx.sym("x");

        let e = ctx.rational(3, 2).unwrap() * x.clone() + ctx.int(2);
        let v = ctx.eval_rational(&e, &[("x", Rational::from_integer(&Integer::from_i64(4)))]);
        assert_eq!(v.map(|r| r.to_string()), Some("8".to_string()));

        // 负幂
        let e = x.clone().pow(-2);
        let v = ctx.eval_rational(&e, &[("x", Rational::from_integer(&Integer::from_i64(2)))]);
        assert_eq!(v.map(|r| r.to_string()), Some("1/4".to_string()));

        // Float 与函数不进精确求值
        assert_eq!(ctx.eval_rational(&ctx.float(1.5), &[]), None);
        let s = ctx.call("sin", std::slice::from_ref(&x));
        assert_eq!(ctx.eval_rational(&s, &[("x", Rational::zero())]), None);

        // 数值求值：初等函数 + 与展开一致性
        let sv = ctx.eval_float(&s, &[("x", 0.5)]);
        assert!((sv.unwrap() - 0.5f64.sin()).abs() < 1e-15);
        let e = (x.clone() + ctx.int(1)).pow(10);
        let g = ctx.expand(&e);
        let pts = [("x", 1.7f64)];
        let a = ctx.eval_float(&e, &pts).unwrap();
        let b = ctx.eval_float(&g, &pts).unwrap();
        assert!((a - b).abs() < 1e-9 * a.abs().max(1.0));

        // log 非正数 → None
        let lg = ctx.call("log", std::slice::from_ref(&x));
        assert_eq!(ctx.eval_float(&lg, &[("x", -1.0)]), None);
    }
}
