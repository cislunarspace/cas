//! 假设系统（D5，对标 MATLAB `assume`，P1 后期）。
//!
//! - 闭谓词集（9 个）：real / rational / integer / even / odd /
//!   positive / negative / nonzero / finite，位集表示；
//! - 构造时闭包补全（`even⇒integer⇒rational⇒real`、`positive⇒
//!   real∧nonzero∧¬negative` 等）并检冲突（`positive∧negative` 为
//!   编程错误）；
//! - 绑定：同一 Context 内同名符号**一次性设定**，重设不同值报错；
//! - 查询：`query(e, P) → True | False | Unknown`，自底向上按头函数
//!   事实表传播（`exp(x)>0`、`x^2≥0`（x real）等），遇 Unknown 截断；
//! - 与化简的耦合只有一条路：规则 guard 调查询，仅 `True` 触发
//!   （宁可少化简，不可错化简）。不做 SAT / 一阶逻辑引擎。

use crate::Inner;
use crate::node::Node;

/// 闭谓词集（D5：P0 固定 9 个）。
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Predicate {
    Real,
    Rational,
    Integer,
    Even,
    Odd,
    Positive,
    Negative,
    NonZero,
    Finite,
}

/// 三值查询结果。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Trinary {
    True,
    False,
    Unknown,
}

impl From<bool> for Trinary {
    fn from(b: bool) -> Self {
        if b { Trinary::True } else { Trinary::False }
    }
}

/// 假设位集（u16；位序与 `Predicate::bit` 一致）。恒存闭包后的形态。
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Assumptions(u16);

impl Predicate {
    pub fn bit(self) -> u16 {
        match self {
            Predicate::Real => 1 << 0,
            Predicate::Rational => 1 << 1,
            Predicate::Integer => 1 << 2,
            Predicate::Even => 1 << 3,
            Predicate::Odd => 1 << 4,
            Predicate::Positive => 1 << 5,
            Predicate::Negative => 1 << 6,
            Predicate::NonZero => 1 << 7,
            Predicate::Finite => 1 << 8,
        }
    }
}

/// "已绑定"标记位（区分未绑定与绑定空集：前者查询返回 Unknown）。
const BOUND: u16 = 1 << 15;

impl Assumptions {
    pub fn none() -> Self {
        Assumptions(0)
    }

    /// 单谓词构造（经闭包补全）。
    pub fn with(p: Predicate) -> Self {
        Assumptions(p.bit()).close()
    }

    /// 多谓词并集构造（经闭包补全 + 绑定标记）。
    pub fn union(preds: &[Predicate]) -> Self {
        let mut bits = 0u16;
        for &p in preds {
            bits |= p.bit();
        }
        let mut a = Assumptions(bits).close();
        a.0 |= BOUND;
        a
    }

    /// 是否有绑定（无绑定 ⇒ 一切查询 Unknown）。
    pub fn bound(&self) -> bool {
        self.0 & BOUND != 0
    }

    pub fn has(&self, p: Predicate) -> bool {
        self.0 & p.bit() != 0
    }

    /// 闭包补全 + 冲突检测。返回补全后的位集；冲突（positive∧negative）
    /// 属编程错误，直接 panic（与符号名校验同级的契约）。
    pub fn close(mut self) -> Self {
        // even/odd ⇒ integer ⇒ rational ⇒ real
        if self.has(Predicate::Even) || self.has(Predicate::Odd) {
            self.0 |= Predicate::Integer.bit();
        }
        if self.has(Predicate::Integer) {
            self.0 |= Predicate::Rational.bit();
        }
        if self.has(Predicate::Rational) {
            self.0 |= Predicate::Real.bit();
        }
        // positive/negative ⇒ real ∧ nonzero ∧ finite
        for p in [Predicate::Positive, Predicate::Negative] {
            if self.has(p) {
                self.0 |=
                    Predicate::Real.bit() | Predicate::NonZero.bit() | Predicate::Finite.bit();
            }
        }
        // nonzero ∧ real：与零假设无冲突（零不是 nonzero）
        // 冲突检测
        assert!(
            !(self.has(Predicate::Positive) && self.has(Predicate::Negative)),
            "假设冲突：positive ∧ negative"
        );
        // even ∧ odd 冲突
        assert!(
            !(self.has(Predicate::Even) && self.has(Predicate::Odd)),
            "假设冲突：even ∧ odd"
        );
        self
    }

    /// 非负（positive ∨ zero）——sqrt(x^2) 规则的判据。
    pub fn nonnegative(&self) -> bool {
        self.has(Predicate::Positive)
    }
}

impl Inner {
    /// 符号当前假设（无绑定返回空）。
    pub(crate) fn assumptions_of(&self, sym_id: u32) -> Assumptions {
        *self
            .sym_assumptions
            .get(&sym_id)
            .unwrap_or(&Assumptions::none())
    }

    /// 三值查询：自底向上传播，遇 Unknown 截断。
    pub(crate) fn query_at(&self, id: u32, p: Predicate, depth: u32) -> Trinary {
        assert!(depth <= 10_000, "表达式嵌套过深");
        match &self.nodes[id as usize] {
            Node::Int(v) => {
                // 字面量谓词精确判定（符号/奇偶/整数性）
                match p {
                    Predicate::Integer
                    | Predicate::Rational
                    | Predicate::Real
                    | Predicate::Finite => Trinary::True,
                    Predicate::Positive => (!v.is_zero() && v.sign() > 0).into(),
                    Predicate::Negative => (v.sign() < 0).into(),
                    Predicate::NonZero => (!v.is_zero()).into(),
                    Predicate::Even | Predicate::Odd => match v.to_i64() {
                        Some(k) => {
                            let want_even = matches!(p, Predicate::Even);
                            (k.rem_euclid(2) == i64::from(!want_even)).into()
                        }
                        None => Trinary::Unknown, // 大数奇偶可判但 P1 保守
                    },
                }
            }
            Node::Rat(_) => {
                // 有理数：rational/real/finite 恒真；符号看分子
                match p {
                    Predicate::Rational | Predicate::Real | Predicate::Finite => Trinary::True,
                    _ => Trinary::Unknown, // 精确符号判定 P1 不做（避免大数分解）
                }
            }
            Node::Float { .. } => match p {
                Predicate::Real | Predicate::Finite => Trinary::True,
                _ => Trinary::Unknown,
            },
            Node::Sym(s) => {
                let a = self.assumptions_of(*s);
                if !a.bound() {
                    return Trinary::Unknown;
                }
                if a.has(p) {
                    return Trinary::True;
                }
                // 闭包蕴含的否定：positive⇒¬negative、negative⇒¬positive、
                // even⇒¬odd、odd⇒¬even
                match p {
                    Predicate::Negative if a.has(Predicate::Positive) => Trinary::False,
                    Predicate::Positive if a.has(Predicate::Negative) => Trinary::False,
                    Predicate::Odd if a.has(Predicate::Even) => Trinary::False,
                    Predicate::Even if a.has(Predicate::Odd) => Trinary::False,
                    _ => Trinary::Unknown,
                }
            }
            Node::Add { args: sp } => {
                // 正性：全部 positive ⇒ positive；实性：全部 real ⇒ real
                self.query_seq(self.node_args(*sp), p, depth, QueryOp::All)
            }
            Node::Mul { args: sp } => {
                match p {
                    // 实性：全部 real ⇒ real
                    Predicate::Real
                    | Predicate::Rational
                    | Predicate::Integer
                    | Predicate::Finite => {
                        self.query_seq(self.node_args(*sp), p, depth, QueryOp::All)
                    }
                    Predicate::Positive | Predicate::Negative => {
                        // 符号 = 各因子符号的乘积：奇数个 negative 翻转
                        let mut neg = false;
                        for &a in self.node_args(*sp) {
                            match self.query_at(a, Predicate::Positive, depth + 1) {
                                Trinary::True => {}
                                Trinary::False => neg = !neg,
                                Trinary::Unknown => return Trinary::Unknown,
                            }
                        }
                        if p == Predicate::Positive {
                            (!neg).into()
                        } else {
                            neg.into()
                        }
                    }
                    _ => Trinary::Unknown,
                }
            }
            Node::Pow { base, exp } => {
                // x^2（x real）非负；x^k（k 偶、x real）非负——sqrt(x^2)→x 判据
                if p == Predicate::Positive {
                    // exp(x) 型此处不进；x^k：k 正且 base positive ⇒ positive
                    match self.query_at(*base, Predicate::Positive, depth + 1) {
                        Trinary::True => {
                            match self.query_at(*exp, Predicate::Positive, depth + 1) {
                                Trinary::True => Trinary::True,
                                Trinary::Unknown => Trinary::Unknown,
                                Trinary::False => Trinary::Unknown, // x^0 = 1 仍 positive——保守 Unknown
                            }
                        }
                        _ => Trinary::Unknown,
                    }
                } else if p == Predicate::Real {
                    let b = self.query_at(*base, Predicate::Real, depth + 1);
                    let e = self.query_at(*exp, Predicate::Real, depth + 1);
                    match (b, e) {
                        (Trinary::True, Trinary::True) => {
                            // base real 且 exp 整数 ⇒ real；exp 实数 ⇒ 需 base>0
                            match self.query_at(*exp, Predicate::Integer, depth + 1) {
                                Trinary::True => Trinary::True,
                                _ => Trinary::Unknown,
                            }
                        }
                        _ => Trinary::Unknown,
                    }
                } else {
                    Trinary::Unknown
                }
            }
            Node::Fn { head, args: sp } => {
                let name = self.fn_names[*head as usize].as_ref();
                let args = self.node_args(*sp);
                match (name, p) {
                    // exp(x) > 0 恒成立；exp(x) real 当 x real
                    ("exp", Predicate::Positive) => Trinary::True,
                    ("exp", Predicate::Real | Predicate::Finite | Predicate::NonZero) => {
                        Trinary::True
                    }
                    ("exp", _) => Trinary::Unknown,
                    // sqrt(x)：x≥0 ⇒ real；x>0 ⇒ positive
                    ("sqrt", Predicate::Real) => {
                        if args.len() == 1 {
                            match self.query_at(args[0], Predicate::Positive, depth + 1) {
                                Trinary::True => Trinary::True,
                                _ => Trinary::Unknown,
                            }
                        } else {
                            Trinary::Unknown
                        }
                    }
                    ("sqrt", Predicate::Positive) => {
                        if args.len() == 1 {
                            self.query_at(args[0], Predicate::Positive, depth + 1)
                        } else {
                            Trinary::Unknown
                        }
                    }
                    // sin/cos/tan：real ⇒ real
                    ("sin" | "cos" | "tan", Predicate::Real) => {
                        if args.len() == 1 {
                            self.query_at(args[0], Predicate::Real, depth + 1)
                        } else {
                            Trinary::Unknown
                        }
                    }
                    // log(x)：x>0 ⇒ real
                    ("log", Predicate::Real) => {
                        if args.len() == 1 {
                            match self.query_at(args[0], Predicate::Positive, depth + 1) {
                                Trinary::True => Trinary::True,
                                _ => Trinary::Unknown,
                            }
                        } else {
                            Trinary::Unknown
                        }
                    }
                    _ => Trinary::Unknown,
                }
            }
        }
    }

    fn query_seq(&self, ids: &[u32], p: Predicate, depth: u32, _op: QueryOp) -> Trinary {
        let mut all = true;
        for &a in ids {
            match self.query_at(a, p, depth + 1) {
                Trinary::True => {}
                Trinary::False => return Trinary::False,
                Trinary::Unknown => all = false,
            }
        }
        if all { Trinary::True } else { Trinary::Unknown }
    }
}

#[derive(Clone, Copy)]
enum QueryOp {
    All,
}
