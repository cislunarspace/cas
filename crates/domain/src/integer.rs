//! 任意精度整数：i64 内联 + 大数后备。
//!
//! 规范形：凡能放进 i64 的值一律是 `Small`，`Big` 只存越界值。
//! 这条不变量由 [`Integer::from_big`] 在所有构造路径上保证，
//! 因此相等、序、哈希都与具体表示无关，跨表示比较无需分配。

use num_bigint::BigInt;
use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};

/// `BigInt::sign()`（`Sign` 枚举）到 `Ordering` 的换算。
pub(crate) fn big_sign(b: &BigInt) -> Ordering {
    match b.sign() {
        num_bigint::Sign::Minus => Ordering::Less,
        num_bigint::Sign::NoSign => Ordering::Equal,
        num_bigint::Sign::Plus => Ordering::Greater,
    }
}

/// 大数 gcd（非负）。用 num-integer 的实现（num-bigint 对 BigUint 有
/// 优化算法）；手写欧几里得在几千比特的操作数上是平方级灾难——大幂折叠
/// 产生的有理系数会直接把它引爆（perf_probe 实测）。
pub(crate) fn big_gcd(a: &BigInt, b: &BigInt) -> BigInt {
    use num_integer::Integer as _;
    a.gcd(b)
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum Repr {
    Small(i64),
    Big(Box<BigInt>),
}

/// 任意精度整数，规范形见模块文档。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Integer(Repr);

impl Integer {
    pub const fn zero() -> Self {
        Integer(Repr::Small(0))
    }

    pub const fn one() -> Self {
        Integer(Repr::Small(1))
    }

    pub const fn from_i64(v: i64) -> Self {
        Integer(Repr::Small(v))
    }

    /// 解析十进制字面量（可带负号）；超出 i64 自动落大数表示。
    pub fn parse(s: &str) -> Option<Self> {
        if let Ok(v) = s.parse::<i64>() {
            return Some(Integer(Repr::Small(v)));
        }
        BigInt::parse_bytes(s.as_bytes(), 10).map(Self::from_big)
    }

    pub fn is_zero(&self) -> bool {
        matches!(self.0, Repr::Small(0))
    }

    pub fn is_one(&self) -> bool {
        matches!(self.0, Repr::Small(1))
    }

    pub fn is_negative(&self) -> bool {
        self.sign() < 0
    }

    pub fn sign(&self) -> i32 {
        match &self.0 {
            Repr::Small(v) => v.signum() as i32,
            Repr::Big(b) => match big_sign(b) {
                Ordering::Less => -1,
                Ordering::Equal => 0,
                Ordering::Greater => 1,
            },
        }
    }

    pub fn to_i64(&self) -> Option<i64> {
        match &self.0 {
            Repr::Small(v) => Some(*v),
            Repr::Big(_) => None,
        }
    }

    /// 非负值且不超过 u64 时返回（Rational 的 Small 分母可到 u64 上界；
    /// Big 表示按值判断——(i64::MAX, u64::MAX] 区间的值在 Integer 是 Big，
    /// 但作为 Rational 分母仍应落 Small，故此处不能只看表示）。
    pub fn to_u64(&self) -> Option<u64> {
        match &self.0 {
            Repr::Small(v) if *v >= 0 => Some(*v as u64),
            Repr::Small(_) => None,
            Repr::Big(b) => u64::try_from(b.as_ref()).ok(),
        }
    }

    /// u64 落位：能进 i64 用 Small，否则 Big（规范形要求）。
    pub fn from_u64(v: u64) -> Self {
        match i64::try_from(v) {
            Ok(s) => Integer(Repr::Small(s)),
            Err(_) => Integer(Repr::Big(Box::new(BigInt::from(v)))),
        }
    }

    /// 数值的二进制位数（0 返回 0）。expr 层用它做幂折叠的规模守卫。
    pub fn bit_len(&self) -> u64 {
        match &self.0 {
            Repr::Small(v) => {
                let m = v.unsigned_abs();
                (u64::BITS - m.leading_zeros()) as u64
            }
            Repr::Big(b) => b.bits(),
        }
    }

    pub fn abs(&self) -> Self {
        match &self.0 {
            Repr::Small(v) => match v.checked_abs() {
                Some(a) => Integer(Repr::Small(a)),
                // i64::MIN：|v| = 2^63 越界，落大数
                None => Self::from_big(-self.as_big()),
            },
            Repr::Big(b) => match big_sign(b) {
                Ordering::Less => Self::from_big(-(b.as_ref().clone())),
                _ => self.clone(),
            },
        }
    }

    pub fn neg(&self) -> Self {
        match &self.0 {
            Repr::Small(v) => match v.checked_neg() {
                Some(n) => Integer(Repr::Small(n)),
                None => Self::from_big(-self.as_big()),
            },
            Repr::Big(b) => Self::from_big(-(**b).clone()),
        }
    }

    pub fn add(&self, other: &Self) -> Self {
        if let (Repr::Small(a), Repr::Small(b)) = (&self.0, &other.0) {
            if let Some(s) = a.checked_add(*b) {
                return Integer(Repr::Small(s));
            }
        }
        Self::from_big(self.as_big() + other.as_big())
    }

    pub fn sub(&self, other: &Self) -> Self {
        if let (Repr::Small(a), Repr::Small(b)) = (&self.0, &other.0) {
            if let Some(s) = a.checked_sub(*b) {
                return Integer(Repr::Small(s));
            }
        }
        Self::from_big(self.as_big() - other.as_big())
    }

    pub fn mul(&self, other: &Self) -> Self {
        if let (Repr::Small(a), Repr::Small(b)) = (&self.0, &other.0) {
            if let Some(p) = a.checked_mul(*b) {
                return Integer(Repr::Small(p));
            }
        }
        Self::from_big(self.as_big() * other.as_big())
    }

    /// 非负整数幂。规模守卫由调用方（expr 层）负责，本方法不做上限检查。
    pub fn pow(&self, exp: u32) -> Self {
        if let Repr::Small(v) = &self.0 {
            if let Some(p) = v.checked_pow(exp) {
                return Integer(Repr::Small(p));
            }
        }
        Self::from_big(self.as_big().pow(exp))
    }

    /// 辗转相除，结果非负；`gcd(0, 0) = 0`。
    pub fn gcd(&self, other: &Self) -> Self {
        if let (Repr::Small(a), Repr::Small(b)) = (&self.0, &other.0) {
            Integer(Repr::Small(
                gcd_u64(a.unsigned_abs(), b.unsigned_abs()) as i64
            ))
        } else {
            Self::from_big(big_gcd(&self.as_big(), &other.as_big()))
        }
    }

    /// 精确除法；调用方保证整除（Rational 规范化路径）。
    pub fn div_exact(&self, divisor: &Self) -> Self {
        Self::from_big(self.as_big() / divisor.as_big())
    }

    fn from_big(b: BigInt) -> Self {
        match i64::try_from(&b) {
            Ok(v) => Integer(Repr::Small(v)),
            Err(_) => Integer(Repr::Big(Box::new(b))),
        }
    }

    /// 由大数构造（自动落回 Small 表示，若值在 i64 内）。
    pub fn from_bigint(b: BigInt) -> Self {
        Self::from_big(b)
    }

    fn as_big(&self) -> BigInt {
        match &self.0 {
            Repr::Small(v) => BigInt::from(*v),
            Repr::Big(b) => (**b).clone(),
        }
    }

    /// 导出为大数表示（Rational 跨表示运算的提升通道）。
    pub fn to_bigint(&self) -> BigInt {
        self.as_big()
    }
}

pub(crate) fn gcd_u64(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

impl Ord for Integer {
    fn cmp(&self, other: &Self) -> Ordering {
        match (&self.0, &other.0) {
            (Repr::Small(a), Repr::Small(b)) => a.cmp(b),
            (Repr::Big(a), Repr::Big(b)) => a.cmp(b),
            // Big 一定在 i64 范围之外：正大数比一切 i64 大，负大数比一切 i64 小
            (_, Repr::Big(b)) => big_sign(b).reverse(),
            (Repr::Big(a), Repr::Small(_)) => big_sign(a),
        }
    }
}

impl PartialOrd for Integer {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Hash for Integer {
    fn hash<H: Hasher>(&self, state: &mut H) {
        match &self.0 {
            Repr::Small(v) => state.write_i64(*v),
            Repr::Big(b) => b.hash(state),
        }
    }
}

impl fmt::Display for Integer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.0 {
            Repr::Small(v) => write!(f, "{v}"),
            Repr::Big(b) => write!(f, "{b}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn 小数路径与规范形() {
        assert_eq!(
            Integer::from_i64(2).add(&Integer::from_i64(3)),
            Integer::from_i64(5)
        );
        assert!(!Integer::from_i64(5).is_one());
        // 越界即落 Big，且 from_big 回收可表示值
        let big = Integer::from_i64(i64::MAX).add(&Integer::one());
        assert_eq!(big.to_string(), "9223372036854775808");
        assert_eq!(big.sub(&Integer::one()), Integer::from_i64(i64::MAX));
    }

    #[test]
    fn i64_min_边界() {
        let m = Integer::from_i64(i64::MIN);
        // |MIN| = 2^63 只能落大数
        assert_eq!(m.abs().to_string(), "9223372036854775808");
        assert_eq!(m.neg().to_string(), "9223372036854775808");
        assert_eq!(m.mul(&Integer::from_i64(1)), m);
        // MIN 是负数，且小于一切 Small
        assert!(m.is_negative());
        assert!(m.cmp(&Integer::from_i64(i64::MAX)) == Ordering::Less);
    }

    #[test]
    fn 跨表示序与相等() {
        let big = Integer::parse("9223372036854775808").unwrap();
        assert_eq!(big, Integer::parse("9223372036854775808").unwrap());
        assert!(Integer::from_i64(-1).cmp(&big) == Ordering::Less);
        assert!(big.cmp(&Integer::from_i64(i64::MAX)) == Ordering::Greater);
        let neg_big = big.neg();
        // -(2^63) 的值恰为 i64::MIN：归一化后两者相等
        assert_eq!(neg_big.cmp(&Integer::from_i64(i64::MIN)), Ordering::Equal);
        let below_min = Integer::parse("-9223372036854775809").unwrap();
        assert_eq!(below_min.cmp(&Integer::from_i64(i64::MIN)), Ordering::Less);
    }

    #[test]
    fn 幂与位数() {
        assert_eq!(Integer::from_i64(2).pow(10).to_string(), "1024");
        assert_eq!(
            Integer::from_i64(2).pow(100).to_string(),
            "1267650600228229401496703205376"
        );
        assert_eq!(Integer::from_i64(0).bit_len(), 0);
        assert_eq!(Integer::from_i64(255).bit_len(), 8);
        assert_eq!(
            Integer::parse("1267650600228229401496703205376")
                .unwrap()
                .bit_len(),
            101
        );
    }

    #[test]
    fn 最大公约数() {
        assert_eq!(
            Integer::from_i64(12).gcd(&Integer::from_i64(18)),
            Integer::from_i64(6)
        );
        assert_eq!(
            Integer::from_i64(0).gcd(&Integer::from_i64(5)),
            Integer::from_i64(5)
        );
        let big = Integer::parse("1000000007").unwrap();
        assert_eq!(big.gcd(&Integer::from_i64(3)), Integer::one());
    }
}
