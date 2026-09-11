//! 环与域的代数 trait 层（D2）。
//!
//! P0-M2 引入 Ring/Field，poly 内核泛型于此。EuclideanDomain（gcd，
//! M3 的 PRS/因式分解用）随后引入。实现类型自带规范形（相等 ⇔ 值相等），
//! trait 运算保持规范形。

use crate::{Integer, Rational};

/// 交换环（含单位元）。
pub trait Ring: Sized + Clone + PartialEq {
    fn zero() -> Self;
    fn one() -> Self;
    /// 由 i64 构造系数（poly 偏导的指数乘子等用途）。
    fn from_i64_coeff(v: i64) -> Self;
    fn add(&self, other: &Self) -> Self;
    fn sub(&self, other: &Self) -> Self;
    fn mul(&self, other: &Self) -> Self;
    fn neg(&self) -> Self;
    fn is_zero(&self) -> bool;
    fn is_one(&self) -> bool;

    /// 非负整数幂；默认重复乘，实现方可覆写快速路径。
    fn pow_u32(&self, e: u32) -> Self {
        let mut acc = Self::one();
        for _ in 0..e {
            acc = acc.mul(self);
        }
        acc
    }
}

/// 域：环 + 乘法逆（零无逆）。
pub trait Field: Ring {
    fn inv(&self) -> Option<Self>;
}

impl Ring for Integer {
    fn zero() -> Self {
        Integer::zero()
    }
    fn one() -> Self {
        Integer::one()
    }
    fn from_i64_coeff(v: i64) -> Self {
        Integer::from_i64(v)
    }
    fn add(&self, other: &Self) -> Self {
        Integer::add(self, other)
    }
    fn sub(&self, other: &Self) -> Self {
        Integer::sub(self, other)
    }
    fn mul(&self, other: &Self) -> Self {
        Integer::mul(self, other)
    }
    fn neg(&self) -> Self {
        Integer::neg(self)
    }
    fn is_zero(&self) -> bool {
        Integer::is_zero(self)
    }
    fn is_one(&self) -> bool {
        Integer::is_one(self)
    }
    fn pow_u32(&self, e: u32) -> Self {
        Integer::pow(self, e)
    }
}

impl Ring for Rational {
    fn zero() -> Self {
        Rational::zero()
    }
    fn one() -> Self {
        Rational::one()
    }
    fn from_i64_coeff(v: i64) -> Self {
        Rational::from_integer(&Integer::from_i64(v))
    }
    fn add(&self, other: &Self) -> Self {
        Rational::add(self, other)
    }
    fn sub(&self, other: &Self) -> Self {
        Rational::sub(self, other)
    }
    fn mul(&self, other: &Self) -> Self {
        Rational::mul(self, other)
    }
    fn neg(&self) -> Self {
        Rational::neg(self)
    }
    fn is_zero(&self) -> bool {
        Rational::is_zero(self)
    }
    fn is_one(&self) -> bool {
        Rational::is_one(self)
    }
    fn pow_u32(&self, e: u32) -> Self {
        // Rational 规范形保证既约，幂仍既约（跳过 gcd 的快速路径）
        Rational::pow_reduced(self, e)
    }
}

impl Field for Rational {
    fn inv(&self) -> Option<Self> {
        Rational::inv_reduced(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn 整数环() {
        let (a, b) = (Integer::from_i64(7), Integer::from_i64(-6));
        assert!(a.add(&b).is_one());
        let c = Integer::from_i64(7);
        assert!(!c.mul(&Integer::from_i64(1)).is_one() || c.is_one());
        assert_eq!(a.sub(&a), <Integer as Ring>::zero());
        assert_eq!(Integer::from_i64(-3).pow_u32(3), Integer::from_i64(-27));
        assert_eq!(<Integer as Ring>::from_i64_coeff(5), Integer::from_i64(5));
    }

    #[test]
    fn 有理数域() {
        let r = Rational::from_ints(&Integer::from_i64(3), &Integer::from_i64(4)).unwrap();
        let inv = <Rational as Field>::inv(&r).unwrap();
        assert!(inv.mul(&r).is_one());
        assert!(Rational::zero().inv().is_none());
        assert!(r.pow_u32(2).to_string() == "9/16");
    }
}
