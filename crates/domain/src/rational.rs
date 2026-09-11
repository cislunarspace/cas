//! 任意精度有理数。
//!
//! 规范形：分母为正、分子分母既约、零表示为 `0/1`。
//! `Small{ n: i64, d: u64 }` 是常用路径；分子越界（分母恒可放 u64）时落 `Big`。
//! 既约性使相等与哈希和表示无关；`d = 1` 合法（expr 层在节点上把它归一为 Int）。

use super::integer::{Integer, big_gcd, big_sign, gcd_u64};
use num_bigint::BigInt;
use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};

#[derive(Clone, Debug, PartialEq, Eq)]
enum Repr {
    Small { n: i64, d: u64 },
    Big(Box<(BigInt, BigInt)>),
}

/// 任意精度有理数，规范形见模块文档。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Rational(Repr);

impl Rational {
    pub fn zero() -> Self {
        Rational(Repr::Small { n: 0, d: 1 })
    }

    pub fn one() -> Self {
        Rational(Repr::Small { n: 1, d: 1 })
    }

    pub fn from_integer(v: &Integer) -> Self {
        match v.to_i64() {
            Some(n) => Rational(Repr::Small { n, d: 1 }),
            None => Rational(Repr::Big(Box::new((v.to_bigint(), BigInt::from(1))))),
        }
    }

    /// 由整数分子分母构造并规范化；`d = 0` 返回 `None`。
    pub fn from_ints(n: &Integer, d: &Integer) -> Option<Self> {
        if d.is_zero() {
            return None;
        }
        // 符号挪到分子，分母取正
        let (n, d) = if d.is_negative() {
            (n.neg(), d.abs())
        } else {
            (n.clone(), d.clone())
        };
        if n.is_zero() {
            return Some(Rational::zero());
        }
        let normalized = match (n.to_i64(), d.to_i64()) {
            (Some(sn), Some(sd)) => {
                let sd = sd as u64;
                let g = gcd_u64(sn.unsigned_abs(), sd);
                let nn = (sn as i128 / g as i128) as i64;
                Rational(Repr::Small { n: nn, d: sd / g })
            }
            _ => Self::from_big_parts(n.to_bigint(), d.to_bigint()),
        };
        Some(normalized)
    }

    pub fn is_zero(&self) -> bool {
        matches!(self.0, Repr::Small { n: 0, .. })
    }

    pub fn is_one(&self) -> bool {
        matches!(self.0, Repr::Small { n: 1, d: 1 })
    }

    pub fn is_negative(&self) -> bool {
        self.sign() < 0
    }

    pub fn sign(&self) -> i32 {
        match &self.0 {
            Repr::Small { n, .. } => n.signum() as i32,
            Repr::Big(p) => match big_sign(&p.0) {
                Ordering::Less => -1,
                Ordering::Equal => 0,
                Ordering::Greater => 1,
            },
        }
    }

    pub fn num(&self) -> Integer {
        match &self.0 {
            Repr::Small { n, .. } => Integer::from_i64(*n),
            Repr::Big(p) => Integer::from_bigint(p.0.clone()),
        }
    }

    pub fn den(&self) -> Integer {
        match &self.0 {
            // Small 表示的分母是 u64，可超过 i64::MAX——必须经 from_u64 落位
            Repr::Small { d, .. } => Integer::from_u64(*d),
            Repr::Big(p) => Integer::from_bigint(p.1.clone()),
        }
    }

    pub fn neg(&self) -> Self {
        match &self.0 {
            Repr::Small { n, d } => Rational(Repr::Small { n: -n, d: *d }),
            Repr::Big(p) => Self::from_big_parts(-p.0.clone(), p.1.clone()),
        }
    }

    pub fn abs(&self) -> Self {
        if self.is_negative() {
            self.neg()
        } else {
            self.clone()
        }
    }

    /// 乘法逆元；零无逆元。
    pub fn inv(&self) -> Option<Self> {
        Self::from_ints(&self.den(), &self.num())
    }

    pub fn add(&self, other: &Self) -> Self {
        match (&self.0, &other.0) {
            (Repr::Small { n: n1, d: d1 }, Repr::Small { n: n2, d: d2 }) => {
                let n = (*n1 as i128) * (*d2 as i128) + (*n2 as i128) * (*d1 as i128);
                let d = (*d1 as u128) * (*d2 as u128);
                Self::from_parts_wide(n, d)
            }
            _ => Self::from_big_parts(
                self.n_big() * other.d_big() + other.n_big() * self.d_big(),
                self.d_big() * other.d_big(),
            ),
        }
    }

    pub fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    pub fn mul(&self, other: &Self) -> Self {
        match (&self.0, &other.0) {
            (Repr::Small { n: n1, d: d1 }, Repr::Small { n: n2, d: d2 }) => {
                let n = (*n1 as i128) * (*n2 as i128);
                let d = (*d1 as u128) * (*d2 as u128);
                Self::from_parts_wide(n, d)
            }
            _ => Self::from_big_parts(self.n_big() * other.n_big(), self.d_big() * other.d_big()),
        }
    }

    /// 有理数相除；除零返回 `None`。
    pub fn div(&self, other: &Self) -> Option<Self> {
        Some(self.mul(&other.inv()?))
    }

    /// 非负整数幂。规模守卫由调用方（expr 层）负责。
    pub fn pow(&self, exp: u32) -> Self {
        Self::from_ints(&self.num().pow(exp), &self.den().pow(exp)).expect("分母非零")
    }

    /// 非负整数幂（快速路径）：既约有理数的幂仍既约——`gcd(num,den)^m = 1`，
    /// 无需再跑 gcd。大幂折叠产生几万比特常数时 gcd 是平方级瓶颈
    /// （perf_probe case 1431 实测），此路径跳过它。
    ///
    /// 落位判定必须与 from_big_parts 一致：分母按 u64 判（不是 i64），
    /// 否则同值出现两种表示，破坏"规范形与表示无关"（29^13 案例实测）。
    pub fn pow_reduced(&self, exp: u32) -> Self {
        let n = self.num().pow(exp);
        let d = self.den().pow(exp);
        match (n.to_i64(), d.to_u64()) {
            (Some(sn), Some(sd)) => Rational(Repr::Small { n: sn, d: sd }),
            _ => Rational(Repr::Big(Box::new((n.to_bigint(), d.to_bigint())))),
        }
    }

    /// f64 近似（数值交叉验证路径；大数经十进制字符串中转，精度为 f64 语义）。
    pub fn to_f64(&self) -> f64 {
        match &self.0 {
            Repr::Small { n, d } => (*n as f64) / (*d as f64),
            Repr::Big(p) => {
                let n: f64 = p.0.to_string().parse().unwrap_or(f64::INFINITY);
                let d: f64 = p.1.to_string().parse().unwrap_or(f64::INFINITY);
                n / d
            }
        }
    }

    /// 乘法逆元（快速路径，前提：自身既约——Rational 规范形保证）；零无逆元。
    /// 逆元的既约性继承自原数，跳过 gcd。
    pub fn inv_reduced(&self) -> Option<Self> {
        match &self.0 {
            Repr::Small { n: 0, .. } => None,
            Repr::Small { n, d } => {
                // 逆元 = d/n（符号归分子）；Small 的 d 是 u64，可超 i64::MAX
                if *d <= i64::MAX as u64 {
                    let num = *d as i64;
                    Some(Rational(Repr::Small {
                        n: if *n < 0 { -num } else { num },
                        d: n.unsigned_abs(),
                    }))
                } else {
                    let num = Integer::from_u64(*d);
                    let num = if *n < 0 { num.neg() } else { num };
                    Some(Rational(Repr::Big(Box::new((
                        num.to_bigint(),
                        BigInt::from(n.unsigned_abs()),
                    )))))
                }
            }
            Repr::Big(p) => {
                let (num, den) = if p.0.sign() == num_bigint::Sign::Minus {
                    (-p.1.clone(), -p.0.clone())
                } else {
                    (p.1.clone(), p.0.clone())
                };
                Some(Rational(Repr::Big(Box::new((num, den)))))
            }
        }
    }

    /// 与整数比较：`self - o` 的符号。
    pub fn cmp_integer(&self, o: &Integer) -> Ordering {
        match (&self.0, o.to_i64()) {
            (Repr::Small { n, d }, Some(i)) => (*n as i128).cmp(&((i as i128) * (*d as i128))),
            _ => self.n_big().cmp(&(o.to_bigint() * self.d_big())),
        }
    }

    fn from_parts_wide(n: i128, d: u128) -> Self {
        debug_assert!(d > 0);
        if n == 0 {
            return Rational::zero();
        }
        let g = gcd_u128(n.unsigned_abs(), d);
        let (n, d) = (n / g as i128, d / g);
        match (i64::try_from(n), u64::try_from(d)) {
            (Ok(sn), Ok(sd)) => Rational(Repr::Small { n: sn, d: sd }),
            _ => Self::from_big_parts(BigInt::from(n), BigInt::from(d)),
        }
    }

    fn from_big_parts(n: BigInt, d: BigInt) -> Self {
        debug_assert!(d.sign() != num_bigint::Sign::Minus);
        let g = big_gcd(&n, &d);
        let (n, d) = if g == BigInt::from(1) {
            (n, d)
        } else {
            (n / &g, d / &g)
        };
        // 归约后能放回 Small 的必须放回：规范形与表示无关的前提
        match (i64::try_from(&n), u64::try_from(&d)) {
            (Ok(sn), Ok(sd)) => Rational(Repr::Small { n: sn, d: sd }),
            _ => Rational(Repr::Big(Box::new((n, d)))),
        }
    }

    fn n_big(&self) -> BigInt {
        match &self.0 {
            Repr::Small { n, .. } => BigInt::from(*n),
            Repr::Big(p) => p.0.clone(),
        }
    }

    fn d_big(&self) -> BigInt {
        match &self.0 {
            Repr::Small { d, .. } => BigInt::from(*d),
            Repr::Big(p) => p.1.clone(),
        }
    }
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

impl Ord for Rational {
    fn cmp(&self, other: &Self) -> Ordering {
        match (&self.0, &other.0) {
            (Repr::Small { n: n1, d: d1 }, Repr::Small { n: n2, d: d2 }) => {
                let l = (*n1 as i128) * (*d2 as i128);
                let r = (*n2 as i128) * (*d1 as i128);
                l.cmp(&r)
            }
            _ => (self.n_big() * other.d_big()).cmp(&(other.n_big() * self.d_big())),
        }
    }
}

impl PartialOrd for Rational {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Hash for Rational {
    fn hash<H: Hasher>(&self, state: &mut H) {
        match &self.0 {
            Repr::Small { n, d } => {
                n.hash(state);
                d.hash(state);
            }
            Repr::Big(p) => {
                p.0.hash(state);
                p.1.hash(state);
            }
        }
    }
}

impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.0 {
            Repr::Small { n, d: 1 } => write!(f, "{n}"),
            Repr::Small { n, d } => write!(f, "{n}/{d}"),
            Repr::Big(p) if p.1 == BigInt::from(1) => write!(f, "{}", p.0),
            Repr::Big(p) => write!(f, "{}/{}", p.0, p.1),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(n: i64, d: i64) -> Rational {
        Rational::from_ints(&Integer::from_i64(n), &Integer::from_i64(d)).unwrap()
    }

    #[test]
    fn 规范化() {
        assert_eq!(r(4, 8), r(1, 2));
        assert_eq!(r(-6, 4), r(-3, 2));
        assert_eq!(r(6, -4), r(-3, 2));
        assert_eq!(r(0, 7), Rational::zero());
        assert_eq!(r(9, 3).to_string(), "3");
        assert_eq!(r(7, 3).to_string(), "7/3");
        assert_eq!(
            Rational::from_ints(&Integer::from_i64(1), &Integer::zero()),
            None
        );
    }

    #[test]
    fn 四则() {
        assert_eq!(r(1, 2).add(&r(1, 3)), r(5, 6));
        assert_eq!(r(1, 2).sub(&r(1, 3)), r(1, 6));
        assert_eq!(r(2, 3).mul(&r(3, 4)), r(1, 2));
        assert_eq!(r(1, 2).div(&r(3, 4)), Some(r(2, 3)));
        assert_eq!(r(1, 2).div(&Rational::zero()), None);
        // 1/3 * 3 归一到整数分母
        assert_eq!(r(1, 3).mul(&r(3, 1)), Rational::one());
    }

    #[test]
    fn 大数路径() {
        let big = Integer::parse("9223372036854775808").unwrap(); // 2^63
        let q = Rational::from_ints(&big, &Integer::from_i64(2)).unwrap();
        assert_eq!(q.to_string(), "4611686018427387904");
        // 分子越界时 Big 表示与 Small 表示的序一致
        assert!(q.cmp(&r(1, 1)) == Ordering::Greater);
        // 2^62 × 1/2 = 2^61，跨表示乘法保持精确
        assert_eq!(q.mul(&r(1, 2)).to_string(), "2305843009213693952");
        // i64::MIN 作分子
        let m = Integer::from_i64(i64::MIN);
        let qm = Rational::from_ints(&m, &Integer::from_i64(2)).unwrap();
        assert_eq!(qm.to_string(), "-4611686018427387904");
    }

    #[test]
    fn 幂与逆元() {
        assert_eq!(r(2, 3).pow(3), r(8, 27));
        assert_eq!(r(3, 2).inv(), Some(r(2, 3)));
        assert_eq!(r(-3, 2).inv(), Some(r(-2, 3)));
        assert_eq!(Rational::zero().inv(), None);
    }

    #[test]
    fn 与整数比较() {
        assert_eq!(
            r(5, 2).cmp_integer(&Integer::from_i64(2)),
            Ordering::Greater
        );
        assert_eq!(r(5, 2).cmp_integer(&Integer::from_i64(3)), Ordering::Less);
        assert_eq!(r(4, 2).cmp_integer(&Integer::from_i64(2)), Ordering::Equal);
        let big = Integer::parse("99999999999999999999").unwrap();
        assert_eq!(r(1, 1).cmp_integer(&big), Ordering::Less);
    }

    #[test]
    fn 快速路径落位一致性() {
        // 29^13 超过 i64::MAX 但在 u64 内：pow_reduced/inv_reduced 的落位判定
        // 必须与 from_ints/from_big_parts 一致，否则同值出现两种表示，
        // 破坏"规范形与表示无关"（proptide 回归案例：(-29/2)^-13 往返失败）。
        let d13 = Integer::parse("10260628712958602189").unwrap();
        let r = Rational::from_ints(&Integer::from_i64(-2), &Integer::from_i64(29)).unwrap();
        let p = r.pow_reduced(13);
        let expect = Rational::from_ints(&Integer::from_i64(-8192), &d13).unwrap();
        assert_eq!(p, expect); // PartialEq 比较表示：此处即断言同表示
        assert_eq!(p.to_string(), "-8192/10260628712958602189");

        // 分母超 i64 的逆元：分子落 Big，与通用路径同表示
        let q = Rational::from_ints(&Integer::from_i64(3), &d13).unwrap();
        let inv = q.inv_reduced().unwrap();
        let expect_inv = Rational::from_ints(&d13, &Integer::from_i64(3)).unwrap();
        assert_eq!(inv, expect_inv);
    }
}
