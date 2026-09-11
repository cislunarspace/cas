//! cas 数值层。
//!
//! P0-M1 范围：精确数值类型 [`Integer`] 与 [`Rational`]，i64 内联 + 大数后备的
//! 混合表示，规范形由构造保证（见各类型文档）。Ring/Field trait 层随 M2
//! 多项式内核一起引入（见 `cas/DESIGN.md` D2）。
//!
//! 数字后端：M1 采用保守备选 num-bigint（MIT/Apache）。默认后端的最终
//! 决策（numerica vs num-bigint vs rug 性能档）在 M2 基准门控点完成，
//! 本模块的公共接口不随后端切换变化。

pub mod integer;
pub mod rational;
pub mod ring;

pub use integer::Integer;
pub use rational::Rational;
pub use ring::{Field, Ring};
