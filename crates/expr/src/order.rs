//! 确定性全序（D6）与跨上下文结构相等。
//!
//! 表达式全序（`ORDER_VERSION = 1`）：
//! 1. 精确数（Int/Rat 统一）按有理数值比较，先于其余一切；
//! 2. **符号类桶**：Sym 与"底为 Sym 的 Pow"同桶（使 Mul 因子输出形如
//!    `x^2*y`、`x*y^2`，按底名排）——桶内按（底名字节序, Sym 先于幂, 指数）；
//! 3. 之后依次：Float < 符号类 < Fn < 复合底 Pow < Mul < Add；
//! 4. Fn 按（head 名, 参数字典序）；Pow 按（底, 指数）；Mul/Add 按参数
//!    字典序（参数已升序），平局比长度。
//!
//! 符号按名字节序而非创建序号——跨进程、跨构造路径确定性的关键。
//! 桶化保证传递性：符号底的幂在符号桶内线性化，不与类别秩交叉。

use crate::Inner;
use crate::node::Node;
use std::cmp::Ordering;

const MAX_DEPTH: u32 = 10_000;

/// 比较桶。精确数在进入前已单独处理。需要访问 arena：符号底的幂归符号桶。
fn family_of(inner: &Inner, n: &Node) -> u8 {
    match n {
        Node::Int(_) | Node::Rat(_) => 0,
        Node::Float { .. } => 1,
        Node::Sym(_) => 2,
        Node::Pow { base, .. } => match &inner.nodes[*base as usize] {
            Node::Sym(_) => 2,
            _ => 4,
        },
        Node::Fn { .. } => 3,
        Node::Mul { .. } => 5,
        Node::Add { .. } => 6,
    }
}

fn is_exact(n: &Node) -> bool {
    matches!(n, Node::Int(_) | Node::Rat(_))
}

fn cmp_exact(a: &Node, b: &Node) -> Ordering {
    match (a, b) {
        (Node::Int(x), Node::Int(y)) => x.cmp(y),
        (Node::Rat(x), Node::Rat(y)) => x.cmp(y),
        (Node::Int(x), Node::Rat(y)) => y.cmp_integer(x).reverse(),
        (Node::Rat(x), Node::Int(y)) => x.cmp_integer(y),
        _ => unreachable!("精确数配对"),
    }
}

impl Inner {
    pub(crate) fn cmp_ids(&self, a: u32, b: u32) -> Ordering {
        self.cmp_at(a, b, 0)
    }

    fn cmp_at(&self, a: u32, b: u32, depth: u32) -> Ordering {
        assert!(depth <= MAX_DEPTH, "表达式嵌套超过 {MAX_DEPTH} 层");
        let na = &self.nodes[a as usize];
        let nb = &self.nodes[b as usize];
        if is_exact(na) && is_exact(nb) {
            return cmp_exact(na, nb);
        }
        if is_exact(na) {
            return Ordering::Less;
        }
        if is_exact(nb) {
            return Ordering::Greater;
        }
        let (fa, fb) = (family_of(self, na), family_of(self, nb));
        if fa != fb {
            return fa.cmp(&fb);
        }
        match (na, nb) {
            (Node::Float { bits: x, .. }, Node::Float { bits: y, .. }) => {
                f64::from_bits(*x).total_cmp(&f64::from_bits(*y))
            }
            (Node::Sym(x), Node::Sym(y)) => {
                self.sym_names[*x as usize].cmp(&self.sym_names[*y as usize])
            }
            // 符号桶内的 Sym 与 符号底幂：按底名，同底 Sym 在前
            (Node::Sym(x), Node::Pow { base, .. }) => match &self.nodes[*base as usize] {
                Node::Sym(t) => self.sym_names[*x as usize]
                    .cmp(&self.sym_names[*t as usize])
                    .then(Ordering::Less),
                _ => unreachable!("复合底幂不在符号桶"),
            },
            (Node::Pow { base, .. }, Node::Sym(y)) => match &self.nodes[*base as usize] {
                Node::Sym(t) => self.sym_names[*t as usize]
                    .cmp(&self.sym_names[*y as usize])
                    .then(Ordering::Greater),
                _ => unreachable!("复合底幂不在符号桶"),
            },
            (Node::Pow { base: b1, exp: e1 }, Node::Pow { base: b2, exp: e2 }) => {
                // 同为符号底（桶内）或同为复合底（Pow 桶）：按（底, 指数）
                self.cmp_at(*b1, *b2, depth + 1)
                    .then_with(|| self.cmp_at(*e1, *e2, depth + 1))
            }
            (Node::Fn { head: h, args: sa }, Node::Fn { head: k, args: sb }) => self.fn_names
                [*h as usize]
                .cmp(&self.fn_names[*k as usize])
                .then_with(|| self.cmp_seq(self.node_args(*sa), self.node_args(*sb), depth)),
            (Node::Mul { args: sa }, Node::Mul { args: sb })
            | (Node::Add { args: sa }, Node::Add { args: sb }) => {
                self.cmp_seq(self.node_args(*sa), self.node_args(*sb), depth)
            }
            _ => unreachable!("同桶必有同类配对"),
        }
    }

    fn cmp_seq(&self, a: &[u32], b: &[u32], depth: u32) -> Ordering {
        for (&x, &y) in a.iter().zip(b.iter()) {
            let o = self.cmp_at(x, y, depth + 1);
            if o != Ordering::Equal {
                return o;
            }
        }
        a.len().cmp(&b.len())
    }
}

/// 跨上下文结构相等（规范形保证：结构相等 ⇔ 序列化相等）。
pub(crate) fn deep_eq(x: &Inner, a: u32, y: &Inner, b: u32) -> bool {
    deep_eq_at(x, a, y, b, 0)
}

fn deep_eq_at(x: &Inner, a: u32, y: &Inner, b: u32, depth: u32) -> bool {
    assert!(depth <= MAX_DEPTH, "表达式嵌套超过 {MAX_DEPTH} 层");
    match (&x.nodes[a as usize], &y.nodes[b as usize]) {
        (Node::Int(u), Node::Int(v)) => u == v,
        (Node::Rat(u), Node::Rat(v)) => u == v,
        (Node::Float { bits: u, .. }, Node::Float { bits: v, .. }) => u == v,
        (Node::Sym(u), Node::Sym(v)) => x.sym_names[*u as usize] == y.sym_names[*v as usize],
        (Node::Fn { head: h, args: sa }, Node::Fn { head: k, args: sb }) => {
            x.fn_names[*h as usize] == y.fn_names[*k as usize]
                && seq_eq(x, y, x.node_args(*sa), y.node_args(*sb), depth)
        }
        (Node::Pow { base: b1, exp: e1 }, Node::Pow { base: b2, exp: e2 }) => {
            deep_eq_at(x, *b1, y, *b2, depth + 1) && deep_eq_at(x, *e1, y, *e2, depth + 1)
        }
        (Node::Mul { args: sa }, Node::Mul { args: sb })
        | (Node::Add { args: sa }, Node::Add { args: sb }) => {
            seq_eq(x, y, x.node_args(*sa), y.node_args(*sb), depth)
        }
        _ => false,
    }
}

fn seq_eq(x: &Inner, y: &Inner, a: &[u32], b: &[u32], depth: u32) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b.iter())
            .all(|(&p, &q)| deep_eq_at(x, p, y, q, depth + 1))
}
