//! 显式变换（L1）：expand 与 subst（DESIGN.md D4）。
//!
//! M2 的 expand 用 DAG 分配实现：每一步都经规范形构造器，同类项即时合并，
//! 中间规模被规范形钳制在真实支撑集大小（(x+y+z+w)^20 全程 ≤ 1771 项）。
//! 多项式内核承接的快速 expand（D3 桥接）属 M3。

use crate::Inner;
use crate::node::Node;
use std::collections::HashMap;

/// 幂展开的指数上限。
const MAX_POW_EXPAND: u64 = 10_000;
/// 展开的项数规模守卫：底式项数 × 指数超过此值则保留未展开。
const MAX_EXPAND_TERMS: u64 = 1_000_000;

/// 节点的 owned 视图（避免边读边写 arena 的借用冲突）。
enum Info {
    Add(Vec<u32>),
    Mul(Vec<u32>),
    Pow { base: u32, exp: u32 },
    Fn { head: u32, args: Vec<u32> },
    Atom(u32),
}

fn info_of(inner: &Inner, id: u32) -> Info {
    match &inner.nodes[id as usize] {
        Node::Add { args: sp } => Info::Add(inner.node_args(*sp).to_vec()),
        Node::Mul { args: sp } => Info::Mul(inner.node_args(*sp).to_vec()),
        Node::Pow { base, exp } => Info::Pow {
            base: *base,
            exp: *exp,
        },
        Node::Fn { head, args: sp } => Info::Fn {
            head: *head,
            args: inner.node_args(*sp).to_vec(),
        },
        _ => Info::Atom(id),
    }
}

impl Inner {
    /// 展开：先尝试多项式快路径（纯多项式子树整体进 poly 域，M3），
    /// 不适用则 DAG 分配。两条路径产出同一规范形节点（同构性有测试钉死）。
    pub(crate) fn expand_at(&mut self, id: u32, depth: u32) -> u32 {
        self.expand_impl(id, depth, true)
    }

    /// `fast = false` 强制纯 DAG 分配（同构性测试的对照路径）。
    pub(crate) fn expand_impl(&mut self, id: u32, depth: u32, fast: bool) -> u32 {
        assert!(depth <= 10_000, "表达式嵌套过深");
        // 快路径：整个子树是纯多项式（含乘积分配与整数幂）→ poly 域一次算完
        if fast {
            if let Info::Mul(_) | Info::Pow { .. } = info_of(self, id) {
                if let Some((ring, var_ids, p)) = crate::poly_bridge::to_poly(self, id) {
                    return crate::poly_bridge::from_poly(self, &ring, &var_ids, &p);
                }
            }
        }
        match info_of(self, id) {
            Info::Atom(_) => id,
            Info::Add(args) => {
                let v: Vec<u32> = args
                    .iter()
                    .map(|&a| self.expand_impl(a, depth + 1, fast))
                    .collect();
                self.make_add(&v)
            }
            Info::Mul(args) => {
                let v: Vec<u32> = args
                    .iter()
                    .map(|&a| self.expand_impl(a, depth + 1, fast))
                    .collect();
                // 逐因子分配；每步 make_add 合并同类项，钳制中间规模
                let mut acc = self.sum_terms(v[0]);
                for &f in &v[1..] {
                    let ft = self.sum_terms(f);
                    let mut next = Vec::with_capacity(acc.len() * ft.len());
                    for &a in &acc {
                        for &b in &ft {
                            next.push(self.make_mul(&[a, b]));
                        }
                    }
                    let s = self.make_add(&next);
                    acc = self.sum_terms(s);
                }
                self.make_add(&acc)
            }
            Info::Pow { base, exp } => {
                let k = match &self.nodes[exp as usize] {
                    Node::Int(v) => v.to_i64(),
                    _ => None,
                };
                if let Some(k) = k {
                    if (2..=MAX_POW_EXPAND as i64).contains(&k) {
                        let b = self.expand_impl(base, depth + 1, fast);
                        let bt = self.sum_terms(b);
                        if bt.len() as u64 * k as u64 <= MAX_EXPAND_TERMS {
                            let one = self.lit_int(1);
                            let mut acc = vec![one];
                            for _ in 0..k {
                                let mut next = Vec::with_capacity(acc.len() * bt.len());
                                for &a in &acc {
                                    for &b in &bt {
                                        next.push(self.make_mul(&[a, b]));
                                    }
                                }
                                let s = self.make_add(&next);
                                acc = self.sum_terms(s);
                            }
                            return self.make_add(&acc);
                        }
                    }
                }
                // 负幂/超大指数/非整数指数：只展开内部
                let b = self.expand_impl(base, depth + 1, fast);
                let e = self.expand_impl(exp, depth + 1, fast);
                self.make_pow(b, e)
            }
            Info::Fn { head, args } => {
                let v: Vec<u32> = args
                    .iter()
                    .map(|&a| self.expand_impl(a, depth + 1, fast))
                    .collect();
                self.fn_node_by_id(head, &v)
            }
        }
    }

    /// `id` 作为和式的项列表（非 Add 即单项）。
    fn sum_terms(&self, id: u32) -> Vec<u32> {
        match &self.nodes[id as usize] {
            Node::Add { args: sp } => self.node_args(*sp).to_vec(),
            _ => vec![id],
        }
    }

    /// 代换：`map` 为符号表 id → 替换节点 id。底向上重建，重建经规范形
    /// 构造器——未命中的子树经 intern 去重回原节点，故 `subst(e, x→x) == e`。
    pub(crate) fn subst_at(&mut self, id: u32, map: &HashMap<u32, u32>, depth: u32) -> u32 {
        assert!(depth <= 10_000, "表达式嵌套过深");
        match info_of(self, id) {
            Info::Atom(aid) => match &self.nodes[aid as usize] {
                Node::Sym(s) => map.get(s).copied().unwrap_or(aid),
                _ => aid,
            },
            Info::Add(args) => {
                let v: Vec<u32> = args
                    .iter()
                    .map(|&a| self.subst_at(a, map, depth + 1))
                    .collect();
                self.make_add(&v)
            }
            Info::Mul(args) => {
                let v: Vec<u32> = args
                    .iter()
                    .map(|&a| self.subst_at(a, map, depth + 1))
                    .collect();
                self.make_mul(&v)
            }
            Info::Pow { base, exp } => {
                let b = self.subst_at(base, map, depth + 1);
                let e = self.subst_at(exp, map, depth + 1);
                self.make_pow(b, e)
            }
            Info::Fn { head, args } => {
                let v: Vec<u32> = args
                    .iter()
                    .map(|&a| self.subst_at(a, map, depth + 1))
                    .collect();
                self.fn_node_by_id(head, &v)
            }
        }
    }
}
