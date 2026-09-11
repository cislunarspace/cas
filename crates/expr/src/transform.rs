//! 显式变换（L1）：expand 与 subst（DESIGN.md D4）。
//!
//! M2 的 expand 用 DAG 分配实现：每一步都经规范形构造器，同类项即时合并，
//! 中间规模被规范形钳制在真实支撑集大小（(x+y+z+w)^20 全程 ≤ 1771 项）。
//! 多项式内核承接的快速 expand（D3 桥接）属 M3。

use crate::Inner;
use cas_domain::Rational;
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

    /// 有理函数约化：分解为 数值系数 · 额外因子 · num/den（全部因子须为
    /// 纯多项式），gcd(num, den) 约去后经规范形构造器重建。非多项式成分
    /// （函数、非整数指数）或无公因子时返回原节点（同 id）。
    /// 退化点约定与 make_mul 一致（x/x → 1）。
    pub(crate) fn cancel_at(&mut self, id: u32) -> u32 {
        let args: Vec<u32> = match &self.nodes[id as usize] {
            Node::Mul { args: sp } => self.node_args(*sp).to_vec(),
            _ => vec![id],
        };
        let mut coeff = Rational::one();
        let mut extras: Vec<u32> = Vec::new(); // 浮点因子等，原样保留
        let mut num: Vec<u32> = Vec::new();
        let mut den: Vec<u32> = Vec::new();
        // 先取 owned 判定（避免边读 nodes 边构造的借用冲突）
        enum Piece {
            Coeff(Rational),
            Extra,
            Num,
            Den(u32, i64),
        }
        for &a in &args {
            let piece = match &self.nodes[a as usize] {
                Node::Int(v) => Piece::Coeff(Rational::from_integer(v)),
                Node::Rat(r) => Piece::Coeff(r.clone()),
                Node::Float { .. } => Piece::Extra,
                Node::Pow { base, exp } => {
                    let k = match &self.nodes[*exp as usize] {
                        Node::Int(v) => v.to_i64(),
                        _ => None,
                    };
                    match k {
                        Some(k) if k > 0 => Piece::Num,
                        Some(k) => Piece::Den(*base, -k),
                        None => return id, // 非整数指数：不强行约化
                    }
                }
                _ => Piece::Num,
            };
            match piece {
                Piece::Coeff(v) => coeff = coeff.mul(&v),
                Piece::Extra => extras.push(a),
                Piece::Num => num.push(a),
                Piece::Den(base, m) => {
                    let lit = self.lit_int(m);
                    den.push(self.make_pow(base, lit));
                }
            }
        }
        if den.is_empty() {
            return id; // 无分母：构造器已做能做的合并
        }
        let n_id = if num.is_empty() {
            self.lit_int(1)
        } else {
            self.make_mul(&num)
        };
        let d_id = self.make_mul(&den);
        // 公共环：分子分母变元的并集（按名字节序）
        let mut syms: Vec<u32> = Vec::new();
        crate::poly_bridge::collect_syms(self, n_id, &mut syms, 0);
        crate::poly_bridge::collect_syms(self, d_id, &mut syms, 0);
        let (ring, var_ids, vi) = crate::poly_bridge::ring_for(self, &syms);
        let np = match crate::poly_bridge::to_poly_with(self, n_id, &ring, &vi) {
            Some(p) => p,
            None => return id,
        };
        let dp = match crate::poly_bridge::to_poly_with(self, d_id, &ring, &vi) {
            Some(p) => p,
            None => return id,
        };
        if dp.is_constant() {
            return id; // 常分母已由构造器折叠
        }
        if np.is_zero() {
            return self.lit_int(0);
        }
        let g = np.gcd(&dp);
        if g.is_constant() {
            return id; // 互素：无可约化
        }
        let n2 = np.exact_div(&g).expect("gcd 整除分子");
        let d2 = dp.exact_div(&g).expect("gcd 整除分母");
        let mut n2 = n2;
        // 分母约成常数时并入系数（gcd 本原规范化的符号差异在此吸收：
        // 如 (x^2-y^2)/(x-y) 的 gcd 归一为 y-x，d2 = -1）；再把负号沉入
        // 多项式，使输出首系数为正——得到与手写一致的规范形态。
        if d2.is_constant() {
            let cd = d2
                .terms()
                .next()
                .map(|(_, c)| c.clone())
                .unwrap_or_else(Rational::one);
            coeff = coeff.mul(&cd.inv_reduced().expect("gcd 非零"));
            if coeff.is_negative() && !n2.is_constant() {
                coeff = coeff.neg();
                n2 = n2.neg();
            }
        }
        let mut factors: Vec<u32> = Vec::with_capacity(4);
        factors.push(self.lit_rational(&coeff));
        factors.extend(extras.iter().copied());
        factors.push(crate::poly_bridge::from_poly(self, &ring, &var_ids, &n2));
        if !d2.is_constant() {
            let d_expr = crate::poly_bridge::from_poly(self, &ring, &var_ids, &d2);
            let neg1 = self.lit_int(-1);
            factors.push(self.make_pow(d_expr, neg1));
        }
        self.make_mul(&factors)
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
