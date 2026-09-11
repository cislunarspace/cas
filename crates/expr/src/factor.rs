//! expr 层因式分解桥（M4/M5）：纯多项式表达式 → `Poly::factor`
//! （1/2 变元完全，≥3 变元部分分解）→ 规范形重建 `cont · Π 因子^重数`。
//!
//! 非多项式（函数/浮点/负幂）成分时返回原节点。

use crate::poly_bridge;

impl crate::Inner {
    pub(crate) fn factor_at(&mut self, id: u32) -> u32 {
        let (ring, var_ids, p) = match poly_bridge::to_poly(self, id) {
            Some(x) => x,
            None => return id,
        };
        let (cont, facs) = p.factor();
        if facs.is_empty() {
            return id; // 常数
        }
        let mut factors: Vec<u32> = Vec::with_capacity(facs.len() + 1);
        factors.push(self.lit_rational(&cont));
        for (g, m) in facs {
            let ge = poly_bridge::from_poly(self, &ring, &var_ids, &g);
            if m == 1 {
                factors.push(ge);
            } else {
                let mi = self.lit_int(m as i64);
                factors.push(self.make_pow(ge, mi));
            }
        }
        self.make_mul(&factors)
    }
}
