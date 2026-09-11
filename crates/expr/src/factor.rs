//! expr 层因式分解（M4 桥）：纯一元多项式表达式 → `factor_univariate`
//! → 规范形重建 `cont · Π 因子^重数`。
//!
//! 非一元（含多符号）、非多项式（函数/浮点/负幂）成分时返回原节点。

use crate::poly_bridge;

impl crate::Inner {
    pub(crate) fn factor_at(&mut self, id: u32) -> u32 {
        let (ring, var_ids, p) = match poly_bridge::to_poly(self, id) {
            Some(x) => x,
            None => return id,
        };
        if ring.nvars() != 1 {
            return id; // 多变元因式分解属 M5
        }
        let (cont, facs) = p.factor_univariate();
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
