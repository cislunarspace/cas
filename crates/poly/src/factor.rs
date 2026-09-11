//! 一元因式分解（M4，设计 §5）：ℚ[x] → 整数本原部分 → Yun 平方自由分解
//! → 模 p Cantor–Zassenhaus（DDF+EDF）→ 线性 Hensel 提升（Mignotte 界，
//! 二叉树两因子逐层）→ 确定性子集组合。
//!
//! **确定性**：CZ 的随机多项式取自固定种子的 xorshift，质数表固定，
//! 子集组合按二进制序枚举——同输入恒得同输出（D6 纪律）。
//!
//! 内部表示：稠密系数向量（index = 次数）。`IPoly` 为 i128 整系数；
//! `Vec<u64>` 配模数 p 或 p^k（系数 ∈ [0, 模数)，u128 中转）。

use crate::Poly;
use cas_domain::{Integer, Rational};
use std::sync::Arc;

const PRIMES: &[u64] = &[
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
];

/// CZ 随机源（固定种子）。
struct Det(u64);
impl Det {
    fn new() -> Self {
        Det(0xC0FF_EE01)
    }
    fn next(&mut self) -> u64 {
        let mut x = self.0 | 1;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
}

// ── 稠密整系数多项式 ──────────────────────────────────────────

type IPoly = Vec<i128>;

fn ip_trim(v: &mut IPoly) {
    while v.last() == Some(&0) {
        v.pop();
    }
}

fn ip_deg(v: &IPoly) -> i32 {
    v.len() as i32 - 1
}

/// 精确除法（首项系数须逐项整除，余数须为零），否则 None。
fn ip_exact_div(a: &IPoly, b: &IPoly) -> Option<IPoly> {
    if b.is_empty() {
        return None;
    }
    if a.is_empty() {
        return Some(vec![]);
    }
    if ip_deg(a) < ip_deg(b) {
        return None;
    }
    let db = ip_deg(b) as usize;
    let lb = *b.last()?;
    let mut r = a.clone();
    let mut q = vec![0i128; a.len() - b.len() + 1];
    while ip_deg(&r) >= ip_deg(b) && !r.is_empty() {
        let dr = ip_deg(&r) as usize;
        let lr = *r.last()?;
        if lr % lb != 0 {
            return None;
        }
        let t = lr / lb;
        q[dr - db] = t;
        for (i, &c) in b.iter().enumerate() {
            r[dr - db + i] -= t * c;
        }
        ip_trim(&mut r);
    }
    if r.is_empty() { Some(q) } else { None }
}

/// 本原部分与 ℤ 内容：f = cont · pp，pp 整系数且系数 gcd 为 1、首项为正。
fn ip_primitive(f: &IPoly) -> (i128, IPoly) {
    if f.is_empty() {
        return (0, vec![]);
    }
    let mut g: i128 = 0;
    for &c in f {
        g = gcd_i128(g, c.abs());
    }
    if g == 0 {
        g = 1;
    }
    let mut pp: IPoly = f.iter().map(|&c| c / g).collect();
    if *pp.last().unwrap() < 0 {
        for c in &mut pp {
            *c = -*c;
        }
        g = -g;
    }
    (g, pp)
}

fn gcd_i128(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

// ── 𝔽_p 稠密多项式 ────────────────────────────────────────────

type Fp = Vec<u64>;

fn fp_trim(v: &mut Fp) {
    while v.last() == Some(&0) {
        v.pop();
    }
}

fn fp_deg(v: &Fp) -> i32 {
    v.len() as i32 - 1
}

fn fp_monic(mut v: Fp, p: u64) -> Fp {
    fp_trim(&mut v);
    if let Some(&l) = v.last() {
        let inv = fp_pow(l, p - 2, p);
        for c in &mut v {
            *c = *c % p * inv % p;
        }
    }
    v
}

fn fp_pow(mut b: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1u64;
    b %= p;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % p;
        }
        b = b * b % p;
        e >>= 1;
    }
    r
}

fn fp_sub(a: &Fp, b: &Fp, p: u64) -> Fp {
    let mut out = vec![0u64; a.len().max(b.len())];
    for (i, &x) in a.iter().enumerate() {
        out[i] = x;
    }
    for (i, &y) in b.iter().enumerate() {
        out[i] = (out[i] + p - y % p) % p;
    }
    fp_trim(&mut out);
    out
}

fn fp_mul(a: &Fp, b: &Fp, p: u64) -> Fp {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut out = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if x == 0 {
            continue;
        }
        for (j, &y) in b.iter().enumerate() {
            out[i + j] =
                ((out[i + j] as u128 + x as u128 * y as u128 % p as u128) % p as u128) as u64;
        }
    }
    fp_trim(&mut out);
    out
}

fn fp_divrem(a: &Fp, b: &Fp, p: u64) -> (Fp, Fp) {
    assert!(!b.is_empty(), "fp 除式为零");
    let mut r = a.clone();
    let mut q = vec![0u64; (a.len() as i32 - fp_deg(b)).max(0) as usize + 1];
    let db = fp_deg(b) as usize;
    let inv = fp_pow(*b.last().unwrap(), p - 2, p);
    while !r.is_empty() && fp_deg(&r) >= db as i32 {
        let dr = r.len() - 1;
        let t = *r.last().unwrap() * inv % p;
        q[dr - db] = t;
        for (i, &c) in b.iter().enumerate() {
            r[dr - db + i] = (r[dr - db + i] + p - t * c % p) % p;
        }
        fp_trim(&mut r);
    }
    fp_trim(&mut q);
    (q, r)
}

/// 扩展 Euclid：返回 (monic g, s, t)，s·a + t·b = g。
fn fp_xgcd(a: &Fp, b: &Fp, p: u64) -> (Fp, Fp, Fp) {
    let (mut r0, mut r1) = (a.to_vec(), b.to_vec());
    let (mut s0, mut s1) = (vec![1u64], vec![]);
    let (mut t0, mut t1) = (vec![], vec![1u64]);
    while !r1.is_empty() {
        let (q, r) = fp_divrem(&r0, &r1, p);
        r0 = r1;
        r1 = r;
        let s2 = fp_sub(&s0, &fp_mul(&q, &s1, p), p);
        s0 = s1;
        s1 = s2;
        let t2 = fp_sub(&t0, &fp_mul(&q, &t1, p), p);
        t0 = t1;
        t1 = t2;
    }
    if r0.is_empty() {
        return (vec![], vec![], vec![]);
    }
    let inv = fp_pow(*r0.last().unwrap(), p - 2, p);
    let mul = |v: &Fp| -> Fp { v.iter().map(|&c| c * inv % p).collect() };
    (mul(&r0), mul(&s0), mul(&t0))
}

fn fp_gcd(a: &Fp, b: &Fp, p: u64) -> Fp {
    fp_xgcd(a, b, p).0
}

fn fp_deriv(a: &Fp, p: u64) -> Fp {
    let mut out = vec![0u64; a.len().saturating_sub(1)];
    for i in 1..a.len() {
        out[i - 1] = a[i] * (i as u64) % p;
    }
    fp_trim(&mut out);
    out
}

/// x^e mod f（平方乘，指数 u128——p^d 类指数会超 u64）。
fn fp_powmod(mut b: Fp, mut e: u128, f: &Fp, p: u64) -> Fp {
    let mut r: Fp = vec![1];
    fn rem(v: Fp, f: &Fp, p: u64) -> Fp {
        let (_, r) = fp_divrem(&v, f, p);
        r
    }
    b = rem(b, f, p);
    while e > 0 {
        if e & 1 == 1 {
            r = rem(fp_mul(&r, &b, p), f, p);
        }
        b = rem(fp_mul(&b, &b, p), f, p);
        e >>= 1;
    }
    r
}

/// 确定性随机 𝔽_p 多项式，次数 < n 且非零。
fn fp_rand(n: usize, p: u64, det: &mut Det) -> Fp {
    let mut v: Fp = (0..n).map(|_| det.next() % p).collect();
    fp_trim(&mut v);
    v
}

/// 不同次数分解（f monic squarefree）：返回 (次数 d 的因子之积, d) 列表。
fn ddf(f: &Fp, p: u64) -> Vec<(Fp, u32)> {
    let mut out = vec![];
    let mut fcur = f.clone();
    let mut h: Fp = vec![0, 1]; // x
    let mut d = 1u32;
    loop {
        h = fp_powmod(h, p as u128, &fcur, p); // x^{p^d}
        // DDF 判据：gcd(f, x^{p^d} − x)——减 x，不是减 1
        let xx: Fp = vec![0, 1];
        let g = fp_gcd(&fcur, &fp_sub(&h, &xx, p), p);
        if fp_deg(&g) > 0 {
            let q = fp_monic(g, p);
            // 除掉该因子后应取商（余式必为零）
            let (quot, rem) = fp_divrem(&fcur, &q, p);
            debug_assert!(rem.is_empty(), "ddf 因子必整除");
            out.push((q.clone(), d));
            fcur = quot;
            // 关键：h 须对"剩余部分"取模（新 fcur），对提取因子取模会断
            // 同余链 x^{p^d} ≡ h (mod fcur)，使后续轮次的 d 归属整体慢一拍
            //（x^34-1 实测：16 次因子被标到 d=17）
            h = {
                let (_, r) = fp_divrem(&h, &fcur, p);
                r
            };
        }
        let df = fp_deg(&fcur);
        if df <= 0 {
            break;
        }
        if df <= d as i32 {
            if df > 0 {
                out.push((fcur.clone(), df as u32));
            }
            break;
        }
        d += 1;
    }
    out
}

/// 等次数分解（Cantor–Zassenhaus，f 的全部不可约因子次数 = d）。
fn edf(f: &Fp, d: u32, p: u64, det: &mut Det) -> Vec<Fp> {
    let n = fp_deg(f);
    if n <= d as i32 {
        return vec![f.clone()];
    }
    loop {
        let q = fp_rand(f.len(), p, det);
        if q.is_empty() {
            continue;
        }
        // t = q^{(p^d − 1)/2} mod f
        let e = ((p as u128).pow(d) - 1) / 2;
        let t = fp_powmod(q, e, f, p);
        if t.is_empty() {
            continue;
        }
        let one: Fp = vec![1];
        let g = fp_gcd(f, &fp_sub(&t, &one, p), p);
        let dg = fp_deg(&g);
        if dg > 0 && dg < n {
            // 剩余因子 = 商（g 整除 f，余式恒为零）
            let (quot, rem) = fp_divrem(f, &g, p);
            debug_assert!(rem.is_empty(), "edf 因子必整除");
            let mut out = edf(&fp_monic(g, p), d, p, det);
            out.extend(edf(&fp_monic(quot, p), d, p, det));
            return out;
        }
    }
}

/// 模 p 完全分解（monic squarefree 输入）→ monic 不可约因子列表。
fn cz_factor(f_monic_sqfree: &Fp, p: u64, det: &mut Det) -> Vec<Fp> {
    let mut out = vec![];
    for (g, d) in ddf(f_monic_sqfree, p) {
        out.extend(edf(&g, d, p, det));
    }
    out
}

// ── 线性两因子 Hensel（模数 p^step，s·a + t·b ≡ 1 mod p 恒有效）──

/// 把 A ≡ a·b (mod p) 提升到 A ≡ a'·b' (mod p^k)。A、a、b 为非负稠密
/// 系数（mod p^k 视域）。返回 (a', b')。
#[allow(clippy::too_many_arguments)]
/// 提升不变量破坏时返回 None（保守回退：整体不分解，绝不产生错误因子）。
/// 已知在高次复合因子上偶发（见 DESIGN §17 挂账），修复前以保守行为兜底。
fn hensel2(
    target: &[u64],
    a0: &Fp,
    b0: &Fp,
    _s: &Fp,
    t: &Fp,
    p: u64,
    k: u32,
    pk: u64,
) -> Option<(Vec<u64>, Vec<u64>)> {
    let mut a: Vec<u64> = a0.to_vec();
    let mut b: Vec<u64> = b0.to_vec();
    let mut step_mod = p; // p^step
    for _step in 1..k {
        // prod = a·b mod p^{step+1}
        let next_mod = step_mod * p; // p^{step+1}
        let prod = umul(&a, &b, next_mod);
        // err = (A − prod)/p^step mod p
        let mut err: Fp = vec![];
        for i in 0..prod.len().max(target.len()) {
            let ai = target.get(i).copied().unwrap_or(0) % next_mod;
            let pi = prod.get(i).copied().unwrap_or(0);
            let diff = (ai + next_mod - pi) % next_mod;
            if diff % step_mod != 0 {
                return None; // 提升不变量破坏：保守放弃
            }
            err.push(diff / step_mod % p);
        }
        fp_trim(&mut err);
        if !err.is_empty() {
            // 解 σ·b + τ·a ≡ err (mod p)：σ = (err·t) mod a 后，τ 由精确商
            // (err − σ·b)/a 给出——两者不能各自独立取模（否则差 ab 的倍数）
            let sig = {
                let prod_e = fp_mul(&err, t, p);
                let (_, r) = fp_divrem(&prod_e, a0, p);
                r
            };
            let tau = {
                let pr = fp_sub(&err, &fp_mul(&sig, b0, p), p);
                let (q, r) = fp_divrem(&pr, a0, p);
                if !r.is_empty() {
                    return None; // 校正整除性破坏：保守放弃
                }
                q
            };
            for (i, &c) in sig.iter().enumerate() {
                a[i] = (a[i] + step_mod % pk * c % pk) % pk;
            }
            for (i, &c) in tau.iter().enumerate() {
                b[i] = (b[i] + step_mod % pk * c % pk) % pk;
            }
        }
        step_mod = next_mod;
    }
    trim_u(&mut a);
    trim_u(&mut b);
    Some((a, b))
}

fn trim_u(v: &mut Vec<u64>) {
    while v.last() == Some(&0) {
        v.pop();
    }
}

fn umul(a: &[u64], b: &[u64], m: u64) -> Vec<u64> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut out = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if x == 0 {
            continue;
        }
        for (j, &y) in b.iter().enumerate() {
            out[i + j] = (out[i + j] as u128 + x as u128 * y as u128 % m as u128) as u64 % m;
        }
    }
    trim_u(&mut out);
    out
}

/// 二叉树多重 Hensel：f̄（monic，模 p^k 表示的目标）≡ Π G_i (mod p^k)，
/// G_i ≡ g_i (mod p)。叶子顺序与 g_i 一致。
fn hensel_all(target: &[u64], facs: &[Fp], p: u64, k: u32, pk: u64) -> Vec<Vec<u64>> {
    if facs.len() == 1 {
        return vec![target.to_vec()];
    }
    let mid = facs.len() / 2;
    let (left, right) = facs.split_at(mid);
    let a0 = left
        .iter()
        .cloned()
        .reduce(|x, y| fp_mul(&x, &y, p))
        .unwrap();
    let b0 = right
        .iter()
        .cloned()
        .reduce(|x, y| fp_mul(&x, &y, p))
        .unwrap();
    let (_, s, t) = fp_xgcd(&a0, &b0, p);
    let Some((al, br)) = hensel2(target, &a0, &b0, &s, &t, p, k, pk) else {
        // 保守回退：所有因子之积按 target 返回（等同整体不分解）
        return vec![target.to_vec(); facs.len()];
    };
    let mut out = hensel_all(&al, left, p, k, pk);
    out.extend(hensel_all(&br, right, p, k, pk));
    out
}

// ── Zassenhaus 主流程 ─────────────────────────────────────────

/// 本原 squarefree S（首项正）的不可约因子（ℤ[x]，本原、首项正）。
fn zassenhaus(s: &IPoly) -> Vec<IPoly> {
    let d = ip_deg(s);
    if d <= 1 {
        return vec![s.to_vec()];
    }
    let lc = *s.last().unwrap() as i64; // 界内安全
    let amax: i128 = s.iter().map(|c| c.abs()).max().unwrap();
    // 选质数：不除尽首项、模 p 后仍 squarefree
    let mut chosen: Option<(u64, Fp)> = None;
    for &p in PRIMES {
        if lc % p as i64 == 0 {
            continue;
        }
        let fmod: Fp = s
            .iter()
            .map(|&c| (c.rem_euclid(p as i128)) as u64)
            .collect();
        let fmod = fp_monic(fmod, p);
        let df = fp_deriv(&fmod, p);
        if fp_gcd(&fmod, &df, p).len() <= 1 {
            chosen = Some((p, fmod));
            break;
        }
    }
    let (p, fmod) = match chosen {
        Some(x) => x,
        None => {
            return vec![s.to_vec()]; // 找不到好质数：视为不可约（测试范围外）
        }
    };
    let mut det = Det::new();
    let facs_p = cz_factor(&fmod, p, &mut det);
    if facs_p.len() <= 1 {
        return vec![s.to_vec()]; // 模 p 不可约 ⇒ ℚ 上不可约
    }
    // Mignotte 界与提升模数；上限 2^60（u64 乘法中间量安全）——
    // 超界则保守不分解（多精度 Hensel 属后续工作，见 DESIGN）
    let b_bound: u128 = (2u128 << d.min(120)) * amax as u128 * lc.unsigned_abs() as u128;
    let mut pk = p as u128;
    let mut k = 1u32;
    while pk <= 2 * b_bound {
        pk *= p as u128;
        k += 1;
        if pk > 1 << 60 {
            return vec![s.to_vec()];
        }
    }
    let pk = pk as u64;
    // 目标：f̄ = lc^{-1}·s mod p^k（monic）
    let lc_inv = inv_mod(lc.rem_euclid(pk as i64) as u64, pk);
    // 注意 u128 中转：c 与 lc_inv 都可近 pk（≈2^60），u64 乘法静默溢出——
    // 曾致"整体符号翻转后 Hensel 第一步不变量破坏"（定点复现 seed=2）
    let a: Vec<u64> = s
        .iter()
        .map(|&c| ((c.rem_euclid(pk as i128) as u64 as u128 * lc_inv as u128) % pk as u128) as u64)
        .collect();
    let lifted = hensel_all(&a, &facs_p, p, k, pk);
    // 确定性子集组合
    combine(s, &lifted, pk, lc)
}

fn inv_mod(a: u64, m: u64) -> u64 {
    // 扩展 Euclid（m 不必素数）
    let (mut r0, mut r1) = (m as i128, a as i128 % m as i128);
    let (mut s0, mut s1) = (0i128, 1i128);
    while r1 != 0 {
        let q = r0 / r1;
        let r2 = r0 - q * r1;
        r0 = r1;
        r1 = r2;
        let s2 = s0 - q * s1;
        s0 = s1;
        s1 = s2;
    }
    debug_assert_eq!(r0, 1, "逆元不存在");
    s0.rem_euclid(m as i128) as u64
}

/// 子集组合：候选 = Π_{i∈S} G_i mod pk → 对称化 → 本原部分 → 试除。
fn combine(f: &IPoly, lifted: &[Vec<u64>], pk: u64, lc: i64) -> Vec<IPoly> {
    let r = lifted.len();
    if r > 16 {
        return vec![f.to_vec()]; // 组合爆炸保护（测试范围外）
    }
    let half = pk / 2;
    let sym = |v: Vec<u64>| -> IPoly {
        v.iter()
            .map(|&c| {
                if c > half {
                    c as i128 - pk as i128
                } else {
                    c as i128
                }
            })
            .collect()
    };
    let mut factors: Vec<IPoly> = vec![];
    let mut used = vec![false; r];
    let mut f_cur = f.clone();
    // 子集按基数升序枚举：避免把多个真因子的积当成单因子（漏拆）
    let mut masks: Vec<u64> = (1..(1u64 << (r - 1))).collect();
    masks.sort_by_key(|m| m.count_ones());
    loop {
        let mut found = false;
        for &mask in &masks {
            // 跳过已用因子
            if (0..r - 1).any(|i| mask >> i & 1 == 1 && used[i]) {
                continue;
            }
            if mask.count_ones() as usize > ip_deg(&f_cur).max(0) as usize {
                continue;
            }
            // 候选乘积（模 pk）；提升因子是首一的，真因子须乘 lc(f) 后
            // 再对称化取本原部分（经典 Zassenhaus 组合），同时试 ±两号
            let mut cand: Vec<u64> = vec![1];
            for (i, g) in lifted.iter().enumerate().take(r - 1) {
                if mask >> i & 1 == 1 {
                    cand = umul(&cand, g, pk);
                }
            }
            // 首项比例自由度：试 lc^j（j = 0..=|S|；|S| 个首一因子的积与
            // 真因子的首项差 lc 的某次幂——模 p 下因子内部降次时 j 可 >1）
            let lc_u = (lc as i128).rem_euclid(pk as i128) as u64;
            let k = mask.count_ones();
            let mut hit = None;
            let mut scaled = cand.clone();
            for j in 0..=k {
                if j > 0 {
                    scaled = umul(&scaled, &[lc_u], pk);
                }
                let cand_ip = ip_trim_mut(sym(scaled.clone()));
                if cand_ip.is_empty() {
                    continue;
                }
                let (_, pp) = ip_primitive(&cand_ip);
                if pp.is_empty() || ip_deg(&pp) == 0 || ip_deg(&pp) > ip_deg(&f_cur) {
                    continue;
                }
                if let Some(q) = ip_exact_div(&f_cur, &pp) {
                    hit = Some((pp, q));
                    break;
                }
            }
            if let Some((pp, q)) = hit {
                factors.push(pp);
                f_cur = q;
                for (i, u) in used.iter_mut().enumerate().take(r - 1) {
                    if mask >> i & 1 == 1 {
                        *u = true;
                    }
                }
                found = true;
                break;
            }
        }
        if !found {
            break;
        }
    }
    if !f_cur.is_empty() && ip_deg(&f_cur) > 0 {
        factors.push(f_cur);
    }
    factors
}

fn ip_trim_mut(mut v: IPoly) -> IPoly {
    ip_trim(&mut v);
    v
}

// ── Yun 平方自由分解（复用 Poly<Rational> 的 gcd/exact_div/deriv）──

fn squarefree_parts(f: &Poly<Rational>) -> Vec<(Poly<Rational>, u32)> {
    let df = f.deriv(0);
    if df.is_zero() {
        return vec![(f.clone(), 1)];
    }
    let c = f.gcd(&df);
    let mut w = match f.exact_div(&c) {
        Some(v) => v,
        None => return vec![(f.clone(), 1)],
    };
    let mut z = match df.exact_div(&c) {
        Some(v) => v.sub(&w.deriv(0)),
        None => return vec![(f.clone(), 1)],
    };
    let mut out = vec![];
    let mut i = 1u32;
    while !w.is_constant() {
        let g = w.gcd(&z);
        if !g.is_constant() {
            out.push((g.clone(), i));
            w = w.exact_div(&g).expect("Yun 整除");
            z = z.exact_div(&g).expect("Yun 整除").sub(&w.deriv(0));
        } else {
            z = z.sub(&w.deriv(0));
        }
        i += 1;
    }
    out
}

// ── 公共 API ─────────────────────────────────────────────────

impl Poly<Rational> {
    /// 一元 ℚ[x] 因式分解：返回（常数内容, [(本原因子, 重数)]），
    /// 满足 cont · Π 因子^重数 == self。因子本原、首项正；输出确定性。
    pub fn factor_univariate(&self) -> (Rational, Vec<(Poly<Rational>, u32)>) {
        assert_eq!(self.ring().nvars(), 1, "factor_univariate 仅支持一元");
        if self.is_zero() {
            return (Rational::zero(), vec![]);
        }
        if self.is_constant() {
            return (self.terms().next().unwrap().1.clone(), vec![]);
        }
        // 整数本原化
        let mut den_lcm = Integer::one();
        for (_, c) in self.terms() {
            den_lcm = den_lcm.mul(&c.den());
        }
        let mut num_gcd = Integer::zero();
        for (_, c) in self.terms() {
            let scaled = c.num().mul(&den_lcm.div_exact(&c.den()));
            num_gcd = if num_gcd.is_zero() {
                scaled
            } else {
                num_gcd.gcd(&scaled)
            };
        }
        if num_gcd.is_zero() {
            num_gcd = Integer::one();
        }
        let cont = Rational::from_ints(&num_gcd, &den_lcm).unwrap();
        let mut dense: IPoly = {
            let mut v = vec![0i128; self.degree() as usize + 1];
            for (e, c) in self.terms() {
                let scaled = c.mul(&Rational::from_integer(&den_lcm));
                let r = scaled.div(&Rational::from_integer(&num_gcd)).unwrap();
                assert!(r.den().is_one(), "本原化后必为整系数");
                v[e[0] as usize] = r.num().to_i64().expect("本原化后系数落 i64") as i128;
            }
            v
        };
        // 前置条件：zassenhaus 要求首项正；负号并入内容
        let mut cont = cont;
        if *dense.last().unwrap() < 0 {
            for c in &mut dense {
                *c = -*c;
            }
            cont = cont.neg();
        }
        // 重数分解 + Zassenhaus
        let ring = self.ring().clone();
        let mut factors: Vec<(Poly<Rational>, u32)> = vec![];
        for (s, m) in squarefree_parts(&from_dense(&dense, &ring)) {
            if s.is_constant() {
                continue;
            }
            let dense_s = to_dense(&s);
            for h in zassenhaus(&dense_s) {
                factors.push((from_dense(&h, &ring), m));
            }
        }
        // 确定性排序：次数升序，同次按系数字典序
        factors.sort_by(|a, b| cmp_dense(&to_dense(&a.0), &to_dense(&b.0)));
        (cont, factors)
    }
}

fn to_dense(p: &Poly<Rational>) -> IPoly {
    let mut v = vec![0i128; p.degree() as usize + 1];
    for (e, c) in p.terms() {
        assert!(c.den().is_one(), "to_dense 需整系数");
        v[e[0] as usize] = c.num().to_i64().expect("整系数落 i64") as i128;
    }
    v
}

fn from_dense(v: &IPoly, ring: &Arc<crate::PolyRing>) -> Poly<Rational> {
    let items: Vec<(Vec<u32>, Rational)> = v
        .iter()
        .enumerate()
        .filter(|(_, c)| **c != 0)
        .map(|(i, &c)| {
            (
                vec![i as u32],
                Rational::from_integer(&Integer::from_i64(c as i64)),
            )
        })
        .collect();
    Poly::from_terms(ring.clone(), items)
}

fn cmp_dense(a: &IPoly, b: &IPoly) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    match ip_deg(a).cmp(&ip_deg(b)) {
        Ordering::Equal => {}
        o => return o,
    }
    for i in (0..a.len()).rev() {
        match a[i].cmp(&b[i]) {
            Ordering::Equal => {}
            o => return o,
        }
    }
    std::cmp::Ordering::Equal
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{MonOrder, PolyRing};

    fn ring() -> Arc<PolyRing> {
        PolyRing::new(["x"], MonOrder::DegRevLex)
    }

    fn from_coeffs(cs: &[i64]) -> Poly<Rational> {
        let items: Vec<(Vec<u32>, Rational)> = cs
            .iter()
            .enumerate()
            .filter(|(_, c)| **c != 0)
            .map(|(i, &c)| {
                (
                    vec![i as u32],
                    Rational::from_ints(&Integer::from_i64(c), &Integer::from_i64(1)).unwrap(),
                )
            })
            .collect();
        Poly::from_terms(ring(), items)
    }

    fn check_refold(f: &Poly<Rational>) {
        let (cont, facs) = f.factor_univariate();
        let mut prod = Poly::constant(ring(), cont);
        for (g, m) in &facs {
            prod = prod.mul(&g.pow(*m));
        }
        assert_eq!(prod, *f, "重展开不等于原式: f={f:?} factors={facs:?}");
    }

    #[test]
    fn 二项式与已知分解() {
        let f = from_coeffs(&[-1, 0, 1]); // x^2 - 1
        let (c, facs) = f.factor_univariate();
        assert_eq!(facs.len(), 2);
        let mut prod = Poly::constant(ring(), c);
        for (g, m) in &facs {
            prod = prod.mul(&g.pow(*m));
        }
        assert_eq!(prod, f);

        // x^4 + 4 = (x² − 2x + 2)(x² + 2x + 2)
        let f = from_coeffs(&[4, 0, 0, 0, 1]);
        let (_, facs) = f.factor_univariate();
        assert_eq!(facs.len(), 2, "x^4+4 应分解为两个二次因子: {facs:?}");
        check_refold(&f);

        // x^4 + 1 在 ℚ 上不可约
        let f = from_coeffs(&[1, 0, 0, 0, 1]);
        let (_, facs) = f.factor_univariate();
        assert_eq!(facs.len(), 1);
    }

    #[test]
    fn x_n_minus_1_全族() {
        for n in 1..=64u32 {
            let mut cs = vec![0i64; n as usize + 1];
            cs[0] = -1;
            cs[n as usize] = 1;
            let f = from_coeffs(&cs);
            check_refold(&f);
        }
        // x^12 − 1 有 6 个不可约因子（分圆多项式个数 = d|12 的 φ(d)>0 的因子）
        let mut cs = vec![0i64; 13];
        cs[0] = -1;
        cs[12] = 1;
        let (_, facs) = from_coeffs(&cs).factor_univariate();
        assert_eq!(facs.len(), 6);
    }

    #[test]
    fn 随机积重展开() {
        let mut det = Det::new();
        for _ in 0..600 {
            let nf = 2 + det.next() % 2; // 2–3 个因子
            let mut f = Poly::constant(ring(), Rational::one());
            for _ in 0..nf {
                let deg = 1 + det.next() % 6;
                let cs: Vec<i64> = (0..=deg).map(|_| (det.next() % 13) as i64 - 6).collect();
                let g = from_coeffs(&cs);
                if g.is_zero() {
                    continue;
                }
                f = f.mul(&g);
            }
            // 有理内容
            let scale = Rational::from_ints(
                &Integer::from_i64((det.next() % 7) as i64 - 3),
                &Integer::from_i64(1 + (det.next() % 5) as i64),
            )
            .unwrap();
            let f = f.mul(&Poly::constant(ring(), scale));
            if !f.is_zero() && !f.is_constant() {
                check_refold(&f);
            }
        }
    }

    #[test]
    fn 重数与内容() {
        // (2x − 2)^3 = 2^3 (x − 1)^3
        let f = from_coeffs(&[-8, 24, -24, 8]);
        let (c, facs) = f.factor_univariate();
        assert_eq!(facs.len(), 1);
        assert_eq!(facs[0].1, 3);
        assert_eq!(
            c,
            Rational::from_ints(&Integer::from_i64(8), &Integer::from_i64(1)).unwrap()
        );
        check_refold(&f);
    }
}

#[cfg(test)]
mod hensel_regression {
    use super::*;

    /// u64 溢出回归：整体系数符号翻转后 Hensel 第一步不变量破坏
    /// （lc⁻¹·c 在 u64 下静默回绕；修复为 u128 中转）。
    #[test]
    fn 溢出回归_符号翻转() {
        let f1: IPoly = vec![-1, 2, 5, 2, -6];
        let f2: IPoly = vec![6, -3, 1, 1, 6, 3, -4];
        let f3: IPoly = vec![4, -6, 0, -1, -5];
        let mut f = vec![0i128; 15];
        for (i, &x) in f1.iter().enumerate() {
            for (j, &y) in f2.iter().enumerate() {
                for (k, &z) in f3.iter().enumerate() {
                    f[i + j + k] += x * y * z;
                }
            }
        }
        ip_trim(&mut f);
        assert_eq!(zassenhaus(&f).len(), 3, "原始版本");
        let mut neg = f.clone();
        for c in &mut neg {
            *c = -*c;
        }
        assert_eq!(zassenhaus(&neg).len(), 3, "符号翻转版本（曾触发 u64 溢出）");
    }
}
