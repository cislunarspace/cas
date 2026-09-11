//! 固定键哈希（D6 实现纪律：全库禁 `RandomState`）。
//!
//! 哈希值只用于 intern 桶分配，从不直接进入输出；即使如此也统一用
//! 固定种子的 FxHash，保证任何哈希相关状态跨进程可复现。

use crate::node::Node;
use std::hash::{BuildHasherDefault, Hash, Hasher};

pub(crate) type FxBuild = BuildHasherDefault<FxHasher>;

const K: u64 = 0x517c_c1b7_2722_0a95;

#[derive(Default)]
pub(crate) struct FxHasher {
    hash: u64,
}

impl FxHasher {
    fn add(&mut self, v: u64) {
        self.hash = (self.hash.rotate_left(5) ^ v).wrapping_mul(K);
    }
}

impl Hasher for FxHasher {
    fn write(&mut self, bytes: &[u8]) {
        let mut rest = bytes;
        while rest.len() >= 8 {
            let mut buf = [0u8; 8];
            buf.copy_from_slice(&rest[..8]);
            self.add(u64::from_le_bytes(buf));
            rest = &rest[8..];
        }
        if !rest.is_empty() {
            let mut buf = [0u8; 8];
            buf[..rest.len()].copy_from_slice(rest);
            self.add(u64::from_le_bytes(buf));
        }
    }

    fn write_u8(&mut self, v: u8) {
        self.add(v as u64);
    }

    fn write_u32(&mut self, v: u32) {
        self.add(v as u64);
    }

    fn write_u64(&mut self, v: u64) {
        self.add(v);
    }

    fn write_i64(&mut self, v: i64) {
        self.add(v as u64);
    }

    fn write_usize(&mut self, v: usize) {
        self.add(v as u64);
    }

    fn finish(&self) -> u64 {
        self.hash
    }
}

/// 节点的结构哈希：类别标记 + 字段 + 子节点哈希（子节点已是规范形，
/// 其哈希存于 `hashes`，按索引取用，不重算）。
pub(crate) fn hash_node(node: &Node, args: &[u32], hashes: &[u64]) -> u64 {
    let mut h = FxHasher::default();
    let tag = match node {
        Node::Int(_) => 0,
        Node::Rat(_) => 1,
        Node::Float { .. } => 2,
        Node::Sym(_) => 3,
        Node::Fn { .. } => 4,
        Node::Pow { .. } => 5,
        Node::Mul { .. } => 6,
        Node::Add { .. } => 7,
    };
    h.write_u8(tag);
    match node {
        Node::Int(v) => v.hash(&mut h),
        Node::Rat(r) => r.hash(&mut h),
        Node::Float { bits, .. } => h.write_u64(*bits),
        Node::Sym(s) => h.write_u32(*s),
        Node::Fn { head, .. } => {
            h.write_u32(*head);
            for &a in args {
                h.write_u64(hashes[a as usize]);
            }
        }
        Node::Pow { base, exp } => {
            h.write_u64(hashes[*base as usize]);
            h.write_u64(hashes[*exp as usize]);
        }
        Node::Mul { .. } | Node::Add { .. } => {
            h.write_usize(args.len());
            for &a in args {
                h.write_u64(hashes[a as usize]);
            }
        }
    }
    h.finish()
}
