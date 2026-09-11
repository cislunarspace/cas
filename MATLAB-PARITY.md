# cas 与 MATLAB Symbolic Math Toolbox 对标清单

日期：2026-09-12。对标口径：MATLAB 官方文档列出的符号工具箱核心操作，
逐项给出 cas 的对应 API、验证状态与已知差距。验证以 sympy 为可执行
oracle（语义对拍 = 随机点交叉精确/浮点求值）；MATLAB 本体因无许可
不在本机执行，语义以公开文档为准。

## 核心操作对标

| MATLAB 操作 | cas API | 状态 | 验证 |
|---|---|---|---|
| `sym`/`syms` 声明 | `ctx.sym` / `sym!` 宏 | ✅ | M1：往返 1000 例 |
| 算术与规范形 | 运算符重载（`+ - * / ^`） | ✅ | M1：构造路径无关 500 例 |
| `expand` | `ctx.expand` | ✅ | oracle 5000 例 100% 一致；(1+x+y+z+w)^20 98 ms（sympy 同款 959 ms） |
| `factor` | `ctx.factor` | ✅（有限制） | 单变元 1000 例 100%；双变元 300 例 265 一致/35 漏拆/0 错拆 |
| `simplify` | `ctx.simplify` | ✅（保守集） | oracle 300 例 100%（常量折叠 + sin²+cos²→1；非全自动最简） |
| `cancel` | `ctx.cancel` | ✅ | oracle 1000 例 100% |
| `diff` | `ctx.diff` | ✅ | oracle 300 例 100%；含链式/积商/幂三档/初等函数表 |
| `taylor` | `ctx.taylor` | ✅（限闭式域） | oracle 300 例 100%（展开点各阶导数可精确求值；非闭式点如 cos(4+x²)@0 暂返回 0） |
| `subs` | `ctx.subst` | ✅ | M2 单测（恒等代换回原节点） |
| 数值求值 `subs`/`vpa` | `ctx.eval_rational` / `ctx.eval_float` | ✅ | M1/M2 单测 + 全部 oracle 的求值通道 |
| 表达式化简 `collect` | （未单独提供） | ➖ | collect 语义可由 expand+规范形覆盖 |
| `assume` | （D5 已设计未实现） | ➖ | P1 后期：闭谓词集 |
| `solve` | 非目标 | — | 设计文档明确排除（P2 起按需评估） |
| `int` 积分 | 非目标 | — | 同上 |
| 符号矩阵 `det`/`rank` | 非目标（P2） | — | — |

## 已知差距（按影响排序）

1. **taylor 闭式域限制**：实现基于"各阶导数在展开点的精确求值"，
   非闭式点返回 0。正解为级数算术（对截断级数直接运算），属后续工作。
2. **simplify 保守**：常量折叠 + 无条件恒等式的定向规则集，不做
   假设推理与启发式搜索（确定性优先，设计 D4 的明示取舍）。
3. **双变元 factor 的 35/300 漏拆**与 PRS 高次膨胀（M5b，子结果式）。
4. **假设系统未实现**：`assume`/`assumeAlso` 对应能力待 P1 后期。

## 验证方式汇总

oracle（sympy 1.14，版本锁定）六种对拍模式：expand 5000 例 / cancel
1000 例 / factor 单变元 1000 + 双变元 300 例 / diff 300 例 / taylor
300 例 / simplify 300 例——除 factor 双变元为"通过-保守（0 错拆）"外
全部 100% 语义一致。语料种子固定、判据分档（错拆硬失败、漏拆记录、
次数守恒校验）。运行方式：

```bash
cd cas && cargo run --release -p xtask -- oracle --op <expand|cancel|factor|diff|taylor|simplify> [--cases N]
```
