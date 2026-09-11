# cas：Rust 通用符号代数引擎 · 设计文档 v0.1

- 状态：草案，待评审。
- 日期：2026-09-11。生态调研结论截至本日（§2 附来源）。
- 定位：可当依赖用的泛用 CAS 库——精确算术、表达式树、化简、代换、展开、求导、
  多项式代数（含 Gröbner）、代码生成。性能对标 C++ 系内核（FLINT / Singular 量级），
  oracle 用 sympy；不是 sympy 的 Rust 复刻。
- 本文档只定设计与验收，不含实现。评审通过后：D1–D6 各拆一篇 ADR（`docs/adr/016+`）；
  cas 子工作区获得自己的 CONTEXT.md（届时根目录引入 CONTEXT-MAP.md），术语先见 §12。

---

## 1. 目标与非目标

目标（按优先级排序，冲突时按序取舍）：

1. **正确**：精确符号计算——ℤ/ℚ 多项式环全套、初等函数求导与级数、Gröbner 基、
   符号线性代数。正确性有差分测试与属性测试兜底（§8），不靠人眼。
2. **确定**：表达式与项的确定性全序进入公共 API 承诺（D6）。打印、序列化、代码生成
   三方同序。动机是实际踩过的坑：符号推导产物与数值生成器之间行序失同步会静默产出
   错误数值——确定性排序是防这类事故的基础设施，不是锦上添花。
3. **快**：单核显著快于 sympy（核心负载目标 10× 量级）；中期在 expand / factor /
   Gröbner 上接近 Singular / Symbolica。性能分级承诺见 §7：正确性是硬门槛，
   性能是带目标的跟踪指标。
4. **库形态**：部分 crate no_std、无 Python 运行时、依赖最小、可离线生成数值代码
   （Rust / C，含公共子表达式消除）。

非目标：

- 不做 notebook、绘图、物理单位等应用层。
- 不承诺"全自动最简化"：simplify 是定向规则 + 代价贪心（D4），确定性优先于"最简"。
- P0–P2 不做多项式系统之外的通用方程求解；积分不在任何已排期内。
- 不自造大整数：数字后端复用现有库（D2），且后端可换。

## 2. 生态调研结论（2026-09-11 检索）

| 项目 | 许可 | 定位 | 对本项目 |
|---|---|---|---|
| Symbolica | 非标准"source available"：个人免费、非商业单核免费、商业付费 | 最成熟的 Rust 原生 CAS，性能标杆 | 不可当依赖；架构参考与性能对照 |
| numerica | MIT | Symbolica 1.0 拆出的数值层：自动升降精度整数、有理数、有限域、误差追踪浮点、双数；后端可插 rug/malachite | **数字后端首选候选**（D2） |
| malachite | LGPL-3.0-only | 纯 Rust 任意精度算术 | 不做默认依赖（许可），可经 numerica 透传 |
| Algebraeon | GPL-3.0 | 纯 Rust 计算代数：Berlekamp–Zassenhaus / van Hoeij 因式分解、精确线代、代数数 | 不可链接；作算法参考与第二 oracle |
| Savage | AGPL-3.0 | 自述"primitive"的玩具级 CAS | 无关 |
| mathcore | MIT | 0.3.x，8 commits，很年轻 | 观察，不依赖 |
| egg | MIT | egraph 重写引擎，不是 CAS 底座 | 核心不用；P3+ 可作离线实验工具（D4） |
| symengine | MIT 系 C++ 库的 Rust 绑定 | 绑定层，构建重、非 Rust 原生 | 不采用 |
| sympy | BSD | 功能最全的参照系 | oracle 与功能对标；性能不对标 |

结论成立：**许可干净（MIT/Apache 兼容）、Rust 原生、覆盖 表达式→多项式→Gröbner→代码生成
全链的依赖不存在**。numerica 的出现把"数字层"从自建清单里划掉了，但表达式层、化简/假设、
多项式算法编排、确定性契约、代码生成仍然全缺。

## 3. crate 划分与依赖图

```
依赖方向自上而下，无环：

        ┌──────────────────────────────────────────────┐
        │                     cas                      │  facade：simplify 编排、prelude
        └───────┬──────────────┬───────────────┬───────┘
        ┌───────┴──────┐ ┌─────┴──────┐ ┌──────┴───────┐
        │    parse     │ │    print   │ │   codegen    │
        └───────┬──────┘ └─────┬──────┘ └──────┬───────┘
                └──────────────┼───────────────┘
                        ┌──────▼───────┐
                        │     expr     │   arena + hash-consing、假设、规范化、poly 桥接
                        └──┬────────┬──┘
                  ┌────────▼──┐  ┌──▼────────┐
                  │   poly    │  │  domain   │
                  └────────┬──┘  └───────────┘
                           └─► domain
```

| crate | 对外接口（宽度刻意控制） | 藏住的复杂度（深度来源） |
|---|---|---|
| `domain` | 4 个 trait（Semiring/Ring/EuclideanDomain/Field）+ ℤ、ℚ 两个现成域 | 混合内联小数/大数、后端切换（numerica/num-bigint/rug）、ℚ 规范化 |
| `poly` | `PolyRing` 身份 + `Poly<D>` 的约 12 个代数方法 + `terms()` 稳定迭代器 | 项序、稀疏算法（堆乘、PRS、Hensel）、稠密切换、Gröbner 全套 |
| `expr` | `Context` + `Expr`（运算符重载、expand/factor/subst/diff…） | arena/hash-consing、构造即规范化、假设闭包、poly 往返桥 |
| `parse` | `parse(&mut Context, &str) -> Result<Expr>` | 文法、错误位置 |
| `print` | `plain(&Context, &Expr) -> String`、`latex(...)` | 布局、括号极小化、LaTeX 转义 |
| `codegen` | `codegen_rust / codegen_c（target, 输入集）-> 源码字符串` | CSE、发射调度、数值一致性 |
| `cas` | prelude + `simplify` | 跨 crate 规则编排（规则表本身） |

两点说明：

- **新增 `domain`**（不在原清单）：`poly` 必须泛型于系数域，且不得反向依赖表达式层；
  把 Domain trait 放进 `expr` 会让 poly→expr 纠缠。代价是多一个只有几百行接口的 crate，
  换来 poly 可独立复用与测试。
- **`cas` facade 的浅是角色不是缺陷**，用假想删除法检验：删掉它，simplify 的规则编排会在
  每个调用方重现——它在扛活，保留。

## 4. 六项设计决策

### D1 表达式表示：arena 索引 DAG + 结构化哈希（hash-consing）

| 选项 | 结论 |
|---|---|
| A. `Box<Expr>` 递归树 | 否。无结构共享（嵌套幂、重复代换链内存指数膨胀）、逐节点堆分配与指针追踪、相等性需整树遍历、深递归爆栈 |
| B. **arena 索引 + hash-consing（推荐）** | 节点存 `Vec<Node>`，`Expr` 是 32 位索引；构造时经 intern 表去重 |
| C. 持久化容器（`im` 等） | 否。仍是逐节点分配，无 O(1) 相等性，收益与 B 重叠但代价更高 |

选 B 的理由：

- **规范形即身份**。hash-consing 保证结构相同的子表达式只有一个节点：
  相等性 = 索引相等，O(1)；去重、缓存、记忆化天然可用。
- **结构共享**。`(1+x+y)^n` 的中间量、反复代换产生的重复子式只存一份，
  这是正规化推导（Lie 级数、Poisson 括号链）的常态负载。
- **内存与 cache**。节点定长（目标 24 字节 + 子节点连续存于参数槽），
  遍历是索引算术而非指针追踪；`Expr` 是 12 字节 `Clone` 句柄（共享引用 + u32 索引），
  传参与存储廉价。
- **并行安全（预留）**。arena 只追加不修改，P4 可演进为原子追加的共享 arena
  （Symbolica 同路线）或线程本地子 arena + 规范合并；P0–P3 单线程所有权即可。

明确的代价与对策：

- **句柄无生命周期参数**（M1 实测定稿）：`Expr` 携带 `Rc<RefCell<arena>>` 共享
  引用 + `u32` 索引。运算符重载（`+`/`*`/`pow`）经共享引用在自身上下文里构造，
  `ctx.expand(&e)` 这类 API 无借用纠缠；arena 随最后一个句柄释放（无需 reset，
  stale 句柄在 M1 不存在）。代价：跨 Context 混用句柄属未定义行为契约（debug
  断言随需要引入），RefCell 使构造期不可并发——P4 并行方案落地时一并处理。
- **arena 只增不减**：节点在 Context 存续期内不回收。对策：Context 按问题作用域创建；
  P4 再评估压缩/epoch 回收，P0 明确不解决、写入文档。
- **构造开销**：每次建节点过一次哈希表。换来的是消除重复与 O(1) 相等；
  大批量构造路径（poly→expr 重建）提供批量接口绕过逐节点 intern。
- **递归深度**：遍历 API（`post_order` / `children`）用显式栈的迭代器实现；
  M1 的打印与比较暂用带 10 000 层上限的递归（proptest 深度内无碍），
  迭代化为硬化票——深链是本仓库负载的常态。

节点草图（示意，字段允许实现期微调）：

```rust
pub enum Node {
    Int(IntVal),                 // i64 内联；大数存独立 slab 的索引
    Rat(RatVal),                 // (i64,i64) 内联；大数索引；分母恒正、既约
    Float { bits: u64, prec: u16 },   // 仅字面量；-0.0 归一为 0.0；禁止 NaN
    Sym(u32),                    // 符号表索引；符号按名唯一
    Fn  { head: u32, args: Span },    // sin/exp/log/用户函数；head 按名 intern
    Pow { base: u32, exp: u32 },
    Mul { args: Span },          // n 元、已扁平化、按全序升序、数值系数居首
    Add { args: Span },          // 同上
}
// Span = (offset: u32, len: u16) 指向共享参数槽 Vec<u32>
```

**n 元 Add/Mul 是规范形的不变量核心**：`(a+b)+c`、`a+(b+c)`、`a+c+b` 构造出同一个
索引。扁平化、排序、同类合并、数值折叠全部发生在构造器里，任何时刻取出的
`Expr` 都已是规范形——这是 D4"构造即规范"层的基础，也是 hash-consing 成立的前提。

### D2 系数域：泛型 Domain trait + 可换数字后端

| 选项 | 结论 |
|---|---|
| A. enum 闭集（Int/Rat/… 变体） | 否。加 𝔽_p、扩域要改 enum 并波及全部 match；类型系统帮不上忙 |
| B. **泛型 trait（推荐）** | `Poly<D: Domain>` 单态化零开销；expr 层节点内嵌 ℤ/ℚ 快路径 |
| C. `dyn` 动态分发 | 否。热路径多一次虚调用，且丧失由类型区分算法可用性的能力 |

trait 阶层（放 `domain` crate）：

```rust
pub trait Semiring          { type Elem: Clone; /* zero one add mul is_zero */ }
pub trait Ring: Semiring    { /* neg */ }
pub trait EuclideanDomain: Ring { /* div_rem gcd norm */ }   // ℤ 有；域上退化为带余除
pub trait Field: Ring       { /* inv */ }                    // ℚ 𝔽_p 有
```

实现排期：

- **ℤ**：i64 内联 + 大数 fallback 的混合表示（numerica 的"自动升降精度整数"正是这个形态）。
- **ℚ**：既约有理数，分母恒正；与 ℤ 的互相转换显式。
- **𝔽_p**（P2）：u64 + Montgomery 乘法，奇素数 p < 2^63。F4 / 模 p Gröbner 的前提。
- **ℚ(α) 代数扩域**（P2+）：模极小多项式的余式表示。根式、单位根化简需要。
- **Float 字面量不进 Ring**：浮点运算不可结合，与"构造即折叠"的规范形不变量冲突。
  浮点只存在于字面量与求值出口（codegen / `eval_float`）；转精确域必须显式
  `rationalize`。这条边界是刻意画的，防止 sympy 那类 Float/Exact 混算的语义泥潭。

数字后端选择规则：**默认后端必须 MIT/Apache-2.0 兼容；LGPL 后端只做 opt-in feature**。

| 后端 | 许可 | 角色 |
|---|---|---|
| numerica | MIT | 默认后端首选（M1 spike 验证 API 稳定性与性能；自动升降精度、𝔽_p、双数都在里面） |
| num-bigint + num-rational | MIT/Apache | 保守备选：成熟但慢，spike 不通过时回退 |
| rug（GMP） | LGPL | `feature = "gmp"` 性能档，文档明示静态链接合规注意 |
| malachite | LGPL-3.0-only | 不直接依赖；可经 numerica 的后端 feature 透传给用户 |

M2 末设基准决策点：若系数算术在 profile 中占比 > 30% 且 rug 档显著加速，
重新评估默认档（仍受许可规则约束）。

### D3 多项式内核：有序稀疏分布式 + 基准门控的密度切换

**表示**：`Poly<D> = Arc<PolyRing> + Vec<Term>`，`Term = ExpVec + D::Elem`。
指数向量 `SmallVec<[u32; 4]>`；项按 ring 的单项式序严格升序——有序是除法、
Gröbner、`terms()` 稳定性的共同前提。ℚ 多项式内部保持本原 ℤ 表示 + 容度
（去分母的标准技巧，系数域小一半、gcd 稳定）。

`PolyRing` 的身份 =（变元表，项序）。变元表缺省按符号全序（D6）排布，可显式覆盖；
项序缺省 degrevlex（Gröbner 友好），lex / grlex / neglex 可选。ring 身份参与
serde 序列化——多项式的字节形态由 ring 完全决定。

**乘法与除法**（P0 只做稀疏路径，先正确后快）：

- 稀疏 × 稀疏：Johnson 堆合并（Maple 系 Monagan–Pearce 路线），按目标项序惰性产出。
- 单变元大次数：Karatsuba（P1）。
- **稠密切换与 Kronecker 代入留作基准门控决策点**（P1/P2）：阈值 η（项密度、变元数、
  次数的函数）由基准标定，不写拍脑袋常数。切换逻辑藏在 `mul` 内部，接口不变。
- 多除子归约：首单项式匹配的堆归约，Gröbner 归约复用同一原语。

**GCD 与有理函数规范化**：子结果式 PRS + 内容驱动的本原化；ℚ 上先转 ℤ。
`cancel`（分子分母去公因子）建立在其上。

**因式分解**（P0 最重的一项，拆里程碑）：

- 单变元 ℚ[x]：平方自由分解 → Cantor–Zassenhaus 模 p 分解 → Hensel 提升 → 因子重组。
- 多变元 ℤ[x₁..xₙ]：本原部分 + 主变元递归 + 多变元 Hensel。
- van Hoeij 重组（因子数大时）延后 P2；不可分解不报错，返回原式。
- GPL 的 Algebraeon 是现成参考实现，对拍时当第二 oracle（它自己用 malachite）。

**Gröbner（P2）**：Buchberger + Gebauer–Möller 消冗准则 + sugar 度；单项式哈希表归约；
F4 矩阵化作为演进路线（接口上把"符号消元"写成对归约原语的调用，不在 P2 承诺 F4/F5）；
模 p 计算 + 有理重构用于中间膨胀控制。

### D4 化简：分层流水线，确定性优先于"最简"

```
L0  构造即规范（永远在线，免费）
    扁平化、全序排序、数值折叠、同类合并 —— D1 的构造器不变量。
L1  显式定向变换（用户指挥）
    expand / factor / cancel / together / trigsimp(子集) / series。
    每个都是纯函数式的单入口，行为可独立对拍。
L2  simplify()（按需）
    有序规则表自底向上应用至不动点，燃料上限（缺省 8 轮防振荡）；
    代价函数 c(e) = Σ 节点加权和，贪心接受；同代价时按规则表次序定先。
```

- **规则表有序即确定性承诺**：规则的应用次序、代价函数、平局裁决全部固定并文档化；
  版本间变更规则集视为行为变更，写入 CHANGELOG。这是 L2 与 sympy `simplify`
  的本质差异——我们要同一输入恒得同一输出。
- **规则条件消费假设查询（D5），仅 True 触发**，Unknown 一律不动。
  例：`log(a)+log(b) → log(a·b)` 仅当 `a>0 ∧ b>0` 查询为 True。
- **规则实现 P0/P1 手写每个 head 的匹配器**，不引入通用模式匹配引擎——
  避免搜索不确定性进入核心路径，也避免 egg 类系统的重写语义渗入公共 API。
  P3+ 若规则数膨胀（>100），再评估声明式宏或把 egraph 作为**离线工具**
  （离线搜出来的规则人工固化进规则表），不进库核心。
- 对拍口径见 §8："语义一致"是硬指标，"谁更短"只做分类统计，不做门槛。

### D5 假设系统：闭谓词集 + 闭包表 + 三值查询

- **谓词集 P0 固定 9 个**：real / rational / integer / even / odd / positive /
  negative / nonzero / finite。位集表示。
- **闭包表**在符号构造时补全并检冲突：`even ⇒ integer ⇒ rational ⇒ real`；
  `positive ⇒ real ∧ nonzero ∧ ¬negative`；`positive ∧ negative` 是构造错误。
- **绑定方式：Context 内同名符号一次性设定**，重设不同值报错。取舍说明：
  sympy 的 `Symbol('x', positive=True) ≠ Symbol('x')`（同名异设双对象）是匹配、
  缓存与用户心智的长期痛点；Context 按问题作用域创建后，单绑定语义足够且对拍更清晰。
- **查询**：`query(e, P) -> True | False | Unknown`，自底向上按头函数事实表传播，
  遇 Unknown 截断。事实表示例：`exp(x)` 给 real ⇒ positive ∧ nonzero；
  `x^2` 给 x real ⇒ nonnegative；`|x|` ⇒ nonnegative；`sqrt(x)` 给 x real ∧
  nonnegative ⇒ real ∧ nonnegative。
- **与规则的耦合只有一条路**：规则 guard 调查询，仅 True 触发（D4）。
  可靠性优先于激进——宁可少化简，不可错化简。
- **不做** SAT / 一阶逻辑引擎（sympy 新假设系统的复杂度来源，与确定性和性能目标冲突）；
  用户自定义谓词 P3+ 再议，trait 留扩展位。

### D6 规范化与确定性全序（公共 API 承诺）

两层全序，分开承诺：

**（1）表达式全序**——Add/Mul 参数排序、打印次序、序列化次序的基础
（M1 实现定稿，桶序保证传递性）：

1. 精确数（Int/Rat 统一按有理数值比较）先于其余一切；Float 按值排在精确数后。
2. **符号类桶**：Sym 与"底为 Sym 的 Pow"同桶，按底名字节序排，同底 Sym 在
   幂前（x < x^2），同为幂再比指数。这使 Mul 因子输出形如 `x^2*y`、`x*y^2`
   （sympy 习惯，利于对拍比对）。**符号按名字节序而非创建序号**——跨进程、
   跨构造路径确定性的关键。注意：若按类别秩把 Sym 与 Pow 分桶再"按底比较"，
   会构造出 `a^2 < x < sin(a) < a^2` 的环（M1 实现踩过，排序静默失效），
   同桶线性化是传递性的必要条件。
3. 之后依次：符号类桶 < Fn（head 名，参数字典序）< 复合底 Pow（底，指数）
   < Mul/Add（参数字典序，平局比参数个数）。
4. 数值在 Add/Mul 中天然居首（秩最低），即 `3*x*y`、`3 + x + y` 的形态。

**（2）单项式序**——多项式内部项序，由 `PolyRing` 显式携带（D3），与（1）独立。

**承诺机制**（三层防线，针对"行序失同步"这一事故类）：

- **版本化**：`ORDER_VERSION = 1` 写入文档与 serde 头。破坏全序 = 破坏性变更，
  须升主版本（或按 semver 约定的次版本 + 显式迁移说明），CHANGELOG 列出受影响的
  输出形态。
- **golden 快照进 CI**：同一表达式经不同构造路径（proptest 生成两条路径）→
  序列化字节必须相等；跨平台矩阵（linux/mac/windows）对固定语料的打印字节必须相等。
- **三方同序契约**：`poly.terms()` 迭代次序 = serde 序列化次序 = 打印次序 =
  codegen 发射次序。codegen 产物头部注释携带输入序列化的内容哈希，
  数值侧据此对账——符号侧任何重排都能被显式发现，而不是静默错位。
- **实现纪律**：全部哈希用固定键混合器，禁 `RandomState`；任何哈希表的迭代
  在影响输出前必须经全序排序消费。`-0.0` 在构造时归一为 `0.0`；NaN 禁止入
  Float 字面量（构造断言）。

## 5. 分期与验收基准

硬门槛（不达标不算完成）与跟踪指标（记录、开票、不阻塞）分开标注。
性能数字以基准参考机校准后写入 `benches/BASELINE.md` 固化。

### P0：ℤ/ℚ 多项式全链（M1–M5）

| 里程碑 | 内容 | 验收 |
|---|---|---|
| M1 | workspace + expr 规范形 + print/parse | print∘parse == id（proptest 1000 例）**硬**；构造路径无关的序列化一致（500 例）**硬**；docs 首例编译 **硬** |
| M2 | poly ℤ/ℚ 内核 + expand/subst/eval + 差分 harness 上线 | 5000 例随机 expand 语义一致（sympy 交叉求值，§8）**硬**；expand (1+x+y+z)^20 单核 < 1 s **跟踪** |
| M3 | gcd / cancel / 因式分解 API 骨架 | 随机互质与含公因子用例 gcd 正确 **硬**；cancel 语义对拍 **硬** |
| M4 | 单变元因式分解 | x^n−1 全族（n ≤ 64）+ 1000 例随机积：因子重展开与原式恒等 **硬**；单变元 factor ≥ 5× sympy **跟踪** |
| M5 | 多变元因式分解 + 基准报告 v1 | 500 例随机积同 M4 判据 **硬**；expand ≥ 10× sympy、多变元 factor ≥ 2× sympy **跟踪**；未达标项开改进票 |

### P1：初等函数层

- 求导：`diff` 对拍 sympy 2000 例 **硬**；嵌套深度 20 的链式求导正确 **硬**。
- 级数：Taylor 展开至 O(x^n)，对拍 **硬**；`Derivative`、`O()` 未求值节点入树。
- 三角/指对数规则集（首批约 50 条）：语义一致率报告（分类统计，不设"更简"门槛）**跟踪**。
- serde 往返 100% **硬**；打印-解析往返在新节点类型上扩展 **硬**。

### P2：Gröbner 与线代

- Gröbner：生成子对基自归约 100% **硬**；随机理想成员判定（成员化简为零）**硬**；
  cyclic-6、katsura-8 可解 **硬**；cyclic-7 耗时报告、katsura-8 ≥ 10× sympy **跟踪**。
- 𝔽_p：域运算与模 p Gröbner，与 Singular（容器，可选档）对拍 **跟踪**。
- 符号线代：随机矩阵 det/rank 与 sympy 一致（≤ 8×8，ℚ）**硬**。

### P3：代码生成（实际出口）

- 生成的 Rust/C 源在 CI 中真实编译并运行 **硬**。
- 数值一致性：f64 随机点 ≤ 4 ulp（禁 fast-math，对照库内 `eval_float`）**硬**。
- CSE：泊松括号负载上共享节点发射数 −30% **跟踪**。
- 变体：f32/f64 目标；双数后端（numerica 已有双数类型，衔接自动微分出口）**跟踪**。

### P4：并行与增量

- rayon：批量 expand 与 Gröbner 对集合归约 4 核 ≥ 2.5× **跟踪**。
- arena 压缩设计验证（回收 > 70% 死节点）**跟踪**；
  差分缓存（serde 存档 + 输入哈希命中跳算）在正规化流水线试点 **跟踪**。

## 6. P0 API 草图

签名允许实现期微调；语义承诺（全序稳定、构造即规范、错误路径）不随微调变。

```rust
// ── 构造与运算（expr）────────────────────────────────────────
use cas::prelude::*;

let mut ctx = Context::new();
let x = ctx.sym("x");
let y = ctx.sym_with("y", &[Real, Positive]);   // 假设在建符号时一次性设定

let e = (x.clone() + y.clone()).pow(3) - x.pow(3) - y.pow(3);
// 运算符重载构造规范形：扁平化、排序、数值折叠即刻发生

// ── 变换（L1，纯函数式单入口）──────────────────────────────
let g = ctx.expand(&e);            // 3*x*y^2 + 3*x^2*y（全序见 D6）
let f = ctx.factor(&e)?;           // 3*x*y*(x + y)
let c = ctx.cancel(&(g / (x.clone() * y.clone())));   // (x + y) * 3 → 规范有理函数
let h = ctx.subst(&g, &[(y, x.clone())]);

// ── 求值（精确）────────────────────────────────────────────
let v = ctx.eval_rational(&h, &[(x, Rational::from(7))])?;   // ℚ 上精确值

// ── 假设查询（三值）────────────────────────────────────────
use cas::Predicate::*;
assert_eq!(ctx.query(&y, Positive), True);
assert_eq!(ctx.query(&(x.clone() - x), NonZero), False);

// ── 多项式直接使用（绕过表达式层）──────────────────────────
use cas_poly::{PolyRing, Poly, order::DegRevLex};
let ring = PolyRing::new([x, y], DegRevLex);          // ring 身份 =（变元表, 项序）
let p = Poly::from_expr(&ctx, &ring, &e)?;
let q = &p * &p;
let (quot, rem) = p.div_rem(&[&q]);
for (exp, coeff) in q.terms() { /* 稳定次序：= serde = 打印 = codegen */ }

// ── 解析与输出（parse / print）─────────────────────────────
let e2 = cas::parse(&mut ctx, "3*x*y^2 + 3*x^2*y")?;
let tex = cas::latex(&ctx, &g);
let txt = cas::plain(&ctx, &g);     // plain 输出可被 parse 读回（往返承诺）

// ── sym! 宏：批量建符号 ────────────────────────────────────
let (a, b, c_) = sym!(&mut ctx, a, b, c);
```

错误处理：`CasError` 单一枚举（Parse、NotPolynomial、AssumptionConflict、
UnsupportedDomain…），无 panic 出现在公共变换入口（越界句柄属编程错误，除外）。

## 7. 基准计划

**套件**（criterion，全部种子固定、语料入库）：

| 套件 | 用例 |
|---|---|
| bench-expand | (1+x+y+z)^16/20/24；随机稀疏积（3 变元、deg 12、单因子 30/60/120 项）；**泊松括号**（6 变元、截断 8 阶——qiao 正规化的真实负载） |
| bench-factor | 单变元 deg 32/64/128 随机 2–3 因子积；x^n−1（n = 30..64）；3 变元 deg 8/12 |
| bench-groebner（P2） | cyclic-4..7、katsura-6..10、eco-10，grevlex |
| bench-simplify（P1） | 100 例固定种子三角/指对数混合 |
| bench-diff（P1） | 嵌套深度 5/10/20 链式 |
| bench-codegen（P3） | 泊松括号负载：生成 + 编译 + 运行端到端 |

**对照系**：sympy（版本锁定进 `oracle/requirements.txt`，CI 容器内跑，永远在档）；
Singular（容器，可选档）；Symbolica（许可限手动档，不进 CI）。

**CI 策略**：正确性差分每 PR 跑（子集，缓存 venv）；性能基准每周定时专用 runner，
阈值 ±20% 报警不阻塞；`cargo bench` 本地随时。内存用 dhat 针对性测峰值与分配次数，
`/usr/bin/time -v` 做粗粒度交叉验证。

**决策点**（基准驱动，不拍脑袋）：M2 末决定默认数字后端（D2）；
P1 末决定稠密切换阈值 η（D3）；P3 末评估 van Hoeij 与 F4 的优先级。

## 8. 正确性策略

**差分测试（sympy oracle）**：

- 语料：种子化语法生成器（Rust 侧，`xtask oracle` 产 JSONL：表达式 + 指定操作）。
- python 侧计算 sympy 结果，回写 JSONL。
- **一致性判据以语义为主、字符串为辅**：随机 ℚ 点（小整数格点）交叉精确求值，
  两侧相等才算一致；字符串归一化相等仅用于分类。
- 分类报告：agree / cosmetic（仅形态差）/ ours-longer / theirs-longer /
  **semantic-mismatch（= bug，硬门槛零容忍）**。
- 一致率与分类进 CI 产物归档。

**属性测试（proptest）**：

- print∘parse == id；构造路径无关（同一表达式两条构造路径序列化字节相等）。
- expand 幂等；`expand(a)*expand(b)` 展开后 == `expand(a*b)`。
- factor → 重展开 == 原式（往返恒等）；gcd Bézout 校验（ℤ[x] 上组合出 gcd）。
- 代换后求值 == 求值后代换（ℚ 点精确）；除法不变量 `p = q·d + r` 且 r 次数规约。
- ℚ 构造既约性；假设闭包无冲突泄漏。

**golden 快照**：D6 的跨平台、跨路径字节级钉死，`xtask snapshot` 显式刷新并走 review。

## 9. 工程约束

- **目录**：`cas/` 子工作区（留在 qiao 仓；自包含，可整体迁出为独立仓——开放问题 #4）。

```
cas/
├── Cargo.toml            # workspace
├── crates/{domain,expr,poly,parse,print,codegen,cas}/
├── xtask/                # oracle / snapshot / bench / mem 任务入口
├── oracle/               # python：sympy 对拍脚本 + requirements 锁定
├── benches/              # criterion 套件 + BASELINE.md
└── DESIGN.md
```

- **MSRV**：1.85，edition 2024；CI 三轨（msrv / stable / nightly 的 clippy+fmt）；
  提升 MSRV 视为破坏性变更，提前一个 minor 公告。
- **no_std**：`domain/poly/print` 接口与数据结构按 no_std+alloc 设计（hashbrown）；
  `expr` 视数字后端的 std 依赖决定；`parse/codegen/cas` 仅 std。
  **尽力目标而非硬承诺**——后端限制如实写进各 crate 文档，CI 不设 no_std 门槛。
- **依赖清单**（全量）：num-bigint（M1 已用）+ num-integer（M1 已用，大数 gcd）
  + numerica（M2 决策点）+ hashbrown（no_std 时）+ serde（P1）
  + rayon（P4）；dev：proptest、criterion、dhat。 rug 走 feature。除此之外不加。
- **CI 矩阵**：linux 必跑；mac/windows 一档（目的只有一个：抓平台相关的排序/哈希泄漏）。
- **文档**：docs.rs 首页三例（展开 / 求导 / 生成代码），P0 期先落前两个占位；
  每个 crate 的稳定性承诺（全序、terms() 次序、规则表次序）在类型级 doc 里重申。
- **许可**：推荐 `MIT OR Apache-2.0` 双许可（开放问题 #1）；CONTRIBUTING 与
  基准复现说明随 M1 一并入库。

## 10. 风险登记

| 风险 | 影响 | 缓解 |
|---|---|---|
| 多变元因式分解实现复杂（Hensel 边界条件多） | P0 最可能的延期项 | 拆 M4/M5 两级；Algebraeon（GPL，只读不链）作参考与第二 oracle；不可分解返回原式不算失败 |
| arena 只增不减 | 长会话内存增长 | Context 作用域化 + 文档明示；P4 压缩 |
| 无 GMP 默认档，系数密集负载可能慢 | 性能目标 | D2 决策点 + rug feature；numerica 本身有 rug 后端可透传 |
| 追平 sympy 长尾的功能蔓延 | 范围失控 | §1 非目标清单 + 差分报告"theirs-longer"只统计不立项 |
| numerica 年轻（2026-07 v2.2），API 可能漂移 | 后端迁移成本 | Domain trait 隔离 + num-bigint 备选档；锁定版本 |
| 单人维护的 bus factor | 项目连续性 | 决策全部 ADR 化；接口窄、文档齐，降低接手成本 |

## 11. 开放问题（评审时拍板，不阻塞骨架启动）

1. **许可证**：推荐 `MIT OR Apache-2.0`；若你倾向纯 MIT 也可，只影响依赖合规表。
2. **默认数字后端**：推荐 numerica（MIT、自动升降精度、与 Symbolica 同源），
   M1 spike 后确认；保守备选 num-bigint。
3. **命名**：目录暂用 `cas/`，crate 前缀 `cas-*`。发布名（crates.io 占位）待定，
   机械重命名不影响设计。
4. **仓库位置**：暂留 qiao 仓子工作区（自包含可迁出）；若确定长期独立发布，
   M1 前直接建独立仓更省事。

## 12. 术语表

- **规范形（canonical form）**：构造器保证的不变量形态——n 元 Add/Mul 已扁平化、
  按全序排序、数值已折叠、同类已合并。规范形即身份（hash-consing）。
- **hash-consing / intern**：构造时查表去重，结构相同的子表达式只存一个节点，
  索引相等即结构相等。
- **arena**：节点只追加不修改的存储区；`Expr` 是其中的 32 位索引句柄。
- **表达式全序**：D6 (1)，Add/Mul 排序与打印次序的基础；版本化承诺。
- **单项式序 / 项序**：多项式内部项的序（degrevlex/lex/…），由 PolyRing 携带，
  与表达式全序独立。
- **首单项式（leading monomial）**：项序下最大的单项式；除法与 Gröbner 的锚点。
- **容度 / 本原部分（content / primitive part）**：多项式系数的 gcd 与去除它后的
  本原多项式；ℚ 计算转 ℤ 的桥梁。
- **PRS（多项式余式序列）**：Euclid 辗转的多项式版；子结果式变体防系数爆炸。GCD 的基础。
- **Hensel 提升**：模 p 因子分解提升回 ℤ 的算法；因式分解的核心步骤。
- **Gröbner 基**：理想的一组生成元，使首单项式归约唯一；多项式系统消元的基础。
- **假设（assumption）/ 谓词**：符号携带的性质（real/positive/…），闭包表补全，
  三值查询（True/False/Unknown）。
- **差分测试 / oracle**：以 sympy 为参照系的批量对拍；语义一致 = 随机 ℚ 点交叉求值相等。
- **CSE（公共子表达式消除）**：代码生成时把重复子式提为临时变量；
  hash-consing 的共享集是它的直接输入。
- **golden 快照**：字节级钉死输出的确定性测试；显式刷新并走 review。

## 13. M1 落地记录（2026-09-11）

范围：`cas/` workspace 五 crate（domain/expr/print/parse/cas）全部落地并全绿；
poly 与 codegen 目录随 M2 / P3 创建。M1 验收全部达成：

- print∘parse 往返 1000 例（含跨上下文字节一致）**通过**；构造路径无关 500 例
  （同上下文节点恒等 + 跨上下文结构相等）**通过**；docs 首例 doctest **通过**；
  clippy `-D warnings` 与 rustfmt 干净；黄金快照钉入 `cas/tests/roundtrip.rs`。

实现期确定的细则（均为 D1/D2/D6 的具体化，不改变决策）：

- **句柄形态**：`Expr = Rc<RefCell<arena>> + u32`（见 D1 修订）。Float 节点
  M1 无 prec 字段，P1 多精度时加；**负浮点入口即拆为 `-1*|v|`**，arena 内
  Float 恒非负、`-0.0` 归一 `0.0`。
- **幂折叠规模守卫**：`MAX_FOLD_BITS = 2^16`（约 8 KB 常数）。除挡字面量炸弹
  外，还截断"折叠结果的再折叠"链——((p/q)^-253)^-256 类链条会滚出几十万比特
  常数，打印/解析/归约全为平方级负担（perf_probe 实测单例 14 s）。
- **数值性能两条硬教训**（属性测试随机语料逼出来的）：
  1. 大数 gcd 必须用 num-integer 的实现；手写欧几里得在几千比特操作数上
     是平方级灾难。
  2. 既约有理数的幂与逆元仍既约（gcd(num,den)^m = 1），构造时可跳过 gcd——
     这是快速路径，但**落位判定（Small/Big）必须与通用路径逐字节一致**，
     否则同值两表示，往返测试当场失败（29^13 案例：分母在 (i64, u64] 区间）。
- **全序的传递性陷阱**：Sym 与 Pow 分桶后"按底比较"会构造出
  `a^2 < x < sin(a) < a^2` 的环，`sort_by` 静默失效。符号类必须同桶线性化
  （D6 已按此修订）。
- 数字后端按 D2 保守档（num-bigint）落地；numerica 评测在 M2 多项式内核
  动工前完成（依赖 num-integer 已新增，见 §9 依赖清单）。

**真实语料测评（第 0 层，2026-09-11）**：`xtask corpus` 以 qiao 正规化流水线
的 `hpoly_generated.rs`（源 Hamilton_poly15.mat，828 项 × 36 参数符号系数，
只读）为语料——表达式族与 plain 文法逐字重合（`.powi(n)`→`^n` 机械转写）。
结果：解析 828/828、往返同节点 828/828、跨上下文字节一致 828/828；
规范形 1.50 MB（原文的 83.5%，数值折叠 + 合并同类项的效果）；arena 唯一
节点 88 284 个（跨表达式 hash-consing 共享）；确定性摘要 FNV64 跨运行恒定
（393538e72ab33202）。吞吐：release 下解析 559 条/s（均值 2.2 KB/条）、
打印 6.9 万条/s——解析侧是 M2 的首要 profile 对象（往返再解析反而更慢，
疑与 intern 表增长有关，待基准定位）。泊松括号/Lie 变换的对拍按第 1 层
计划在 M2 落地。

已知限制（开票跟踪）：打印/比较递归 + 万层深度上限；构造期 RefCell 独占
（单线程）；无 serde（P1）；无 `Context::reset`（句柄与 arena 同生命周期）。

## 14. M2 落地记录（2026-09-11）

交付：`cas-domain` 增 Ring/Field trait 层；新 crate `cas-poly`（grevlex/lex
有序稀疏多项式：加/乘/幂/偏导/求值，构造即规范化）；`cas-expr` 增
expand（DAG 分配 + 规模守卫）/subst/eval_rational/eval_float；`xtask` 增
poisson/oracle/bench 子命令与 corpus --eval 档；`oracle/sympy_expand.py`。

验收（对照 §5 M2 行）：

- **硬门槛「5000 例随机 expand 语义一致」通过**：sympy 交叉精确求值
  15000 点值全部相等，0 不一致、0 跳过；我方 expand+求值 25.3 万条/s，
  sympy 1029 条/s（同口径，含其 expand+subs）。
- **跟踪指标 expand (1+x+y+z+w)^20 < 1 s 达标**：10626 项（二项式校验✓）
  830 ms；sympy 同式 971 ms——当前持平而非 10×，原因明确：M2 的 expand
  走 DAG 分配（朴素算法），poly 内核承接的快速 expand 属 M3 桥接。
- 泊松括号（第 1 层负载，q1..q3 p1..p3 布局与 poly_poisson 一致）：
  反对称/双线性/Jacobi ×10 组随机多项式**精确**通过；与中心差分交叉验证
  最大相对差 1.17e-9；deg6×94 项 与 deg5×50 项 → 2256 项，1.6 ms。
- 真实语料数值对拍（corpus --eval）：828 条 × 5 点 × 双通道
  （parse→eval 与 parse→print→parse→eval，独立 f64 求值器为基准）
  0 违例，最大相对差 2.4e-11。

实现偏差（均为 M3 排期调整，不改变决策）：expand 的 poly 承接版、
gcd/cancel、numerica 后端决策点统一挪至 M3——M2 的系数负载未触
num-bigint 瓶颈，且 Ring trait 已隔离后端切换。

工程坑记录（对拍 harness 类）：sympy `parse_expr` 默认 transformations
**不含** convert_xor（`^` 被解析为异或）——oracle 协议显式替换 `^`→`**`；
子进程管道协议在请求量超过管道缓冲时必须独立线程写 stdin（5000 例
请求约 2.5 MB，同步写与响应读取互锁）。

## 15. M3 进展（2026-09-11）：expand 多项式快路径

expr↔poly 桥（`poly_bridge`）：纯多项式子树（Int/Rat/Sym/非负整幂/加乘）
整体提取进 `Poly<Rational>`（变元序 = 名字节序），poly 域完成分配/求幂
（项数上界 C(t+k-1, t-1) 饱和估计守卫，上限 100 万项），重建经规范形
构造器。含浮点/函数/负幂的子树自动让位给 DAG 分配。

同构性由测试钉死：快慢两条路径对随机多项式与混合表达式产出**同一
arena 节点**（60 组 + 混合用例）。

基准变化（release）：(1+x+y+z+w)^20 830 ms → **98 ms**（sympy 0.99 s，
**10.1×**，达成 §7 跟踪指标）；^16 243 ms → 42 ms。M3 剩余：gcd/cancel
（PRS + 本原 ℚ 表示）、numerica 后端决策点。
