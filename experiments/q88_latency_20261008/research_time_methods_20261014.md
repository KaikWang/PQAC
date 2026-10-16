# q88 完整 Show 时间优化复盘与下一轮调研

> 历史调研快照：下文的“最新”及“下一轮”均指 2026-10-14 调研时点。
> 后续第 1 项及第 2 项的调用内转换共享已实施；当前状态和剩余任务以
> [优化待办](OPTIMIZATION_BACKLOG.md) 为准，避免重复实施已完成方向。

日期：2026-10-14。研究依据为现有源码、已保存实验日志和下列论文/作者实现；本轮没有新增性能测量，也没有修改生产算法。论文 PDF 被 ePrint 的 403 拒绝，因此论文侧仅确认可读的官方摘要/项目说明；具体可移植细节仍须通过完整论文和作者代码核对，不能声称已经完成安全论证。

## 最新已保存测量

- 初始完整 Show 生成/验证 3.246313 / 1.657460 秒。
- 更新后后台准备 0.810717 秒，随后完整 Show 3.123185 / 1.719312 秒。
- 准备+生成 3.933902 秒；UAV 更新+准备+生成 3.971126 秒。
- 初始完整 packet 平均 212.489 KiB，更新后 213.481 KiB；accumulator envelope 约 188.56 KiB。
- 三次 fresh 生命周期样本，非最坏情况保证，仍未达到 3秒/1.5秒。

来源：results/background_bulk_fixed_20261014/REPORT.md、lifecycle_measurement.json。

## 已经做过的事

保留：full-ring compact CRT 布局、公开 top-digit/加权范数界、auxiliary radix=3779 与 scaled Z binding、原始 witness 范数组合、批量 split/fused phi、六内部 prime 点积并行、限定大小 allocation pool、Gaussian chunk 并行、immutable public bank、resident statement/Pack 参数、动态 RHS、checked exact-lift cache、批量整数多项式折叠、更新后 background prepare。原约束/随机性/参数门保持。

没有稳定收益或未作为默认：group32（nn/nmax 无变化）、LDD addcheck 四线程、fused/linear 八线程、projection products 并行、sparse phi blocks、CPU affinity、pseudo-Mersenne folded center。GPU 小点积规模下无已保留收益。Pack folding 阈值已经试过，不将重复扫阈值列为优先项。

LNP folded relation/相关 retry-invariant 内容已有复用，不把 fresh masks/challenges/z1/z2 跨重试缓存当作优化方法。

## 瓶颈依据（历史 profile，不是当前新测量）

show_opt_20261014/q88_show_detail_current.log：
- Prove：LNP 前 folding 1.770 s，LNP 0.882 s，LNP 后 folding 0.614 s；outer attempts=18。
- Verify：Dachshund reduction 0.278 s，首轮 LDD 0.543 s，第二轮 0.243 s；LNP Verify 0.104 s。
- 首轮 Rq sparse phi 聚合约 0.166 s；其输入 sparse constraints 为 8946。

这些数字来自旧的单次 detailed profile，不能直接加到最新均值或保证相同比例。下一轮应同时记录参数维度、每轮时间、重试次数，确认瓶颈在当前默认配置仍然成立。

## 候选顺序

| 顺序 | 候选 | 作用阶段 | 与已有尝试的区别 | 边界 |
|---|---|---|---|---|
| 1 | 合成 refresh 矩阵移到 CRT lane 循环之外 | 初始准备、更新后的准备/冷 Show | 不再让三个 lane 重复执行同一 SHAKE 展开与矩阵求和 | 输出应逐系数与旧路径相同，lane 仍分别执行精确乘法和整除 |
| 2 | C 端公共多项式转换缓存；同一调用中的 common witness 转换共享 | 同上 | Python public bank 已缓存；C 的 FLINT 输入仍重复转换 | cache 绑定 PP/上下文；更新失效规则明确；private conversion 仅限同一调用或绑定完整 witness |
| 3 | native 公开系数采用 zero/identity/monomial/sparse block 表示 | Prove/Verify 的首轮 folding | 改计算表示而非对原稠密操作增加线程；已有 scalar-block 特例不重复实现 | 只利用公开系数结构，不能按私有 witness 零值删变量；保留原 transcript/约束 |
| 4 | quotient/carry/auxiliary 的等价消元和共享 | Prove/Verify、可能缩小 proof | 真正减少 backend 维度，而非改分组名称或阈值 | 保留原始 witness 及范数约束的逻辑等价；重新证明整体 beta 界和模数 aliasing，独立 profile |
| 5 | 对 LNP 响应计算和拒绝检查作批量融合 | Prove，尤其尾部波动 | 减少每次 retry 内部的重复转换/扫描，不复用 challenge-dependent 响应 | 同一随机流和完整接受事件；不得减少 soundness repetitions、Gaussian 宽度或放松范数门 |
| 6 | 对照作者的新 succinct toolkit 编译器与 batching 实现 | 编译/内核优化 | 逐项移植，核对公开输出而不是整库替换 | OURS 的 ABI/cache/guard 与 native 修改必须保留，原 proof 可交叉验证 |

### 第一项的具体证据

common/core/src/acc96_exact_lift.c：memver_exact_lift_lane 中，六个 refresh block 的完整 CRT 矩阵展开/选择求和位于 lane 函数内；acc96_memver_exact_lift 又调用该函数三次。令 selector 中 1 的个数为 h，矩阵展开调用由 3*6*(1+h) 次改为 6*(1+h) 次。这仅说明该部分展开减少 2/3，不等于总准备时间减少 2/3。合成结果包含三条 lane，可共享只读六个多项式；不会新增 witness/proof 字段，也不更换模数。

### 第三项的具体证据

constraints.c 的 rotation aggregation 已有 coordinate-zero scalar 特例和 lazy scratch。因此新的目标是其他公开稀疏/monomial block，而非再次实现现有 scalar 优化。先统计实际命中数、稠密 expansion 数和数据访问量；只在足够覆盖首轮 8946 条关系时编写内核。可将约束编译为固定稀疏工作表，再保持输出顺序/精确约减；不按 secret handle/selector 分支削减证明关系。

### 第四项须注意

当前编译器注册 refresh_z 和 scaled refresh_z，并通过六个 binding 证明 scaled Z 正确；handle 也有 full-ring 与 q19 shared-handle 表示。不能直接删除它们。必须证明代入后的主关系、二进制/常数约束、范数界及 quotient/carry 重建均等价，并实际减少 nn/nmax。group32 曾经通过 Clear Verify 但 nn/nmax 不变，应避免重复把这种参数变化当作维度收益。

## 相关论文及作者实现

1. LaBRADOR: Compact Proofs for R1CS from Module-SIS：https://eprint.iacr.org/2022/1341
   研究递归 amortization/folding。用于理解 witness/关系结构对每轮工作量的影响，不将非 ZK 的紧凑 proof 指标直接替代完整 ZK Show。
2. A Toolkit for Succinct Lattice-Based Zero Knowledge Proofs：https://eprint.iacr.org/2026/1289
   将 LaBRADOR 与 LNP22 零知识结合，架构与当前方案相关；官方代码/benchmarks：https://github.com/lazer-crypto/lazer 。官方 README 给出对应 commit 与 benchmark_expansion/compression/membership/blind_sign。
   本地 membership benchmark 使用 Merkle path/Falcon key 的知识证明，不是我们的可撤销 lattice accumulator，不能直接拿它的秒数/大小比较。
3. Lether: Practical Post-Quantum Account-Based Private Blockchain Payments：https://eprint.iacr.org/2026/076
   官方摘要明确指出优化 LaZer/LNP22 的多个 proof components batching。移植前区分单个证明内的计算批量化、协议组件合并和多个请求的批量验证；后者提升 throughput 不等于单请求 latency。
4. Compact Lattice Signatures via Iterative Rejection Sampling：https://eprint.iacr.org/2024/2052
   适用于其特定签名构造；不能据此直接减少当前 LNP 重复次数或改变接受事件，作为协议重设计方向而非近期开关。
5. Orbweaver: Succinct Linear Functional Commitments from Lattices：https://eprint.iacr.org/2024/2026
   预处理 linear functional commitments，可做长期参照；有 structured SRS/不同假设，而且需要处理当前 quadratic/range relations，不作为直接替换方案。

另以 /tmp 独立 checkout 比较了 upstream Labrador HEAD c1ae0e6 与本地快照 3525558。两者接口/实现差异很大，其中大量差异是我们已加入的缓存、guard、精确 CRT、序列化与并行修改；未认定 upstream HEAD 对当前关系具有性能优势，未改动生产 checkout 或自动升级。

## 接下来如何判定成功

- 首先获取当前配置的分阶段 profile，区分 native folding、Dachshund、LNP（记录 retry counts）、C prepare 和 serialization。
- 第一项与第二项用 exact output 对照；第三项用所有内部 NTT prime 的逐元素比较、旧 backend 交叉 Verify。
- 表示/维度变化独立 profile：重算 norm/beta/alias bounds，Clear Verify、Pack Prove/Verify、serialize/reload Verify、binding/invalid norm/changed context negatives。
- 最终测初始 Show、每次更新后的首个 Show、warm repeated Show；至少多个 fresh 生命周期，进一步观察 p50/p95/max。记录有/无后台准备的准备+生成总时间。
- 记录真实 serialized accumulator 和完整 Show packet bytes，不用估算大小替代，随机压缩波动单列。
- 如无稳定收益则撤回候选，不放松安全参数、不做 PGO、不要求特定服务器 CPU 编号。
