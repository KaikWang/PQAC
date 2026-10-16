# q88 CRT 有效 Pack 维度优化：group4

日期：2026-10-13。服务器：autodl-pqa11，Intel Xeon Gold 6430。
独立 Python compiler/worker 副本，生产 `OURS` 和 `lazer` 源码保持原状。
使用与基线相同的生产 C 可执行文件和 native proof engine。

## 后续代码核对更正

后续逐项核对发现：生产 compiler 的 `aggregate_idx` 在 compact witness 且未启用
`REVAC_Q96_AGGREGATE_TIGHT_NORM` 时，也设置 `accumulator_witness=True`。
所以本轮 group4 v2 **也调整了这些 aggregate product 数位的 JL requirements**。
下文“product auxiliary 的单独 requirements 不变”表述不准确，以本更正为准。
lift/carry 与 refresh 的 requirements 确实未进入该分组。原始 input norm、关系和实测数据不受此说明更正影响。
`beta0` 范数门已运行通过，但不能据此推断所有辅助见证的既有安全论证完全不受影响。
下一轮 `dimension_v3` 用独立集合明确排除 aggregate auxiliary，并重建、验证参数。

## 完整 show 对比

本次重新生成五份完整匿名凭证，包含匿名签名与完整 CRT 累加器证明。
签名和累加器生成沿用基线并行方式。下表均为预热后在线请求，不含首次构造。

| 项目 | 基线五次 | group4 五次 |
|---|---:|---:|
| 完整 show 生成平均 | 5.312709 s | 3.555605 s |
| 完整 show 生成范围 | 5.156323–5.485187 s | 3.178731–3.958801 s |
| 完整 show wire 验证平均 | 2.821800 s | 1.891000 s |
| 完整 show wire 验证范围 | 2.731–2.870 s | 1.835–1.975 s |
| 累加器 worker 生成平均 | 5.243889 s | 3.482310 s |
| 累加器 worker 验证平均 | 2.810638 s | 1.879319 s |
| 累加器序列化 envelope 平均 | 194.37 KiB | 188.49 KiB |
| 累加器 envelope 范围 | 194.00–194.50 KiB | 188.12–188.77 KiB |
| 压缩完整 show packet 平均 | 214.60 KiB | 213.70 KiB |
| 压缩完整 show packet 范围 | 212.90–219.11 KiB | 212.08–216.28 KiB |
| 最后一份 envelope 实际字节数 | 199170 | 193305 |
| UAV/GCS 预热 | 5.916 / 2.111 s | 5.913 / 2.066 s |

完整 show 生成和验证平均分别降低 33.07%、32.99%。累加器 envelope 降低约 3.03%。
匿名签名拒绝采样会影响完整 packet 大小，应比较范围及均值，不把单次波动归因于累加器。
每次 envelope 的 KiB 来自 wire 实际序列化长度（日志四舍五入），最后一份另保存真实二进制文件。

| 次数 | 完整生成 s | wire 验证 s | 累加器 envelope KiB |
|---|---:|---:|---:|
| 1 | 3.853016 | 1.911 | 188.34 |
| 2 | 3.433692 | 1.867 | 188.12 |
| 3 | 3.178731 | 1.867 | 188.46 |
| 4 | 3.353783 | 1.975 | 188.74 |
| 5 | 3.958801 | 1.835 | 188.77 |

五次在线样本全部满足 4 s / 2 s。额外新进程、反例后恢复的有效证明验证分别出现
2.060079 s、2.100669 s。因此这是已测在线性能结果，尚不能保证任何启动状态或负载下都低于 2 s。

追加 group8 的三次 fixture 完整证明试验：生成 3.354043、3.068888、3.119671 s；
验证 1.933709、1.811946、1.979171 s，序列化 192601、192633、193273 bytes，均验证通过。
其验证平均为 1.908275 s，与此前 group4 fixture 的 1.907877 s 基本一致，
没有增加验证时间余量。因此当前保留拥有五次新完整 show 证据、范数预算更小的 group4。

## 改动及收益来源

参考 Lazer-Acc 的 witness/范数组分方式，并检查 native `dachshund.c` 的近似向量合并条件。
该条件累计输入 norm，并以候选向量中最小的 JL requirement 判断是否允许合并。
原编译器许多独立数位向量的 requirement 过细，产生大量独立 Dachshund 输出组。

候选只为最多四个 accumulator 数位向量设置共同公开 requirement：

`R_group = next_power_of_two(sum(original_requirement_i))`。

这样为 native 合并提供空间。实际 witness 系数与原始输入 norm 不变，原逻辑向量数也不变。
真正减少的是 Dachshund 输出及后续 Pack 折叠维度：

| 参数 | 基线 | group4 |
|---|---:|---:|
| compiler 逻辑 witness vectors | 141 | 141 |
| degree-256 等价 witness polynomials | 16352 | 16352 |
| Dachshund approximate merged groups | 108 | 43 |
| Dachshund 输出 r | 118 | 57 |
| 第一次 Pack 输入向量数 | 120 | 59 |
| 第一次 Pack nn / nmax | 36232 / 4320 | 35492 / 2160 |
| 第二次 Pack nn / nmax | 83450 / 4320 | 25182 / 1750 |

稀疏零项、scratch pool、向量排序、fast CRT 原生实验没有获得同样稳定的整体收益，
本候选未加载这些替代库。

## 约束、范数和参数绑定

- 原始三个 CRT 素数、degree-4096、24-bit hidden handle、128-bit refresh selector 不变。
- hidden handle binary/product、refresh binding、三条 CRT lift、q19 四行 link 与 shared-handle 约束均沿用完整生产 compiler。
- lift/carry、refresh、product auxiliary 的单独范数 requirements 不变。
- accumulator 数位的中间 JL requirements 改变，原始 input norm 不变。
- 重新执行 `_audit_norm_profile`：所有 requirement 均为二次幂且覆盖公开输入 norm 与 JL slack。
- 保守 Cauchy–Schwarz 重组界 `log2(norm) <= 75.059 < 80`，仍通过原有 `beta0 < 2^80` 门。
- native 重新生成 Dachshund/Pack 参数，仍执行库中的 commitment rank/SIS 参数门。
- group4 编译策略绑定进 static profile ID 及每次 show context；生产 group1 verifier 拒绝 group4 envelope。

上述是实现中的范数、参数和功能检查；它不是对整个密码构造的独立形式化安全证明。
当前保留为独立候选，没有替换生产 profile，也没有减少新挑战所需的响应采样。

## 验证证据

1. 捕获的完整生产 witness 三次：Clear Verify、Pack Prove、序列化、反序列化 Pack Verify 全通过。
   实际 envelope：193273、192601、192889 bytes；这是 fixture 测量，区别于上面的五次新 show。
2. 五次新完整 show 全部 Pack Verify 通过，生产实体集成测试通过。
3. 新进程仅使用公开 statement 重建 verifier 参数后，保存的 compact proof 验证通过。
4. 改 accumulator coefficient 被公开 digest 校验拒绝。
5. 改仍为 canonical 的 q19 link RHS 被 native Pack Verify 拒绝（system check 失败）。
6. 截断 envelope、添加尾部数据被反序列化拒绝；恢复有效公开输入后原证明再次验证通过。
7. baseline group1 对 group4 的证明因 profile ID 不一致而拒绝。

## 复现与文件

运行 `../../native_sparse/run_norm_group4.sh`，或 README 中的服务器绝对路径。
`measurement.json` 为五次新 show 及检查的机器可读汇总。
`q88_norm_g4_online.show.log` / `q88_norm_g4_online.worker.log` 为原始计时。
`q88_norm_g4_envelope.bin` 为最后一份新 show 的公开累加器证明。
`q88_norm_acc_g4_proof.*` 为独立 fixture 证明实测。
该次原始 JSON 的 `profile` 字段来自旧 benchmark 模板，误写为 `production CRT unchanged`；
同名原始日志明确记录 group size 4。benchmark 已修复此标签，当前 `measurement.json` 标识准确。
`q88_norm_acc_g4.log` 为详细参数生成审计；g8/g16 文件是其他组大小的参数筛选。
`q88_norm_acc_g8_proof.*` 为追加 group8 的完整 fixture 证明测量。
`q88_norm_g4_independent_verify.log`、`q88_norm_g4_negative_checks.log`、
`q88_norm_g4_baseline_rejection.log` 分别记录独立验证、反例检查和 profile 绑定拒绝。
私有 witness fixture 仅留在服务器，没有归档到仓库。

基线两份 compiler/worker 源文件 SHA256：

```text
136cfc7df30157f55dda2af74e3cf69d45f2a2d2fad6cba4d0b5b91e73f9d7ed  revac_acc96_full_ring.py
cc2d16f1d4fed2b595db04464fb2d8d7796c585325673a4fceeabe3916729d74  revac_acc96_memver_worker.py
```
