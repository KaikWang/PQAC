# 完整 q88 CRT v3：维度、算子与分配优化

2026-10-13；服务器 autodl-pqa11；独立 profile；基线源码未改。

## 五次完整新 show

| 次数 | 完整生成 s | wire 验证 s | 累加器 worker 生成 s | 累加器 worker 验证 s | envelope KiB |
|---|---:|---:|---:|---:|---:|
| 1 | 3.378518 | 1.665 | 3.305721 | 1.645463 | 188.68 |
| 2 | 3.303526 | 1.660 | 3.235388 | 1.648863 | 188.68 |
| 3 | 2.882511 | 1.627 | 2.810483 | 1.615944 | 188.71 |
| 4 | 3.110630 | 1.676 | 3.035525 | 1.665337 | 188.71 |
| 5 | 3.519819 | 1.715 | 3.450743 | 1.697814 | 188.68 |

| 项目 | 原生产基线五次平均 | v3 五次平均 |
|---|---:|---:|
| 完整 show 生成 | 5.312709 s | 3.239001 s |
| 完整 wire 验证 | 2.821800 s | 1.668600 s |
| 累加器 worker 生成 | 5.243889 s | 3.167572 s |
| 累加器 worker 验证 | 2.810638 s | 1.654684 s |
| accumulator envelope | 194.37 KiB | 188.692 KiB |

完整生成平均降低 39.03%，验证降低 40.87%。
最后一份实际 envelope 为 193209 bytes。
完整 packet 212.29–212.39 KiB；包含匿名签名、公开实例及其他包装。
计时排除首次 UAV/GCS 预热 5.829/2.005 s。
数据为不同随机新 show；不能把固定 fixture 的最小值当成全流程保证。
完整 show 和 accumulator worker 都未在五次测试中全部达到 3 s / 1.5 s。

## 固定 fixture 与选择过程

- 较严格的 v3 初始布局约 3.77/2.14 s，未保留作为选定配置。
- 将 auxiliary radix 与原始 accumulator radix 解耦，保持 auxiliary requirement，减少近似分组。
- bounded unbalanced top + scaled Z binding：first Pack nn=31164，second nn=22169。
- 临时分配缓存 1024 MiB：三次生成 3.081/3.409/3.093 s，验证 1.690/1.706/1.575 s。
- 六 prime 分量点积：三次生成 2.905/2.677/2.799 s；验证 1.553/1.510/1.591 s。
- 三线程点积的三次验证 1.752/1.660/1.609 s，未优于六线程，选定六线程。
- LDD Add-check 4 线程验证 1.930/1.879/1.943 s，明显更慢，保留原单线程。
- group32 参数门与 Clear Verify 通过，但所有 Pack nn/nmax 与 group16 一致，无维度收益；恢复 group16。
- 完整新 show 的结果比固定 fixture 差，见上表；没有把最佳三次结果替代五次新 show。

pool_online5 是早期分配缓存测量，期间短暂重叠独立交叉验证工作，不用于最终性能对比。
sprod_online5 没有其他编译或验证实验并发。

## 检查范围

- 完整主 relation、quotient/carry、CRT lift、hidden-handle/refresh/product 和 q19 binding 保留。
- 原 native 参数生成安全门、公开加权 beta 界检查通过；log2 bound=73.453438 < 80。
- Clear Verify、Pack Verify、serialize→deserialize→verify 通过。
- 独立新进程原 liblabrador38.so 成功验证新引擎生成的 proof，见 q88_v3_crossbackend.log。
- 同一 public profile 的 negative checks：shared handle 翻转、非二进制 handle、非零高次 handle、
  scaled Z 绑定篡改、动态 accumulator RHS 和 q19 RHS 篡改、截断与尾随 proof 均拒绝；恢复后验证通过。
- mixed radix 边界与随机值 101611 个精确 reconstruction 检查通过；公开 top bound=9624。
- batch/fused audit 对旧原生算子逐元素对照通过（64/16 次），见 kernel_audit2.log。

没有独立完整密码学安全认证。新范数布局属于实验，生产参数未替换；
原分组 v2 对辅助 requirement 的说明错误已在其 REPORT/README 顶部更正，v3 显式排除这些辅助分组。
缓存只复用 allocations；Gaussian masks、challenges、LNP responses 均不跨重试复用。

## 源码与复现

../dimension_v3/profile.sh 为选定环境；run_online.sh 为完整新 show runner。
../dimension_v3/README.md 说明布局与构建。所有 source 以本地为准，单向同步远端。
私有 fixtures 未同步回本地。保留完整 logs、measurement.json 和实际 public proof/envelope。

下一步仍需减少约 0.15–0.2 s 在线验证耗时，并降低 LNP 重试造成的生成时间波动。
本次未更改 folding threshold、生产参数、重复次数或 Gaussian 分布，也未使用 PGO。
