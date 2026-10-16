# 完整 Show checked preflight 并行实验：2026-10-14

## 状态

已修改 OURS，但异步 preflight 仅作为显式 opt-in 候选，生产默认继续串行 preflight。保留 q88 三素数 CRT、完整约束、原始范数界、随机采样和协议字段。未使用 PGO、CPU 绑定或 GPU 专属优化。服务器 autodl-pqa11，最终生产 native 算子与上轮相同。

不能宣称本轮稳定达到 4 秒 / 2 秒或 3 秒 / 1.5 秒。初步三次更新后样本均小于 4 秒，但最终复测出现 5.935 秒样本，全部纳入统计，未删除。

## 实现

- `revac_uav_show_prove` 在 online parallel q88 路径将 MemVer cache 检查、整体范数检查和 exact CRT lift 构造委托已有 accumulator proof 线程，与签名证明并行。
- 证明线程完成原有 checked preflight 后才能向 Pack worker 发送请求。任何检查失败，整个 Show 返回失败；线程创建失败和异步任务失败仍经过 checked serial fallback。
- 签名线程读取凭证字段，accumulator proof 线程独占修改 MemVer/lift cache；finish join 后才返回完整 Show。
- `REVAC_Q88_ASYNC_PRECHECK=0` 在新进程恢复原串行 preflight，用于对照。
- 新增完整 Show 负面测试：无效范数界、损坏 witness 均拒绝；恢复 witness 后 checked cache 恢复通过。

## 最终配置测量

这三轮测试时显式等价于开启异步候选（当时默认开启）；最终因不稳定收益改回默认关闭。三次独立生命周期，每次初始 Show 3 次，因此初始共 9 个样本，更新后共 3 个样本。生成包含序列化，验证包含完整 wire 解码、签名及 accumulator 验证。

| 场景 | 上轮生成均值 | 本轮生成均值 | 上轮验证均值 | 本轮验证均值 |
|---|---:|---:|---:|---:|
| 初始预热 Show | 3.250626 | 3.683521 | 1.756955 | 1.747198 |
| 更新后首次 Show | 4.153042 | 4.416124 | 1.592983 | 1.724296 |

| 更新后样本 | 生成秒 | 验证秒 | packet bytes |
|---|---:|---:|---:|
| run1 | 3.670200 | 1.630987 | 217371 |
| run2 | 5.935276 | 1.946318 | 217272 |
| run3 | 3.642898 | 1.595584 | 217294 |

run2 首个初始 Pack Prove 为 5.685662 秒，更新后 Pack Prove 为 4.920745 秒；均为 worker attempt 1/8 成功，不是完整 Pack 外层失败重试。其验证也较慢。该轮没有开启 detailed LNP profiling，无法确定耗时增加是否来自 LNP 拒绝采样或调度/负载；不能将原因断言为其中任何一种。

其他两轮更新后的异步 checked-lift 耗时分别为 0.578300 / 0.542354 秒，已与签名重叠；总生成 3.670200 / 3.642898 秒。

## 大小

- 初始完整 Show：平均 213.231879 KiB，范围 212.146484–215.722656 KiB。
- 更新后完整 Show：平均 212.219076 KiB，范围 212.179688–212.276367 KiB。

新增 cache/调度状态不写入 wire，证明字段和参数不变；序列化长度的样本差异不被解释为协议压缩收益。

## 完整阶段（秒）

本轮全测试总时间包括三次初始 Show 和新增两个 Show 负面测试，上轮只有一次初始 Show，所以墙钟总时长不能直接比较。同步没有真实网络 RTT。

| 阶段 | 上轮均值 | 本轮均值 |
|---|---:|---:|
| runtime_initialization | 0.001592 | 0.001636 |
| ta_setup | 4.155744 | 3.983595 |
| keygen_3uav_1gcs_and_test_identity_selection | 0.008515 | 0.007728 |
| registration_3uav_1gcs | 0.000111 | 0.000108 |
| issue_uav1 | 3.342763 | 3.234556 |
| issue_uav2 | 3.350000 | 3.226988 |
| issue_uav_local | 3.408270 | 3.194104 |
| prover_prewarm | 2.795224 | 2.687770 |
| verifier_prewarm | 1.449310 | 1.444997 |
| ta_tail4_prefetch | 0.503437 | 0.483009 |
| uav_tail4_prefetch_auth_cache | 0.005576 | 0.006034 |
| ta_delete_delta6_and_compatible_tail4 | 1.220191 | 1.141910 |
| ta_update_manifest_export | 0.000109 | 0.000076 |
| uncached_local_update_with_clear_verify | 1.312026 | 1.279413 |
| uav_delta6_update | 0.038181 | 0.038400 |
| post_update_show_generation_including_serialization | 4.153042 | 4.416124 |
| post_update_show_wire_verification | 1.592983 | 1.724296 |
| ta_memwitsync_export | 3.407960 | 3.323591 |
| memwitsync_recovery_with_clear_checks | 4.568612 | 4.593346 |
| full_test_wall_including_checks_and_cleanup | 44.689126 | 57.750251 |
| issue_three_credentials_total | 10.101032 | 9.655648 |
| initial_show_generation_including_serialization | 3.250626 | 3.683521 |
| initial_show_wire_verification | 1.756955 | 1.747198 |

## 验证与候选实验

最终三轮完整生命周期成功，含 Clear MemVer、三条 exact lift、Pack Verify、完整 packet 反序列化验证、缓存失效/恢复、完整异步 Show 负面测试、撤销/同步篡改拒绝。

稀疏 phi block 聚合候选为 8946 条约束创建 26826 个区块任务，974 个输出块。四次 audit 的所有六条内部 NTT prime 输出与串行结果逐元素一致，Clear Verify、Pack reload、binding negative checks 通过。但固定输入三次 profile 测得关闭候选 verify 均值 1.643663 秒，开启 1.634261 秒，首轮折叠均值没有下降，不视为稳定收益。已从生产 native engine 撤回，代码/日志保存在 ../show_overlap_20261014/candidates 与 fixed_input。

## 重现

```bash
REVAC_Q88_ASYNC_PRECHECK=1 \
PQAC_REPETITIONS=3 \
PQAC_DEMO=/root/autodl-tmp/pq_ac/OURS/code/entities/build-baseline-20261013/revocable_ac_demo \
PQAC_RESULTS_DIR=/root/autodl-tmp/pq_ac/OURS/experiments/q88_latency_20261008/results/recheck/run1 \
bash /root/autodl-tmp/pq_ac/OURS/code/entities/tools/measure_q88_show.sh
```

本轮前的源码备份：/home/wkk/pq_ac/production_backups/q88_pre_overlap_20261014/sources.tar.gz。

## 源码 SHA-256

- `OURS/code/entities/uav/revac_uav.c`: `d1f49b8b4ee6621dc33e2a2affe6d55f499ddf614d21e19d6fcff35080888e02`
- `OURS/code/entities/common/revac_labrador_bridge.c`: `0d05ee6280d6b7bfe2e9de547bb62564ae0b05c16a3e4c5cbbfa5444371e4adf`
- `OURS/code/entities/test/revocable_ac_demo.c`: `62320433be5cb4e948a9734ece2e8c5c814557385f60d2b0d8aada07cdcb8ecb`
- `lazer/src/labrados/constraints.c`: `d19cf83f10b35fbb704483246ad453df1a1fb0e406e5db2f5fe7839c26200091`

## 同一二进制显式开关对照

每个配置一个 fresh 生命周期，初始三次 Show，更新后一次 Show；不同随机性、顺序采样，无配对显著性结论。

| 场景 | 串行生成秒 | 异步生成秒 | 串行验证秒 | 异步验证秒 |
|---|---:|---:|---:|---:|
| 初始 | 3.552808 | 3.359766 | 1.918544 | 1.854090 |
| 更新后 | 4.301943 | 4.695898 | 1.650789 | 1.677414 |

初始候选有下降，更新后候选没有收益。保留作为显式 opt-in 实验，不作为默认生产优化。验证端本轮没有保留任何算法改动，不将验证差异声称为优化效果。生产算法为上轮 checked-lift cache + 串行 preflight。
