# CRT 共享 exact-lift 准备优化：2026-10-14

## 本轮结论

保留 `acc96_exact_lift.c` 的 per-call 共享 refresh 矩阵与公共 witness 整数转换。三个 CRT lane 读取本次调用内的相同不可变数据，用完清理；未引入跨调用私有 witness 缓存。三素数、范数/窗口/tag gate、逐 lane exact divisibility、证明关系与序列化字段均未修改。该改动减少准备工作；warm native Pack Prove/Verify 算法没有保留新增改动。

稀疏旋转与公开零 RHS 两个 native 候选均已撤回，生产 `constraints.c` 与本轮修改前备份逐字相同。保留候选源码及日志仅供复查，不参与生产编译。

## 固定输入 exact-lift 对照

同一 witness/state、同一二进制，三种模式交错轮换，各五次；此诊断组额外 15 次调用不属于下面正常生命周期性能样本。

| 模式 | 均值秒 |
|---|---:|
| 原始三 lane | 0.363406822 |
| 共享 refresh | 0.292523700 |
| 共享 refresh + witness 转换 | 0.284149203 |

完整共享比参考减少 21.81%。样本量较小，公共转换在共享 refresh 之上的独立收益不作稳定百分比承诺。

模式 4 完整生命周期将共享结果与原实现比较：返回状态、成功输出的全部 lift 系数及最大位数完全相等，成功和拒绝路径均有 PASS。prepare 审计 Clear Verify 通过，后续新鲜 Show 的 native guard、Pack Verify、序列化及反序列化验证通过；负面测试拒绝损坏 lift cache/witness/范数并可恢复。详见 `../shared_lift_20261014/common_audit_bench/`。

## 三轮 fresh 完整生命周期

服务器 autodl-pqa11，每轮 fresh setup/三凭证签发/初始 Show/撤销/更新/准备/更新后 Show/同步恢复。生成时间包含完整 packet 序列化，验证包含 packet 解码、签名和累加器验证。TA、网络和准备未混入在线 Show；准备与更新到生成的总和另列。

| 项目 | 上轮均值秒 | 本轮均值秒 | 解释 |
|---|---:|---:|---|
| 初始完整 Show 生成 | 3.246313 | 3.561744 | 本轮变慢，不能宣称初始 Show 提速 |
| 初始完整 Show 验证 | 1.657460 | 1.682232 | 本轮变慢，无新增 Verify 优化 |
| 更新后准备 wall | 0.810717 | 0.690663 | 共享准备有固定输入对照支持 |
| 更新后完整 Show 生成 | 3.123185 | 3.028749 | 均值降低，但 native 算法未变，含随机波动 |
| 更新后完整 Show 验证 | 1.719312 | 1.700903 | 本轮均值，无稳定加速结论 |

| 样本 | 初始生成 | 初始验证 | 更新后准备 | 更新后生成 | 更新后验证 | 准备+生成 | UAV 更新+准备+生成 |
|---|---:|---:|---:|---:|---:|---:|---:|
| run1 | 4.058819 | 1.798762 | 0.739119 | 3.125926 | 1.860126 | 3.865044 | 3.904282 |
| run2 | 3.545075 | 1.671849 | 0.679733 | 2.734031 | 1.644241 | 3.413764 | 3.450219 |
| run3 | 3.081339 | 1.576085 | 0.653138 | 3.226291 | 1.598342 | 3.879429 | 3.915439 |

准备+更新后生成平均 **3.719412 s**（上轮 3.933902 s）；UAV 更新+准备+生成平均 **3.756647 s**（上轮 3.971126 s）。两组 fresh 样本不配对，均值差不能全部归因于共享准备。

初始生成均值 3.561744 s，比上轮 3.246313 s 更慢。最慢初始完整 Show 为 4.058819 s，其中 native Pack Prove 3.931640 s；worker 请求总耗时 3.944333 s。该请求使用缓存并通过 native guard，且应用层为 attempt 1/8；没有内部 LNP 重试次数记录，不能确定是内部重试还是系统负载造成。此时 exact-lift 已在先前 prewarm 完成，改动不在该在线路径中，但仍保留整个慢样本，不宣称所有 Show 已低于 4 s。3 s / 1.5 s 尚未稳定达成。

## 实际序列化大小

| 样本 | 初始 accumulator bytes | 初始完整 packet bytes | 更新后 accumulator bytes | 更新后完整 packet bytes |
|---|---:|---:|---:|---:|
| 1 | 193305 | 217247 | 193081 | 215303 |
| 2 | 193081 | 215198 | 193273 | 217254 |
| 3 | 193305 | 217276 | 193305 | 217282 |

完整 packet 均值：初始 **211.497721 KiB**，更新后 **211.536133 KiB**；accumulator envelope 分别 **188.701497 / 188.691081 KiB**。上轮完整 packet 均值 212.489258 / 213.481120 KiB。本轮未增加字段或 proof 维度；大小仍随随机证明及压缩变化，不能将均值下降当作编码优化。

## 未保留的 native 候选

1. 公开稀疏旋转：与原实现逐项比较 NTT projection 和 width metadata，并通过固定完整 proof/负面绑定/序列化重载验证，但实际 q88 optimized_blocks=0，没有命中，撤回。
2. 公开零 RHS：使用公开 width 确认零项，audit 复查原始 NTT 系数；保留 canonical reduction，并将 getenv 提出高频循环。功能检查通过，但同一 proof 五轮交错 Verify：参考平均 1.685671 s，候选平均 1.690708 s，无稳定收益，撤回。此处重复验证同一 proof 仅是纯 worker 运算微基准，不是完整 Show 的防重放验证。

证据：`../shared_lift_20261014/candidates/`、`../shared_lift_20261014/candidate_logs/`。没有把候选阶段数据混入最终三轮。

## 全流程均值秒

本表列出两轮均值，未改动的阶段差异仅作为观测，不认定为优化效果。MemWitSync 回调为本地认证回调，不含真实网络 RTT；全测试 wall 包含负面测试、LOCAL 回退、恢复和清理。

| 阶段 | 上轮 | 本轮 |
|---|---:|---:|
| runtime_initialization | 0.001501 | 0.001512 |
| ta_setup | 3.957226 | 3.911883 |
| keygen_3uav_1gcs_and_test_identity_selection | 0.011927 | 0.008898 |
| registration_3uav_1gcs | 0.000120 | 0.000108 |
| issue_uav1 | 3.322325 | 3.290773 |
| issue_uav2 | 3.261959 | 3.273088 |
| issue_uav_local | 3.270243 | 3.317137 |
| prover_prewarm | 2.766609 | 2.699944 |
| verifier_prewarm | 1.440459 | 1.517794 |
| ta_tail4_prefetch | 0.501628 | 0.493703 |
| uav_tail4_prefetch_auth_cache | 0.005416 | 0.005421 |
| ta_delete_delta6_and_compatible_tail4 | 1.203822 | 1.245614 |
| ta_update_manifest_export | 0.000147 | 0.000126 |
| uncached_local_update_with_clear_verify | 1.306346 | 1.289877 |
| uav_delta6_update | 0.037224 | 0.037234 |
| post_update_prepare_wall | 0.810717 | 0.690663 |
| post_update_prepare_task | 0.766716 | 0.647256 |
| post_update_show_generation_including_serialization | 3.123185 | 3.028749 |
| post_update_show_wire_verification | 1.719312 | 1.700903 |
| ta_memwitsync_export | 3.425189 | 3.336159 |
| memwitsync_recovery_with_clear_checks | 4.636371 | 4.682660 |
| full_test_wall_including_checks_and_cleanup | 44.408152 | 44.327069 |
| issue_three_credentials_total | 9.854526 | 9.880998 |
| initial_show_generation_including_serialization | 3.246313 | 3.561744 |
| initial_show_wire_verification | 1.657460 | 1.682232 |

## 开关、重现与备份

- `REVAC_Q88_LIFT_SHARED_REFRESH=3`：共享 refresh 与公共转换，现为默认（环境变量未设置或空字符串）。
- `=0`：恢复参考三 lane；`=1`：仅共享 refresh；`=2`：共享 refresh 并逐项审计；`=4`：完整共享并逐项审计。
- `REVAC_Q88_LIFT_BENCHMARK=1`：demo 启用上面的三模式交错诊断，默认关闭。
- `REVAC_Q88_POST_UPDATE_PREWARM=1`：demo 单列更新后准备；应用接口仍遵守借用对象在 finish 前不变的生命周期要求。
- 修改前备份 `/home/wkk/pq_ac/production_backups/q88_pre_shared_lift_20261014/`。

最终三轮使用显式 mode 3；测完仅将未设置环境变量时的默认改为 3，没有改运算逻辑，并在本地/服务器重新构建成功。

```bash
REVAC_Q88_LIFT_SHARED_REFRESH=3 REVAC_Q88_POST_UPDATE_PREWARM=1 \
PQAC_REPETITIONS=1 \
PQAC_DEMO=/root/autodl-tmp/pq_ac/OURS/code/entities/build-baseline-20261013/revocable_ac_demo \
PQAC_RESULTS_DIR=/root/autodl-tmp/pq_ac/OURS/experiments/q88_latency_20261008/results/shared_lift_recheck/run1 \
bash /root/autodl-tmp/pq_ac/OURS/code/entities/tools/measure_q88_show.sh
```

审计另用 mode 4 与 REVAC_Q88_PREPARE_AUDIT=1，不将审计计入正常性能。

## 最终源码 SHA-256

- `code/entities/common/core/src/acc96_exact_lift.c`: `71bba8e64b0cde6bc93561c33e43f014dff2ce2c08b6ed0c7ba002d39bc9a6c8`
- `code/entities/test/revocable_ac_demo.c`: `d2e94f731d69c55ddc35158ea202b2747ab99ddab23c8b3cc22bd5123c29a59e`
- `../lazer/src/labrados/constraints.c`: `d19cf83f10b35fbb704483246ad453df1a1fb0e406e5db2f5fe7839c26200091`
