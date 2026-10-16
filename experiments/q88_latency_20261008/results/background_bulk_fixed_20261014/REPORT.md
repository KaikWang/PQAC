# 更新后后台预热与批量 exact-lift：2026-10-14

## 保留的实现

1. `revac_q88_prepare_async_start/finish`：前台生成并编码 q19 临时 opening，后台仅运行独占的整数/CRT scratch、检查 witness 和安装 worker witness。借用的 state/UAV 在 finish 前不能修改、销毁或同时执行 Show；finish join 后验证 digest-bound lift cache。新 Show 正常生成新 nonce、opening、mask/challenge，不复用证明或随机响应。
2. `acc96_exact_lift.c`：保留原整数乘法，使用 shift/truncate/sub/add 在 X^4096+1 中做准确折叠；scalar contributions 使用本次已转换的相同 rhs 多项式，不再逐系数重新做 CRT 转换。所有整除检查、window/tag 检查和范数界均保留。默认 bulk=1；REVAC_Q88_BULK_LIFT=0 可恢复参考实现，=2 做每个乘加与参考逐系数算法的严格等价审计。
3. 修复 warm prover prewarm 的两个旧缺陷：verifier-only 提前返回不再吞掉 prepare_prover 请求；full-ring wrapper 使用关键字传递 prepare_prover，避免落入两个闲置的位置参数。修复前，更新后 U 请求只 rebind RHS 却未加载新 witness，导致后续 compact R 请求被完整关系 guard 拒绝并重试。修复后完整新 witness 在发布 last_prover 之前被安装。

后台接口已在 production demo 的更新后路径接入（REVAC_Q88_POST_UPDATE_PREWARM=1）。默认业务时序不强制等待这项额外工作；应用可在更新后 idle 期 start，并在 Show/再次更新前 finish。实验立即 finish，以如实报告全部准备开销，没有用 sleep 隐藏工作。

## 完整 Show 实测

autodl-pqa11，三次 fresh 生命周期，初始/更新后各一个完整 Show。生成包括 compact serialization，验证包括完整 packet 解码、签名及 Pack verification。生产 native 引擎恢复原算子，没有保留 folded-center 候选。

| 项目 | 上轮均值秒 | 本轮均值秒 |
|---|---:|---:|
| 初始 Show 生成 | 3.250626 | 3.246313 |
| 初始 Show 验证 | 1.756955 | 1.657460 |
| 更新后准备 | 无单列，首个 Show 内计算 | 0.810717 |
| 更新后 Show 生成 | 4.153042 | 3.123185 |
| 更新后 Show 验证 | 1.592983 | 1.719312 |

| 样本 | 准备 | Show 生成 | 准备 + 生成 | UAV 更新 + 准备 + 生成 | Show 验证 |
|---|---:|---:|---:|---:|---:|
| 1 | 0.809731 | 3.116689 | 3.926420 | 3.963519 | 1.786460 |
| 2 | 0.794328 | 3.122224 | 3.916552 | 3.954035 | 1.600032 |
| 3 | 0.828093 | 3.130641 | 3.958734 | 3.995823 | 1.771444 |

准备+生成平均 3.933902 秒，最大 3.958734 秒。UAV 更新+准备+生成平均 3.971126 秒，最大 3.995823 秒。该总和不含 TA 撤销生成、网络传输，也不含最后 verifier 时间。

三次样本的初始与准备后 Show 均达到 4 秒 / 2 秒；并非最坏情况保证，仍未达到 3 秒 / 1.5 秒。上下两组三样本不配对，随机重试/负载存在，不能将全部均值差归因于算法优化。

## 相同二进制 bulk 对照

control0/control1 各一 fresh 生命周期，未开启更新后预热。首次更新后 MemVer/checked lift preflight：参考 0.566038 秒，bulk 0.436204 秒，减少 0.129834 秒。此单样本对照与逐项整数审计支持减少重复转换的判断，不作为整体 Show 稳定速度保证。

## 大小

- 初始完整 packet 平均 212.489258 KiB，范围 212.359375–212.739258 KiB；accumulator envelope 平均 188.586914 KiB。
- 更新后完整 packet 平均 213.481120 KiB，范围 212.314453–215.789062 KiB；accumulator envelope 平均 188.555664 KiB。

未增加 packet 字段或证明维度，proof modulus/CRT/security parameters 不变。更新后 packet 均值比上轮约增加 0.61%，这里记录实际随机证明及压缩结果，不宣称严格零字节增长。

## 验证与失败候选

- bulk 审计的所有乘加、scalar contribution 与原 fmpz coefficient loop 严格相等；三个 CRT lane 的原 exact divisibility、窗口和 tag 检查通过。
- 初始/更新后 resident witness Clear Verify 审计均通过，后续 compact R 请求成功；三轮无 rejected-invalid-native 或 proof-generation-guard-failed 日志。
- 完整生命周期、packet 反序列化后 Verify、lift cache 篡改拒绝、完整 Show 损坏 witness/无效范数界拒绝及恢复均通过。
- 首次预热尝试失败的数据保留在 ../background_bulk_20261014/run1 与 run2，首次修复遗漏 wrapper 位置参数的失败保留在 refresh_audit。修复后的审计是 refresh_audit_fixed，不能将失败/审计组纳入正常性能均值。
- native folded-center 候选：1,000,007 个 standalone INT128 输入及 33,901,824 个实际 native 系数与原整数余数严格相同，Clear Verify/Pack reload/negative binding tests 均通过。但固定输入三次 Verify 平均参考约 1.671 秒、候选约 1.933 秒，未显示稳定收益，撤回生产代码与库。候选代码/日志保存在 native_candidate。
- 本轮验证算法未保留变化，验证均值差不能当作 native Verify 加速成果。

## 全流程统计（秒）

| 阶段 | 本轮均值 |
|---|---:|
| runtime_initialization | 0.001501 |
| ta_setup | 3.957226 |
| keygen_3uav_1gcs_and_test_identity_selection | 0.011927 |
| registration_3uav_1gcs | 0.000120 |
| issue_uav1 | 3.322325 |
| issue_uav2 | 3.261959 |
| issue_uav_local | 3.270243 |
| prover_prewarm | 2.766609 |
| verifier_prewarm | 1.440459 |
| ta_tail4_prefetch | 0.501628 |
| uav_tail4_prefetch_auth_cache | 0.005416 |
| ta_delete_delta6_and_compatible_tail4 | 1.203822 |
| ta_update_manifest_export | 0.000147 |
| uncached_local_update_with_clear_verify | 1.306346 |
| uav_delta6_update | 0.037224 |
| post_update_prepare_wall | 0.810717 |
| post_update_prepare_task | 0.766716 |
| post_update_show_generation_including_serialization | 3.123185 |
| post_update_show_wire_verification | 1.719312 |
| ta_memwitsync_export | 3.425189 |
| memwitsync_recovery_with_clear_checks | 4.636371 |
| full_test_wall_including_checks_and_cleanup | 44.408152 |
| issue_three_credentials_total | 9.854526 |
| initial_show_generation_including_serialization | 3.246313 |
| initial_show_wire_verification | 1.657460 |

全测试总时长包括负面测试、LOCAL 回退、同步/恢复和清理；同步回调不含真实网络 RTT。

## 使用接口

```c
revac_q88_prepare_async_t *job = NULL;
if (!revac_q88_prepare_async_start(&job, &state, &uav)) { /* handle failure */ }
/* idle work that does not mutate state/UAV or run Show */
double prepare_seconds = 0;
int ready = revac_q88_prepare_async_finish(job, &prepare_seconds);
/* only if ready: generate Show with fresh nonce */
```

## 重现

```bash
REVAC_Q88_POST_UPDATE_PREWARM=1 REVAC_Q88_BULK_LIFT=1 \
PQAC_REPETITIONS=1 \
PQAC_DEMO=/root/autodl-tmp/pq_ac/OURS/code/entities/build-baseline-20261013/revocable_ac_demo \
PQAC_RESULTS_DIR=/root/autodl-tmp/pq_ac/OURS/experiments/q88_latency_20261008/results/bg_recheck/run1 \
bash /root/autodl-tmp/pq_ac/OURS/code/entities/tools/measure_q88_show.sh
```

等价/准备审计另加 REVAC_Q88_BULK_LIFT=2、REVAC_Q88_PREPARE_AUDIT=1，不将审计时间混入正常测量。

## 源码 SHA-256

- `OURS/code/entities/common/core/src/acc96_exact_lift.c`: `8d12a8be2fb7193b565ffc7f3e3559deb061779a6bd038d94844b9ddbf3b39f0`
- `OURS/code/entities/common/revac_labrador_bridge.c`: `f3ef59f3aecb50cedcbc9fdcb2f10df6205977699ed797e77dca6b0536fa4163`
- `OURS/code/entities/common/revac_labrador_bridge.h`: `2d1b107a15a0f265c1d8db0cd11678528386548804be46724740ab4968d80e2b`
- `OURS/code/entities/test/revocable_ac_demo.c`: `3fa1741ab679d53adc6bdca33f431da7a0bde9d2ca2a9a509001b668822b3a8a`
- `OURS/code/entities/zk_lazer/revac_acc96_full_ring.py`: `1bb169a2778b32378dea00709664d8b804cd461276f39e1d56b040a20527cabc`
- `OURS/code/entities/zk_lazer/revac_acc96_memver_labrador.py`: `e19848e3500e082916f68b0819b568f1479678eede26709517db4f65fa937651`
- `OURS/code/entities/zk_lazer/revac_acc96_memver_worker.py`: `f3d3079ef957642a1fe474ad9fc821c48576967d6dbdf845f617fd52ee3027c0`
- `lazer/src/labrados/polx.c`: `9049cf9855dd07893b4db5237815d3348e743d8bc63ddc33996fefafbcb522c3`
