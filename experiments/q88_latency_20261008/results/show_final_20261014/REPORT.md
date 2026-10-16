# 完整 Show 优化复测：2026-10-14

服务器：autodl-pqa11，Xeon Gold 6430；每组为三次独立全流程，每次新 worker、Setup、三个凭证签发、撤销、更新、Show 和同步恢复。最终配置：fused=4、linear=4、sprod=6；无 CPU 绑定、无投影并行分支、无 PGO。

## 原因与保留的改动

更新后的第一次 Show 会发生 MemVer 缓存未命中，旧路径先做模关系检查，再在异步证明线程重复生成 exact lift。新路径保留整体范数检查，并用三条精确 CRT 整除关系、tag/window 检查完成预检查；缓存其 exact lift，bridge 按状态摘要校验后复用。不把此计算挪到更新阶段。缓存绑定原有所有 witness 块、handle、累加器、序号、PP 摘要，并新增范数上界绑定。改变 witness 或范数界的负面测试及恢复测试三轮全通过。

相同固定输入的 profiling 曾出现 LNP outer attempts 4 与 18，说明随机拒绝采样能明显改变生成时间；不能据三次非配对样本断言所有时间差都来自实现。验证路径没有因这个 C 缓存改变，验证改善/变慢不能归因于缓存。

## 全流程时间

单位为秒；前值是上轮 lifecycle_opt_20261014，后值是本报告 show_final_20261014。实现未变的行不声称获得算法优化。

| 项目 | 上轮均值 | 本轮均值 | 变化 | 本轮实现变化 |
|---|---:|---:|---:|---|
| 运行时初始化 | 0.001688 | 0.001592 | -5.69% | 无 |
| TA Setup | 4.045091 | 4.155744 | +2.74% | 无 |
| 3 UAV + 1 GCS 密钥生成（含测试身份选择） | 0.008250 | 0.008515 | +3.21% | 无 |
| 3 UAV + 1 GCS 注册 | 0.000109 | 0.000111 | +2.17% | 无 |
| 签发 UAV1 | 3.285408 | 3.342763 | +1.75% | 无 |
| 签发 UAV2 | 3.401617 | 3.350000 | -1.52% | 无 |
| 签发 LOCAL UAV | 3.531164 | 3.408270 | -3.48% | 无 |
| 三个凭证签发合计 | 10.218190 | 10.101032 | -1.15% | 无 |
| Prover 首次准备 | 3.286395 | 2.795224 | -14.95% | 复用已检查 exact lift |
| Verifier 首次准备 | 1.444547 | 1.449310 | +0.33% | 无 |
| TA T4 预生成 | 0.511594 | 0.503437 | -1.59% | 无 |
| UAV T4 验签缓存 | 0.005408 | 0.005576 | +3.10% | 无 |
| TA 撤销 + DELTA6 + 兼容 T4 | 1.219896 | 1.220191 | +0.02% | 无 |
| 更新 manifest 导出 | 0.000114 | 0.000109 | -4.90% | 无 |
| 无缓存 LOCAL 更新及明文检查 | 1.311312 | 1.312026 | +0.05% | 无 |
| UAV DELTA6 更新 | 0.037055 | 0.038181 | +3.04% | 无 |
| 更新后完整 Show 生成（含序列化） | 4.673874 | 4.153042 | -11.14% | 复用已检查 exact lift |
| 更新后完整 Show wire 验证 | 1.780982 | 1.592983 | -10.56% | 无 |
| TA MemWitSync 生成 | 3.438360 | 3.407960 | -0.88% | 无 |
| MemWitSync 恢复及检查 | 4.591347 | 4.568612 | -0.50% | 无 |
| 全测试墙钟（含负面测试及清理） | 45.588072 | 44.689126 | -1.97% | 汇总项 |
| 预热后完整 Show 生成（含序列化） | 3.162552 | 3.250626 | +2.78% | 复用已检查 exact lift |
| 预热后完整 Show wire 验证 | 1.664477 | 1.756955 | +5.56% | 无 |

完整 Show 更新后生成范围 3.827660–4.417067 秒，验证范围 1.555628–1.645093 秒。未更新的预热 Show 生成 3.148498–3.377073 秒，验证 1.697739–1.811327 秒。更新后平均生成仍超过 4 秒，更未达到 3 秒 / 1.5 秒。

签发平均每凭证：上轮 3.406063 秒，本轮 3.367011 秒；本轮没有再次更改签发算法。

## 同一当前二进制关闭缓存的对照

另跑三次新 worker 完整生命周期，仅设置 REVAC_Q88_CHECKED_LIFT_CACHE=0。顺序采样、不同随机性，非配对统计检验。

| 项目 | 缓存关闭秒 | 缓存开启秒 | 变化 |
|---|---:|---:|---:|
| Prover 首次准备 | 3.258286 | 2.795224 | -14.21% |
| 预热 Show 生成 | 3.123730 | 3.250626 | +4.06% |
| 预热 Show 验证 | 1.704367 | 1.756955 | +3.09% |
| 更新后 Show 生成 | 4.293197 | 4.153042 | -3.26% |
| 更新后 Show 验证 | 1.725219 | 1.592983 | -7.66% |

这个更近的对照显示：准备时间减少约 0.46 秒；更新后生成仅减少约 0.14 秒（3.26%），不能把与上轮比较的 11.14% 全算作缓存收益。初始 Show 生成/验证没有稳定改善，验证时间差不是 C lift 缓存的直接收益。下一步应在保留全部校验的前提下，使首次更新后 lift/预检查与匿名签名更充分重叠，并检查 Pack Verify 中的主要折叠算子。

## 序列化大小

| 项目 | 上轮平均 KiB | 本轮平均 KiB |
|---|---:|---:|
| 预热 Show packet | 212.328776 | 214.304362 |
| 更新后 Show packet | 212.418620 | 212.183919 |

新缓存只驻留本地内存，不写入 packet；未增加协议字段或改变 proof 参数。本轮初始 packet 实测约增加 1%，更新后约减少 0.11%；不能把随机证明和压缩长度差异宣称为压缩优化。

本轮初始 accumulator envelope 为 192761–192825 bytes；更新后为 193177–193305 bytes（188.649–188.774 KiB）。这些是完整 accumulator envelope，与完整 Show packet 分开统计。

## 验证与未保留实验

- 三轮完整生命周期成功，包含独立原始 Clear MemVer、三条 exact CRT lift、Pack Verify、完整 packet 反序列化验证、撤销/同步篡改拒绝。
- 八线程候选未在完整 Show 中获得稳定收益，撤回至 4/4。
- 投影并行先因阈值未触发，修正后四次 NTT 结果逐元素对照通过，但普通测量没有验证收益，撤回生产修改。候选代码保存在 ../show_opt_20261014/parallel_projection_labradoodle.c。
- 通用物理核选择候选未显示稳定收益，撤回。
- 已通过的功能/等价检查不构成新的密码学安全认证。
- 全测试总时长包含负面测试、额外 LOCAL 与同步检查，不等于最小业务流程；同步回调不含真实网络延迟。

## 重现

```bash
PQAC_DEMO=/root/autodl-tmp/pq_ac/OURS/code/entities/build-baseline-20261013/revocable_ac_demo \
PQAC_REPETITIONS=1 \
PQAC_RESULTS_DIR=/root/autodl-tmp/pq_ac/OURS/experiments/q88_latency_20261008/results/recheck/run1 \
bash /root/autodl-tmp/pq_ac/OURS/code/entities/tools/measure_q88_show.sh
```

设置 REVAC_Q88_CHECKED_LIFT_CACHE=0 可关闭新 exact-lift 复用；用新进程启动。备份位于 production_backups/q88_pre_show_optimization_20261014/sources.tar.gz。

## 源码 SHA-256

- `code/entities/uav/revac_uav.c`: `de71d069e10e024d3e695d707a090f6e33e66b7f2870dc64678084574bfc8e7a`
- `code/entities/uav/revac_uav.h`: `fd52b05705af1ba537ee391d8c766747c7e2eac1c88b21b815e607a1be8f6a08`
- `code/entities/common/revac_labrador_bridge.c`: `f8c781fba3469e0cce742b9eac7ee55580dfc5f41d105de0148671aa21efae52`
- `code/entities/common/revocable_ac.h`: `37d78106e01dc00bb24a46e4f4e26a8e1e6a8e198a7e430f5b153129453b524a`
- `code/entities/test/revocable_ac_demo.c`: `a7f9a1133962b653da0d31280c5cd3e7d62010f41fdf2357629f4c5c11320573`
- `code/entities/zk_lazer/revac_q88_profile.py`: `dfbc32259df2003e21b40fb4ed2265b65adb34980cd7d92e0bf589221c755624`
- `lazer/src/labrados/labradoodle.c`: `ad87a1906792a5df3947cd2fc24a4b32c13d195dcf4308c2091b4bfa21159dec`
