# 连续三次撤销更新与完整 Show：2026-10-14

同一 fresh 生命周期、同一 survivor UAV，同一 worker、同一生产 CRT profile。生产路径为 checked-lift cache + 串行 preflight，ASYNC_PRECHECK 默认关闭。独立实验 binary 从当前 demo 生成，未修改生产源文件。

初始状态先按生产方式预热；第一次撤销使用已有 DELTA6 路径。第二、三次为真实 TA 撤销，经签名 manifest/T4 或 LOCAL 更新，均通过 Clear MemVer。第二次撤销在一次正常 MemWitSync 恢复后进行，因此这是包含同步的业务链，不是三个完全相同更新路径的配对基准。Show 内不会额外预热；首个更新后 Show 含状态变化引起的缓存未命中成本。所有 Show 包含序列化并经完整 guarded wire verification。

| 状态 | 第一次生成 / 验证（秒） | 后续生成 / 验证均值（秒） |
|---|---:|---:|
| 初始已预热状态 | — | 3.213675 / 1.803847（3次） |
| 更新 1 次 | 3.901841 / 1.775424 | 3.364892 / 1.799525（3次） |
| 更新 2 次 | 4.104993 / 1.966831 | 3.452501 / 1.826400（2次） |
| 更新 3 次 | 3.691997 / 1.787329 | 3.082020 / 1.767171（2次） |

## 证据与解释

- 布局审计在所有状态相同：vectors=142，backend_polynomials=16448，approx_vectors=109。这里是进入 Dachshund 前的 compiler witness 布局，不是折叠后 nn。生产缓存还要求 dimensions/constraint_count 完全匹配并复用原 Pack 参数。
- 生成方面：初始样本已预热；改变 accumulator/sequence/witness 后，旧私有 witness 与 checked-lift cache 不能再命中。首次 Show 多做 exact lift、payload 编码、动态 witness 安装。重复 Show 可复用这部分，仍使用新 nonce、随机性、opening 和 challenge。
- 验证方面：更新只改变公开 accumulator RHS/context，同一 profile 的维度和 Pack 参数保持不变，不存在更新越多算法工作量必然越小的机制。冷热缓冲、随机证明差异、进程调度和主机负载均可能影响墙钟值，本实验不能分别量化其贡献。第二次更新后 verifier 更慢，直接否定“更新后必然更快”。
- 这仅是一条业务链的三次更新测试，不能外推认证 50000 次更新性能、安全性或容量。当前允许窗口内不会将所有历史更新证明串接在 Show 中，但 witness 范数和窗口边界仍须检查。

## 两条不同的优化路线

1. 在更新后 idle/background prepare 新状态，可以使后续 Show 使用 warm path；必须同时报告 update + prepare 的就绪时间，不能声称这部分计算消失。
2. 真正减少总耗时：减少 exact-lift 重复公共多项式转换/乘法、复用状态无关编译内容，并优化 native Pack folding/refresh。验证初始和更新后使用同一个优化；不能通过更新次数获得约定的固定 verifier 时间。

## 重现

```bash
python make_demo.py
python build_demo.py /root/autodl-tmp/pq_ac/OURS/code/entities/build-baseline-20261013
PQAC_DEMO=/root/autodl-tmp/pq_ac/OURS/experiments/q88_latency_20261008/sequence_show/sequence_show_demo \
PQAC_REPETITIONS=3 REVAC_Q88_LAYOUT_AUDIT=1 \
PQAC_RESULTS_DIR=/root/autodl-tmp/pq_ac/OURS/experiments/q88_latency_20261008/results/sequence_recheck \
bash /root/autodl-tmp/pq_ac/OURS/code/entities/tools/measure_q88_show.sh
```

JSON: sequence_measurement.json；完整 stdout/worker 时间：show.log、worker.log。
