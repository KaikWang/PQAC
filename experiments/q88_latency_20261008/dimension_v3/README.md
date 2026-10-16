# q88 完整 CRT 独立优化 profile v3

生产 Python/C/native 源码保持不变。此目录使用生产完整 degree-4096 CRT 关系、
hidden-handle 与 refresh 绑定、CRT lift、四行 q19 link 及 shared handle。
保留主 witness radix 32767、原证明域、LNP 重复次数、Gaussian 参数和刷新窗口。

## 实际改动

- 明确按原始 accumulator limb 分组，排除 aggregate auxiliary；按 limb 加权重建 beta 界。
- 从公开 50000-delete 上界推导原 witness 最高位界；handle scalar 精确范数由 q19 binary
  shared-handle 和 degree-zero 绑定保证。生产安全参数门重新运行。
- 三个 auxiliary 数位：两个 balanced radix-3779 low limbs 和一个未饱和的 bounded top limb。
  最高位公开系数界 9624；准确覆盖整个 proof field。原 auxiliary JL requirement 保留。
- 引入 scaled refresh Z = 4095 Z，追加六个 full-ring 线性 binding，保留原 Z exact norm。
  降低 quadratic 系数的最大值，native quadratic left decomposition 从三位降至两位。
- 批量 split-multiply 延迟 refresh，融合 challenge 线性聚合；audit 模式逐元素对照旧实现。
- 临时分配 exact-size pool（选定缓存上限 1024 MiB）；六个内部 CRT prime 的点积分别计算。
  每个 prime 内部的累加顺序不变；避免嵌套 OpenMP。默认原生函数仍保持串行，选定 profile 显式开启。

实际首次 Pack `nn` 35492 → 31164；二次 `nn` 为 22169，`nmax` 为 1631。
这不是原始逻辑 witness 数量同比减少。新增 scaled Z 会增加 96 个 base-ring 多项式。
参数选择使用原 native 安全性门；不能把 Clear/Pack 验证通过称为独立密码学安全认证。

## 复现（autodl-pqa11）

```bash
cd /root/autodl-tmp/pq_ac/OURS/experiments/q88_latency_20261008/dimension_v3
make LIBRARY=liblabrador38_sprod.so -j2
bash run_online.sh
```

`PQAC_RESULTS_DIR` 指定结果目录；`PQAC_REPETITIONS` 默认 5。
`PQAC_NATIVE_LIBRARY` 可切换独立引擎。基线库与生产 C demo 均不改动。
构建使用现有 native 的 AVX512/OpenMP flags，没有 PGO 或新服务器专属参数。
主要算法改动不依赖这台 CPU，但本构建继承的指令集需求仍需在迁移时检查。

公开边界检查：source profile.sh 后用 Sage Python 执行 check_mixed_radix.py。
完整 proof fixtures 只存服务器；禁止将 private prover payload 纳入结果或仓库。
交叉验证：同一 Python profile 下 unset LD_PRELOAD，再用 bench_worker.py --verify-artifact。

完整测量与局限见 ../results/dimension_v3_20261013/REPORT.md。
在线时间不含 statement/参数初次构建；实际 serialized accumulator envelope 与完整 show packet
单独记录。当前目标 3s/1.5s 尚未达成。
