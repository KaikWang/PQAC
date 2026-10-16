# q88 独立编译器优化实验

当前保留候选为 **accumulator norm grouping v2，group size 4**。
生产源文件位于 `OURS/code/entities/zk_lazer`；此目录中的两个同名文件为独立副本。

后续核对更正：compact 模式下部分 aggregate product 数位同样被原 compiler 标为
accumulator witness，因此 v2 也扩大了它们的 JL requirements。
下文对 product auxiliary requirements 不变的描述不准确，见结果报告的更正节。
后续 `../dimension_v3` 使用独立分类，明确排除这些 auxiliary。

## 已测有效的方法

`revac_acc96_full_ring.py` 仅对 accumulator 数位向量设置分组公共 JL requirement。
每组最多四个已有向量，公共 requirement 为原 requirement 总和向上取二次幂。
原始输入 norm、witness 系数、所有线性/二次约束保持一致。
lift/carry、refresh 和 product auxiliary 的单独 requirements 保持一致。
编译后重新运行原有 `_audit_norm_profile` 和 native 参数生成。

`revac_acc96_memver_worker.py` 将分组配置绑定到 compiler profile ID，
使旧参数、新参数以及不同组大小不能在常驻缓存或 compact proof 中混用。
默认 `REVAC_Q88_NORM_GROUP_SIZE=1`；候选须显式设置为 4。

## 运行完整 show

在已同步并完成基线编译的服务器上运行：

```bash
bash /root/autodl-tmp/pq_ac/OURS/experiments/q88_latency_20261008/native_sparse/run_norm_group4.sh
```

默认使用 Sage Python 和 `build-baseline-20261013/revocable_ac_demo`。
可设置 `PQAC_PYTHON`、`PQAC_DEMO`、`PQAC_RESULTS_DIR`、`PQAC_REPETITIONS`。
脚本启动专用 worker，等待 socket 就绪，并在退出时回收该 worker。
使用现有生产 proof engine；候选不依赖本目录的原生库。

## 独立证明与验证

`bench_worker.py` 使用服务器私有 fixture，执行 Clear Verify、Pack Prove、
compact proof 序列化、独立 verifier statement 反序列化验证。
`--verify-artifact` 可以在新进程重建公开 verifier 后验证已保存的证明。
`--negative-checks` 检查公开 accumulator、q19 link RHS 以及证明封装篡改；
`--expect-rejection` 用于确认错误 profile 被拒绝。
私有 fixture 不纳入本目录或结果归档。

## 未采用的原生优化

`polx.c`、`scratch_pool.h`、`Makefile` 为稀疏零项及 scratch pool 的独立实验。
它们尚未给出稳定的完整证明时间收益，当前 group4 候选无需加载这些实现。
不同近似范数向量排序也未减少实际 Pack 维度，当前使用 baseline 排序。

## 结果与安全范围

五次新完整 show 均低于 4 秒/2 秒，见
`../results/norm_group4_20261013/REPORT.md`。
该配置改变中间 JL requirement，不能描述为“所有范数参数不变”。
现有全局 accumulator 范数门仍为 `beta0 < 2^80`；当前保守重组界为
`log2(norm) <= 75.059`。这些是代码范数审计、参数生成与功能验证结果，
不等同于对整个密码构造的独立形式化安全证明。
