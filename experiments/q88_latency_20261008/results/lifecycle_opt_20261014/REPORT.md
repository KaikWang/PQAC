# 签发、首次准备与 MemWitSync 优化

2026-10-14；autodl-pqa11；Intel Xeon Gold 6430；当前生产 q88 compact CRT v3。
对照为同日 lifecycle_v3_20261014 的三次独立全流程；优化后另跑三次独立全流程。
每轮重新 Setup、注册和签发三个凭证，撤销前后各生成、验证一次完整 Show。
双方计时范围相同；正常测量关闭参考对照审计和原子 profiling counters，没有并发构建任务。

## 实测前后对比（秒）

| 项目 | 修改前平均 | 修改后平均 | 缩短 | 修改后范围 |
|---|---:|---:|---:|---:|
| 三个凭证签发合计 | 17.355199 | 10.218190 | 41.12% | 9.872333–10.396463 |
| Prover 首次准备 | 5.880202 | 3.286395 | 44.11% | 3.069814–3.529742 |
| Verifier 首次准备 | 1.997021 | 1.444547 | 27.66% | 1.411001–1.485043 |
| TA MemWitSync 生成 | 5.971942 | 3.438360 | 42.42% | 3.306314–3.518424 |
| MemWitSync 恢复及检查 | 7.127618 | 4.591347 | 35.58% | 4.472991–4.739953 |
| 整个回归测试程序 | 60.202880 | 45.588072 | 24.28% | 44.717450–46.026305 |

单凭证整体平均：5.785066 →
3.406063 s，缩短
41.12%。
这里是每轮三个 Issue 的合计再除以三个，不把注册/初始化计入单个 Issue。

| 轮次 | 单凭证平均 | Prover 准备 | Verifier 准备 | TA Sync 生成 | Sync 恢复及检查 |
|---|---:|---:|---:|---:|---:|
| run1 | 3.461924 | 3.529742 | 1.485043 | 3.490341 | 4.739953 |
| run2 | 3.290778 | 3.259628 | 1.411001 | 3.306314 | 4.472991 |
| run3 | 3.465488 | 3.069814 | 1.437598 | 3.518424 | 4.561096 |

全部源数据见 lifecycle_measurement.json 和每轮 measurement.json。

## 实现

1. **带证明界的 Gaussian 候选提前判定**：单独的 64-bit MPFR workspace 对同一
   candidate 的同一 256-bit uniform 做 RNDD/RNDU 区间判定。只有 uniform 严格落在
   acceptance interval 外才作提前接受/拒绝；重叠或 coarse interval 无效时进入原
   高精度 acceptance loop。原高精度 workspace、proposal support/radius、输入 center/
   variance、随机流和有理 Klein sampler 不变。原 workspace 不被较低精度覆盖。
2. **memoize 常量**：仅在输入 variance 的精确 dyadic 值与工作精度一致时复用 sigma/alpha；
   center 每次重新设置。精度变化或 rational sampler 切换会使 cache 失效。
   `cached_variance` 按 source 的 MPFR precision 精确复制，避免误把邻近 variance 当成相等。
3. **批量公参 CRT 展开**：11-byte SHAKE candidate 分为 32/32/24-bit words，以原始 q88
   模数作完整 lexicographic rejection，然后精确计算三个 CRT projection。相同 seed/domain
   和 stream 长度递增规则产生相同候选顺序。中间乘积小于 2^62，全程 uint64 无截断。
4. **只共享 public constants**：PP entry 内 prover/verifier 共享不可写 refresh-bank 和
   centered Atrap/Ba arrays。没有共享私有 witness/native statement ownership、Gaussian samples
   或 challenge-dependent 数据。proof profile ID、系数、native 参数和序列化 ABI 不变。

点 Gaussian 的支持范围和 Bernoulli 随机数均保留。使用 MPFR 有方向的正确舍入
构造区间的依据见 [MPFR 官方手册](https://mpfr.org/mpfr-current/mpfr.html)。上述“同一随机数
判定不变”是本实现的区间包围推导和参考对照结果，不是手册对本项目的安全认证。
旧的有限支持/Arb midpoint 等安全论证范围不会因这次加速自动扩大。

## 验证

- audit/show.log：12297422 次 point Gaussian 的同 candidate、同 uniform 参考高精度
  对照全部一致；该完整审计在后续 variance memoization 加入前运行，不能称为对最终代码
  所有 cache 条件的 12297422 次审计。
- gaussian_unit.log：最终 cache/screen 实现通过原 Gaussian + exact rational Klein 回归；
  point Gaussian 另做 3050 次 high-precision reference 对照，全一致；8-thread Klein 通过。
- audit/worker.log：完整 refresh-bank 和 Atrap/Ba CRT coefficient arrays 与原 scalar 路径
  逐项对照通过。新数组也用于正常三轮全部完整证明的生成、验证。
- run1/run2/run3：全部 Setup、Issue、明文 MemVer、完整 Pack Prove/Verify、wire 编码/解码、
  撤销拒绝、DELTA6 更新后重新完整 Show、整数 CRT lift 与 MemWitSync 回归通过。
- reload_old_proof.log：新 frontend 成功 reload/verify 修改前的生产 proof；证实相同 fixture 下
  public statement 和 profile binding 保持兼容。

## 其他时间及局限

TA Setup 本轮 4.045 s（对照 3.908 s）。
DELTA6 update 本轮 0.0371 s。
撤销前完整 packet 生成/验证为
3.163/
1.664 s。
更新后完整 packet 为 4.674/
1.781 s，较前一组三次
4.319/
1.618 s 慢；本次未单独定位差异原因，
不把本次 Issue/preparation/sync 的收益推广为每个 Show 都变快。
完整 Show 的 3 s/1.5 s 目标仍未达成。

完整回归约 45.59 s 包含三个 Issue、首次预热、两份 Show、LOCAL 回退、两次 TA Sync 生成、
全部反例及 cleanup，不是最短业务路径。Sync 恢复走本地认证 callback，不包含网络 RTT；
该阶段包含再次 TA 生成和 client 验签/安装/明文检查，不能与单独 TA Sync 生成相加解释为一次恢复。
每轮 public packet 的实际字节数见 measurement.json；这次算法没有增添 proof fields。

## 开关与恢复

`REVAC_GAUSSIAN_CERTIFIED_SCREEN=0` 关闭区间筛选；`=2` 开启同随机数参考对照。
`REVAC_GAUSSIAN_SCREEN_REPORT=1` 开启计数；正常无需开启。
`REVAC_Q88_PUBLIC_BANK_AUDIT=1` 开启原 public-expansion 比较。

旧源码备份：/home/wkk/pq_ac/production_backups/q88_pre_lifecycle_optimization_20261014/sources.tar.gz。
变动 Gaussian 代码后需要重建生产 C demo；Python 改动需重启 resident worker。
本地源码为准，单向同步 server；没有 PGO，没有 SDE，私有 fixtures 没有下载。
