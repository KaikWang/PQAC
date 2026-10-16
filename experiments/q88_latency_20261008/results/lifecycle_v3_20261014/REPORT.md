# q88 生产优化版全流程复测

日期：2026-10-14（Asia/Shanghai）；服务器 autodl-pqa11；Intel Xeon Gold 6430。
按照服务器恢复后的新请求，重新完成三次独立运行，不复用上一日测量数值。
每次重新创建 TA、3 个 UAV、1 个 GCS，并新签发三个凭证；使用当前生产 compact CRT v3
worker/compiler 与 native liblabrador38.so，未使用实验 PYTHONPATH 或 LD_PRELOAD。
C 测试源码重新编译，并明确指定 build-baseline-20261013/revocable_ac_demo，
该目录里的 binary 已更新为最新生产代码；目录名不代表运行旧生产方案。

## 三轮计时（单位：秒）

| 阶段及范围 | 平均 | 最小–最大 |
|---|---:|---:|
| 运行时初始化（arith/对象/限流器） | 0.001732 | 0.001448–0.002262 |
| TA Setup（包含累加器公参及 trapdoor） | 3.907941 | 3.804714–3.984844 |
| 3 UAV + 1 GCS KeyGen（含确定性测试句柄选择） | 0.007837 | 0.006041–0.010470 |
| 3 UAV + 1 GCS 注册 | 0.000131 | 0.000108–0.000143 |
| UAV1 凭证签发 | 5.885059 | 5.798415–6.002665 |
| UAV2 凭证签发 | 5.801065 | 5.696680–5.936135 |
| LOCAL 测试 UAV 凭证签发 | 5.669076 | 5.634328–5.691864 |
| 三个凭证签发合计 | 17.355199 | 17.242640–17.454081 |
| UAV prover 初次准备/预热 | 5.880202 | 5.649175–6.047403 |
| GCS verifier 初次准备/预热 | 1.997021 | 1.970453–2.041240 |
| 撤销前完整 Show packet 生成（含压缩序列化） | 3.290934 | 3.271961–3.311827 |
| 撤销前完整 Show wire 验证 | 1.717721 | 1.685338–1.758726 |
| TA T4 预分发数据生成 | 0.469127 | 0.458827–0.475439 |
| UAV T4 验签及缓存 | 0.005468 | 0.005031–0.005754 |
| TA Delete + DELTA6 构造 + 兼容 T4 | 1.154831 | 1.138088–1.174872 |
| TA 更新 packet/manifest 导出（主要复用已缓存数据） | 0.000082 | 0.000069–0.000104 |
| UAV DELTA6 验签与选择性见证更新 | 0.037369 | 0.036699–0.038590 |
| 无 T4 缓存 LOCAL 更新 + 明文 MemVer | 1.256142 | 1.203111–1.301233 |
| 撤销更新后完整 Show packet 生成（含压缩序列化） | 4.318946 | 4.054742–4.816634 |
| 撤销更新后完整 Show wire 验证 | 1.617646 | 1.570868–1.706326 |
| TA MemWitSync 重新采样、签名与导出 | 5.971942 | 5.915843–6.000548 |
| UAV MemWitSync 恢复（含 TA 再采样、验签、安装及明文检查） | 7.127618 | 6.936232–7.410581 |
| 完整测试程序 wall time（含全部反例、额外路径及 cleanup） | 60.202880 | 59.668326–61.162375 |

三个凭证总签发平均 17.355199 s，单凭证总体平均
5.785066 s。

签发后完整 Show 为 3.291/
1.718 s；撤销更新后为
4.319/
1.618 s。
两者包含匿名签名、完整 accumulator proof、same-handle link、压缩编码与 wire 验证。
更新后生成明显慢于本轮撤销前生成；本次没有剖析确定其原因，不能直接归因于见证大小或 LNP 重试。
当前数据没有达成 3 s / 1.5 s 目标。

## 完整 Show packet 实际字节数

| 轮次 | 撤销前 bytes | 撤销前 KiB | 更新后 bytes | 更新后 KiB |
|---|---:|---:|---:|---:|
| run1 | 215240 | 210.195 | 217298 | 212.205 |
| run2 | 215377 | 210.329 | 217335 | 212.241 |
| run3 | 217163 | 212.073 | 217253 | 212.161 |

撤销前完整 packet 平均 210.866 KiB；
更新后平均 212.202 KiB。
各轮原始 measurement.json 另列实际 accumulator envelope 字节数。
last_show_packet.bin 是撤销前的完整公开 wire artifact；更新后的 packet 大小来自 C 对
完整序列化缓冲区 show_wire_len 的直接计数，并通过完整 wire verification。

## 计时解释

- 主流程为真实 q88；代码中的 acc96/q96 是保留的兼容命名。
- TA Setup 含默认初始化/公参及 trapdoor 认证；运行时初始化、实体 KeyGen/注册另列。
- Issue 使用生产的区间认证扰动与精确有理 Klein SampleLeft；表中的单次 Issue
  是一个凭证调用，三个凭证合计行是前三个 Issue 的合计，不重复加入总时间。
- prover/verifier 预热独立记录：不能把预热当成每次 Show 在线延迟，也不能漏掉首次部署成本。
- 删除时间包含融合 DELTA6 层与兼容 T4 的构造；T4 prefetch 是撤销前生成的数据，单独计时。
- DELTA6 命中后见证更新与无缓存 LOCAL 回退是不同分支，本程序各测一次。
- MemWitSync 通过同进程认证 TA callback 模拟请求；这些时间包含计算与验签，没有真实网络 RTT。
- 测试为篡改拒绝先直接导出一次 sync，再经 recovery callback 再生成一次 sync；
  因而完整程序 wall time 包含两次同步生成，不能当成仅一次业务恢复的端到端延迟。
- 完整程序约 60.20 s 还包含全部负例、公开输入往返、LOCAL 回退、精确 lift 与 cleanup；
  它不是单个匿名凭证的在线生成/验证时间，也不是最短业务路径的时间。

## 验证状态

三轮 C demo 和 measurement summarizer 均返回成功。
每轮撤销前、更新后两份完整 Show 都通过 signature/accumulator/shared-handle wire verification。
已撤销用户拒绝、窗口/参数摘要绑定、signed manifest/DELTA6 的缺层/错序/跨序列/篡改拒绝、
精确 CRT lift、MemWitSync 的签名和 nonce 重放拒绝、修复损坏 witness 后的明文 MemVer 均通过。
这些是功能与参数门证据，不是新增的独立密码学安全认证。

## 证据和复现

- lifecycle_measurement.json：三轮阶段均值/范围及每轮原始阶段数值。
- run1/run2/run3 中的 show.log、worker.log、measurement.json：阶段与证明详细计时。
- 每轮 last_show_packet.bin / last_accumulator_envelope.bin：公开序列化结果。
- source_sha256.txt：本次计时相关源码与 production profile 的 SHA256。

```bash
export PQAC_DEMO=/root/autodl-tmp/pq_ac/OURS/code/entities/build-baseline-20261013/revocable_ac_demo
export PQAC_REPETITIONS=1
export PQAC_RESULTS_DIR=/path/to/fresh-lifecycle-run
bash /root/autodl-tmp/pq_ac/OURS/code/entities/tools/measure_q88_show.sh
```

三个独立 results/runN 目录完成后运行 tools/summarize_q88_lifecycle.py 生成统计。
没有从远端覆盖本地源码；未下载私有 prover fixture；没有使用 SDE 或 PGO。
