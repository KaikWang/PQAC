# 生产 q88 compact CRT v3 替换与完整 Show packet 实测

日期：2026-10-13；服务器：autodl-pqa11；五次新 Show。
生产 Python compiler/worker 和 lowercase lazer native engine 已更新；
生产 C demo 已重编译。本次测量的 Python path 只有生产 zk_lazer 与 lazer/python，
没有实验 compiler path，也没有 LD_PRELOAD。新的公开配置由生产
`zk_lazer/revac_q88_profile.py` 默认启用。

## 完整 packet 时间与大小

生成时间包含 `revac_uav_show_prove`、匿名签名、完整 accumulator proof，以及
`revac_show_proof_pack_compact(..., deflate=1)` 的压缩序列化。
验证从收到完整 packet 到完成 `revac_gcs_verify_show_wire_guarded_parallel`；
包含 packet 解码、匿名签名、累加器证明、same-handle binding 和 guard 检查。
初始化、首次 statement/参数构造、Issue、撤销更新不计入在线 Show 时间。

| 次数 | 完整 packet 生成 s | wire 验证 s | 压缩序列化 s | packet bytes | packet KiB | accumulator bytes |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 3.792042 | 1.985664 | 0.039279 | 215253 | 210.208 | 193081 |
| 2 | 3.438231 | 1.854092 | 0.039820 | 222346 | 217.135 | 192889 |
| 3 | 4.002032 | 1.916744 | 0.033185 | 217382 | 212.287 | 193049 |
| 4 | 3.062824 | 1.913723 | 0.043597 | 217522 | 212.424 | 193081 |
| 5 | 3.028476 | 1.724339 | 0.042683 | 217555 | 212.456 | 193113 |

| 指标 | 平均 | 范围 |
|---|---:|---:|
| 完整 packet 生成 | 3.464721 s | 3.028476–4.002032 s |
| 完整 wire 验证 | 1.878912 s | 1.724339–1.985664 s |
| packet 序列化/压缩 | 0.039713 s | 0.033185–0.043597 s |
| 完整 packet | 212.902 KiB | 210.208–217.135 KiB |
| 仅 accumulator envelope | 188.518 KiB | 188.368–188.587 KiB |

累加器 worker 平均生成 3.345228 s、
验证 1.866084 s；不要与完整 packet 时间混用。
最后的 `last_show_packet.bin` 实际大小 217555 bytes；
其 `last_accumulator_envelope.bin` 大小 193113 bytes，
均与 C 日志逐字节计数一致。

此前独立优化实验的 3.239/1.669 s 是上一组随机新 Show 的测量；它的生成时间未包含
packet 压缩序列化。本轮新测量为 3.465/1.879 s，
不能用旧的较小结果替代生产接入后的结果。两组验证耗时也存在差异；本轮没有测定差异原因。
旧生产基线完整 Show 平均 5.312709/2.821800 s、完整压缩 packet 平均 214.598 KiB；
其生成计时也不含 packet 序列化，比较时需注意口径。当前五次仍未全部达成 3 s/1.5 s。

## 替换内容与安全边界

- 新默认 full-ring public norm/layout 取自已验证 dimension_v3，保留原三个 CRT 素数。
- weighted original witness beta0 implication、public window bound、auxiliary requirement
  preservation、scaled refresh Z 的六个 full-ring bindings 均保留；native 参数门通过。
- 完整 accumulator 主关系、整数 quotient/carry、CRT lift、hidden handle/refresh 产品绑定、
  q19 四行 link 与逐 bit shared handle binding 保留；没有修改 Gaussian 分布或 LNP 重复次数。
- native 使用 validated batch split、fused linear aggregation、六 prime 点积并行和
  bounded exact-size allocation pool；Pack reducer 保持串行。没有启用未选中的 zero-refresh
  或 Dachshund cost 选择实验。没有 PGO。
- 原 native build flags 与指令集需求保留；部署到其他 CPU 时须能执行对应 AVX512 指令。
- 新旧 statement layout 使用不同 static_profile_id；升级双方须一致并重启旧 resident workers。
  原凭证、CRT accumulator 参数和 private wire ABI 保留，但旧 layout proof 不应直接按新 statement 验证。
- 功能和参数门检查不是独立完整密码学安全认证；原系统的安全论证缺口不会因验证通过自动消失。

## 验证证据

`show.log`：五次完整新 Show 全部成功，完成签名/accumulator/same-handle wire 验证。
完整流程还运行了旧窗口、credential/revocation/sync 的现有回归检查。
`validation.log`：生产 compiler Clear Verify、native Pack Prove、序列化、反序列化、Pack Verify；
shared handle 翻转、非二进制 full handle、非 constant handle、scaled refresh Z 绑定篡改、
accumulator RHS/q19 RHS 篡改、截断/尾随 proof 均拒绝；恢复后有效 proof 通过。
`crossbackend.log`：新生产 proof 在升级前原生 engine 上验证。
`profile_mismatch.log`：旧 statement layout 拒绝新 layout 的 proof。
源码 hashes 保存于 source_sha256.txt。

## 复现与恢复

```bash
export PQAC_RESULTS_DIR=/path/to/new-results
bash /root/autodl-tmp/pq_ac/OURS/code/entities/tools/measure_q88_show.sh
```

runner 使用生产 worker/native；`PQAC_PYTHON`、`PQAC_DEMO` 可指定部署解释器和 C binary。
measurement.json、last_show_packet.bin、last_accumulator_envelope.bin 均是实际序列化证据。
私有 fixture 未下载到本地。

源码及旧本地库备份：`/home/wkk/pq_ac/production_backups/q88_pre_v3_20261013/baseline.tar.gz`。
服务器旧库：`/root/autodl-tmp/pq_ac/production_backups/q88_pre_v3_20261013/liblabrador38.so`。
备份含原文件 SHA256、Git HEAD 和恢复说明。清洁环境下设置
`REVAC_Q88_OPTIMIZED_PROFILE=0` 可选旧 compiler layout；完全回滚可恢复归档并重建 C demo。
