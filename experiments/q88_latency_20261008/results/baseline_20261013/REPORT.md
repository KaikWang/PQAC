# 生产 q88 CRT 基线实测（2026-10-13）

服务器：autodl-pqa11，Intel Xeon Gold 6430。当前基线源码在独立 build-baseline-20261013 目录重新编译，未修改基线代码。5 次完整 show 均通过验证。

| 项目 | 平均 | 最小–最大 |
|---|---:|---:|
| 完整 show 生成 | 5.313 s | 5.156–5.485 s |
| 完整 show wire 验证 | 2.822 s | 2.731–2.870 s |
| 累加器 worker 生成请求 | 5.244 s | 5.075–5.417 s |
| 累加器 worker 验证请求 | 2.811 s | 2.720–2.861 s |
| native Pack Prove | 5.231 s | 5.063–5.404 s |
| native Pack Verify | 2.797 s | 2.706–2.846 s |
| 累加器证明 envelope | 194.370 KiB | 194.000–194.500 KiB |
| 压缩完整 show packet | 214.598 KiB | 212.900–219.110 KiB |

时间为预热后在线时间；首次 UAV/GCS 预热分别 5.916 / 2.111 s，另计。

匿名签名和累加器生成并行。C 日志的“累加器 Labrador MemVer 证明”阶段是签名完成后的剩余等待时间，不能当作累加器完整生成耗时；此处使用 worker 请求完整计时。完整 show 验证包含反序列化及完整凭证验证；累加器 worker 验证包含实例重绑定、证明反序列化和 Pack Verify。

每轮 KiB 来自程序日志（保留两位小数）；最后一次实际序列化 envelope 为 199170 bytes。完整 show packet 包含匿名签名、累加器证明及其他协议字段。

原始日志：q88baseline_measure.log、q88baseline_measure.worker.log。机器可读结果：measurement.json。

## 复现入口

构建目标：OURS/code/entities 的 Release revocable_ac_demo，REVAC_LAZER_PYTHON=/root/miniconda3/envs/sage/bin/python。

使用 q88 专用 revac_acc96_memver_worker.py，并指定独立 Unix socket；运行 revocable_ac_demo 时设置 REVAC_LABRADOR_ONLINE=1、REVAC_LABRADOR_USE_SDE=0、REVAC_LABRADOR_Q96_SOCKET=/tmp/q88baseline_measure_20261013.sock、REVAC_Q96_SHOW_REPETITIONS=5。worker 使用生产默认参数，未启用单素数实验布局。
