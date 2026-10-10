# 可撤销匿名凭证参数 profile

本目录把跨 C、Python、LaZer/LaBRADOR 和同步协议使用的参数集中为一个可审计的 manifest。手工修改的唯一来源是 JSON；C 头文件、Python 模块和审计报告均由脚本生成。

当前活动参数为 `revac128-candidate-v1.json`。它冻结了目前实验使用的精确数值和 profile ID，但仍是 **candidate**，不是已经完成安全证明的 128-bit 发布参数。未完成项会显示在生成的审计报告中，且 `--require-release` 会拒绝它。

当前接入实体流程的累加器参数覆盖为
`revac-q88-window50000-candidate-v2.json`（候选 D v2）：88-bit CRT
模数、24-bit handle、balanced base-26007 的 6-digit gadget，以及
50000 次撤销窗口。历史 `acc96_*` C/Python 符号仅作为兼容 API
保留；`RV88*` wire magic 保证旧 q96 载荷不会被误解码。运行

```bash
python3 scripts/params/check_q88_candidate.py
```

可同时检查 C 宏、Python wire ABI、窗口算术上限和结构尺寸。
该 profile 已经实现，但仍是参数候选；`beta_1 < q` 的条件算术
上限 56245 不等价于已完成组合安全证明。

`revac112pq-q96-candidate-v1.json` 是论文最终参数的独立候选清单，目标为 128-bit classical / 112-bit quantum、每窗口至少 500 次删除；它取代安全余量和关系规模均较差的 q80 搜索稿。它尚未完整接入生产算术，不能替代当前活动 manifest。运行

```bash
python3 scripts/params/check_integrated_security_candidate.py
python3 scripts/params/check_integrated_security_candidate.py --require-release
```

第一条只检查已有筛选证据是否自洽；第二条在所有证明、实现、实验和独立复核完成前必须失败。

## 检查与重新生成

在仓库根目录运行：

```bash
python3 scripts/params/generate_revac_profile.py
python3 scripts/params/generate_revac_profile.py --check
```

生成文件为：

- `code/entities/common/core/include/revac_profile_generated.h`
- `code/entities/zk_lazer/revac_profile_generated.py`
- 与 manifest 同名的 `*.generated.md` 审计报告

CMake 构建依赖 `revac_params_check`。如果 manifest 与生成文件不一致，构建会失败，而不是继续使用一部分旧参数。

## 后续修改规则

1. 不覆盖已经用于论文正式实验的 manifest；复制为新文件并修改 `profile_name`，例如 `revac128-candidate-v2.json`。
2. 用新 manifest 重新生成绑定：

   ```bash
   python3 scripts/params/generate_revac_profile.py \
     --manifest params/revac128-candidate-v2.json
   ```

3. 重新配置构建目录：

   ```bash
   cmake -S code/entities -B code/entities/build \
     -DREVAC_PARAMETER_MANIFEST="$PWD/params/revac128-candidate-v2.json"
   ```

4. 参数改变会产生新的 SHA3-256 profile ID。该 ID已绑定到撤销 context 公参哈希和 LaBRADOR proof profile，因此旧 context/proof 会被拒绝。
5. 凡是改变环、签名/证明、累加器 gadget、截断采样器或公共矩阵维度的 profile，必须重新运行 Setup 和 Issue，并重新生成凭证、累加器见证、context、证明以及论文实验数据。不要尝试把旧持久化对象解释成新 profile。

## 发布门槛

只有满足以下条件后，才可把 manifest 改为 `release_status: frozen`：

- Issue/Show 与标签友好采样器参数完成重新推导并能由后端实现；
- 平衡 base-3 累加器的正确性、范数增长、删除容量和安全参数完成推导；
- 同一句柄承诺的绑定/隐藏及组合安全证明完成；
- LaBRADOR grouped-right 范数预算由证明支持，不再是经验策略；
- 组合统计距离、查询界和所有安全估计达到论文声明；
- 完整回归、跨设备互操作和正式实验通过。

最后运行：

```bash
python3 scripts/params/generate_revac_profile.py --require-release --check
```

该命令成功才表示代码层的发布门槛已解除；它不能替代理论证明和人工审计。
