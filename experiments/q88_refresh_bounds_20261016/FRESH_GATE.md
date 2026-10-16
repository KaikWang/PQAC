# 生产新鲜见证检查：2026-10-16

用户要求：补齐刷新正确性推导中新鲜系数界 F 尚未强制检查的前提。

## 已实现的保证

`ACC96_FRESH_COEFF_BOUND = 2^30-1 = 1073741823` 现在定义在关系头文件中；证明编码的 `ACC96_WINDOW_FRESH_COEFF_BOUND` 引用同一个常量，避免两个定义漂移。

新增 `acc96_relation_fresh_witness_check`：

- 检查 left、v12 和 22 个 v2 多项式的每个中心系数都在 [-F,F]；
- 检查全部 right 和 refresh_z 系数为 -1、0、1；
- 检查 tag 是常数 1，且 head 的 CRT 表示规范；三元和 tag 比较也要求精确规范表示；
- 无效指针或未初始化上下文返回失败。

selector 保存为 16 字节，解释为恰好 128 个比特，不存在额外的非二进制编码值。

## 接入位置

| 路径 | 失败行为 |
|---|---|
| acc96_sampleleft_secure 返回前 | 返回 0，清零输出，无内部重采样 |
| acc96_add 的 sampler 回调后 | 返回 0，清零输出，不插入成员数据库 |
| revac_ta_export_memwitsync 回调后 | 返回失败，不生成响应 |
| revac_q96_sync_response_pack | 拒绝生成不满足新鲜界的签名响应，wire=NULL、len=0 |
| revac_q96_sync_response_unpack | 认证解码后、原 MemVer 之前检查；失败不覆盖接收方输出或 sequence |

没有裁剪系数，没有把超界候选改成别的见证，没有默认重试开关。即使调用者配置其他 sampler 回调，Add 与 TA 同步路径也不会绕过新鲜检查。

这道检查仅用于新签发/重新采样的同步见证。普通 MemVer、Show、MemWitUp/DELTA 和已更新见证的 window gate 不受此新鲜限制，否则会错误拒绝合法更新后的 right。

参数、proof 字段、序列化格式、50000 Delete 计数限制以及总范数门均未改动。新鲜检查是更严格的安装策略，不是附录 E 长期安全归约的替代品。

## 验证结果

本地 Release 全量构建通过。以下五项定向 CTest 通过：

1. `revac_acc96_fresh_gate_demo`：全部 24 个 head 多项式的 +F、-F 接受，±(F+1) 拒绝；right/z 三元边界、tag、非规范 CRT、空输入拒绝。一个原成员等式和 L2 门均接受的超 F 候选仍被 Add 拒绝，输出清零、成员表/sequence/accumulator 不变；相同结构的 F 边界候选可正常入库。大于 F 的候选仍可通过普通较大 window gate，证明两种界没有混用。
2. `revac_acc96_sampleleft_demo`：真实截断 sampler 产生的新见证可成功签发并通过成员验证。
3. `revac_acc96_lifecycle_demo`：更新后的成员验证正常；原窗口计数边界保持。该旧测试用元数据定位 50000/50001 边界，不是 50000 次真实撤销实验。
4. `revac_acc96_proof_layout_demo`：编码布局检查通过。
5. `revac_acc96_fresh_sync_demo`：真实 TA Setup、凭证签发、签名 MemWitSync 响应往返和 UAV 恢复安装成功，篡改签名被拒绝。还使用真实 TrapGen 关系构造一个原成员等式与 L2 门均接受、但超 F 的等价预像：直接打包和 TA 的自定义 sampler 回调均拒绝它；回调只调用一次，没有静默重采样；恢复原 sampler 后仍正常成功。

第五项的等价预像测试利用 trunc_b[0]=R0+A*R1：令 a=F+1，把 v2[0] 加 a、v12 加 R1*a、left 加 (R0-13^2)*a，原成员等式保持不变。这使测试确实区分新鲜门与原等式/L2 门，而非只用一个本来就无效的见证测试拒绝。

证据：`fresh_gate_validation/ctest.log`、`fresh_gate_validation/signed_sync.log`。

## 未完成的运行及环境边界

- `autodl-pqa11` 的现有 SSH 端口拒绝连接；没有声称代码已部署到服务器。
- 尝试本地完整 Show 流程，三次签发成功，但 native Labrador worker 在建立 statement 时触发 Illegal instruction。本地 CPU flags 缺少 AVX-512；没有运行 SDE 或更改 native 后端。失败日志保留在 `fresh_gate_validation/full_flow/`，本次不能声称 Pack/完整 Show 通过。
- 定向签名同步测试不依赖 native Labrador，已在上述硬件限制下完成真实签发、同步和安装路径验证。

## 对论文命题的影响

现在可把以下内容作为程序不变量：**凡由生产 sampler 返回成功、或经生产签发/同步安装路径接受的新鲜见证，均满足新鲜界检查。**因此 RESEARCH.md 第 3.1 节的确定性 B_fresh 对这些成功输出成立。

仍需分析输出失败概率与可用性。若理想采样分布为 P、接受事件为 E，成功输出分布为 P|E；不能说增加检查后成功输出分布完全不变。完整接口在失败时输出 abort，其与原接口的差异受 Pr[not E] 控制，未来若外部重试还需纳入查询与统计损失。当前小样本没有超界事件，不等于已证明其概率可忽略。

多次刷新后攻击者保留所有旧见证的安全归约仍是 RESEARCH.md 第 5–6 节的问题，本次没有把它标记为解决。

## 复现

```bash
cmake -S code/entities -B code/entities/build-ours
cmake --build code/entities/build-ours -j 4
ctest --test-dir code/entities/build-ours --output-on-failure \
  -R '^revac_acc96_(fresh_gate|fresh_sync|lifecycle|sampleleft|proof_layout)_demo$'
python experiments/q88_refresh_bounds_20261016/audit_bounds.py \
  --output experiments/q88_refresh_bounds_20261016/bounds_with_fresh_gate.json
```

修改前快照：`/home/wkk/pq_ac/production_backups/q88_pre_fresh_gate_20261016/`。
旧 `bounds.json` / `evidence.json` 保留研究初始状态；新数值状态见 `bounds_with_fresh_gate.json`，数值界没有变化，sampler 检查标记变为 true。
