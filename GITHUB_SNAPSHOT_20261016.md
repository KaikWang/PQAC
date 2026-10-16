# q88 工作区快照（2026-10-16）

本地代码已分仓记录在同一个 GitHub 仓库的不同分支中：
https://github.com/KaikWang/PQAC

| 目录 | 分支 | 固定版本 |
| --- | --- | --- |
| OURS | experiment/q88-proof-size-bridge | 本文件所在提交 |
| lazer | backup/pqac-20261010-lazer | `674238831c0eaf8e77b20644746ea9a89819b69a` |
| lazer/src/labrados | backup/pqac-20261010-labrados | `e5a8a8a2e367dae5bcb6afde704df30ac3f8c43f` |

## 恢复到新的空目录

```bash
git clone --branch experiment/q88-proof-size-bridge git@github.com:KaikWang/PQAC.git OURS
git clone --branch backup/pqac-20261010-lazer git@github.com:KaikWang/PQAC.git lazer
git -C lazer checkout 674238831c0eaf8e77b20644746ea9a89819b69a
git -C lazer submodule init
git clone --no-checkout git@github.com:KaikWang/PQAC.git lazer/src/labrados
git -C lazer/src/labrados checkout e5a8a8a2e367dae5bcb6afde704df30ac3f8c43f
git -C lazer submodule status
```

修改过的 labrados 提交保存在 PQAC，而非上游仓库；直接执行上游 submodule update 可能无法下载该提交，故此处显式从 PQAC 恢复。源码、已跟踪配置、当前未忽略的实验代码和结果包含在快照中；被既有 .gitignore 排除的本地构建产物与缓存不属于该 Git 快照。记录本快照不代表新增实验或安全验证已经完成。
