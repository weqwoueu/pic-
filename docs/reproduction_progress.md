# 复现进度记录

更新时间：2026-05-11

## 本地仓库

- 当前分支：`main`
- 本轮目标：准备在电光云上复现师兄论文里的 collisionless plasma expansion 算例。
- 重点算例：Maxwellian 初始速度分布，`mi/me = 400`，目标时刻 `omega_pi t = 50`。
- 对应脚本：`scripts/setup_maxwell_mi400_case.sh`
- 画图脚本：`scripts/plot_maxwell_mi400.py`

## 云端状态

电光云 SSH 已配置成本地别名：

```bash
ssh dghpc
```

云端工作目录：

```bash
/data/home/dg001947/pic-
```

已完成：

- 当前仓库代码已上传到云端。
- 已执行 Maxwell `mi/me=400` 证明算例配置。
- 已加载 `intel/2022.1` 和 `cmake/3.23.5`。
- `PIC-IFE_GEC/1DPIC` 已在云端编译成功。

云端可执行文件：

```bash
/data/home/dg001947/pic-/PIC-IFE_GEC/1DPIC
```

云端 Slurm 脚本：

```bash
/data/home/dg001947/pic-/run_maxwell_mi400.slurm
```

## 当前卡点

代码已经能编译，当前卡在电光云账号的 Slurm 队列权限。

测试过这些分区：

```bash
comp
dg83
dg9a
zq0009
```

提交时都报：

```text
User's group not permitted to use this partition
```

结论：这不是代码问题，是账号 `dg001947` 还没有被绑定到可提交作业的计算分区/计费策略。

需要联系电光云管理员：

```text
账号 dg001947 可以登录，也能看到 comp/dg83/dg9a/zq0009 队列，
但是 sbatch 提交时报：
User's group not permitted to use this partition。
请帮我开通可提交的 Slurm 分区/绑定计费策略。
```

## 权限开通后的运行命令

```bash
ssh dghpc
cd /data/home/dg001947/pic-
sbatch -p comp run_maxwell_mi400.slurm
squeue -u dg001947
tail -f PIC-IFE_GEC/run.log
```

如果管理员指定其他可用分区，把 `comp` 换成对应分区，例如：

```bash
sbatch -p dg83 run_maxwell_mi400.slurm
```

## 跑完后检查

```bash
ls PIC-IFE_GEC/OUTPUT/Field/Average_x_020000.dat
ls PIC-IFE_GEC/OUTPUT/Velocity/velocity_IJ_3020000.dat
```

画图：

```bash
cd /data/home/dg001947/pic-
python3 scripts/plot_maxwell_mi400.py PIC-IFE_GEC 20000
```

预期图片：

```bash
PIC-IFE_GEC/figures/maxwell_mi400_t020000.png
```

## 本轮修正

云端编译时发现两个生成算例后的 Fortran 问题，本轮已固化到 `scripts/setup_maxwell_mi400_case.sh`：

- `OUTPUT_velocity.f90` 使用 `nt` 时需要 `Use TimeControl, Only: nt`。
- `Output_Energy.f90` 的旧 `WRITE(543, 17)`/`17 FORMAT` 组合在两物种输出改写后需要清理。

后续换电脑时，重新拉仓库后优先验证：

```bash
bash scripts/setup_maxwell_mi400_case.sh 1000 20000
cd PIC-IFE_GEC
rm -rf build
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(which ifort)"
cmake --build build -j"$(nproc)"
```
