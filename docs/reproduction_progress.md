# 复现进度记录

更新时间：2026-05-11

## 本地仓库

- 当前分支：`main`
- 本轮目标：准备在电光云上复现师兄论文里的 collisionless plasma expansion 算例。
- 重点算例：Maxwellian 初始速度分布，`mi/me = 400`，目标时刻 `omega_pi t = 50`。
- 对应脚本：`scripts/setup_maxwell_mi400_case.sh`
- 画图脚本：`scripts/draw_png.py`
- 后处理脚本：`scripts/postprocess_maxwell_mi400.py`

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

保存到 `verification_runs` 后画图：

```bash
cd /data/home/dg001947/pic-
python3 scripts/draw_png.py verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000 20000
```

预期图片：

```bash
verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000/postprocessed/maxwell_mi400_t020000.png
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

## 2026-07-08 当前进度

已经完成标准 Maxwellian 小粒子数基准测试：

```text
case = maxwellian_mi400_thermal_ppc1000_nt20000
distribution = Maxwellian
left boundary = thermal-reservoir
species = electron + ion
mi/me = 400
particles per cell per species = 1000
nt = 20000
dt = 0.05 / omega_pe
target time = omega_pi t = 50
```

服务器保存目录：

```bash
~/pic-/verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000
```

该目录已保存：

```text
case_config.txt
pic.inp
controlflow.txt
global_diagnostics.csv
physics_parameter.inp
normalize.inp
Average_x_012000.dat
Average_x_020000.dat
velocity_IJ_3012000.dat
velocity_IJ_3020000.dat
run.log
```

诊断结论：

- Slurm 任务正常结束，`run finish at after it=20000`。
- 日志中只出现 `Electron` 和一个离子物种 `He(+)`，没有 `Ar(+)`，说明已经是两物种标准测试。
- `He(+)` 只是 MCC/gas 数据库显示名；标准 case 中质量和电荷由 `INPUT/pic.inp` 覆盖为 `mi/me=400, Z=1`。
- 粒子数守恒正常：`Ne_domain + Ne_lost` 基本回到初始电子数。
- `Ebalance_error` 在 `ppc=1000` 小粒子数测试下最终约为 6%，可作为流程验证；正式论文精度需要提高 ppc 并做收敛性。

已新增本地后处理入口：

```bash
python3 scripts/draw_png.py verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000 12000
python3 scripts/draw_png.py verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000 20000
python3 scripts/postprocess_maxwell_mi400.py verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000
```

这些脚本会读取 `verification_runs/<case_name>/` 中的扁平归档文件，并把生成物写入：

```text
verification_runs/<case_name>/postprocessed/
```

下一步：

1. 把本地脚本同步到服务器，在已完成的 `maxwellian_mi400_thermal_ppc1000_nt20000` 目录上运行画图和后处理。
2. 检查 `postprocessed/profiles_t30.csv`、`profiles_t50.csv`、`gamma_fit_t30.csv`、`gamma_fit_t50.csv` 是否正常。
3. 确认后处理无误后，再提高粒子数或做 `dx/dt/ppc` 收敛性。

## 2026-07-14 后处理与收敛准备

- 已从 `pdf材料/工作建议.pdf` 原始页面核对统一拟合流程，固定区间为 `10^-3 <= ni/n0 < 1`，区间由离子密度 `ni` 判定，拟合关系仍为 `Te(ne)`。
- 已修复速度文件七列记录被 Intel Fortran 自动换行造成的读取失败，并修正热速度、粒子计数和 cell 坐标列映射。
- 后处理现在输出 `gamma_e` 标准误差、残差标准误差、`R^2`、拟合空间范围、排除点数、逐点纳入标记和拟合散点图。
- 当前 `ppc=1000` Maxwellian 预验证结果为：t30 `gamma_e=1.0387 +/- 0.0052`，t50 `gamma_e=1.0269 +/- 0.0052`。该结果接近论文值，但不能替代正式粒子数收敛。
- 已为程序增加统一的 `PIC_RANDOM_SEED` 入口，同时控制 intrinsic `RANDOM_NUMBER` 与历史 `DRandom`。
- `setup_maxwell_mi400_case.sh` 已支持参数 `ppc nt boundary seed dt dx`，其中 `dx=1.0/0.5`、`dt=0.05/0.025`。
- 已增加 `archive_verification_case.sh`，正式作业成功后按 case 名自动归档并拒绝覆盖已有结果。
- 第一阶段 Maxwellian 收敛顺序和资源规模见 `docs/numerical_validation_plan.md`。
