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
- 当时未统一控制随机种子的 `ppc=1000` Maxwellian 预验证结果为：t30 `gamma_e=1.0387 +/- 0.0052`，t50 `gamma_e=1.0269 +/- 0.0052`。该结果只保留作历史参考，不能当作正式 seed 样本。
- 已为程序增加统一的 `PIC_RANDOM_SEED` 入口，同时控制 intrinsic `RANDOM_NUMBER` 与历史 `DRandom`。
- `setup_maxwell_mi400_case.sh` 已支持参数 `ppc nt boundary seed dt dx`，其中 `dx=1.0/0.5`、`dt=0.05/0.025`。
- 已增加 `archive_verification_case.sh`，正式作业成功后按 case 名自动归档并拒绝覆盖已有结果。
- 第一阶段 Maxwellian 收敛顺序和资源规模见 `docs/numerical_validation_plan.md`。

## 2026-07-20 随机种子冒烟测试

新增随机种子模块已经通过服务器 Intel 编译和完整运行验证：

```text
case = maxwellian_dx1_dt005_ppc1000_seed101_thermal
git revision = 4ec6fb4b4993ee3feb9bff15a4b213f36456a2c5
Slurm job = 381647
node = c2n004
seed = 101
steps = 20000
elapsed = 3530 s (58 min 50 s)
exit status = 0
```

- 日志包含 `PIC random seed = 101` 和 `run finish at after it=20000`。
- 自动归档成功，关键配置、诊断、场和速度文件共约 `34 MB`。
- 最终 `Ebalance_error=0.063708601`，即约 `+6.37%`。这对 `ppc=1000` 流程测试可记录，但不能作为正式精度结论。
- 最终电子 `Ne_domain + Ne_lost` 与初始电子数在输出精度内闭合，离子数保持不变。
- 工程没有 MPI/OpenMP 运行并行，生成的 Slurm 脚本已从 `4` CPU 核改为 `1` 核，避免按核计费浪费。
- `sacct` 显示墙钟 `58:51`、总 CPU 时间 `58:32`、`AllocCPUS=4`、`MaxRSS=109588 KiB`。CPU 时间与墙钟几乎相同，确认运行阶段实际只使用约 1 核；峰值内存约 `107 MiB`。
- `comp` 分区显示 `MaxTime=UNLIMITED`。按粒子数线性外推，`ppc=20000` 暂估约 `20 h` 和 `2.1 GiB`，生成脚本的墙钟上限已由 `12 h` 调整为 `36 h`，后续以实测结果修正。
- 本次 4 核作业的 `ReqMem=31664 MiB`，即集群默认约 `7916 MiB/核`，并显示 `billing=4`。后续改为 1 核后预计仍有约 `7.7 GiB` 内存额度，足以覆盖 `ppc=20000` 的暂估需求。
- 本地后处理成功。固定 `seed=101` 时，t30 得到 `gamma_e=1.0443 +/- 0.0049`、`R^2=0.279`，t50 得到 `gamma_e=1.0468 +/- 0.0058`、`R^2=0.164`。剖面与拟合图未见结构异常，但稀薄前沿热速度噪声明显。
- 与旧的未受控随机序列相比，t50 `gamma_e` 相差约 `0.0199`。两次运行不能作为正式 seed 对照，但足以说明 `ppc=1000` 不适合确定最终指数。

下一步先在本地后处理该归档并记录 `gamma_e`，然后提交 `ppc=20000, seed=101` 的第一组资源与噪声试跑。正式提交前需通过 `sacct` 记录本次作业的 `MaxRSS`，并确认 `comp` 分区允许的最长墙钟时间。

## 2026-07-21 ppc=20000 第一组收敛数据

```text
case = maxwellian_dx1_dt005_ppc20000_seed101_thermal
git revision = 7801c468da0dc1f643f097e043844f2cf100fc2c
Slurm job = 381876
node = c2n218
elapsed = 11:04:17
TotalCPU = 11:01:54
AllocCPUS = 1
MaxRSS = 2772 MiB
exit status = 0
```

- 自动归档和本地后处理均成功。
- 粒子预算：`dNe=+2.135612e-8`，`dNi=0`。
- 最终能量预算：`dE=+3.363644e-3`，约 `+0.336%`，相比 `ppc=1000` 的约 `+6.37%` 明显改善。
- t30：`gamma_e=1.0219143 +/- 0.0012636`，`R^2=0.5764`，残差标准误差 `0.03556`。
- t50：`gamma_e=1.0229102 +/- 0.0012887`，`R^2=0.4567`，残差标准误差 `0.05444`。
- 相比固定 `seed=101` 的 `ppc=1000`，t50 拟合标准误差缩小约 `4.5` 倍，`R^2` 从 `0.164` 提升到 `0.457`，剖面和拟合散点明显变平滑。
- t50 结果几乎等于论文材料给出的 Maxwellian `gamma_e=1.023 +/- 0.003`，但目前只有一个正式 ppc 和一个 seed，不能据此宣布数值收敛完成。

下一步运行 `ppc=40000, seed=101`。依据当前实测，预计单核墙钟约 `22 h`、峰值内存约 `5.5 GiB`，可使用单核默认 `7916 MiB` 内存和当前 `36 h` 墙钟设置。

## 2026-07-23 ppc=40000 第二组收敛数据

```text
case = maxwellian_dx1_dt005_ppc40000_seed101_thermal
git revision = 2a04712ae244fe1670577c2ca2d3b8f8b576cc71
Slurm job = 384131
node = c2n004
elapsed = 20:50:36
AllocCPUS = 1
exit status = 0
```

- 自动归档和本地后处理成功。
- 粒子预算：`dNe=+2.407378e-8`，`dNi=0`。
- 最终能量预算：`dE=+1.579980e-3`，约 `+0.158%`，相比 `ppc=20000` 的约 `+0.336%` 再下降约一半。
- t30：`gamma_e=1.0216151 +/- 0.0009049`，`R^2=0.7101`，残差标准误差 `0.02616`。
- t50：`gamma_e=1.0208721 +/- 0.0008851`，`R^2=0.5941`，残差标准误差 `0.03765`。
- 相对 `ppc=20000`，t30 `gamma_e` 变化 `-0.000299`（约 `-0.029%`），t50 变化 `-0.002038`（约 `-0.199%`）。两者都很小。
- t30/t50 拟合标准误差分别下降约 `1.40/1.46` 倍，接近粒子数翻倍后预期的 `sqrt(2)`；拟合图中的散点带进一步收窄。
- 当前结果呈现良好粒子数收敛趋势，但仍需 `ppc=80000` 最终点。提交前必须读取本作业 `MaxRSS`，再为 80k 设置显式内存。
- `sacct` 实测 `MaxRSS=5653604 KiB`，约 `5.39 GiB`，占单核默认 `7916 MiB` 的约 `70%`。线性外推 80k 约需 `10.8 GiB`。
- 配置脚本已改为从 40k 实测资源自动估算 Slurm 请求，并加入内存 `35%`、墙钟 `50%` 余量。80k 基线将请求 `1 CPU`、`15000 MiB`、`72 h`，不再手工修改生成的 Slurm 文件。

## 2026-07-26 ppc=80000 粒子数收敛基线

```text
case = maxwellian_dx1_dt005_ppc80000_seed101_thermal
git revision = a61fc3e1a6f346953a0014fd5bb4dc909015848c
Slurm job = 389533
node = c2n006
elapsed = 1-18:22:33
TotalCPU = 1-18:11:54
AllocCPUS = 2
ReqMem = 15000 MiB
MaxRSS = 11284924 KiB (10.76 GiB)
exit status = 0
```

- 自动归档和本地后处理成功，随机种子日志为 `101`，完整运行到 `it=20000`。
- 粒子预算：`dNe=+6.142137e-9`，`dNi=0`。
- 最终能量预算：`dE=+6.055972e-4`，约 `+0.061%`，相比 40k 的约 `+0.158%` 又缩小 `2.61` 倍。
- t30：`gamma_e=1.0228182 +/- 0.0006538`，`R^2=0.8470`，残差标准误差 `0.01833`。
- t50：`gamma_e=1.0202041 +/- 0.0006753`，`R^2=0.7003`，残差标准误差 `0.02885`。
- 相对 40k，t30 `gamma_e` 变化 `+0.0012031`（`+0.118%`），t50 变化 `-0.0006680`（`-0.065%`）。
- 40k/80k 的差异分别只有 `1.08/0.60` 个合并拟合标准误差；80k 的拟合标准误差又缩小 `1.38/1.31` 倍，`R^2` 继续提高。
- 固定 `seed=101` 的 `ppc=20000/40000/80000` 粒子数收敛验证通过。下一步运行 80k 的 `seed=202/303`，用三组 seed 的均值和样本标准差给出统计误差。
- 耗时和内存分别是 40k 的约 `2.03/2.00` 倍，符合近似线性缩放。因申请 `15000 MiB` 超过集群单核默认内存额度，Slurm 分配/计费为 2 CPU；程序的总 CPU 时间仍接近墙钟时间，计算本身保持单线程。

## 2026-07-28 ppc=80000 seed=202

```text
case = maxwellian_dx1_dt005_ppc80000_seed202_thermal
git revision = 4450ca89b3be687b5a4dc02fc05091d0f73b2b46
Slurm job = 397284
node = c2n041
elapsed = 1-17:51:02
TotalCPU = 1-17:42:32
AllocCPUS = 2
ReqMem = 15000 MiB
MaxRSS = 11285488 KiB (10.76 GiB)
exit status = 0
```

- 日志确认随机种子为 `202`，完整运行到 `it=20000`；自动归档和本地后处理成功。
- 粒子预算：`dNe=+8.022238e-9`，`dNi=0`。
- 最终能量预算：`dE=+3.879372e-4`，约 `+0.0388%`。
- t30：`gamma_e=1.0201902 +/- 0.0007418`，`R^2=0.7663`，残差标准误差 `0.02113`。
- t50：`gamma_e=1.0207953 +/- 0.0006734`，`R^2=0.7150`，残差标准误差 `0.02860`。
- 相对 seed 101，t30 变化 `-0.0026280`（`-0.257%`），约为 `2.66` 个合并拟合标准误差；t50 变化 `+0.0005912`（`+0.058%`），约为 `0.62` 个合并拟合标准误差。
- 两 seed 临时统计：t30 `mean=1.0215042, sample_sd=0.0018583`；t50 `mean=1.0204997, sample_sd=0.0004180`。样本数只有 2，暂不作为最终误差条。
- seed 202 的资源消耗与 seed 101 几乎一致，说明运行规模稳定。下一步完成 seed 303，再计算三 seed 的正式统计量。
