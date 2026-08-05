# PIC-IFE_GEC 使用说明

这个仓库的核心程序在 `PIC-IFE_GEC/`，是一个 Fortran 2D IFE-PIC 程序，编译后生成 `PIC-IFE_GEC/1DPIC`。现在文档入口只保留本文件，换 IDE 或上服务器时照这里走。

## 目录和入口

| 路径 | 作用 |
|------|------|
| `PIC-IFE_GEC/` | CMake 工程根目录，编译和运行都在这里 |
| `PIC-IFE_GEC/INPUT/` | 程序运行输入文件 |
| `PIC-IFE_GEC/OUTPUT/` | 场、速度、能谱等输出 |
| `PIC-IFE_GEC/DUMP/` | 重启 dump 文件 |
| `PIC-IFE_GEC/run_min_case.sh` | 云端已有的最小算例脚本 |
| `scripts/setup_maxwell_mi400_case.sh` | 配置论文 Maxwell `mi/me=400` 证明算例 |
| `scripts/draw_png.py` | 从 `verification_runs/<case_name>/` 读取归档结果并画密度、电子热速度图 |
| `scripts/postprocess_maxwell_mi400.py` | 从归档结果生成论文后处理 CSV 和 `gamma_e` 拟合结果 |
| `scripts/summarize_maxwell_seed_ensemble.py` | 汇总多个 Maxwellian seed 的 `gamma_e` 均值、样本标准差和统计图 |
| `scripts/archive_verification_case.sh` | 按 `case_config.txt` 自动归档正式算例的配置、日志和关键输出 |
| `docs/numerical_validation_plan.md` | 按 PDF 工作建议整理的收敛性测试顺序和判定标准 |
| `docs/术语解释.md` | 物理符号、数值术语、诊断量、工作建议和汇报话术的通俗说明 |

程序用相对路径读写文件，所以运行 `1DPIC` 时工作目录必须是 `PIC-IFE_GEC/`。

## 换 IDE 怎么设置

打开仓库根目录 `/path/to/pic-`。IDE 里的关键设置如下：

| 设置项 | 值 |
|--------|----|
| CMake source directory | `/path/to/pic-/PIC-IFE_GEC` |
| CMake build directory | `/path/to/pic-/PIC-IFE_GEC/build` |
| Fortran compiler | 优先 `ifort`，也可用服务器上的 `ifx` |
| Executable | `/path/to/pic-/PIC-IFE_GEC/1DPIC` |
| Working directory | `/path/to/pic-/PIC-IFE_GEC` |
| Program arguments | 空 |
| Run before launch | `mkdir -p OUTPUT/Field OUTPUT/Velocity OUTPUT/Particle OUTPUT/Global OUTPUT/Phase OUTPUT/Energy OUTPUT/History OUTPUT/Average DUMP` |

远程 IDE 找不到 `ifort` 时，在 IDE 终端先执行：

```bash
module purge
module load intel/2022.1
module load cmake/3.23.5
which ifort
which cmake
```

没有 module 的机器可以试：

```bash
source /opt/oneapi/setvars.sh
```

## 先拉云端

每次改之前先同步云端：

```bash
cd /path/to/pic-
git fetch origin
git status --short --branch
git pull --ff-only
```

如果有本地未提交改动，先提交或 `git stash push -u`，不要直接覆盖。

## 最小算例

先用云端已有脚本确认代码能编译、能运行、能输出：

```bash
module purge
module load intel/2022.1
module load cmake/3.23.5

cd /path/to/pic-/PIC-IFE_GEC
bash ./run_min_case.sh
```

可改小算例步数：

```bash
NT=2000 DT=0.05 bash ./run_min_case.sh
```

成功标志：

- 日志出现 `run finish at after it= 1000`，或你设置的 `NT`
- `OUTPUT/Field/Average_x_001000.dat`
- `OUTPUT/Velocity/velocity_IJ_3001000.dat`
- `DUMP/var0001000dump`

## Maxwell mi/me=400 证明算例

论文里是 `mi/me=400`，不是 `me/mi=400`。脚本里离子质量设置为：

```text
mi = 400 * me = 3.6438D-28 kg
```

这个证明算例设置为：

| 项 | 值 |
|----|----|
| 初始分布 | Maxwellian |
| 计算域 | `[0, 1024 lambda_D0] x [0, 4 lambda_D0]` |
| 初始 plasma slab | `[0, 128 lambda_D0] x [0, 4 lambda_D0]` |
| 网格 | `1025 x 5`，即 `dx=dy=lambda_D0` |
| 时间步 | `dt=0.05 omega_pe^-1` |
| 目标时刻 | `omega_pi t = 50` |
| 左 x 边界 | thermal reservoir（默认）或 specular |
| 右 x 边界 | outflow / delete |
| y 边界 | periodic |
| 温度 | `Te0=1 eV`, `Ti0=0.01 eV` |
| 密度 | `1.0D21 m^-3` |

配置证明算例：

```bash
cd /path/to/pic-
bash scripts/setup_maxwell_mi400_case.sh 1000 20000 thermal 101 0.05 1.0
```

这里 `nt=20000, dt=0.05`，所以：

```text
omega_pi t = nt * dt / sqrt(mi/me) = 20000 * 0.05 / 20 = 50
```

编译和运行：

```bash
cd PIC-IFE_GEC
module purge
module load intel/2022.1
module load cmake/3.23.5

rm -rf build
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(which ifort)"
cmake --build build -j4

sbatch run_maxwellian_dx1_dt005_ppc1000_seed101_thermal.slurm
squeue -u dg001947
```

看运行日志：

```bash
tail -f run.log
```

跑完后关键输出应包括：

```bash
OUTPUT/Field/Average_x_020000.dat
OUTPUT/Velocity/velocity_IJ_3020000.dat
OUTPUT/physics_parameter.inp
```

## 画图和后处理

先把跑完的结果保存到 `verification_runs/<case_name>/`，再回到仓库根目录画图。画图脚本不再依赖 `PIC-IFE_GEC/OUTPUT/`：

```bash
cd /path/to/pic-
python3 scripts/draw_png.py verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000 12000
python3 scripts/draw_png.py verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000 20000
```

输出图片：

```bash
verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000/postprocessed/maxwell_mi400_t012000.png
verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000/postprocessed/maxwell_mi400_t020000.png
```

做论文后处理：

```bash
python3 scripts/postprocess_maxwell_mi400.py verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000
```

三个 seed 都完成单 case 后处理后，在本地汇总统计：

```powershell
python scripts\summarize_maxwell_seed_ensemble.py `
  verification_runs\maxwellian_dx1_dt005_ppc80000_seed101_thermal `
  verification_runs\maxwellian_dx1_dt005_ppc80000_seed202_thermal `
  verification_runs\maxwellian_dx1_dt005_ppc80000_seed303_thermal `
  --reference 1.023 --reference-error 0.003
```

脚本按 seed 等权计算 `gamma_e` 均值、样本标准差和均值标准误。论文结果的随机误差优先使用 `gamma_e_sample_standard_deviation`，不要用单次拟合标准误差替代跨 seed 波动。

生成文件统一放在：

```bash
verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000/postprocessed/
```

主要包括：

```text
profiles_t30.csv
profiles_t50.csv
gamma_fit_t30.csv
gamma_fit_t50.csv
gamma_fit_t30.png
gamma_fit_t50.png
postprocess_summary.txt
```

`gamma_e` 统一按照工作建议中的区间 `10^-3 <= ni/n0 < 1` 拟合；CSV 同时记录标准误差、`R^2`、空间范围、排除点数，profile 文件用 `gamma_fit_included` 标记每个点是否参与拟合。

如果服务器没有画图库：

```bash
python3 -m pip install --user matplotlib numpy
```

图里会有两部分：

- `n/n0` 随 `x/lambda_D0` 的电子和离子密度
- `v_th,e/v_th,e0` 随 `x/lambda_D0` 的电子热速度

## 论文级粒子数

论文参数是每个初始 slab 网格里 `80000` 个宏电子和 `80000` 个宏离子。初始 slab 是 `128 x 4` 个 cell，总粒子数约：

```text
128 * 4 * 80000 * 2 = 81.92 million
```

这会非常吃内存和时间。先用 `1000/cell` 跑通流程，再按 [数值验证计划](numerical_validation_plan.md) 依次跑 `20000/40000/80000 ppc`，不要直接提交最大算例。

```bash
cd ~/pic-
bash scripts/setup_maxwell_mi400_case.sh 20000 20000 thermal 101 0.05 1.0
cd PIC-IFE_GEC
sbatch run_maxwellian_dx1_dt005_ppc20000_seed101_thermal.slurm
```

## Slurm 示例

配置脚本会在 `PIC-IFE_GEC/` 自动生成 Slurm 文件，不需要手写。文件名包含所有关键参数，例如：

```text
run_maxwellian_dx1_dt005_ppc1000_seed101_thermal.slurm
```

提交后用 `squeue` 查看状态；任务正常结束会自动归档，不要同时提交两个使用同一 `PIC-IFE_GEC` 工作目录的 case。

```bash
sbatch run_maxwellian_dx1_dt005_ppc1000_seed101_thermal.slurm
squeue -u dg001947
```

## 常见问题

| 现象 | 处理 |
|------|------|
| `cannot open ./INPUT/mesh.inp` | 工作目录错了，先 `cd PIC-IFE_GEC` 再运行 |
| `cmake` 找不到 `ifort` | 先 `module load intel/2022.1`，再检查 `which ifort`；不要用带 Intel 参数的工程强制切到 `gfortran` |
| 输出没有 `020000` | 确认 `INPUT/pic.inp` 里是 `20000, 0.05`，或重新运行配置脚本 |
| 画图脚本找不到文件 | 先确认结果已经保存到 `verification_runs/<case_name>/`，并且里面有 `Average_x_*.dat` 和 `velocity_IJ_3*.dat` |
| 论文级粒子数跑不动 | 先用 `1000/cell` 或 `5000/cell` 做证明图，再申请更多内存和墙钟时间 |

## 最短命令清单

```bash
cd /path/to/pic-
git pull --ff-only
bash scripts/setup_maxwell_mi400_case.sh 1000 20000 thermal 101 0.05 1.0

cd PIC-IFE_GEC
module purge
module load intel/2022.1
module load cmake/3.23.5
rm -rf build
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(which ifort)"
cmake --build build -j4
sbatch run_maxwellian_dx1_dt005_ppc1000_seed101_thermal.slurm
squeue -u dg001947

# 跑完并自动归档后，在仓库根目录后处理：
cd ..
case_dir=verification_runs/maxwellian_dx1_dt005_ppc1000_seed101_thermal
python3 scripts/draw_png.py "$case_dir" 12000
python3 scripts/draw_png.py "$case_dir" 20000
python3 scripts/postprocess_maxwell_mi400.py "$case_dir"
```

## 输出文件和结果管理

核心原则：`PIC-IFE_GEC/OUTPUT/` 只保存当前这一次运行的临时结果；真正要留作论文材料的结果，统一复制到 `PIC-IFE_GEC/../verification_runs/<case_name>/`。

### 常见根目录文件

在 `PIC-IFE_GEC/` 下面经常会看到这些文件：

| 文件或目录 | 含义 | 是否长期保留 |
|---|---|---|
| `1DPIC` | 编译出来的可执行程序 | 可保留，重新编译会覆盖 |
| `build/` | CMake 编译目录 | 可删，重新 `cmake` 会生成 |
| `case_config.txt` | 当前 case 的配置摘要 | 保存到 `verification_runs` |
| `run.log` | 当前运行日志 | 保存一份，空间紧可 `gzip` |
| `run_maxwellian_*.slurm` | 配置脚本按参数生成的 Slurm 提交脚本 | 可保留 |
| `OUTPUT/` | 当前运行输出 | 不直接长期堆在这里 |
| `DUMP/` | 重启 dump 文件 | 不做重启时可删 |
| `MCCB*.dat`, `SN*.dat`, `SR*.dat`, `10000*.dat` | MCC/碰撞模块生成的中间表或统计文件 | 一般可删 |
| `*.maxwell_mi400.bak` | 配置脚本第一次运行时保存的源码/输入备份 | 先保留 |

### OUTPUT 目录里各文件的作用

| 路径 | 作用 |
|---|---|
| `OUTPUT/global_diagnostics.csv` | 全局守恒诊断，包含粒子数、动能、场能、注入/损失能量、`Ebalance_error` |
| `OUTPUT/Field/Average_x_012000.dat` | `omega_pi t = 30` 的 x 方向平均场/密度剖面 |
| `OUTPUT/Field/Average_x_020000.dat` | `omega_pi t = 50` 的 x 方向平均场/密度剖面 |
| `OUTPUT/Velocity/velocity_IJ_3012000.dat` | `omega_pi t = 30` 的速度/热速度剖面 |
| `OUTPUT/Velocity/velocity_IJ_3020000.dat` | `omega_pi t = 50` 的速度/热速度剖面 |
| `OUTPUT/physics_parameter.inp` | 本次运行输出的物理参数 |
| `OUTPUT/normalize.inp` | 本次运行使用的归一化参数 |
| `OUTPUT/CellVolume.dat` | 网格体积/面积相关输出 |
| `OUTPUT/ElectronChange.dat`, `OUTPUT/IonChange.dat` | 电子/离子数量变化记录 |
| `OUTPUT/PartcountReal.dat` | 实际粒子计数记录 |
| `OUTPUT/Energy/` | 能谱或能量相关输出 |
| `OUTPUT/Particle/` | 粒子采样输出，通常较大 |
| `OUTPUT/Phase/` | 相空间输出 |
| `OUTPUT/History/`, `OUTPUT/Global/`, `OUTPUT/Average/` | 历史量、全局量、平均量输出 |

### 保存一个跑完的 case

以标准小粒子数基准 `maxwellian_mi400_thermal_ppc1000_nt20000` 为例：

```bash
cd ~/pic-/PIC-IFE_GEC

case_dir=../verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000
mkdir -p "$case_dir"

cp case_config.txt INPUT/pic.inp MCC_jw/input/controlflow.txt run.log "$case_dir"/
cp OUTPUT/global_diagnostics.csv OUTPUT/physics_parameter.inp OUTPUT/normalize.inp "$case_dir"/
cp OUTPUT/Field/Average_x_012000.dat OUTPUT/Field/Average_x_020000.dat "$case_dir"/
cp OUTPUT/Velocity/velocity_IJ_3012000.dat OUTPUT/Velocity/velocity_IJ_3020000.dat "$case_dir"/

ls -lh "$case_dir"
```

如果 `run.log` 很大，可以压缩：

```bash
gzip "$case_dir/run.log"
```

### 从 verification_runs 画图和生成后处理结果

保存后的 case 目录是扁平结构，例如：

```text
verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000/
  Average_x_012000.dat
  Average_x_020000.dat
  velocity_IJ_3012000.dat
  velocity_IJ_3020000.dat
  global_diagnostics.csv
  physics_parameter.inp
  normalize.inp
  case_config.txt
```

其中 `12000` 对应文件名里的 `012000`，`20000` 对应 `020000`；`velocity_IJ_3******.dat` 里面开头的 `3` 是原程序固定前缀，不是步数的一部分。

画图：

```bash
cd ~/pic-
python3 scripts/draw_png.py verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000 12000
python3 scripts/draw_png.py verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000 20000
```

生成论文后处理 CSV：

```bash
python3 scripts/postprocess_maxwell_mi400.py verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000
```

所有生成物都放在：

```text
verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000/postprocessed/
```

不要把 `postprocessed/` 里的文件再混回 `OUTPUT/`；`OUTPUT/` 跑完保存后可以清理。

### 清理当前运行垃圾

确认结果已经保存到 `verification_runs/` 之后，再清理当前工作目录：

```bash
cd ~/pic-/PIC-IFE_GEC

rm -f run.log slurm-*.out slurm-*.err test_*.out
rm -f MCCB*.dat SN*.dat SR*.dat 10000*.dat
rm -rf OUTPUT DUMP
mkdir -p OUTPUT/{Field,Velocity,Particle,Global,Phase,Energy,History,Average} DUMP
```

不要删这些：

```text
INPUT/
MCC_jw/
code/
CMakeLists.txt
case_config.txt
run_maxwellian_*.slurm
*.f90
```

### 下次重新跑标准 Maxwellian case

```bash
cd ~/pic-
git pull --ff-only

bash -n scripts/setup_maxwell_mi400_case.sh
bash scripts/setup_maxwell_mi400_case.sh 1000 20000 thermal 101 0.05 1.0

cd PIC-IFE_GEC
module purge
module load intel/2022.1
module load cmake/3.23.5

rm -rf build
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(which ifort)"
cmake --build build -j4

sbatch run_maxwellian_dx1_dt005_ppc1000_seed101_thermal.slurm
squeue -u dg001947
```

### 并行运行不同 case

同一个 `PIC-IFE_GEC` 目录中的可执行文件、`INPUT`、`OUTPUT`、`run.log` 和锁文件都是共享的，因此不能在该目录中并行运行两个 case。确需并行时，为新 case 创建独立 Git worktree，并把归档根目录指回主仓库：

```bash
cd ~
git -C ~/pic- fetch origin main
git -C ~/pic- worktree add --detach ~/pic-dx05 origin/main

cd ~/pic-dx05
PIC_ARCHIVE_ROOT=~/pic-/verification_runs \
  bash scripts/setup_maxwell_mi400_case.sh 80000 20000 thermal 101 0.05 0.5
```

随后只在 `~/pic-dx05/PIC-IFE_GEC` 中编译和提交该 case。生成的 Slurm 文件会使用 worktree 的绝对路径，结果仍归档到 `~/pic-/verification_runs/<case_name>/`。不同 worktree 可以并行；同一个 worktree 仍只能运行一个 case。

任务正常结束后，脚本会自动把配置、日志、全局诊断以及 `t=30/50` 的场和速度输出归档到：

```text
verification_runs/maxwellian_dx1_dt005_ppc1000_seed101_thermal/
```

查看运行：

```bash
tail -f run.log
```

退出查看用 `Ctrl+C`，不会停止 Slurm 任务。真正停止任务才用：

```bash
scancel <JOBID>
```

### 跑完后检查

```bash
cd ~/pic-/PIC-IFE_GEC

squeue -u dg001947
grep -E "run finish|Saving|Total Number of Populations|nt, dt" run.log | tail -30
tail -5 OUTPUT/global_diagnostics.csv
ls OUTPUT/Field/Average_x_012000.dat OUTPUT/Field/Average_x_020000.dat
ls OUTPUT/Velocity/velocity_IJ_3012000.dat OUTPUT/Velocity/velocity_IJ_3020000.dat
```

标准两物种 case 的日志里应该看到 `Saving Electron` 和一个离子物种，例如 `Saving He(+)`。这里 `He(+)` 只是 MCC/gas 数据库里的显示名，标准 case 里质量和电荷已经由 `INPUT/pic.inp` 覆盖成 `mi/me=400, Z=1`。不应该再出现 `Ar(+)`。
