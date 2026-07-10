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
module load compiler/intel/2021.3.0
module load compiler/cmake/3.23.3
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
git pull --rebase
```

如果有本地未提交改动，先提交或 `git stash push -u`，不要直接覆盖。

## 最小算例

先用云端已有脚本确认代码能编译、能运行、能输出：

```bash
module load compiler/intel/2021.3.0
module load compiler/cmake/3.23.3

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
| 左 x 边界 | specular / full reflection |
| 右 x 边界 | outflow / delete |
| y 边界 | periodic |
| 温度 | `Te0=1 eV`, `Ti0=0.01 eV` |
| 密度 | `1.0D21 m^-3` |

配置证明算例：

```bash
cd /path/to/pic-
bash scripts/setup_maxwell_mi400_case.sh 1000 20000
```

这里 `nt=20000, dt=0.05`，所以：

```text
omega_pi t = nt * dt / sqrt(mi/me) = 20000 * 0.05 / 20 = 50
```

编译和运行：

```bash
cd PIC-IFE_GEC
module load compiler/intel/2021.3.0
module load compiler/cmake/3.23.3

rm -rf build
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(which ifort)"
cmake --build build -j"$(nproc)"

mkdir -p OUTPUT/Field OUTPUT/Velocity OUTPUT/Energy OUTPUT/Average DUMP
./1DPIC > run.log 2>&1
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
postprocess_summary.txt
```

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

这会非常吃内存和时间。先用 `1000/cell` 跑通证明流程，再跑论文级：

```bash
cd /path/to/pic-
bash scripts/setup_maxwell_mi400_case.sh 80000 20000
cd PIC-IFE_GEC
./1DPIC > run_paper.log 2>&1
```

## Slurm 示例

保存为 `run_maxwell_mi400.slurm`，分区名按账号修改：

```bash
#!/bin/bash
#SBATCH -J maxwell-mi400
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module load compiler/intel/2021.3.0
module load compiler/cmake/3.23.3

cd /path/to/pic-/PIC-IFE_GEC
mkdir -p OUTPUT/Field OUTPUT/Velocity OUTPUT/Energy OUTPUT/Average DUMP
./1DPIC > run.log 2>&1
```

提交：

```bash
sbatch run_maxwell_mi400.slurm
```

## 常见问题

| 现象 | 处理 |
|------|------|
| `cannot open ./INPUT/mesh.inp` | 工作目录错了，先 `cd PIC-IFE_GEC` 再运行 |
| `cmake` 找不到 `ifort` | 先 `module load compiler/intel/...`，或改用 `ifx` |
| 输出没有 `020000` | 确认 `INPUT/pic.inp` 里是 `20000, 0.05`，或重新运行配置脚本 |
| 画图脚本找不到文件 | 先确认结果已经保存到 `verification_runs/<case_name>/`，并且里面有 `Average_x_*.dat` 和 `velocity_IJ_3*.dat` |
| `SyntaxError: future feature annotations is not defined` | 服务器 Python 偏老；拉取最新版脚本，当前 `draw_png.py` 和 `postprocess_maxwell_mi400.py` 已兼容 Python 3.6 |
| 论文级粒子数跑不动 | 先用 `1000/cell` 或 `5000/cell` 做证明图，再申请更多内存和墙钟时间 |

## 最短命令清单

```bash
cd /path/to/pic-
git pull --rebase
bash scripts/setup_maxwell_mi400_case.sh 1000 20000

cd PIC-IFE_GEC
module load compiler/intel/2021.3.0
module load compiler/cmake/3.23.3
rm -rf build
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(which ifort)"
cmake --build build -j"$(nproc)"
mkdir -p OUTPUT/Field OUTPUT/Velocity OUTPUT/Energy OUTPUT/Average DUMP
./1DPIC > run.log 2>&1

cd ..
case_dir=verification_runs/maxwellian_mi400_thermal_ppc1000_nt20000
mkdir -p "$case_dir"
cp PIC-IFE_GEC/case_config.txt PIC-IFE_GEC/INPUT/pic.inp PIC-IFE_GEC/MCC_jw/input/controlflow.txt PIC-IFE_GEC/run.log "$case_dir"/
cp PIC-IFE_GEC/OUTPUT/global_diagnostics.csv PIC-IFE_GEC/OUTPUT/physics_parameter.inp PIC-IFE_GEC/OUTPUT/normalize.inp "$case_dir"/
cp PIC-IFE_GEC/OUTPUT/Field/Average_x_012000.dat PIC-IFE_GEC/OUTPUT/Field/Average_x_020000.dat "$case_dir"/
cp PIC-IFE_GEC/OUTPUT/Velocity/velocity_IJ_3012000.dat PIC-IFE_GEC/OUTPUT/Velocity/velocity_IJ_3020000.dat "$case_dir"/

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
| `run_maxwell_mi400_thermal.slurm` | Slurm 提交脚本 | 可保留 |
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
run_maxwell_mi400_thermal.slurm
*.f90
```

### 下次重新跑标准 Maxwellian case

```bash
cd ~/pic-
git pull --ff-only

bash -n scripts/setup_maxwell_mi400_case.sh
bash scripts/setup_maxwell_mi400_case.sh 1000 20000 thermal

cd PIC-IFE_GEC
module purge
module load intel/2022.1
module load cmake/3.23.5

rm -rf build
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(which ifort)"
cmake --build build -j4

sbatch run_maxwell_mi400_thermal.slurm
squeue -u dg001947
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
