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
| `scripts/plot_maxwell_mi400.py` | 读取输出并画密度、电子热速度图 |

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

## 画图

回到仓库根目录：

```bash
cd /path/to/pic-
python3 scripts/plot_maxwell_mi400.py PIC-IFE_GEC 20000
```

输出图片：

```bash
PIC-IFE_GEC/figures/maxwell_mi400_t020000.png
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
| 画图脚本找不到文件 | 先确认模拟跑完，并且 `OUTPUT/Field` 和 `OUTPUT/Velocity` 里有数据 |
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
python3 scripts/plot_maxwell_mi400.py PIC-IFE_GEC 20000
```
