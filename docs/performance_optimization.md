# PIC 性能优化与服务器测试

## 当前优化

- 使用 OpenMP 并行粒子推进和电荷/能量沉积；论文显式、无磁场算例启用粒子推进并行。
- 热边界随机数改为每线程独立状态，诊断累计量使用原子更新。
- CG/JPCG 工作数组按尺寸复用，避免每个时间步重复申请和释放。
- 速度谱和能谱只在第 1 步、每 1000 步及最后一步计算。
- `OUTPUT/performance_summary.csv` 记录各模块耗时、占比、每步耗时和线程数。
- 作业脚本默认使用 8 个 OpenMP 线程，并绑定到物理核心。可用 `PIC_OMP_THREADS` 修改。

PETSc CG/GAMG 暂不作为默认求解器。先根据性能 CSV 判断场求解占比；如果正式高粒子数算例仍主要耗在粒子模块，优先继续优化粒子路径，盲目替换场求解器不会带来同等幅度的总加速。

## 第一次服务器测试

先跑约 100 万初始粒子、2000 步的 Maxwellian smoke case，核对并行正确性和速度：

```bash
cd ~/pic-
git pull --ff-only

PIC_OMP_THREADS=8 bash scripts/setup_paper_case.sh \
  maxwellian none 1000 2000 thermal 101 0.05 1.0 400

cd PIC-IFE_GEC
module purge
module load intel/2022.1
module load cmake/3.23.5
rm -rf build
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(which ifort)"
cmake --build build -j8
sbatch run_maxwellian_mi400_dx1_dt005_ppc1000_seed101_thermal.slurm
```

查看作业：

```bash
squeue -u "$USER" -o "%.18i %.16j %.2t %.10M %.6D %R"
tail -f run.log
```

`tail -f` 只是在看日志，程序结束后按 `Ctrl+C` 退出查看。

## 完成判据

```bash
grep -E "OpenMP max threads|run finish|NaN|Infinity|failed|error" run.log
cat run_metadata.txt
cat OUTPUT/performance_summary.csv
tail -5 OUTPUT/global_diagnostics.csv
cp OUTPUT/performance_summary.csv performance_summary_omp8.csv
cp OUTPUT/global_diagnostics.csv global_diagnostics_omp8.csv
```

应满足：

- 出现 `run finish at after it=2000`；
- `OpenMP max threads = 8`；
- 没有 `NaN`、`Infinity`、`failed`；
- `performance_summary.csv` 的 `total` 行线程数为 8；
- 粒子和能量诊断处于原有可接受范围。

## 1 核对照

使用同一参数重新生成 1 核脚本并提交：

```bash
cd ~/pic-
PIC_OMP_THREADS=1 bash scripts/setup_paper_case.sh \
  maxwellian none 1000 2000 thermal 101 0.05 1.0 400

cd PIC-IFE_GEC
sbatch run_maxwellian_mi400_dx1_dt005_ppc1000_seed101_thermal.slurm
```

作业结束后保存 1 核汇总：

```bash
cp OUTPUT/performance_summary.csv performance_summary_omp1.csv
cp OUTPUT/global_diagnostics.csv global_diagnostics_omp1.csv
```

分别保留两次的 `performance_summary.csv`。总加速比为：

```text
speedup = 1 核 total seconds / 8 核 total seconds
parallel efficiency = speedup / 8
```

正式长算例开始前，先比较 1 核和 8 核末行诊断。并行浮点归并会产生末位差异，因此不要求逐位相同，应比较粒子守恒、能量误差和主要物理曲线。

## 线程数选择

可用相同方法测试 2、4、8 核：

```bash
PIC_OMP_THREADS=4 bash scripts/setup_paper_case.sh \
  maxwellian none 1000 2000 thermal 101 0.05 1.0 400
```

对当前单节点、内存带宽受限的粒子程序，默认先用 8 核。只有 8 核相对 4 核仍有明显收益时，再测试 16 核。

## 2026-08-26 线程基准

所有测试均为 Maxwellian、`mi/me=400`、`dx=1`、`dt=0.05`、热边界和随机种子 101。高粒子数测试使用 `ppc=80000`、`nt=200`。

低粒子数 `ppc=1000`、`nt=2000` 的 `comp` 对照中，1 线程主循环耗时为 401.238 s，8 线程为 193.367 s，总加速比为 2.075，耗时下降 51.8%。其中粒子推进加速 4.050 倍，电荷沉积加速 6.755 倍；场求解基本没有并行收益。

| 分区 | 线程 | 主循环耗时/s | 墙钟时间/s | 备注 |
| --- | ---: | ---: | ---: | --- |
| comp | 8 | 393.530 | 641 | 高粒子数基准 |
| comp | 16 | 272.292 | 未记录 | 相对 8 线程主循环加速 1.445 倍 |
| comp | 32 | 210.791 | 459 | comp 可用配置中的最快结果 |
| dg97 | 64 | 87.870 | 247 | 当前正式计算推荐配置 |
| dg97 | 128 | 138.269 | 275 | 线程过多，场求解和同步开销显著增加 |

`comp` 和 `dg97` 的 CPU 型号不同，32 到 64 线程的差异不能解释为纯线程加速比，但可以用于选择实际运行配置。`dg97` 的 64 线程结果满足粒子守恒误差为 0，最终能量平衡误差为 `-1.118684e-4`。按短测试分项耗时估算，`nt=20000` 正式算例约需 2.2 小时。

最终选择：`dg97` 分区、64 个 OpenMP 线程、至少 30000 MiB 作业内存。不再测试 256 线程。

## dg97 编译与运行

`CMakeLists.txt` 默认使用 `PIC_CPU_TARGET=host`。`-xHost` 生成的可执行文件依赖编译节点的指令集，因此不能把在 `comp`/登录环境生成的 AVX-512 可执行文件直接拿到不支持 AVX-512 的 `dg97` 节点运行。最高性能做法是在目标计算分区重新编译：

```bash
cd ~/pic-/PIC-IFE_GEC

srun -p dg97 -N 1 -n 1 -c 8 --mem=4000M -t 00:30:00 \
  bash -lc '
    module purge
    module load intel/2022.1
    module load cmake/3.23.5
    cd /data/home/dg001947/pic-/PIC-IFE_GEC
    rm -rf build-dg97
    cmake -S . -B build-dg97 \
      -DCMAKE_Fortran_COMPILER="$(which ifort)" \
      -DPIC_CPU_TARGET=host
    cmake --build build-dg97 -j8
  '
```

若必须在其他型号的 x86-64 节点编译并拿到 `dg97` 运行，可使用 `-DPIC_CPU_TARGET=avx2`，性能可能略低但兼容性更好。

生成正式 64 线程作业时直接指定分区和内存，不再手工修改 Slurm 文件：

```bash
cd ~/pic-

PIC_OMP_THREADS=64 \
PIC_SLURM_PARTITION=dg97 \
PIC_SLURM_MEMORY_MIB=30000 \
bash scripts/setup_paper_case.sh \
  maxwellian none 80000 20000 thermal 101 0.05 1.0 400
```
