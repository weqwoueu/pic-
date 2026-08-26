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
