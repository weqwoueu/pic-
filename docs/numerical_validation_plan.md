# 论文数值验证执行计划

本计划以 `pdf材料/工作建议.pdf` 和 `pdf材料/3_plasma_expansion (12).pdf` 为准。当前只把已经验证过的 Maxwellian thermal-reservoir 路线投入计算；Kappa 和 Polytropic 必须先完成与初始分布一致的 thermal 边界采样，不能直接复用 Maxwellian 边界后当作正式结果。

## 当前基线

已完成的 `maxwellian_dx1_dt005_ppc1000_seed101_thermal` 是当前可复现的流程验证算例：

- `mi/me=400`，`Z=1`，`Te0/Ti0=100`
- `dx=dy=lambda_D0`，`dt=0.05/omega_pe`
- `omega_pi t=30/50` 对应 `it=12000/20000`
- `ppc=1000`，总初始宏粒子数 `1,024,000`
- `seed=101`
- `gamma_e(t30)=1.0443 +/- 0.0049`，`R^2=0.279`
- `gamma_e(t50)=1.0468 +/- 0.0058`，`R^2=0.164`

这组结果接近论文 Maxwellian 的 `gamma_e=1.023 +/- 0.003`，但粒子数远低于正式设置，只能证明计算、诊断、归档和后处理链路可用。

## 统一拟合规则

所有 case 使用同一规则：

```text
ln(Te/Te0) = (gamma_e - 1) ln(ne/ne0) + b
10^-3 <= ni/n0 < 1
```

每次拟合必须保存：

- 拟合时刻和空间范围
- 纳入及排除的数据点
- 斜率、截距、`gamma_e`
- `gamma_e` 标准误差、残差标准误差和 `R^2`
- 是否排除 Kappa 极端高能尾

`scripts/postprocess_maxwell_mi400.py` 已按该规则输出 CSV、散点拟合图和 `gamma_fit_included` 标记。

## 第一阶段：Maxwellian 验证

先串行执行，不要同时提交到同一个 `PIC-IFE_GEC` 工作目录。生成的 Slurm 脚本带有运行锁和自动归档。

当前 `1DPIC` 没有启用 MPI 或 OpenMP，运行阶段只申请 `1` 个 CPU 核。`cmake --build build -j4` 仅用于加速编译，不代表计算程序会使用 4 核。服务器按核计费时不要为单个 `1DPIC` 作业申请多余核心。

`ppc=1000` 冒烟测试的墙钟时间为 `58:51`，总 CPU 时间为 `58:32`，峰值内存为 `109588 KiB`。集群默认每核约分配 `7916 MiB`，1 核额度足以容纳下一组。按粒子数线性外推只是提交前的保守资源估算，不是性能结论：`ppc=20000` 约为 `20 h` 和 `2.1 GiB` 量级。下一组暂设 `36 h` 墙钟；实际完成后必须用 `sacct` 更新估算，再决定 `ppc=40000/80000` 的资源请求和是否需要显式增加内存。

| 顺序 | ppc | seed | dx | dt | nt | 初始总粒子数 | 目的 |
|---|---:|---:|---:|---:|---:|---:|---|
| 1 | 20,000 | 101 | 1.0 | 0.05 | 20,000 | 20,480,000 | 资源和噪声试跑 |
| 2 | 40,000 | 101 | 1.0 | 0.05 | 20,000 | 40,960,000 | 粒子数收敛 |
| 3 | 80,000 | 101 | 1.0 | 0.05 | 20,000 | 81,920,000 | 论文基线 |
| 4 | 80,000 | 202 | 1.0 | 0.05 | 20,000 | 81,920,000 | 随机种子误差 |
| 5 | 80,000 | 303 | 1.0 | 0.05 | 20,000 | 81,920,000 | 随机种子误差 |
| 6 | 80,000 | 101 | 1.0 | 0.025 | 40,000 | 81,920,000 | 时间步收敛 |
| 7 | 80,000 | 101 | 0.5 | 0.05 | 20,000 | 327,680,000 | 网格收敛，先确认节点内存 |

第 7 组规模约 3.28 亿宏粒子，必须根据第 1-3 组实际峰值内存估算后再提交。若节点资源不够，应先向老师说明资源约束，不能悄悄降低 ppc 后仍称作 PDF 规定的网格收敛。

## 单个 case 的配置命令

参数顺序为：`ppc nt boundary seed dt dx`。

```bash
cd ~/pic-
bash scripts/setup_maxwell_mi400_case.sh 20000 20000 thermal 101 0.05 1.0

cd PIC-IFE_GEC
module purge
module load intel/2022.1
module load cmake/3.23.5

rm -rf build
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(which ifort)"
cmake --build build -j4

sbatch run_maxwellian_dx1_dt005_ppc20000_seed101_thermal.slurm
```

作业成功后会自动归档到：

```text
verification_runs/maxwellian_dx1_dt005_ppc20000_seed101_thermal/
```

下载到本地后运行：

```powershell
python scripts\draw_png.py verification_runs\<case_name> 12000
python scripts\draw_png.py verification_runs\<case_name> 20000
python scripts\postprocess_maxwell_mi400.py verification_runs\<case_name>
```

时间步为 `0.025` 时，脚本会从配置中自动使用 `24000/40000`，后处理命令不需要手工指定步数。

## 第一阶段判定标准

- 同一 seed 重复运行时关键输出应一致；不同 seed 应产生不同统计实现。
- `gamma_e` 随 ppc 增加趋于稳定，噪声大致随 `1/sqrt(Nppc)` 下降。
- `ppc=40000` 与 `80000` 的差异应接近或小于随机种子样本标准差。
- `dx=1.0` 与 `0.5` 的 `gamma_e` 差异若超过约 2%-3%，需要继续检查网格依赖。
- `dt=0.05` 与 `0.025` 应比较 `gamma_e`、前沿位置、高能尾和能量预算。
- 所有正式 case 都必须具有能量预算和电子、离子粒子数预算。

## 后续阶段

1. 完成 Kappa `kappa=2` 与 Polytropic `gamma_e=2` 的初始分布和同分布 thermal 边界参数化。
2. 对这三种代表分布重复网格、时间步、ppc 和 seed 验证。
3. 补充 Polytropic `gamma_e=3` 的 step-like 离散敏感性测试。
4. 数值验证通过后，再开展质量比和 thermal/specular 边界的正式物理对比。
