# pic- · PIC-IFE_GEC

**二维 IFE–PIC**（浸没有限元 + 粒子–网格）Fortran 程序，用于 **低气压放电 / 工艺腔等离子体** 等 **2D** 数值实验。单机可执行 **`1DPIC`**，默认 **无 MPI、无 GPU**。

---

## 本 README 定位（项目总览）

本文件主打三件事：**项目结构**、**功能边界**、**输出结果怎么读**。  
部署与代码使用细节请以 **[docs/DEPLOY.md](docs/DEPLOY.md)** 为主（本文件仍保留完整历史流程，便于追溯）。

### 你先看哪一段

- 想知道项目是做什么的：看「**功能与能力边界**」。
- 想知道目录怎么分工：看「**仓库结构总览**」。
- 想知道跑完会产出什么：看「**输出与重启总览**」。
- 想直接在西北一区跑最小算例：看「**西北一区最小算例（推荐，直接复制）**」。

---

## 功能与能力边界

- **问题类型：** 2D IFE-PIC 等离子体数值模拟（平面/轴对称）。  
- **核心流程：** 粒子推进（PIC）+ 场求解（IFE/SIDG/PPR）+ 碰撞（MCC_jw）。  
- **工程形态：** 单机可执行 `1DPIC`，默认无 MPI、无 GPU。  
- **当前阶段：** 探索/排坑中，优先保证“可编译、可运行、可复现”。

---

## 仓库结构总览

| 路径 | 作用 |
|------|------|
| `PIC-IFE_GEC/` | 主工程根目录（`CMakeLists.txt`、`1DPIC`、`INPUT/`、`OUTPUT/`、`DUMP/`） |
| `docs/` | 文档目录：`README.md` 讲文档关系，`DEPLOY.md` 讲部署与使用 |
| `readme` | 历史环境笔记（oneAPI 等） |
| `lzj_doc` | 旧文档指针 |

---

## 当前阶段（探索 / 排坑中）

本仓库目前处于**探索与排坑阶段**，文档和构建脚本会继续迭代。阅读和使用时请按下面原则：

- 以“**先能编译、再能运行、再改算例**”为顺序，不要一次改太多参数。  
- 若文档和实际行为冲突，以**当前代码与编译日志**为准，并把差异记录到 `docs/DEPLOY.md`。  
- 集群无法直连 GitHub 时，优先用 `git archive + scp` 同步文件内容。  
- 旧目录可能混有历史文件；遇到异常建议在**新目录解包重编译**。

---

**最短路径：** 先看 **[docs/README.md](docs/README.md)** 的文档关系，再按下面 **[第一次运行：复制粘贴](#第一次运行复制粘贴)** 或“西北一区最小算例”执行。

---

## 依赖与工程位置

- **编译器：** Intel **`ifort`**（或 **`ifx`**）；**CMake ≥ 3.10**。  
- **工程目录：** 本仓库里的 **`PIC-IFE_GEC/`**（该目录下有 **`CMakeLists.txt`**）。  
- **生成文件：** **`PIC-IFE_GEC/1DPIC`**（必须在 **`PIC-IFE_GEC`** 目录下运行）。

---

## 运行前快速自检（建议先做）

进入 **`PIC-IFE_GEC/`** 后先检查：

```bash
ls -la CMakeLists.txt INPUT/pic.inp INPUT/object.inp INPUT/ife.inp
```

说明：
- `INPUT/pic.inp`、`INPUT/object.inp`、`INPUT/ife.inp` 是当前仓库可见的基础输入文件。  
- 部分历史分支/算例会额外使用 `INPUT/mesh.inp`；若你的程序日志提示缺少它，请按课题组样例补齐。

---

## 第一次运行：复制粘贴

以下以 **西北一区（西安）** 为例；其它机器只要有 **`ifort` + `cmake`**，删掉 `module` 两行即可。

> 提示：当前为排坑阶段，若编译失败请优先看报错末尾文件名；常见是历史源码被误编译，按 `CMakeLists.txt` 与 `docs/DEPLOY.md` 的最新说明处理。

```bash
# 0) 加载编译器（每开一个新终端都要做，或写入 ~/.bashrc）
module load compiler/intel/2021.3.0
module load compiler/cmake/3.23.3

# 1) 进入你的家目录并进入仓库（按你实际路径改 pic-）
cd /work/home/$USER/pic-/PIC-IFE_GEC

# 2) 确认输入文件在「运行目录」下（至少应有 pic/object/ife 三个输入）
#    ls INPUT/pic.inp INPUT/object.inp INPUT/ife.inp

# 3) 编译
rm -rf build
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(which ifort)"
cmake --build build -j"$(nproc)"

# 4) 确认可执行文件已生成
test -x ./1DPIC && echo "OK: 1DPIC exists" || echo "FAIL: 未生成 1DPIC"

# 5) 运行（工作目录必须是 PIC-IFE_GEC，程序读 ./INPUT/）
#    先补齐常用输出子目录，避免首次输出时报 "No such file or directory"
mkdir -p OUTPUT/Field OUTPUT/Velocity OUTPUT/Particle OUTPUT/Global OUTPUT/Phase OUTPUT/Energy OUTPUT/History DUMP
./1DPIC
```

### 西北一区最小算例（推荐，直接复制）

> 目标：克隆后尽快验证“能编译 + 能运行 + 能输出 + 能正常结束”。

```bash
module load compiler/intel/2021.3.0
module load compiler/cmake/3.23.3
cd /work/home/$USER/pic-/PIC-IFE_GEC
bash ./run_min_case.sh
```

说明：
- `run_min_case.sh` 会先备份 `INPUT/pic.inp` 到 `INPUT/pic.inp.bak`。
- 默认把首个有效 `nt,dt` 改成 `1000, 0.05`（小算例快测）。
- 也可临时改参数：`NT=2000 DT=0.05 bash ./run_min_case.sh`。
- 正常结束标志是日志出现 `run finish at after it= 1000`（或你设定的 `NT`）。

跑起来后，标准输出里会出现读取 `INPUT` 文件的提示；结果在 **`OUTPUT/`**，重启数据在 **`DUMP/`**。  
**改时间步长、网格、边界** 等：见 **[docs/DEPLOY.md](docs/DEPLOY.md)** 中的部署与排错说明。

若报错类似 `./OUTPUT/Field/*.dat` 或 `./OUTPUT/Velocity/*.dat` 不存在，说明输出子目录未创建；先执行上面的 `mkdir -p` 再运行。

### 输出文件速查（汇报常用）

| 路径/文件 | 含义（简要） | 用途 |
|------|------|------|
| `OUTPUT/Field/field_IJ_*.dat` | 二维网格场量（电势/电场/密度等） | 主结果可视化（Tecplot 等） |
| `OUTPUT/Field/Average_x_*.dat` | 沿 `x` 的平均/剖面量 | 一维趋势图、对比工况 |
| `OUTPUT/Velocity/velocity_IJ_*.dat` | 速度相关网格输出 | 速度场诊断 |
| `OUTPUT/Particle/` | 粒子统计/分布输出 | 粒子行为分析 |
| `OUTPUT/Global/` | 全局守恒与总量统计 | 稳定性/收敛性检查 |
| `OUTPUT/Energy/` | 能量相关诊断 | 能量闭合检查 |
| `DUMP/` | 重启快照（checkpoint） | 中断续算、复现实验状态 |

> 小算例验证建议：先把 `nt` 调小（如 1000），确认能完整输出并出现 `run finish`，再回到长步数工况。

### 最小算例跑完后的标准文件结构（`nt=1000` 示例）

若你按本文最小算例流程运行并清理过中间步，常见保留结果如下：

```text
PIC-IFE_GEC/
├─ OUTPUT/
│  ├─ Field/
│  │  ├─ field_IJ_001000.dat
│  │  └─ Average_x_001000.dat
│  ├─ Velocity/
│  │  └─ velocity_IJ_3001000.dat
│  ├─ Energy/                  # 可能为空（取决于输出步长）
│  ├─ Particle/                # 可能为空（取决于输出开关）
│  ├─ Global/                  # 可能为空（取决于输出开关）
│  ├─ Phase/                   # 可能为空（取决于输出开关）
│  ├─ History/                 # 可能为空（取决于输出开关）
│  ├─ CellVolume.dat
│  ├─ physics_parameter.inp
│  ├─ normalize.inp
│  ├─ PartcountReal.dat
│  ├─ ElectronChange.dat
│  └─ IonChange.dat
└─ DUMP/
   ├─ var0001000dump
   ├─ phi0001000dump
   └─ par0001000dump
```

说明：
- `OUTPUT/Field/*001000.dat` 是第 1000 步主要场结果，可直接用于汇报绘图。
- `DUMP/*0001000dump` 是第 1000 步重启快照，可用于续算。
- 目录为空不代表错误，通常由输出步长和开关控制。

### 根目录杂项文件清理（可选）

历史运行可能在 `PIC-IFE_GEC/` 根目录留下 `SR *.dat`、`SN *.dat`、`MCCB *.dat`、`Norm_Error_*.txt`、`ife.msh` 等杂项文件。  
建议先备份再删除（在 `PIC-IFE_GEC/` 目录执行）：

```bash
# 1) 备份杂项文件（若没有匹配文件，命令会跳过）
JUNK_LIST="$(ls -1 | grep -E '^(10000 .*\.dat|MCCB +[0-9]+ +\.dat|SN +[0-9]+ +\.dat|SR +[0-9]+ +\.dat|Norm_Error_.*\.txt|ife\.msh)$' || true)"
if [ -n "$JUNK_LIST" ]; then
  echo "$JUNK_LIST" | tar -czf root_junk_backup_$(date +%Y%m%d_%H%M%S).tar.gz -T -
  echo "$JUNK_LIST" | sed 's/.*/"&"/' | xargs rm -f
fi
```

清理后建议再执行一次 `ls -a`，确认仅保留源码、构建目录与 `OUTPUT/`、`DUMP/` 等必要内容。

---

## 部署：西北一区（模块说明）

登录节点默认 **没有** `ifort`，系统自带 **cmake 过旧**，必须先：

```bash
module load compiler/intel/2021.3.0
module load compiler/cmake/3.23.3
which ifort && ifort --version
which cmake && cmake --version
```

**Slurm 作业脚本里也必须写同样的两行 `module load`**，否则计算节点找不到编译器。

### 获取代码

```bash
cd /work/home/$USER
git clone https://github.com/weqwoueu/pic-.git
cd pic-
```

访问 GitHub **超时**时：在本机 **`git archive` 打包**，用 **`scp`** 传到集群再解压；或问管理员要镜像/代理。

### 无 module 的机器（例如本机 Linux）

若已安装 **`gfortran`** 且想用 GNU 编译（与课题组环境不一致时请慎用）：

```bash
cd PIC-IFE_GEC
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(which gfortran)"
cmake --build build -j"$(nproc)"
```

---

## Slurm 作业示例

把 **`#SBATCH -p`** 改成你账号能用的分区；**首次建议先交互登录手动跑通 `./1DPIC`**，再交作业。

```bash
#!/bin/bash
#SBATCH -J pic-ife
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module load compiler/intel/2021.3.0
module load compiler/cmake/3.23.3

cd /work/home/$USER/pic-/PIC-IFE_GEC
mkdir -p OUTPUT/Field OUTPUT/Velocity OUTPUT/Particle OUTPUT/Global OUTPUT/Phase OUTPUT/Energy OUTPUT/History DUMP
./1DPIC
```

---

## 其它机房 / 本地 Intel

根目录 **`readme`** 中有 **`/opt/oneapi/setvars.sh`** 等示例；**西北一区以 `module` 为准**。

---

**一句话：** 进 **`PIC-IFE_GEC`** → **`module load`（若需要）** → **`cmake` + `cmake --build`** → 同目录先建 **`OUTPUT/*`** 与 **`DUMP`** → **`./1DPIC`**；部署与使用只看 **[docs/DEPLOY.md](docs/DEPLOY.md)**。
