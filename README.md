# pic- · PIC-IFE_GEC

**二维 IFE–PIC**（浸没有限元 + 粒子–网格）Fortran 程序，用于 **低气压放电 / 工艺腔等离子体** 等 **2D** 数值实验。单机可执行 **`1DPIC`**，默认 **无 MPI、无 GPU**。

---

## 后人接手：请先读这两份

| 顺序 | 文档 | 你用它做什么 |
|------|------|----------------|
| 1 | **本文（README）** | 克隆、装环境、**编译出 `1DPIC`**、在集群上**运行** |
| 2 | **[docs/USAGE.md](docs/USAGE.md)** | **第一次跑通**的逐步说明、**改算例**先改谁、`INPUT` **参数表**、输出说明、排错 |

**最短路径：** 先按下面 **[第一次运行：复制粘贴](#第一次运行复制粘贴)** 做一遍；跑通后再打开 `docs/USAGE.md` 改参数。

---

## 文档分工

| 文档 | 内容 |
|------|------|
| **[docs/USAGE.md](docs/USAGE.md)** | 上手流程、能力说明、源码树、`INPUT` 与**参数表**、输出、常见问题 |
| **[docs/README.md](docs/README.md)** | 文档索引 |

---

## 依赖与工程位置

- **编译器：** Intel **`ifort`**（或 **`ifx`**）；**CMake ≥ 3.10**。  
- **工程目录：** 本仓库里的 **`PIC-IFE_GEC/`**（该目录下有 **`CMakeLists.txt`**）。  
- **生成文件：** **`PIC-IFE_GEC/1DPIC`**（必须在 **`PIC-IFE_GEC`** 目录下运行）。

---

## 第一次运行：复制粘贴

以下以 **西北一区（西安）** 为例；其它机器只要有 **`ifort` + `cmake`**，删掉 `module` 两行即可。

```bash
# 0) 加载编译器（每开一个新终端都要做，或写入 ~/.bashrc）
module load compiler/intel/2021.3.0
module load compiler/cmake/3.23.3

# 1) 进入你的家目录并进入仓库（按你实际路径改 pic-）
cd /work/home/$USER/pic-/PIC-IFE_GEC

# 2) 确认输入文件在「运行目录」下（本仓库已带 INPUT/；若没有，从 code/INPUT 拷贝一份）
#    ls INPUT/mesh.inp INPUT/pic.inp INPUT/object.inp INPUT/ife.inp

# 3) 编译
rm -rf build
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(which ifort)"
cmake --build build -j"$(nproc)"

# 4) 确认可执行文件已生成
test -x ./1DPIC && echo "OK: 1DPIC exists" || echo "FAIL: 未生成 1DPIC"

# 5) 运行（工作目录必须是 PIC-IFE_GEC，程序读 ./INPUT/）
mkdir -p OUTPUT DUMP
./1DPIC
```

跑起来后，标准输出里会出现读 **`mesh.inp` / `object.inp` / `pic.inp`** 的提示；结果在 **`OUTPUT/`**，重启数据在 **`DUMP/`**。  
**改时间步长、网格、边界** 等：见 **[docs/USAGE.md](docs/USAGE.md)** 中的「第一次改算例」与参数表。

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
mkdir -p OUTPUT DUMP
./1DPIC
```

---

## 其它机房 / 本地 Intel

根目录 **`readme`** 中有 **`/opt/oneapi/setvars.sh`** 等示例；**西北一区以 `module` 为准**。

---

## 仓库根目录

| 路径 | 说明 |
|------|------|
| **`PIC-IFE_GEC/`** | 唯一 CMake 工程与 **`1DPIC`**；内有 **[PIC-IFE_GEC/README.md](PIC-IFE_GEC/README.md)**（说明文档分工与本目录含义） |
| **`docs/`** | **[USAGE.md](docs/USAGE.md)** 使用说明与参数表 |
| **`readme`** | 历史环境笔记 |
| **`lzj_doc`** | 文档指针 |

---

**一句话：** 进 **`PIC-IFE_GEC`** → **`module load`（若需要）** → **`cmake` + `cmake --build`** → 同目录 **`mkdir -p OUTPUT DUMP`** → **`./1DPIC`**；算例怎么改只看 **[docs/USAGE.md](docs/USAGE.md)**。
