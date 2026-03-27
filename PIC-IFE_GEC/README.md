# PIC-IFE_GEC（本目录 = CMake 工程根）

你在 **`PIC-IFE_GEC/`** 里：这里有 **`CMakeLists.txt`**，编译生成的可执行文件 **`1DPIC`** 也放在**本目录**；运行程序时 **当前工作目录必须是这里**（程序读 **`./INPUT/`**、写 **`./OUTPUT/`**、**`./DUMP/`**）。

---

## 文档在哪、为什么要这样分

| 文档 | 路径 | 解释 |
|------|------|------|
| **介绍 + 部署** | [仓库根目录 `../README.md`](../README.md) | 讲「这是什么」、**怎么克隆**、**怎么 `module load`**、**第一次从编译到 `./1DPIC` 的复制粘贴**、Slurm、Git 连不上怎么办。新人**先读这个**。 |
| **部署 + 使用说明** | [`../docs/DEPLOY.md`](../docs/DEPLOY.md) | 讲**西北一区部署**、**最小算例跑通**、**目录别搞错**、结果文件在哪、常见报错。 |
| **文档索引** | [`../docs/README.md`](../docs/README.md) | 只列上面两个文档的链接和推荐阅读顺序。 |

**为什么分两份：** 根目录 README 解决「**项目是什么、目录与输出怎么看**」；`docs/DEPLOY.md` 解决「**怎么部署、怎么跑通**」。只打开本仓库子目录的人，先看本文件再点链接即可。

---

## 本目录里最常用的命令（摘要）

```bash
module load compiler/intel/2021.3.0   # 西北一区示例；其它机器按实际
module load compiler/cmake/3.23.3

rm -rf build
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(which ifort)"
cmake --build build -j"$(nproc)"
mkdir -p OUTPUT/Field OUTPUT/Velocity OUTPUT/Particle OUTPUT/Global OUTPUT/Phase OUTPUT/Energy OUTPUT/History DUMP
./1DPIC
```

## 最小算例一键跑通（推荐）

在本目录直接执行：

```bash
bash ./run_min_case.sh
```

脚本会自动：
- 备份 `INPUT/pic.inp` 到 `INPUT/pic.inp.bak`
- 将首个有效 `nt,dt` 调整为小算例默认值（`1000, 0.05`）
- 创建常用输出目录
- 重新编译并运行 `./1DPIC`

也可覆盖默认值：

```bash
NT=2000 DT=0.05 bash ./run_min_case.sh
```

跑通最小算例后，重点核对：
- `OUTPUT/Field/field_IJ_001000.dat`
- `OUTPUT/Field/Average_x_001000.dat`
- `OUTPUT/Velocity/velocity_IJ_3001000.dat`
- `DUMP/var0001000dump`、`DUMP/phi0001000dump`、`DUMP/par0001000dump`

完整说明与排错仍以 **[../README.md](../README.md)**、**[../docs/DEPLOY.md](../docs/DEPLOY.md)** 为准。
