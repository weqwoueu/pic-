# PIC-IFE_GEC（本目录 = CMake 工程根）

你在 **`PIC-IFE_GEC/`** 里：这里有 **`CMakeLists.txt`**，编译生成的可执行文件 **`1DPIC`** 也放在**本目录**；运行程序时 **当前工作目录必须是这里**（程序读 **`./INPUT/`**、写 **`./OUTPUT/`**、**`./DUMP/`**）。

---

## 文档在哪、为什么要这样分

| 文档 | 路径 | 解释 |
|------|------|------|
| **介绍 + 部署** | [仓库根目录 `../README.md`](../README.md) | 讲「这是什么」、**怎么克隆**、**怎么 `module load`**、**第一次从编译到 `./1DPIC` 的复制粘贴**、Slurm、Git 连不上怎么办。新人**先读这个**。 |
| **使用说明 + 参数表** | [`../docs/USAGE.md`](../docs/USAGE.md) | 讲**跑通之后先改哪个 `INPUT`**、**目录别搞错**、**当前仓库输入文件与参数表**、结果文件在哪、常见报错。**改算例必看**。 |
| **文档索引** | [`../docs/README.md`](../docs/README.md) | 只列上面两个文档的链接和推荐阅读顺序。 |

**为什么分两份：** 根目录 README 解决「**环境装好、二进制有了、能跑起来**」；`docs/USAGE.md` 解决「**算例怎么改、参数别写乱**」。只打开本仓库子目录的人，先看本文件再点链接即可。

---

## 本目录里最常用的命令（摘要）

```bash
module load compiler/intel/2021.3.0   # 西北一区示例；其它机器按实际
module load compiler/cmake/3.23.3

rm -rf build
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(which ifort)"
cmake --build build -j"$(nproc)"
mkdir -p OUTPUT DUMP
./1DPIC
```

完整说明与排错仍以 **[../README.md](../README.md)**、**[../docs/USAGE.md](../docs/USAGE.md)** 为准。
