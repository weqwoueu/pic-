# PIC-IFE_GEC 部署与运行（西北一区）

本文是**部署与代码使用主文档**。  
目标：克隆后在西北一区快速完成“可编译 + 可运行 + 可输出 + 可结束”的最小算例。

---

## 1) 首次部署：复制粘贴

```bash
module load compiler/intel/2021.3.0
module load compiler/cmake/3.23.3

cd /work/home/$USER/pic-/PIC-IFE_GEC

# 一键最小算例（默认 NT=1000, DT=0.05）
bash ./run_min_case.sh
```

可选参数覆盖：

```bash
NT=2000 DT=0.05 bash ./run_min_case.sh
```

成功标志：
- 日志出现 `run finish at after it= 1000`（或你设置的 NT）。

---

## 2) 手动部署（不用脚本）

```bash
module load compiler/intel/2021.3.0
module load compiler/cmake/3.23.3
cd /work/home/$USER/pic-/PIC-IFE_GEC

rm -rf build
cmake -S . -B build -DCMAKE_Fortran_COMPILER="$(which ifort)"
cmake --build build -j"$(nproc)"

mkdir -p OUTPUT/Field OUTPUT/Velocity OUTPUT/Particle OUTPUT/Global OUTPUT/Phase OUTPUT/Energy OUTPUT/History DUMP
./1DPIC
```

---

## 3) 最小算例跑通核对

建议至少确认这些文件存在：

- `OUTPUT/Field/field_IJ_001000.dat`
- `OUTPUT/Field/Average_x_001000.dat`
- `OUTPUT/Velocity/velocity_IJ_3001000.dat`
- `OUTPUT/PartcountReal.dat`
- `DUMP/var0001000dump`
- `DUMP/phi0001000dump`
- `DUMP/par0001000dump`

---

## 4) 常见问题（部署向）

| 现象 | 原因 | 处理 |
|---|---|---|
| `cmake` 找不到 `ifort` | 没有加载模块 | `module load compiler/intel/2021.3.0` |
| `./OUTPUT/Field/*.dat not found` | 输出子目录未创建 | 先执行 `mkdir -p OUTPUT/Field OUTPUT/Velocity OUTPUT/Particle OUTPUT/Global OUTPUT/Phase OUTPUT/Energy OUTPUT/History DUMP` |
| `git pull` 失败（GitHub 不通） | 集群网络限制 | 用 `git archive + scp` 离线同步 |
| `./INPUT/mesh.inp` 缺失 | 算例依赖该文件但仓库未带 | 从课题组样例补齐 |

---

## 5) 根目录杂项清理（可选）

```bash
cd /work/home/$USER/pic-/PIC-IFE_GEC
JUNK_LIST="$(ls -1 | grep -E '^(10000 .*\.dat|MCCB +[0-9]+ +\.dat|SN +[0-9]+ +\.dat|SR +[0-9]+ +\.dat|Norm_Error_.*\.txt|ife\.msh)$' || true)"
if [ -n "$JUNK_LIST" ]; then
  echo "$JUNK_LIST" | tar -czf root_junk_backup_$(date +%Y%m%d_%H%M%S).tar.gz -T -
  echo "$JUNK_LIST" | sed 's/.*/"&"/' | xargs rm -f
fi
```

---

项目结构、功能边界、输出物理含义请看根目录 **`README.md`**。
