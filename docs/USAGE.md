# PIC-IFE_GEC 使用说明

**编译与集群模块** 见仓库根目录 **[README.md](../README.md)**。  
本文回答：**第一次怎么跑通**、**跑通后先改什么**、**每个 `INPUT` 里是什么**、**结果在哪**、**错了怎么查**。

---

## 当前阶段（探索 / 排坑中）

本说明对应的是当前仓库状态，仍在持续排坑。使用时建议：

- 先追求“能编译、能运行”，再逐步改参数与算例。  
- 若文档与实际编译/运行行为冲突，以**最新代码与终端日志**为准。  
- 发现差异时，优先把结论补到本文，保持团队可复现。

---

## 0. 第一次用：按顺序做

1. 读完根目录 **[README.md](../README.md)** 里的 **「第一次运行：复制粘贴」**，在 **`PIC-IFE_GEC/`** 下编出 **`./1DPIC`** 并成功执行一次。
2. 确认当前目录是 **`PIC-IFE_GEC`**，且存在 **`INPUT/pic.inp`**、**`INPUT/object.inp`**、**`INPUT/ife.inp`**。  
   若你的程序日志明确提示还需要 `INPUT/mesh.inp`，请按课题组样例补齐（当前仓库未自带该文件）。
3. 需要 **Tecplot** 或同类软件查看 **`OUTPUT/`** 里生成的 `.plt` / 数据文件（具体扩展名以运行结果为准）。

---

## 0.1 运行目录约定（必看）

程序用 **相对路径** 读 **`./INPUT/...`**、写 **`./OUTPUT/`**、**`./DUMP/`**。

| 正确 | 错误 |
|------|------|
| `cd .../PIC-IFE_GEC` 再 `./1DPIC` | 在仓库根 `pic-/` 或 `code/` 里直接运行可执行文件 |
| **`INPUT/`** 与 **`1DPIC`** 在同一层（都在 `PIC-IFE_GEC/`） | 只有 `code/INPUT/` 没有上层 `INPUT/` |

验证：

```bash
cd /path/to/PIC-IFE_GEC
pwd
ls -la INPUT/pic.inp INPUT/object.inp INPUT/ife.inp ./1DPIC
```

---

## 0.2 第一次改算例：建议顺序

不要一次改光；按优先级试：

| 优先级 | 文件 | 建议 |
|--------|------|------|
| 1 | **`INPUT/pic.inp`** | 先把 **`nt, dt`** 改小做短测试（例如几千步），确认能跑完再看物理 |
| 2 | **`INPUT/mesh.inp`**（若算例使用） | 改计算域和网格分辨率；与 `object.inp` 几何一致 |
| 3 | **`INPUT/object.inp`** | 边界与物体；**`N_Objects=0`** 时也要满足文件格式（见参数表） |
| 4 | **`INPUT/ife.inp`** | IFE 求解与 **`delta`**（平面 / 轴对称） |
| 5 | **`INPUT/IDG_inf*.inp`** | 仅在用到自适应 DG/IFE 时改 |

改前 **`cp -r INPUT INPUT.bak`** 备份。

---

## 0.3 常见问题（跑不通时查）

| 现象 | 可能原因 | 处理 |
|------|----------|------|
| `cannot open ./INPUT/mesh.inp` | 工作目录不对，或算例依赖该文件但仓库未提供 | 先确认在 `PIC-IFE_GEC` 目录，再从课题组样例补齐 |
| `cmake` 报找不到 `ifort` | 未 `module load` Intel | 见根目录 README |
| 链接错误 / 缺 `ModuleMCCInterface` | 用了不完整的老 CMake | 使用本仓库 **`PIC-IFE_GEC/CMakeLists.txt`**（已包含 `MCC_jw`） |
| `pic.inp` 读崩 | 增删了行，**读入顺序**乱了 | 对照本文 **§5 参数表** 逐行核对 |
| `git clone` 超时 | 集群封 GitHub | `git archive` + `scp` 或镜像 |
| 样例 `pic.inp` 在 `Te_Sec` 后面还有外场行 | **当前程序不读**这些行 | 可忽略或删掉，见 **§5.3** |

---

## 1. 能力范围

| 维度 | 说明 |
|------|------|
| **空间维数** | **2D**（平面或轴对称：由 `ife.inp` 中 `delta` 等控制） |
| **场求解** | **IFE** + **SIDG**；**PPR** 做粒子–网格映射 |
| **粒子** | **PIC**；**IMPIC** 为隐式相关分支（`pic.inp`） |
| **碰撞** | **MCC**：`code/MCC/` 与 **`MCC_jw/`**（`ModuleMCCInterface`） |
| **诊断** | Tecplot 类输出、统计；重启 **`./DUMP/*dump`** |
| **并行** | **单机 `1DPIC`**；**无 MPI**；少量 **OpenMP** |

**不包含：** 三维主流程、GPU、开箱多节点 MPI。

---

## 2. 可执行文件与构建关系

| 产物 | 说明 |
|------|------|
| **`PIC-IFE_GEC/1DPIC`** | 主程序可执行文件（由 `CMakeLists.txt` 构建） |
| **`MCC_jw/code/Main.F90`** | 独立程序；**不参与** `1DPIC` |

**链接进 `1DPIC`：** `code/**/*.f90` + `MCC_jw/code/**/*.f90`（排除 `Main.F90`）+ 根目录 **`OUTPUT_velocity.f90`**、**`Output_Energy.f90`**、**`generate_Elementmap.f90`**。见 **`CMakeLists.txt`**。

---

## 3. 源码目录导览（`PIC-IFE_GEC/code/`）

| 目录 | 角色（简要） |
|------|----------------|
| **`PIC/`** | 主循环相关、**`Dump_2D`**、**`Restart_2D`** |
| **`Data/`** | `Domain_2D`、`Field_2D`、`Particle_2D`、`IFE_Data`、`IFE_INTERFACE`、`IMPIC_Data_2D` |
| **`IFE-Rectangular/`**、`IFE-Solver/`、`IFE-error/` | IFE 网格、求解、误差 |
| **`SIDG/`** | 自适应 DG/IFE |
| **`Setup/`** | **`mesh.inp` / `ife.inp`**、**`Setup_IFE_Mesh_2D`**、**`SetupGrids_2D_QLL`** |
| **`IMPIC/`** | PrePush / PostPush / Move |
| **`PPR/`** | **`GetEfield_SIDG_PPR`** |
| **`In-Output/`** | 输入解析与 Tecplot 输出相关代码 |
| **`MCC/`** | 旧版 MCC |
| **`MCC_jw/`**（与 `code/` 并列） | **`data/`** 截面；**`code/`** JW 实现 |

---

## 4. 输入文件一览

工作目录：**`PIC-IFE_GEC/`**，路径 **`./INPUT/...`**。

| 文件 | 作用 |
|------|------|
| **`INPUT/mesh.inp`**（可选） | 域与网格（部分历史算例使用） |
| **`INPUT/object.inp`** | 物体与边界 |
| **`INPUT/pic.inp`** | PIC 主控、restart、IMPIC、步长、物种与注入 |
| **`INPUT/ife.inp`** | IFE 与 `delta`（平面/轴对称） |
| **`INPUT/IDG_inf*.inp`** | 自适应 DG/IFE |
| **`INPUT/B.DAT`** | 外磁场（与 `pic.inp` 磁场开关） |
| **`MCC_jw/input/`、`MCC_jw/data/`** | 碰撞（启用 JW MCC 时） |

---

## 5. 输入文件参数表（按当前实现与样例整理）

增删行会导致读错或崩溃。首行可空以吞掉标题；**整行** `!` 注释在自由格式下通常可忽略。

### 5.1 `INPUT/mesh.inp`（可选，部分算例使用）

| 顺序 | 变量 | 含义 |
|------|------|------|
| 1 | `xmin, ymin, zmin` | 域下角点 |
| 2 | `xmax, ymax, zmax` | 域上角点 |
| 3 | `nnx, nny, nnz` | 各方向点数 |
| 4 | `hx(1), hx(2)` | 网格步长（或控制量） |
| 5 | `hz` | 第三方向标量（2D 常占位） |

更复杂块由 **`SetupGrids_2D_QLL`** 等再读，以源码与样例为准。

### 5.2 `INPUT/object.inp`

| 顺序 | 内容 | 含义 |
|------|------|------|
| 1 | `N_Objects, vacuum` | 物体个数；真空区域阈值 |
| 每物体 | `Shape, Axis, Dimensions` | 形状、轴、尺寸 |
| | `Locations(1–4, :)` | 定位 |
| | `Regions, Direction, Phi, Eps, Erosion` | 区域、内外、电势、介电常数、刻蚀 |
| `Erosion>0` | `Wall` 若干行 | 壁面 |
| `Erosion≤0` | 5 行占位 | 必须跳过 |
| | `N_Boundary` | 边界段数 |
| 每边界 | `bc_index, bc_value, bc_point_1(2), bc_point_2(2)` | 类型与端点 |
| 周期边界 | `ALQ(1,i), ALQ(2,i)` | 周期参数 |

**`Shape==5`（椭圆）** 几何不满足会 **Stop**。

### 5.3 `INPUT/pic.inp`（`OPEN(5)`）

| 顺序 | 变量 | 含义 |
|------|------|------|
| 0 | 空行 | 吞标题 |
| 1 | `irestart, ilap` | 重启 |
| 2 | `IMPIC_index` | 隐式 PIC |
| 3 | `Bfiled_index` | 外磁场 |
| 4 | `Density_ref` | 密度参考 (m⁻³) |
| 5 | `Temperature_ref` | 温度参考 (eV) |
| 6 | `ParticlePerGrid` | 每格粒子相关 |
| 7 | `affp_bjw(1:3)` | 权重 |
| 8 | `N_part_tot` | 粒子总数相关 |
| 9 | `OutputRemovedParticles, OutputGlobalMoment` | 输出开关 |
| 10 | `ispe_tot` | 物种数 |
| 11 | `nt, dt` | 步数与时间步 |
| 12 | `Erostart` | 刻蚀/输出起始 |
| 13 | `n_updatefld` | 场更新间隔 |
| 14 | `n_pdiag` | 粒子诊断间隔 |
| 15 | `n_stride(1:ispe_tot)` | 输出步长 |
| 16–18 | `n_engydiag, n_gdiag, n_dump` | 诊断与 dump |
| 19–21 | `index_*`, `x/y/zcenter`, `c` | 诊断与常数 |
| 22 | 每物种 | `qs, xm` |
| 23–28 | 场/粒子边界 | `f_periodic`, `f_zeroe`, `periodic`, `pabsorb`, `preflect`, `pemit` |
| 29 | `phiouter(1:4)` | 外电势 |
| 30 | `den0_ref, Te_ref, phi0_ref` | 参考量 |
| 31 | `ispe_inject` 及循环 | 注入块 |
| 32 | `Te_Sec` | 二次电子能量等 |

读完 **`Te_Sec`** 后会关闭输入。若样例在 **`Te_Sec` 之后**还有外场、`Lz` 等行，当前版本通常不读取这些尾部行。

### 5.4 `INPUT/ife.inp`

| 顺序 | 含义 |
|------|------|
| 1 | 是否生成新 IFE 网格 |
| 2 | 界面单元占比初值 |
| 3 | 非线性求解器 |
| 4 | 稀疏线性求解器 |
| 5 | 节点固定/浮点 |
| 6–11 | 迭代与容差 |
| 12 | `delta`：0 平面 / 1 轴对称 |
| 13 | `LOG_modify` |

### 5.5 `INPUT/IDG_inf.inp`

惩罚系数、重复次数、多块 `area_flag` 与边界框；见 **`generate_Adaptive_DGT_IFE.f90`**。**`IDG_inf_step_refine.inp`** 可被程序改写。

---

## 6. 输出与重启

| 位置 | 内容 |
|------|------|
| **`OUTPUT/`** | **`physics_parameter.inp`**、**`normalize.inp`**、**`PartcountReal.dat`**、**`ElectronChange.dat`** / **`IonChange.dat`**、Tecplot 等 |
| **`DUMP/`** | 重启用；**`pic.inp`** 里 **`irestart`** |

---

## 7. 维护索引

| 任务 | 位置 |
|------|------|
| 改算例 | **`INPUT/*.inp`** |
| 改代码 | **`PIC-IFE_GEC/code/`**、**`MCC_jw/`** |
| 改构建 | **`PIC-IFE_GEC/CMakeLists.txt`**（**`1DPIC`**） |
| 主程序入口定位 | 在 `code/` 与 `MCC_jw/code/` 中按 `program` 关键字检索（不同历史版本可能不同） |
| JW 接口 | **`MCC_jw/code/Interface_IFE/MCCInterface.f90`** |

---

## 8. 附录

| 路径 | 说明 |
|------|------|
| **`PIC-IFE_GEC/build.sh`** | 简易 `cmake && make` |
| **`readme`**（仓库根） | 历史 oneAPI 示例 |
