# LACT_sim 可信光学追踪修改实施记录

日期：2026-07-25
修改基线：GitHub `main`，提交 `0048369cd06e87322d9ef3611a4a369a1152fe99`

## 1. 实施原则

本轮修改不以复刻 sim_telarray 为目标。sim_telarray 只用于发现我们可能遗漏的物理
步骤、坐标语义和边界条件。只要 LACT_sim 自身的表示满足明确、完整、可验证的内部
不变量，就不强制改成 sim_telarray 的坐标轴、数组编号或实现结构。

具体原则：

- 不改动本来正确、只是与 sim_telarray 表示不同的内部编号；
- 不把用户配置改成 sim_telarray 风格；
- 修复 LACT_sim 自身会造成未定义行为、错误物理量或不可复现结果的问题；
- 对非法输入 fail-fast，不静默截断成看似正常的结果；
- 输出足够的中间量，使光程、随机流和运行来源可以复核。

## 2. 修改后的主执行流程

`run_corsika_trace` 的目标流程如下：

1. 展开主配置和组件配置。
2. 读取 EventIO 或 PhotonCsv，保留 shower、array、telescope 和 bunch 身份。
3. 把输入坐标转换到 telescope-local 光学坐标。
4. 对 EventIO 应用用户给出的单值参考偏移。
5. 从 bunch 身份派生稳定的 represented-photon 随机流。
6. 应用大气和预几何探测概率。
7. 与镜片求交，做镜面反射和粗糙度散射。
8. 检查入射段和反射段遮挡。
9. 与输出面求交，累计主镜光程和到达时间。
10. 映射到相机像素；启用集光器时，在集光器内部继续追迹。
11. 把集光器内部光程加入时间，再做后几何接受。
12. 累计像素、波形、trigger 和阵列结果。
13. 输出展开配置、代码版本和输入文件 SHA-256。

## 3. EventIO 单值 `-16` 的最终语义

用户接口保持为一个标量：

```ini
source.eventio_reference_z_m=-16
```

内部解释为 telescope-local 平移：

```text
(dx, dy, dz) = (0, 0, -16 m)
```

因此：

```text
(x, y, z) -> (x, y, z - 16 m)
```

这里的 `(0,0)` 是偏移量的 x/y 分量，不是把光子的物理 x/y 清零。

二维 EventIO 原始 `z=0`，所以变成 `z=-16 m`；三维 EventIO 保留自己的原始 z，
再加同一个标量偏移。例如 `(1.25,-2.5,0.75)` 变成
`(1.25,-2.5,-15.25)`。

实现位置：

- `SourceRuntimeConfig::eventio_reference_z_m`
- `applyEventIOReferenceZOffset`
- `transformEventIOBunchToTraceFrame`

旧键 `source.eventio_2d_input_plane_z_m` 仅作为弃用兼容别名。新旧键同时存在且值
冲突时直接报错。

## 4. Collector 正确性修改

修改前 Bezier collector 存在两个相互关联的问题：

- 交点函数找到有效根后，法线函数重新求根；
- 重新求根时可能找不到完全相等的根，并解引用 `end()`，ASan 已证实越界。

修改后每个曲面一次返回完整交点记录：

```text
hit
ray_distance
surface_parameter
position
normal
```

Bezier 法线直接使用同一次求交得到的 Bezier 参数及解析导数，不再重新求根。
SquareCone 在所有壁面候选中按正向光程选择第一个交点；缝和角上对同距离法线做
确定性的平均。反射使用该交点记录携带的法线。

同时增加：

- 稳定的线性/二次方程退化处理；
- NaN、退化法线和反向交点拒绝；
- 显式 500 次反射上限状态；
- 集光器总内部光程；
- 集光器光程导致的时间延迟；
- CSV 调试字段 `collector_path_length_m`、
  `collector_time_delay_ns` 和
  `collector_reflection_limit_reached`。

集光器光程单位在内部几何中是毫米，离开 collector 模块时转换为米，再除以配置的
`propagation.speed_of_light_m_per_ns` 加到 `hit.time_ns`。

## 5. 稳定随机流

新增固定的 SplitMix64 派生规则和阶段编号。represented photon 的身份只来自：

```text
shower_event_id
array_id
telescope_id
source_bunch_index
represented_index
```

绝对位置和方向不再作为随机身份。因此对整个光学系统做坐标平移时，镜面粗糙度抽样
不会改变。

各随机阶段使用不同 seed：

- 缺失波长；
- 预几何接受；
- 镜面粗糙度；
- 后几何接受；
- 镜片位置误差；
- 镜片法线误差；
- 曲率半径误差；
- 反射率误差。

镜片误差还加入 facet ID，因而 CSV 行顺序改变不会改变同一 facet 的误差实现。

## 6. 输入和配置契约

### 6.1 PhotonCsv 和 EventIO

现在拒绝：

- 非有限位置、方向和时间；
- 零长度方向；
- 负或非有限 weight、multiplicity；
- PhotonCsv 非正或非有限 wavelength；
- 二维 EventIO 明显满足 `cx²+cy²>1` 的方向余弦；
- 三维 EventIO 零方向。

二维 EventIO 对浮点舍入造成的微小超界保留有限容差，但不再把明显非法值静默改成
`cz=0`。

### 6.2 相机

相机 CSV 和生成式相机现在检查：

- ID 非负且唯一；
- 坐标和尺寸有限，尺寸为正；
- 像素内部不能重叠，边界接触允许。

重叠检查使用圆和凸多边形的分离轴测试，不依赖像素 ID 或输入行顺序。

### 6.3 Square-cone collector

当前 collector 实现是方形入口，因此启用时明确要求：

- 相机像素为方形；
- 所有像素尺寸一致；
- 入口、出口和高度均为有限正数；
- 入口大于出口。

Bezier collector 可使用配置尺寸。Parabolic collector 的曲面参数仍是固定标定值，
所以只有入口 `0.0243992 m`、出口 `0.0150014 m`、高度 `0.0252718 m`
时才允许启用，防止配置显示一个尺寸、实际追迹另一个尺寸。

### 6.4 效率概率

镜面反射率、滤光片透过率、PDE、大气谱因子及其曲线必须在 `[0,1]`。随机 PE
模式遇到超出 `[0,1]` 的最终概率会报错，不再静默 clamp。

Facet reflectivity 也限制在 `[0,1]`；随机反射率误差采用截断到物理区间的实现。

同时启用 `atmosphere.model` 和额外的大气谱因子时默认报错。确实要表达两个独立
因子时必须显式设置：

```ini
atmosphere.allow_additional_spectral_factor=true
```

显式光线追迹 collector 与解析 `cos(theta)` funnel 同时启用时也默认报错。确实
有意叠乘时必须设置：

```ini
camera.allow_analytic_funnel_with_collector=true
```

## 7. 输出 provenance

CMake 配置时执行 `git describe --always --dirty`。二进制保存该字符串作为
`producer_version`，所以未提交源码会带 `-dirty`。

当输出 HDF5 或 lact_event ROOT 时，程序在追迹前对实际 EventIO/PhotonCsv 输入做
流式 SHA-256。两个输出都写入：

- `producer_version`
- `source_path`
- `source_sha256`

SHA-256 实现有空字符串、`abc` 和文件读取的已知答案测试。

## 8. 新增或扩展的回归测试

- 2D/3D 共用标量 z 偏移，且 x/y 保持不变；
- collector 直射和反射光程、时间、反射上限状态；
- Bezier collector 扫描；
- expectation/stochastic 都拒绝负 multiplicity；
- 随机流不依赖坐标；
- 粗糙度散射的整体平移不变量；
- facet error 的 CSV 行顺序不变量；
- 非法 PhotonCsv 行；
- 非法效率曲线和重复大气/collector 配置；
- 相机像素重叠和非方形 square-cone collector；
- SHA-256 已知答案。

## 9. 完整规模运行补出的输出一致性修复

真实 1000-shower 运行还暴露了小样本测试不容易看到的三个问题：

- NSB 曾由图像/波形与 trigger/ROOT 分别抽样。现在每个 event-telescope 只调用一次
  `generateTimeBinnedNsbPe`，图像、保存波形、trigger、HDF5 和 ROOT 共用同一实现；
- ROOT 的动态树曾因 `SetDirectory(nullptr)` 不能真正执行中途 basket flush。现在
  动态树绑定输出 `TFile`，检查 `Fill`、`FlushBaskets`、`AutoSave` 和最终 `Write`
  的返回值，关闭文件前再解绑，避免超大单对象和所有权冲突；
- HDF5 的图像参考时间和 trigger 时间曾是 `float`，而 ROOT 为 `double`，时间 bin
  边缘可能产生极少数触发差异。相关 HDF5 compound 字段现统一为 `double`。

最终完整结果的 ROOT/HDF5 单镜 trigger 已逐条检查 24,706 条，零差异；第一 shower
诊断结果也已从保存的 31-bin 波形重算 12 条 trigger，零差异。

这些修复遵循的是 LACT_sim 自己的内部一致性要求，与是否采用 sim_telarray 的坐标轴
或数组编号无关。

## 10. 当前验证状态

最终源码已复制到现有 Linux 测试机，在彼此独立的构建目录中完成：

- 普通构建：编译成功，CTest `28/28` 通过；
- UBSan：`UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1`，CTest `28/28`
  通过，未报告未定义行为；
- ASan：`ASAN_OPTIONS=halt_on_error=1:detect_leaks=1`，CTest `28/28`
  通过，未报告内存错误或泄漏；修改前可复现的 Bezier collector 越界不再出现；
- 启用 ROOT/HDF5 的完整构建：编译成功，CTest `28/28` 通过；
- 以五行 PhotonCsv 实际运行 `run_corsika_trace`，HDF5 和 ROOT 输出均成功生成。
- 以真实 CORSIKA 7.741 μ 子文件运行 1000 个 shower，生成 1.1 GB HDF5 和
  201 MB ROOT；两文件均重新打开并完成逐条交叉检查；
- 全量机器检查 `sanity.passed=true`，ROOT/HDF5 trigger
  `checked=24706, failure_count=0`。

实际输入文件的独立 `sha256sum` 为：

```text
5a5c70915a3fe22840d7afaf1d6e81923aa5bd6ded0237d46aff51f5b9aca756
```

HDF5 根属性和 ROOT `config` 树都包含 `producer_version`、`source_path` 和
`source_sha256`；ROOT 文件中的实际来源路径和哈希与上述输入一致。

远端验证使用的是不含 `.git` 目录的源码副本，所以该构建的
`producer_version=source-tree`。在真实 Git 工作树中配置 CMake 时会由
`git describe --always --dirty` 写入提交及 dirty 状态。

`git diff --check` 在最终本地差异上再次执行并通过。

这些结果验证了本轮实现和回归测试覆盖的内部契约。它们不是实验数据或独立光学程序
之间的物理标定，因此不能单凭 CTest 证明所有现实光学参数都已完成外部验证。

## 11. 本轮有意不做的事情

- 不把 LACT_sim 的坐标和编号强行改成 sim_telarray；
- 不实现完整模拟电子学链；
- 不用 sim_telarray 的内部数据结构替换现有程序；
- 不把两个程序逐事件数值完全相同作为正确性定义。

后续外部基准应比较物理可观测量，例如焦斑、离轴响应、collector angular
acceptance、光程分布和总探测效率，并把实现差异与物理错误分开解释。
