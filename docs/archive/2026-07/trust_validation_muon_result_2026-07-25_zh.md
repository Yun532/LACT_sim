# LACT_sim 真实 CORSIKA μ 子全链路检验结果

日期：2026-07-25
代码基线：GitHub `main`，
`0048369cd06e87322d9ef3611a4a369a1152fe99` 加本轮可信光学修改

## 1. 目的

这不是只验证“程序能运行”的 smoke test，而是同时检查：

- 真实 CORSIKA/EventIO 输入读取和坐标转换；
- EventIO 统一 `(0,0,-16 m)` telescope-local 参考偏移；
- 实测约束镜面、粗糙度、遮挡和大气传播；
- 相机像素、Bezier collector、滤光片和 SiPM PDE；
- Cherenkov/NSB/final PE 分量；
- proxy PE 波形、单镜 trigger 和阵列 trigger；
- HDF5/ROOT 结构、配置展开和输入 SHA-256；
- 逐光子诊断量及总体统计的一致性。

## 2. 输入

输入是一份真实 CORSIKA 7.741 μ 子 EventIO：

```text
文件：LACT_review_muon_E100_th0_run000001.zst
大小：56 MiB
shower：1000
望远镜：2
IACT array offset：每个 shower 最多 20 组
CWAVLG：260–1000 nm
SHA-256：20250674adb895df40382f212fa99378dc22fd4ea3f3352bac2e9c868af2962b
```

输入卡标称 `PRMPAR 6`、`ERANGE 4–1000 GeV`、天顶角 `0°`、观测高度
`4400 m`。第一 shower 的实际能量为 `6.136686 GeV`。

## 3. 配置

生产配置：

```text
configs/official_tests/trust_full_muon_1000.cfg
```

第一 shower 逐级诊断配置：

```text
configs/official_tests/trust_full_muon_first_event_diagnostic.cfg
```

两者使用相同物理链：

- LACT2 2026-06-22 实测约束镜面；
- 90° 仰角镜面插值；
- `source.eventio_reference_z_m=-16`；
- MODTRAN 4400 m 沙漠大气；
- 最终结构 27 个遮挡 primitive；
- 1616 个方形像素；
- Bezier square-cone collector；
- 实测镜面反射率、滤光片透过率和 SiPM PDE 曲线；
- SkyCalc 暗夜 NSB，`0.076375 PE/ns/pixel`；
- 2 ns 时间 bin，`[-2,60] ns`，共 31 bin；
- 像素阈值 10 PE、相机 multiplicity 3、阵列 multiplicity 2。

这里的波形仍是 time-binned PE proxy，不是完整模拟电子学脉冲。没有加入 SiPM
饱和、串扰、afterpulse、暗计数和硬件触发板响应。NSB 在每个
event-telescope 上只生成一次；图像、保存波形、trigger 和 ROOT/HDF5 共用同一份
离散实现，不再分别抽样。

## 4. 第一 shower 逐级结果

### 4.1 光学计数

```text
输入 bunch                         1360
输入等效光子                      6397.600
镜面求交且可到输出面（遮挡前）     887
被结构遮挡                         203
  入射段                            189
  反射段                             14
遮挡后到达输出面                    684
命中相机像素区域                    533
collector/探测响应后接受            407
Cherenkov PE                      315.467
```

遮挡后输出面通过率为 `684/887 = 0.771139`。

两台望远镜分别为：

```text
tel 0：1195 bunch，5648.250 photons，355 accepted，281.241 PE
tel 1： 165 bunch， 749.350 photons， 52 accepted， 34.225 PE
```

### 4.2 图像、NSB 和波形

第一 shower 产生 7 个 array event、12 个 event-telescope 图像：

```text
Cherenkov PE 总和                  315.466838
NSB PE 总和                     91556
final PE 总和                   91871.466836
PE 分量最大绝对闭合误差              9.54e-7
稀疏波形行数                       85082
波形 Cherenkov PE 积分             315.466828
波形 NSB PE 积分                 91556
```

图像和波形中的 PE、Cherenkov PE、NSB PE 全部有限且非负。波形积分与图像积分
在 float 精度内一致。

NSB 数量级也与配置闭合：

```text
0.076375 PE/ns/pixel × 62 ns × 1616 pixel × 12 image
≈ 9.19e4 PE
```

### 4.3 Trigger

```text
单镜 trigger：2 / 12
阵列 trigger：0 / 7
```

第一 shower 中通过阈值的两台单镜不属于满足同一 array event 的双镜组合，因此没有
阵列 trigger。这是触发逻辑结果，不是输出丢失；配置使用
`output.save_only_triggered=false`，所以未触发事件也全部保存。

### 4.4 Collector

collector 诊断写出 534 行，包含：

- 接受状态；
- 内部反射次数；
- 是否达到反射次数上限；
- collector 内部光程；
- 对应时间延迟；
- 波长、权重、相对效率；
- 出口位置和方向。

示例光子内部光程 `0.0360559 m`，延迟 `0.120269 ns`，满足
`delay = path / 0.299792458`。

## 5. 输出结构与可复现性

第一 shower 诊断 HDF5：

```text
/events/corsika                    7
/events/corsika_showers         1000
/images/index                     12
/images/dense/*             12 × 1616
/trigger/telescope                12
/trigger/array                     7
/waveforms/samples             85082
/waveforms/time_centers_ns        31
/camera/pixels                  1616
/mirrors/facets                   54
/telescopes/table                  2
```

ROOT 中包含：

- `config`
- `corsika_events`
- `observations`
- `trace_summary`
- `waveforms`
- `waveform_config`
- `camera_pixels`
- `telescopes`
- `optics`

HDF5 和 ROOT 都实际写入：

```text
producer_version = source-tree
source_path       = /tmp/LACT_review_muon_E100_th0_run000001.zst
source_sha256     = 20250674adb895df40382f212fa99378dc22fd4ea3f3352bac2e9c868af2962b
```

测试机使用的是不含 `.git` 目录的源码副本，所以版本为 `source-tree`。在真实 Git
工作树构建时会写入 `git describe --always --dirty`。

## 6. 全 1000 shower 结果

完整运行用时 `487.145820 s`，其中光线追迹 `227.640622 s`。输入的 1000 个
shower 均完成读取和处理；其中 708 个 shower 至少产生一个望远镜输出事件。后者不是
丢失计数：HDF5 的 `/events/corsika_showers` 明确保存了全部 1000 个 shower header。

### 6.1 数据规模

```text
CORSIKA shower header                  1000
array event                           13100
event-telescope 图像/trigger           24706
望远镜                                    2
相机像素                               1616
输入 photon bunch                  7923785
输入等效光子                    39181433.587
```

全量 HDF5 的主要数据集为：

```text
/events/corsika_showers              1000
/events/corsika                     13100
/events/table                       13100
/images/index                       24706
/images/dense/*               24706 × 1616
/trigger/telescope                  24706
/trigger/array                      13100
/camera/pixels                       1616
/mirrors/facets                        54
/telescopes/table                       2
```

### 6.2 光学链路

```text
遮挡前到达输出面                    4686708
被结构遮挡                           900876
  入射段                             782613
  反射段                             118263
遮挡后到达输出面                    3785832
命中相机像素区域                    2863103
collector/探测响应后接受            2183277
Cherenkov PE                    1775990.247
```

遮挡后输出面总通过率为
`3785832 / 4686708 = 0.807365`。这里的 `accepted` 是经过 collector 几何和
探测效率抽样后的 represented-photon 数，不等于 PE 权重总和。

### 6.3 图像与 NSB

全量图像共有 `39,924,896 = 24706 × 1616` 个像素值：

```text
Cherenkov PE 总和                 1775990.246613
NSB PE 总和                     189070567
final PE 总和                  190846557.246563
每幅图像 Cherenkov PE 平均值           71.884977
每幅图像 Cherenkov PE 中位数           53.038520
每幅图像 Cherenkov PE 最大值         1363.134529
PE 分量最大绝对闭合误差                 3.814697e-6
```

所有 Cherenkov、NSB 和 final PE 均为有限非负数。全量文件有意不写逐 bin 波形，
以免再增加约 12 亿个稠密样本；它保留稠密的三个图像分量以及最终 trigger。
逐 bin 波形、波形积分闭合和从保存波形重算 trigger 的检查由第一 shower 诊断包完成。

### 6.4 Trigger

```text
单镜 trigger                  10823 / 24706 = 43.8072%
阵列 trigger                   3075 / 13100 = 23.4733%
```

逐条比较 ROOT 和 HDF5 中的 24,706 条单镜 trigger，`failure_count=0`。第一
shower 诊断包还从保存的 31-bin 波形重新计算全部 12 条 trigger，结果同样
`failure_count=0`。

### 6.5 文件、大小与 SHA-256

```text
全量 corsika_trace.h5
  1124582303 byte
  e64255e0882e2c2d0b2277c2a5f089caf8fc445f67e98e0ddf7ac1656a0aeeda

全量 lact_events.root
  209898581 byte
  2f11d1f7fc37719374b3bcf0a8888c3d6b492a575d0be47ede3813f2e4475610

全量 validation_summary.json
  0c3d22a8426f95ef094a3fcfec041074dd4445d5618e86bbf2f14e16bbeb5997

第一 shower corsika_trace.h5
  3d3b52b11d3ab429810518dd27bd2b773552a761f78ee88ddc399358794531d2

第一 shower lact_events.root
  b1e02a512df2e1db882f917d5f68fbc415bffc79fafaabaecdc4bf516a700ca6
```

下载后在本地重新计算的哈希与服务器端完全一致。

## 7. 自动检查

运行：

```bash
python tools/summarize_trust_validation.py \
  results/trust_validation/muon_full_1000_final/corsika_trace.h5 \
  --root results/trust_validation/muon_full_1000_final/lact_events.root \
  --stdout-log results/trust_validation/muon_full_1000_final/stdout.txt \
  --pretty
```

检查内容包括：

- HDF5 各数据集形状；
- 图像、分量及波形的有限性和非负性；
- `final PE = Cherenkov PE + NSB PE`；
- telescope/array trigger 数量；
- ROOT/HDF5 逐条 trigger 一致性；
- 若保存了波形，则从保存波形重新计算 trigger；
- 输入路径和 SHA-256；
- 日志中的逐 shower 光学计数总和。

全量报告的 `sanity.passed=true`，ROOT/HDF5 trigger 比较
`checked=24706, failure_count=0`。第一 shower 报告的保存波形 trigger 比较
`checked=12, failure_count=0`。

## 8. 完整运行反向发现并修复的问题

这次运行不是单纯“出一个文件”。完整规模和双后端比较实际发现了三个实现级问题：

1. NSB 原先在图像/波形与 trigger/ROOT 两侧独立抽样，统计分布相同，但同一事件不是
   同一随机实现，因而无法从保存波形复算 trigger。现在只生成一次并在所有消费者间
   共用。
2. ROOT 动态树原先脱离 `TFile`，中途 `FlushBaskets/AutoSave` 没有真正流式写入，
   全量结束时会形成超过 1 GiB 的单对象并损坏文件。现在树绑定到文件、检查每次写入
   返回值、结束后安全解绑；全量 ROOT 已能在运行中增长并正常重新打开。
3. HDF5 图像时间曾用 `float`，ROOT 用 `double`。在时间 bin 边界附近，这会令极少数
   trigger 发生翻转。现在两端关键时间字段统一保存为 `double`，最终 24,706 条比较
   为零差异。

这些修复不是为了模仿 sim_telarray 的坐标或数据结构，而是为了满足 LACT_sim 自己的
基本可信条件：同一随机事件只能有一个实现、保存结果可复算、长任务输出可恢复读取、
不同输出后端不改变物理判断。

## 9. 如何检验

建议按两个层次检查：

1. 先打开 `results/trust_validation/muon_full_1000_final/validation_summary.json`，
   检查 `sanity.passed` 和 `root_hdf5_trigger.passed`；
2. 用全量 HDF5 做 shower、图像、PE 和 trigger 统计，用 ROOT 做独立后端读取；
3. 用 `muon_first_event_canonical_nsb` 中的 `trace_summary.csv`、
   `mirror_diagnostic.csv`、`collector_photons.csv` 和稀疏波形从输入到输出逐级追查；
4. 修改配置后重跑同一输入，比较机器可读报告，而不是只比较日志文本。

这套结果证明的是当前代码链的内部一致性、确定性和输出完整性。它仍不能替代望远镜
实测焦斑、离轴 PSF、collector angular acceptance、总探测效率和电子学标定数据。
