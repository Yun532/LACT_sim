# S17351 微单元索引与饱和模型

## 规格映射

- 8 个通道，物理排布为 2 × 4。
- 每通道尺寸为 6.6 × 3.2 mm²。
- 每通道包含 264 × 128 = 33792 个、pitch 为 25 μm 的微单元。
- 相邻通道的 x/y 间隔均为 0.2 mm。
- 8 通道包围框为 13.4 × 13.4 mm²。
- 每个相机 pixel 总计 270336 个微单元。

`s17351_tiled_2x4` 遵循规格书正面图：左列从上到下为 B-1 至
B-4，右列从上到下为 A-1 至 A-4。

## 采用的物理定义

厂家 PDE 直接作为器件的有效 PDE 使用，认为其中已经包含内部几何填充
因子的平均影响。程序不会把一个 25 μm 微单元再划分为人为指定的感光区
和死区，也不会执行 `PDE / fill_factor`。

25 μm pitch 只定义微单元索引。一个已经通过 PDE 抽样的 primary p.e.
根据其 SiPM 坐标分配到对应微单元。关闭恢复时，同一
event/telescope/pixel 内每个微单元最多产生一次 avalanche；开启恢复时，
重复命中的电荷等效 p.e. 按指数恢复比例计算。

明确给出的 0.2 mm 通道间隔仍作为宏观无传感器区域单独处理，它不等于
微单元内部的几何填充因子。

```text
collector exit photon
  -> specified 0.2 mm inter-channel gap
  -> datasheet PDE (applied once)
  -> primary p.e.
  -> 25 um-pitch microcell assignment
  -> optional no-recovery occupancy saturation
  -> fired p.e.
  -> optional single-p.e. waveform
  -> 4 ns samples
```

## 8 通道 gap 是否已经包含在 PDE 中

这个问题与微单元内部的 47% 几何填充因子无关。程序提供独立开关：

```ini
microcell.pde_includes_inter_channel_gaps=true
```

- `false`：把输入 PDE 理解为通道有效区域上的 PDE，随后显式的
  0.2 mm 通道 gap 会额外造成损失。
- `true`（默认）：把输入 PDE 理解为已经对整个 13.4 x 13.4 mm 封装面积平均过。
  程序先除以已知通道面积比例，再用真实 gap 几何逐光子拒绝。

已知通道面积比例为：

```text
f_channel = 8 * (6.6 mm * 3.2 mm) / (13.4 mm * 13.4 mm)
          = 0.9409668
```

因此开关打开时，通道内条件 PDE 的修正系数是
`1 / f_channel = 1.0627367`。修正只针对 8 通道之间的宏观 gap；
不会除以 47%，也不会猜测单个微单元内部的感光区形状。

## 配置

```ini
microcell.enabled=true
microcell.saturation_enabled=true
microcell.model=explicit_exponential_recovery
microcell.recovery_enabled=true
microcell.recovery_time_ns=10.0
microcell.layout=s17351_tiled_2x4
microcell.sensor_size_x_m=0.0134
microcell.sensor_size_y_m=0.0134
microcell.channel_size_x_m=0.0066
microcell.channel_size_y_m=0.0032
microcell.channel_gap_x_m=0.0002
microcell.channel_gap_y_m=0.0002
microcell.microcell_columns_per_channel=264
microcell.microcell_rows_per_channel=128
```

`saturation_enabled=false` 时，primary p.e. 不发生微单元占用损失。
`saturation_enabled=true, recovery_enabled=false` 时，重复命中同一微单元
的后续 p.e. 记为 `saturation_rejected`。`recovery_enabled=true` 时采用

```text
f_recovery = 1 - exp(-delta_t / recovery_time_ns)
```

作为该次 avalanche 的电荷等效 `fired_pe`；同一时刻的重复命中仍为零。
当前 10 ns 不是 S17351 实测值：未检索到厂家公开的该型号恢复时间。
Hamamatsu 的 MPPC 技术说明指出恢复常数由结电容与淬灭电阻的 RC 决定，
并给出典型像元恢复约 15 ns 的量级。因此强度干涉模拟默认扫描 1、10、
30 ns，直到获得器件在实际温度和过压下的标定。

测试配置可以打开逐 hit 输出：

```ini
electronics.output.save_primary_sequence=true
electronics.output.save_fired_sequence=true
electronics.output.save_microcell_decisions=true
collector.debug_photon_output=true
```

正式配置默认关闭这些详细量，仅保存相机图、可选波形和触发所需结果。
