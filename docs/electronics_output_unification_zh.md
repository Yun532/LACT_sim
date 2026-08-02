# 电子学结果与 ROOT、HDF5、CSV 输出

## 单次计算原则

`run_corsika_trace` 在每个 `event_id/telescope_id` 上只构造一次
`CameraElectronicsEvent`：

```text
相机接受的 Cherenkov p.e.
  -> 加入可选 NSB
  -> 微单元 gap 与饱和
  -> Fired p.e.
  -> 可选单 p.e. 模板与 4 ns 采样
  -> 相机触发
  -> 阵列时间修正与阵列触发
  -> ROOT / HDF5 / CSV serializers
```

writer 不再各自抽样 NSB或重新计算触发。因此三种格式中的像素图、波形、
触发标志和时间都来自同一个内存对象。

## 时间参考

当配置为：

```ini
waveform.time_reference=image_first
```

每台望远镜使用该相机第一个 Cherenkov p.e. 的 EventIO 时间作为
`reference_time_ns`。整台相机的所有像素共享这个参考时间和同一组采样 bin；
不同像素只有第一个非零样本的位置不同。

```text
global_time_ns = reference_time_ns + time_center_ns
```

ROOT 的 `observations` 和 `waveforms`、HDF5 的 `images/index` 和
`waveforms/reference_time_ns`、CSV 的相应列都显式保存该值。

## 因果触发时间

触发曲线在当前采样时刻 `t` 只使用已经到达的样本：

```text
[t - coincidence_window, t]
```

`trigger_time_ns` 是第一次确认像素多重数达到要求的采样时刻，不再把结果记到
窗口起点。像素图、PE 序列和波形不受此时间定义修正影响。

阵列时间字段为：

```text
coincidence_time_ns = trigger_time_ns + geometric_delay_ns
```

## 三格式验证

配置：

```text
configs/examples/corsika_electronics_v2_three_format_validation.cfg
```

交叉验证：

```bash
python tools/validate_format_unification.py \
  validation/electronics_chain_v2/format_unification
```

脚本逐项比较：

- 逐像素 Primary Cherenkov/NSB/dark、Fired Cherenkov/NSB/dark、gap 损失和饱和损失；
- 稀疏波形坐标与每个采样值；
- Primary p.e. 与 Fired p.e. 序列；
- `reference_time_ns`、相机触发和阵列符合时间。

CSV 因为是表格式，完整电子学输出会拆成 `camera_pixels.csv`、
`waveforms.csv`、`triggers.csv`、`primary_pe.csv`、`fired_pe.csv` 和
`trace_summary.csv`；它们与单个 ROOT/HDF5 文件表达同一组逻辑数据。

HDF5 的 `/images/dense/pe` 与 ROOT 的 `image_pe` 都表示最终 Fired 总量，
而不是 Primary 总量。各分量使用完整名称，例如
`primary_cherenkov_pe`、`fired_nsb_pe`、`saturation_lost_pe`，不再用含义
不明确的 `cherenkov_pe`、`nsb_pe` 或 `image_fired_pe` 兼容字段。
