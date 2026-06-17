# LACT ROOT 输出中的 NSB 与触发流程

本文档记录当前 `run_corsika_trace` 里，面向 LACT ROOT/pylast 输出的 NSB 生成、触发判断和写出顺序。

这套逻辑目前属于 **pre-electronics p.e. proxy**：NSB 被加到 p.e. waveform/image 上，用来生成 ROOT waveform 并供 pylast adapter 读取。它还不是完整电子学链路里的 ADC waveform。后续如果真实电子学模型接入，NSB 应迁移到电子学前端，与 Cherenkov p.e. arrival 一起进入 SiPM/电子学响应。

## 相关配置

LACT ROOT example 当前主要使用：

```text
nsb.config=../nsb/spectral_skycalc_dark_with_obstruction.cfg
trigger.config=../trigger/example_simple_multiplicity.cfg

output.save_only_triggered=true
output.lact_profile=timeseries_pe
output.lact_root_write_components=false
output.lact_root_auto_flush_mb=200
output.lact_root_flush_events=20

waveform.enabled=true
waveform.source=pe
waveform.time_reference=image_first
waveform.time_bin_width_ns=2
waveform.time_window_start_ns=-2
waveform.time_window_end_ns=60

trigger.pixel_threshold_pe=10
trigger.coincidence_window_ns=10
```

其中：

- `waveform.time_window_start_ns/end_ns` 是最终写入 ROOT waveform 的完整时间窗口。
- `trigger.coincidence_window_ns` 是触发判断使用的滑动窗口宽度。
- `trigger.pixel_threshold_pe` 是单像素触发阈值。
- `trigger.camera_multiplicity` 来自 trigger config，表示单台望远镜中至少多少个像素过阈值。
- `trigger.array_multiplicity` 来自 trigger config，表示一个 event 至少多少台望远镜触发。

当前 example 中完整 waveform 窗口是 `-2..60 ns`，bin 宽是 `2 ns`，即 31 个 time bins。触发窗口是 `10 ns`，即 5 个 bins。

## NSB rate 的计算

程序启动时读取并展开 `nsb.config`。对于常用的：

```text
nsb.model=spectral_flux
```

程序调用：

```cpp
resolveNsbSpectralRate(...)
```

它会读取 NSB 光谱文件，并结合：

- 光学效率曲线
- `nsb.effective_area_m2`
- 像素立体角 `nsb.pixel_solid_angle_sr`

得到：

```text
nsb.rate_pe_per_ns_per_pixel
```

这个 rate 是每个像素、每 ns 的平均 p.e. rate。当前模型假设所有 telescope 使用同一个平均 NSB rate，但每个 `event_id + telescope_id + pixel/bin cell` 使用不同随机种子，因此具体涨落不同。

## 程序启动时的纯 NSB 误触发估计

程序开始光追前会打印纯 NSB 误触发概率估计：

```text
Pure NSB trigger estimate
```

输出包括：

```text
rate_pe_per_ns_per_pixel
effective_window_ns
pixel_threshold_pe
camera_multiplicity
array_multiplicity
n_pixels
n_telescopes
pixel_prob_ge_threshold
camera_single_window_prob
camera_sliding_upper_prob
array_sliding_upper_prob
```

计算逻辑是：

```text
lambda = nsb.rate_pe_per_ns_per_pixel * effective_window_ns
pixel_prob_ge_threshold = P(Poisson(lambda) >= pixel_threshold_pe)
camera_single_window_prob = P(至少 camera_multiplicity 个像素过阈值)
camera_sliding_upper_prob ≈ n_windows * camera_single_window_prob
array_sliding_upper_prob = P(至少 array_multiplicity 台 telescope 触发)
```

这个估计用于判断“纯 NSB 触发的望远镜”是否可以忽略。当前 10 ns trigger window 下，纯 NSB 触发概率极低，因此后续流程只在已有 Cherenkov 信号的 telescope 上考虑 NSB。

## 光追阶段

主循环在 `apps/run_corsika_trace.cpp` 中流式读取 EventIO photon bunch，不预加载全文件。

每个 photon bunch 经过：

```text
EventIO bunch
-> 坐标变换
-> 大气透过
-> 镜面光追
-> 遮挡判断
-> 相机/light collector/SIPM p.e. 响应
```

如果命中相机，则累积 Cherenkov 信号到：

```text
summaries
pixels
waveforms
raw_waveform_hits
```

当前这一步只处理 Cherenkov，不加 NSB，也不做 NSB trigger。

## event 结束后的 ROOT writer 流程

当 `event_id` 切换时，程序调用：

```cpp
lact_root_stream_writer->writeEvent(...)
```

然后进入：

```cpp
prepareLactRootObservations(...)
```

该函数对当前 event 中有 Cherenkov 信息的 telescope 逐台处理。

## Cherenkov-only waveform 构建

对于 `timeseries_pe`，先只用 Cherenkov 构建 waveform：

```text
raw_waveform_hits / waveforms
-> 按 reference_time 分配到 time bin
-> 累积到 waveform_pe
-> 同时累积 image_cherenkov_pe 和 image_pe
```

当前 `time_reference=image_first`，所以相对当前 telescope 的第一束 Cherenkov 时间对齐。

此时：

```text
waveform_pe = Cherenkov only
image_pe = Cherenkov only
image_cherenkov_pe = Cherenkov only
```

## Cherenkov-only 预触发

在未加入 NSB 前，程序先对 Cherenkov-only waveform 做一次 trigger 判断：

```text
在完整 waveform 窗口内滑动 trigger.coincidence_window_ns 窗口
对每个窗口统计 pe >= trigger.pixel_threshold_pe 的像素数
取最大像素数作为该 telescope 的 Cherenkov-only trigger multiplicity
```

如果：

```text
Cherenkov-only multiplicity >= trigger.camera_multiplicity
```

则该 telescope 本身已经触发。它后续一定需要完整 waveform，因此直接进入完整 NSB 生成。

## 轻量 NSB trigger

如果 Cherenkov-only 未触发，但 telescope 有 Cherenkov peak，则程序只在 Cherenkov peak 附近做轻量 NSB trigger。

候选窗口范围是所有包含 `camera_peak_bin` 的 trigger windows：

```text
first_window = max(0, camera_peak_bin - trigger_window_bins + 1)
last_window = camera_peak_bin
```

对这些窗口的 union bins，程序逐 pixel/bin 抽 NSB：

```cpp
sampleTimeBinnedNsbPeCell(...)
```

然后对每个候选 trigger window 计算：

```text
Cherenkov pe in window + light NSB pe in window
```

如果窗口内过阈值像素数达到：

```text
trigger.camera_multiplicity
```

则这个 telescope 被认为可以由 Cherenkov + NSB 触发，后续进入完整 NSB 生成。

如果仍然不触发，则：

```text
不生成完整 NSB waveform
不保存该 telescope
```

这一步是主要优化点：避免给大量最后不会保存的 telescope 生成完整 NSB waveform。

## 严格一致的 NSB 抽样

为了避免“轻量 trigger 用一份 NSB，最终 waveform 又是另一份 NSB”的不一致，当前实现使用 cell 级确定性抽样。

新增函数：

```cpp
sampleTimeBinnedNsbPeCell(
    nsb,
    waveform_cfg,
    event_id,
    telescope_id,
    n_pixels,
    n_bins,
    col,
    bin)
```

它使用：

```text
nsb.seed
event_id
telescope_id
cell = col * n_bins + bin
```

混合成随机种子，然后抽：

```text
Poisson(nsb.rate_pe_per_ns_per_pixel * waveform.time_bin_width_ns)
```

因此同一个：

```text
event_id + telescope_id + pixel/bin
```

在轻量 trigger 和完整 waveform 阶段会得到完全相同的 NSB p.e.。

## 完整 NSB waveform 生成

只有以下情况会生成完整 NSB waveform：

1. `output.save_only_triggered=false`
2. trigger 未启用
3. Cherenkov-only 已经触发
4. Cherenkov-only 未触发，但 peak 附近轻量 NSB trigger 后触发

完整 NSB 生成遍历所有 pixel/bin：

```text
for pixel in camera pixels:
  for bin in waveform bins:
    nsb_pe = sampleTimeBinnedNsbPeCell(...)
    waveform_pe += nsb_pe
    image_pe += nsb_pe
    image_nsb_pe += nsb_pe
```

由于使用同一个 cell 级抽样函数，完整 waveform 会包含轻量 trigger 阶段已经判断过的同一份 NSB。

## 最终 trigger 与保存

完整 NSB 加完后，程序重新对最终 waveform 做 trigger 判断：

```text
Cherenkov + NSB waveform
-> 10 ns sliding window
-> 统计最大过阈值像素数
```

得到：

```text
obs.n_pixels_above_threshold
obs.triggered
obs.trigger_time_ns
```

然后统计同一 `event_id` 下有多少台 telescope triggered：

```text
triggered_telescopes_by_event[event_id]
```

如果：

```text
output.save_only_triggered=true
trigger.enabled=true
```

则只保存满足以下条件的 telescope：

```text
该 telescope 自身 triggered
且同一 event_id 下 triggered telescope 数 >= trigger.array_multiplicity
```

否则该 telescope 的 observation/waveform 不写入 ROOT。

## 当前流程总结

当前顺序可以概括为：

```text
光追 Cherenkov
-> 构建 Cherenkov-only waveform
-> Cherenkov-only trigger
-> 如果已触发：生成完整 NSB waveform
-> 如果未触发：只在 Cherenkov peak 附近做轻量 NSB trigger
-> 如果轻量后仍未触发：丢弃 telescope
-> 如果轻量后触发：生成严格一致的完整 NSB waveform
-> 最终 Cherenkov+NSB trigger
-> array trigger / save_only_triggered 过滤
-> 写 ROOT observations/waveforms
```

## 优化效果预期

旧逻辑是：

```text
所有 telescope 都先生成完整 NSB waveform
再判断 trigger
最后丢掉未触发 telescope
```

新逻辑避免了对未触发 telescope 的完整 NSB waveform 生成。对于大多数只有少数 telescope 被 Cherenkov 或 Cherenkov+NSB 触发的 event，完整 NSB 生成量会近似按“保存 telescope 数 / 有 Cherenkov telescope 数”下降。

例如如果一个 event 有 30 台 telescope 有少量 Cherenkov 信息，但最终只有 3 台需要保存，则完整 NSB waveform 生成量约降到原来的 10% 左右。轻量 trigger 只检查 Cherenkov peak 附近的 10 ns 窗口，成本远低于完整 `-2..60 ns` waveform 生成。

## 当前限制

1. 当前 NSB 仍是 p.e. proxy 层级，不是真电子学 ADC waveform。
2. 轻量 trigger 只在 Cherenkov peak 附近检查，不扫描完整 waveform 的所有窗口。
3. 纯 NSB 触发默认通过概率估计忽略，不为完全没有 Cherenkov 的 telescope 生成 NSB trigger。
4. 所有 telescope 当前使用同一个平均 NSB rate，没有 telescope/pixel 位置依赖。
5. 后续完整电子学接入后，NSB generator 应下沉到电子学前端。

