# LACT Event ROOT 输出设计

本文设计 LACT_sim 新增的重建接口输出。这个格式只保存 LACT_sim 真正产生、并且 pylast 后续重建需要的内容；不复刻 sim_telarray EventIO block，也不把 pylast 的 `DL1/DL2` 结果层级提前写进 LACT_sim。

核心原则：

```text
LACT_sim 文件保存 detector response 输入
pylast adapter 把它转成 ArrayEvent
pylast 后续再产生 DL1/DL2
```

因此本文件里不设计 `/events/dl1` 和 `/events/dl2`。它们不是 LACT_sim 的输出职责。

## 设计目标

1. 保存 pylast 重建前真正需要的 LACT_sim detector response。
2. 默认小文件、低内存、可流式写入。
3. 静态信息只写一次，event/telescope 数据按行 sparse 存储。
4. waveform、响应分量、debug 信息全部可选。
5. 现有 HDF5 继续保留，不改变其作为 LACT_sim 当前分析输出的语义。

## 推荐格式

```text
schema: lact_event_root_v1
container: ROOT TFile
storage: columnar TTrees + RVec variable-length branches
default profile: reco_min
```

推荐新增 reader：

```text
LactEventSource
```

它的职责是读取 `lact_event_root_v1`，按需构造 pylast 的 `ArrayEvent`。这个 reader 是接口层，不是 LACT_sim 主程序对 pylast 的依赖。

## 总体层级

新的层级按 LACT_sim 数据生命周期组织，而不是按 pylast 后续结果组织：

```text
/file
/config
/subarray
/events
/observations
/waveforms        optional
/response         optional
/debug            optional
```

含义：

```text
/file
  文件格式、单位、坐标约定、profile。

/config
  本次 LACT_sim 运行配置和关键组件摘要。

/subarray
  静态望远镜、相机、光学信息。整文件只写一次。

/events
  array event / shower 级别信息。每个 event 一行。

/observations
  event-telescope 级别相机图像。默认核心输出。

/waveforms
  event-telescope 级别 p.e. 时间序列。可选。

/response
  LACT_sim 响应分量。可选，只用于检查。

/debug
  追迹统计和输入 provenance 摘要。可选。
```

## 为什么这样更省

不使用 dense cube：

```text
event x telescope x pixel
event x telescope x time_bin x pixel
```

默认使用 event/telescope 行：

```text
一行 = 一个 event 中一个真正有输出的 telescope
```

图像和 waveform 使用 variable-length arrays：

```text
只写触发 telescope
只写有效 pixel
只写非零 waveform bin
静态 camera geometry 只写一次
```

这样批量读取时也可以只读需要的列，例如只读 `image_pe`，不读 waveform/debug/config 大字段。

## 本版相对上一版的变化

上一版设计更接近 pylast ROOT data-level 容器，包含 `/events/dl0`、`/events/dl1`、`/events/dl2` 等层级。当前版本改为 LACT_sim 自己的 detector-response schema：

```text
旧思路:
  /events/dl0
  /events/dl1
  /events/dl2
  /events/r1
  /events/r0

新思路:
  /events/core
  /observations/images
  /waveforms/pe optional
  /response optional
  /debug optional
```

主要变化：

```text
删除 LACT_sim 不产生的 DL1/DL2 结果层级
把图像从 pylast 命名 /events/dl0 改为 LACT 命名 /observations/images
把 waveform 从默认 dense 设想改为默认 sparse COO
把相机几何从 telescope/event 数据中归一化到 /subarray/cameras
把 response/debug/provenance 全部改成 profile 可选
```

## 预期存储和内存改善

以下估算用变量表示：

```text
E = array events 数
T = 每个 event 有输出的 telescope 数
P = camera pixel 数，new_camera 当前约 1616
B = waveform time bins 数，例如 25 或 70
f = sparse 图像非零或保存 pixel 占比
g = sparse waveform 非零 pixel-bin 占比
```

### 图像层

dense 图像保存：

```text
E * T * P * 4 bytes
```

sparse 图像保存近似：

```text
E * T * (f * P) * (4 bytes image_pe + 4 bytes pixel_id)
```

因此 sparse 图像相对 dense 的大小约为：

```text
2f
```

例如：

```text
f = 0.05  -> 约 dense 的 10%
f = 0.10  -> 约 dense 的 20%
f = 0.25  -> 约 dense 的 50%
```

如果图像几乎全相机都有 NSB 非零值，sparse 不再占优；这时 writer 应该自动切换 `storage=dense`，或者只对 Cherenkov/trigger-clean 前后的目标分量使用 sparse。

### waveform 层

dense waveform 保存：

```text
E * T * P * B * 4 bytes
```

sparse COO waveform 保存近似：

```text
E * T * (g * P * B) * (4 bytes pe + 4 bytes pixel_id + 2 bytes time_bin)
```

因此 sparse COO 相对 dense 的大小约为：

```text
2.5g
```

例如：

```text
g = 0.01  -> 约 dense 的 2.5%
g = 0.05  -> 约 dense 的 12.5%
g = 0.10  -> 约 dense 的 25%
```

如果开启 NSB 并且每个 time bin 都产生大量非零 p.e.，waveform 可能接近全满；这时 COO 不一定省，writer 应按占比自动选择 dense 或直接关闭 waveform profile。

### 静态几何归一化

旧式做法如果在每个 event/telescope 行重复 camera geometry，几何体量近似：

```text
E * T * P * geometry_fields
```

当前设计只写：

```text
unique_camera_types * P * geometry_fields
```

对于同一种 camera 被多台望远镜和多事件复用的情况，几何存储从随事件增长变为常数级。运行越大，收益越明显。

### 内存峰值

不推荐的 dense 聚合内存：

```text
image:    E * T * P * 4 bytes
waveform: E * T * P * B * 4 bytes
```

当前流式写入内存：

```text
image:    T_event * saved_pixels * row buffers
waveform: T_event * saved_pixel_bins * row buffers
```

也就是从“整场运行规模”降到“单个 event/telescope 或单个 event”的规模。以 `P=1616, B=25` 为例，一个 telescope 的 dense waveform 约：

```text
1616 * 25 * 4 bytes = 158 KiB
```

如果一次保留 32 台 telescope，单 event dense waveform 约：

```text
32 * 158 KiB = 5.0 MiB
```

如果保留整场 `E=1447` 个 event 的 dense waveform，则约：

```text
1447 * 5.0 MiB = 7.1 GiB
```

当前设计要求 writer 流式逐 event 写入，不应在内存中持有这类整场 dense waveform。

### 读性能预期

默认 `reco_min` profile 读 pylast 重建输入时，只需要：

```text
/events/core
/events/pointing
/observations/images
/subarray
```

不会读取：

```text
/waveforms
/response
/debug
```

因此相对“把所有层级写在一个大 dense HDF5/ROOT 数组里”的设计，批量重建的 I/O 量主要由图像行决定，debug 和 waveform 不会拖慢默认分析。

## 必需层级

### `/file/info`

单行 TTree 或 key-value metadata。

```text
schema_name              string   "lact_event_root_v1"
schema_version           int32
creator                  string
created_unix             int64
profile                  string
unit_system              string   "SI_rad_pe"
coordinate_convention    string
```

### `/config/run`

保存本次 LACT_sim 配置摘要。完整配置文本可以压缩成 string，关键字段作为 branch 便于筛选。

```text
run_id                   int32
source_kind              string   corsika_eventio / photon_csv / artificial
electronics_model        string   current default: ideal_pe_placeholder
waveform_enabled         bool
waveform_source          string   none / pe / photon_count
nsb_enabled              bool
trigger_enabled          bool
config_text              string
```

### `/subarray/telescopes`

每台望远镜一行。

```text
tel_id                   int32
x_m                      float64
y_m                      float64
z_m                      float64
focal_length_m           float64
mirror_area_m2           float64
n_pixels                 int32
camera_id                int32
```

`camera_id` 用于复用 camera geometry。多台望远镜如果相机完全相同，不需要重复写像素几何。

### `/subarray/cameras`

每种 camera 一行，像素数组只写一次。

```text
camera_id                int32
camera_name              string
n_pixels                 int32
pixel_id                 RVec<int32>
pixel_x_m                RVec<float32>
pixel_y_m                RVec<float32>
pixel_area_m2            RVec<float32>
pixel_shape              RVec<int16>
pixel_width_m            RVec<float32>
pixel_x_fov_rad          RVec<float32>
pixel_y_fov_rad          RVec<float32>
```

优化点：

```text
相机几何从每个 event/telescope 中移出
相同 camera 只写一次
图像行只引用 tel_id/camera_id
```

### `/events/core`

每个 array event 一行。

```text
run_id                   int32
event_id                 int64
shower_id                int64
array_id                 int32
primary_id               int32
energy_tev               float64
alt_rad                  float64
az_rad                   float64
core_x_m                 float64
core_y_m                 float64
x_max_g_cm2              float64
h_first_int_m            float64
h_max_m                  float64
triggered                bool
n_tel_observations       int32
tel_ids                  RVec<int32>
```

映射到 pylast：

```text
ArrayEvent.event_id
ArrayEvent.run_id
ArrayEvent.simulation.shower
```

### `/events/pointing`

如果一个 event 的所有望远镜 pointing 相同，可以只写 array pointing；如果有 per-telescope pointing，则写 telescope pointing 行。

```text
run_id                   int32
event_id                 int64
array_az_rad             float64
array_alt_rad            float64
tel_id                   int32     -1 means array-level only
tel_az_rad               float64
tel_alt_rad              float64
```

映射到 pylast：

```text
ArrayEvent.pointing
```

## 默认观测层

### `/observations/images`

这是默认核心输出。每行是一个 event/telescope 图像。

```text
run_id                   int32
event_id                 int64
tel_id                   int32
camera_id                int32
n_pixels_total           int32
storage                  uint8     0=dense, 1=sparse
pixel_id                 RVec<int32>
image_pe                 RVec<float32>
peak_time_ns             RVec<float32>
image_sum_pe             float32
time_reference_ns        float64
image_status             uint16
```

默认推荐：

```text
storage = sparse
pixel_id = 非零或超过阈值的 pixel
image_pe = 对应 pixel 的 integrated p.e.
peak_time_ns = 对应 pixel 的时间；没有 timing 时可为空
```

如果某个任务明确需要全像素向量，可以设置：

```text
storage = dense
pixel_id = 省略或保存完整 pixel axis
image_pe.size == n_pixels_total
```

但是批量默认不建议 dense。

映射到 pylast：

```text
LactEventSource 在读取时根据 camera geometry 展开成 DL0Camera.image
DL0Camera.peak_time
```

注意：文件本身不叫 `/events/dl0`，因为这是 LACT_sim 输出，不是 pylast 内部数据层。adapter 可以把它映射成 pylast `dl0`。

## 可选 waveform 层

### `/waveforms/pe`

当需要保存你现在做的时间 p.e. 序列时启用。每行是一个 event/telescope waveform。

```text
run_id                   int32
event_id                 int64
tel_id                   int32
camera_id                int32
n_pixels_total           int32
n_time_bins              int32
time_bin_width_ns        float32
time_start_ns            float64
time_reference_ns        float64
storage                  uint8     0=dense_pixel_time, 1=sparse_coo
pixel_id                 RVec<int32>
time_bin                 RVec<int16>
pe                       RVec<float32>
```

默认推荐：

```text
storage = sparse_coo
pixel_id/time_bin/pe 只保存非零 bin
```

如果 waveform 近似全满，writer 可以自动切换为 dense：

```text
storage = dense_pixel_time
pixel_id = 完整 pixel axis 或空
time_bin = 空
pe = flattened [pixel, time_bin]
```

映射到 pylast：

```text
proxy p.e. waveform -> ArrayEvent.r1
```

但是必须设置 metadata：

```text
waveform_semantics = proxy_pe
not_raw_adc = true
```

这样不会被误当成 `R0 ADC`。

## 可选 trigger 层

### `/events/trigger_tel`

只有需要 trigger 细节时写。默认 `events/core.tel_ids` 已经足够表达哪些 telescope 有观测。

```text
run_id                   int32
event_id                 int64
tel_id                   int32
triggered                bool
trigger_time_ns          float64
n_trigger_pixels         int32
trigger_pixel_ids        RVec<int32>
```

## 可选 response 层

### `/response/image_components`

只用于响应检查，不作为 pylast 默认输入。

```text
run_id                   int32
event_id                 int64
tel_id                   int32
storage                  uint8
pixel_id                 RVec<int32>
cherenkov_pe             RVec<float32>
nsb_pe                   RVec<float32>
photon_count             RVec<uint16>
```

### `/response/waveform_components`

只用于 waveform 检查。

```text
run_id                   int32
event_id                 int64
tel_id                   int32
storage                  uint8
pixel_id                 RVec<int32>
time_bin                 RVec<int16>
cherenkov_pe             RVec<float32>
nsb_pe                   RVec<float32>
```

## 可选 debug/provenance 层

### `/debug/trace_summary`

每个 event/telescope 一行，保存计数和耗时相关摘要，不保存完整 photon bunch。

```text
run_id                   int32
event_id                 int64
tel_id                   int32
input_bunches            uint32
input_photons            float64
mirror_hits              uint32
camera_hits              uint32
accepted_pe              float64
first_cherenkov_time_ns  float64
```

### `/debug/input_eventio`

如果输入来自 CORSIKA/EventIO，只保存摘要。

```text
run_id                   int32
event_id                 int64
shower_event_id          int64
array_id                 int32
eventio_block_flags      uint64
raw_bunch_count          uint32
raw_tel_count            uint16
```

不保存原始 EventIO block payload。

## Profiles

### `reco_min`

默认生产模式，最小可重建输入。

```text
/file/info
/config/run
/subarray/telescopes
/subarray/cameras
/events/core
/events/pointing
/observations/images
```

特点：

```text
只保存 integrated p.e. image
默认 sparse
低内存流式写
pylast 从 DL0 等价层开始
```

### `reco_waveform`

需要 pylast 从 waveform extraction 开始时使用。

```text
reco_min
+ /waveforms/pe
```

特点：

```text
保存 p.e. 时间序列
adapter 映射到 pylast R1
比 reco_min 大，但仍不保存 raw ADC
```

### `audit`

验证 LACT_sim 响应时使用。

```text
reco_waveform
+ /response/image_components
+ /response/waveform_components
+ /debug/trace_summary
```

特点：

```text
方便检查 Cherenkov/NSB/组件贡献
不建议批量生产默认开启
```

### `provenance`

需要追溯输入 EventIO 摘要时使用。

```text
reco_min
+ /debug/input_eventio
+ /debug/trace_summary
```

特点：

```text
保存来源摘要
不保存完整 EventIO block
```

### `electronics_adc`

未来真实电子学模式。当前不实现为默认 profile。

```text
reco_waveform
+ /electronics/r0_adc
+ /electronics/calibration
```

只有当 LACT_sim 真的生成 ADC waveform、pedestal、gain、digitization 时才启用。

## 内存和写入策略

writer 应该按 event 流式写入：

```text
for each event:
  写 /events/core 一行
  for each telescope with useful output:
    写 /observations/images 一行
    可选写 /waveforms/pe 一行
    可选写 response/debug 一行
```

不要在内存中构造整场运行的：

```text
all_events x all_telescopes x all_pixels
all_events x all_telescopes x all_pixels x all_time_bins
```

推荐压缩：

```text
ROOT compression: ZSTD
image_pe: float32
peak_time_ns: float32
photon_count: uint16 or uint32
time_bin: int16 if n_time_bins < 32768
```

默认 sparse 阈值：

```text
image_pe > 0
或 image_pe >= configured_min_saved_pe
```

是否允许阈值丢弃必须写入 metadata，避免误认为完整图像。

## 配置建议

默认：

```ini
output.format=root
output.root_path=run_logs/my_run/lact_events.root
output.root_schema=lact_event_root_v1
output.root_profile=reco_min
output.root_image_storage=sparse
output.root_min_saved_pe=0
output.root_store_waveform=false
output.root_store_response=false
output.root_store_debug=false
```

保存 p.e. 时间序列：

```ini
output.root_profile=reco_waveform
output.root_store_waveform=true
output.root_waveform_storage=sparse_coo

waveform.enabled=true
waveform.source=pe
waveform.time_reference=image_first
waveform.time_bin_width_ns=1
waveform.time_window_start_ns=-5
waveform.time_window_end_ns=20
```

响应检查：

```ini
output.root_profile=audit
output.root_store_response=true
output.root_store_debug=true
```

## LactEventSource 输入接口

`LactEventSource` 不要求文件里有 pylast 的 `/events/dl0`、`/events/dl1`、`/events/dl2` 目录。它读取 LACT schema，再构造 pylast 对象。

映射规则：

```text
/events/core
  -> ArrayEvent.event_id
  -> ArrayEvent.run_id
  -> ArrayEvent.simulation.shower

/subarray/*
  -> SubarrayDescription

/events/pointing
  -> ArrayEvent.pointing

/observations/images
  -> ArrayEvent.dl0.tels[tel_id].image
  -> ArrayEvent.dl0.tels[tel_id].peak_time

/waveforms/pe
  -> ArrayEvent.r1.tels[tel_id].waveform
  -> ArrayEvent.r1.tels[tel_id].gain_selection = 0
```

当 sparse image 被读入 pylast 时，reader 负责展开成完整相机长度：

```text
full_image = zeros(n_pixels_total)
full_image[pixel_id] = image_pe
```

当 sparse waveform 被读入 pylast R1 时，reader 负责展开成：

```text
waveform[n_pixels_total, n_time_bins]
```

这个展开只发生在当前 event/telescope 读取时，不需要把整文件展开到内存。

## 与 pylast/RootEvent 的关系

本格式不追求文件内层级直接等于 pylast `RootEventSource`。原因是 pylast ROOT 文件是 pylast 自己的中间/结果容器，而 LACT_sim 文件是 detector response 输入容器。

推荐边界：

```text
LACT_sim writes lact_event_root_v1
LactEventSource reads lact_event_root_v1
pylast modules produce DL1/DL2
pylast DataWriter can write pylast RootEvent output if needed
```

这样 LACT_sim 不需要保存它不产生的 `DL1/DL2`，也不会因为迁就 pylast ROOT 文件而牺牲存储结构。

## 与现有 HDF5 的关系

现有 HDF5 继续用于 LACT_sim 自身检查和数组分析：

```text
/images/dense
/waveforms
/trigger
/metadata
```

新的 ROOT 输出用于重建接口：

```text
/events/core
/observations/images
/waveforms/pe optional
/subarray
```

两者可以同时存在。ROOT 输出不要求复用 HDF5 的 dense 布局。

## 结论

推荐默认实现：

```text
lact_event_root_v1 + reco_min
```

默认只写：

```text
run/config 摘要
subarray/camera 静态信息
event/shower truth
pointing
sparse event-telescope image_pe + peak_time
```

可选写：

```text
sparse p.e. waveform
response components
trace/debug summaries
future electronics ADC
```

不写：

```text
pylast DL1
pylast DL2
完整 sim_telarray EventIO block
默认 R0 ADC
```

这个设计的优化点在于：静态信息归一化、事件观测 sparse 化、waveform COO 化、debug/provenance profile 化，并且 reader 只在读取当前 event 时展开成 pylast 需要的完整相机数组。
