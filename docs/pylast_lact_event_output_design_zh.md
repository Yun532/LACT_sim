# LACT_sim 到 pylast 的事件输出设计

本文定义一个面向 pylast 读取的 LACT_sim detector-response 输出格式。它不是
sim_telarray EventIO 的复制格式，也不是 pylast 后续 DL1/DL2 的结果文件。它只保存
LACT_sim 已经产生、且 pylast 后续重建需要读取的事件输入。

## 1. 目标和边界

目标：

```text
LACT_sim optical/camera response
  -> lact_event 文件
  -> pylast LactEventSource
  -> pylast ArrayEvent
  -> pylast image cleaning / Hillas / reconstruction
```

LACT_sim 当前有两类相机输出，需要同时支持：

```text
1. integrated p.e. image
   每个 event/telescope/pixel 一个积分后的 p.e. 值。

2. time-series p.e.
   每个 event/telescope/pixel 有按时间 bin 累积的 p.e. 序列。
```

非目标：

```text
不复刻完整 EventIO block payload
不伪造 raw ADC R0
不在 LACT_sim 中产出 DL1/DL2
不把逐光子 debug 命中作为默认重建输入
```

## 2. 推荐文件格式

推荐新增 ROOT 输出：

```text
schema_name    = lact_event_root
schema_version = 1
container      = ROOT TFile
storage        = TTrees + variable-length branches
default profile = image_pe
```

原因：

```text
pylast 已有 ROOT 读写基础
event/telescope 行可以自然保存 variable-length image/waveform
静态 subarray/camera/optics 信息只写一次
后续可用 uproot/ROOT 快速检查
```

现有 HDF5 继续作为 LACT_sim 分析和调试格式保留。`lact_event_root_v1` 是给
pylast adapter 的正式输入契约。

## 3. Profile

### image_pe

默认最小重建格式。

保存：

```text
文件元信息
CORSIKA/shower truth
subarray / telescope / camera / optics
event-level pointing
event-telescope integrated p.e. image
trigger summary
trace summary
```

不保存：

```text
time-series waveform
逐光子 whiteboard/debug 命中
response component 分解
```

### timeseries_pe

时间序列格式。

在 `image_pe` 基础上增加：

```text
time bin definition
event-telescope-pixel p.e. waveform
reference time
per-pixel time_mean / time_rms / time_peak
```

integrated p.e. image 仍然必须保存，作为 waveform 沿时间积分后的快速入口和
校验量。

### debug_full

调试格式。

在 `timeseries_pe` 基础上可选增加：

```text
Cherenkov / NSB / response components
collector debug rows
whiteboard hit rows
obstruction counters
input photon-bunch counters
```

此 profile 不作为 pylast 常规重建输入。

## 4. ROOT 文件结构

建议 TFile 顶层对象：

```text
TNamed /file/schema_name
TNamed /file/schema_version
TNamed /file/producer
TNamed /file/profile
TNamed /file/config_text

TTree config
TTree corsika_events
TTree telescopes
TTree camera_pixels
TTree optics
TTree observations
TTree waveforms          optional
TTree trace_summary      optional
TTree response_components optional
```

ROOT 中不需要真的建目录层级；上面名字是逻辑层级。实现时可以用扁平 TTree 名称。

## 5. 静态信息

### config

一行，保存文件级元信息。

字段：

```text
schema_name                 string
schema_version              int
profile                     string
producer                    string
producer_version            string
source_kind                 string      # EventIO / PhotonCsv / Synthetic
source_path                 string
source_sha256               string      # optional
run_id                      int
coordinate_convention       string
event_id_mode               string
config_text                 string
expanded_config_text        string
```

### telescopes

每台 telescope 一行。

字段：

```text
telescope_id                int
name                        string
array_x_north_m             double
array_y_west_m              double
array_z_up_m                double
radius_m                    double
pointing_az_deg             double
pointing_el_deg             double
camera_id                   int
optics_id                   int
```

坐标约定：

```text
array_x_north_m : CORSIKA magnetic North positive
array_y_west_m  : West positive
array_z_up_m    : Up positive
pointing_az_deg : North toward East
pointing_el_deg : above horizon
```

### camera_pixels

每个 camera pixel 一行。

字段：

```text
camera_id                   int
pixel_id                    int
x_m                         double
y_m                         double
size_m                      double
shape_code                  int         # 1 square, 2 hexagonal, 3 circular
```

### optics

每个 optics type 一行。

字段：

```text
optics_id                   int
name                        string
num_mirrors                 int
mirror_area_m2              double
equivalent_focal_length_m   double
effective_focal_length_m    double
```

## 6. Event Truth

### corsika_events

每个 array event 一行。

字段：

```text
event_id                    long
shower_event_id             int
array_id                    int
run_id                      int
primary_type                int
energy_gev                  double
theta_deg                   double
phi_deg                     double
azimuth_north_to_east_deg   double
altitude_deg                double      # derived, optional but recommended
core_x_north_m              double
core_y_west_m               double
array_rotation_deg          double
h_first_int_m               double      # from simtel MCShower when present, otherwise NaN
x_max_g_cm2                 double      # from simtel MCShower when present, otherwise NaN
h_max_m                     double      # from simtel MCShower when present, otherwise NaN
starting_grammage_g_cm2     double      # from simtel MCShower when present, otherwise NaN
has_simtel_mc_shower        bool
```

说明：

```text
event_id 是 LACT_sim 和 pylast adapter 的主键。
当 event_id_mode=event_array100 时：
  event_id = shower_event_id * 100 + array_id
```

对于多 CSCAT/array offset，`core_x_north_m/core_y_west_m` 保存该 `array_id`
对应的有效 core 坐标，而不是只保存原 shower header 中的默认 core。

当前实现会合并 CORSIKA `EVTH` 和 simtel `MCShower`。`EVTH` 提供 shower id、
core 和 array rotation；`MCShower` 存在时补全 primary、energy、altitude、azimuth、
`h_first_int_m`、`x_max_g_cm2`、`h_max_m`、`starting_grammage_g_cm2`，
并将 `has_simtel_mc_shower=true`。如果输入 EventIO 文件没有 `MCShower` 块，
这些 shower-profile truth 字段保留为 `NaN`。

## 7. Integrated p.e. Image

### observations

每个 event/telescope 一行。默认只写有输出或触发的 telescope；是否只写 triggered
由配置决定，但文件中必须记录策略。

字段：

```text
event_id                    long
telescope_id                int
triggered                   bool
n_pixels_camera             int
n_pixels_saved              int
pixel_id                    vector<int>
image_pe                    vector<float>
image_time_mean_ns          vector<float>  # optional but recommended
image_time_rms_ns           vector<float>  # optional but recommended
image_time_peak_ns          vector<float>  # required for timeseries_pe, optional for image_pe
total_pe                    double
time_first_ns               double
time_mean_ns                double
time_rms_ns                 double
time_peak_ns                double         # camera-summed p.e. waveform peak time, optional
impact_parameter_m          double         # optional
n_pixels_above_threshold    int
trigger_time_ns             double         # 首次达到相机 multiplicity 的窗口中心
trigger_first_time_ns       double         # 与 trigger_time_ns 相同，显式语义字段
trigger_max_multiplicity_time_ns double     # 最大 multiplicity 窗口中心，仅用于诊断
```

`pixel_id` 和 `image_pe` 使用 sparse 表示，只保存需要保存的 pixel。对于
`image_pe` profile，pylast adapter 读入后应重建成完整长度为 `n_pixels_camera`
的 dense vector/matrix，未出现的 pixel 填 0。

如果配置要求 dense 写出，可以把 `pixel_id` 写成完整 `[0, ..., n_pixels-1]`，
这样 reader 逻辑不变。

`image_time_mean_ns` 和 `image_time_rms_ns` 使用 Cherenkov 信号权重的原始
到达时间矩计算，不用波形峰值代替，也不把 NSB p.e. 加入分母。因此负的相对
到达时间是有效值。`image_time_peak_ns` 则来自包含 NSB 的离散 p.e. 波形，
三者描述的是不同物理量。

`image_time_peak_ns` 的定义：

```text
对于每个 pixel，在该 pixel 的 p.e. time-series 中寻找最大 p.e. 的 time bin。
time_peak_ns = reference_time_ns + time_centers_ns[peak_bin]
```

如果多个 time bin 并列最大，writer 使用最早的 bin。若该 pixel 的总 p.e. 为 0，
`image_time_peak_ns` 写 NaN。`timeseries_pe` profile 必须写该字段；`image_pe`
profile 如果没有 waveform，可以省略或写 NaN。

event/telescope 级别的 `time_peak_ns` 定义为相机总 p.e. time-series
`sum_pixel pe[pixel, time_bin]` 的峰值时间，同样使用最早并列峰值。

## 8. Time-Series p.e.

### waveform metadata

可放在 `config` 或单独 `waveform_config` TTree 中。

字段：

```text
waveform_enabled            bool
waveform_source             string      # pe_proxy / cherenkov_pe / total_pe
time_reference              string      # absolute / image_first / image_mean
time_bin_width_ns           double
time_window_start_ns        double
time_window_end_ns          double
n_time_bins                 int
```

### waveforms

每个 event/telescope 一行，使用 sparse COO 保存非零 pixel-bin。

字段：

```text
event_id                    long
telescope_id                int
n_pixels_camera             int
n_time_bins                 int
pixel_id                    vector<int>
time_bin                    vector<unsigned short>
pe                          vector<float>
```

约束：

```text
pixel_id.size == time_bin.size == pe.size
time_bin in [0, n_time_bins)
同一个 event/telescope/pixel/time_bin 不应重复；如重复，reader 可累加但 writer 应避免。
```

`observations.image_pe` 必须等于对应 waveform 对时间积分后的结果，允许浮点误差：

```text
image_pe[pixel] ~= sum_t waveform_pe[pixel, t]
```

`observations.image_time_peak_ns` 必须由同一份 waveform 计算得到。pylast reader
可以直接使用该缓存，也可以重新从 `waveforms` 计算并在 debug 模式下比较差异。

## 9. Trace Summary

### trace_summary

每个 event/telescope 一行，供质量检查和 debug。

字段：

```text
event_id                    long
telescope_id                int
input_bunches               unsigned long
input_photons               double
blocked_by_obstruction      unsigned long
blocked_incoming            unsigned long
blocked_reflected           unsigned long
hit_mirror                  unsigned long
hit_output_plane            unsigned long
hit_camera                  unsigned long
accepted_camera             unsigned long
lost_between_pixels         unsigned long
unique_hit_pixels           int
signal_pe                   double
time_mean_ns                double
time_rms_ns                 double
```

这部分不是 pylast 重建必需字段，但强烈建议保存，方便验证 LACT 输出是否合理。

## 10. pylast 映射

新增 reader：

```cpp
class LactEventSource : public EventSource
```

它与 `SimtelEventSource` 平级，不修改 `SimtelEventSource` 去兼容 LACT 输出。

初始化时读取：

```text
config            -> simulation_config / metaparam
telescopes        -> subarray telescope positions
camera_pixels     -> CameraGeometry
optics            -> OpticsDescription
corsika_events    -> event truth lookup
```

遍历事件时读取：

```text
corsika_events[event_id]       -> ArrayEvent.simulation.shower
observations[event_id, tel]    -> ArrayEvent.r1
waveforms[event_id, tel]       -> ArrayEvent.r1 when profile=timeseries_pe
observations[event_id, tel]    -> ArrayEvent.dl0 image/peak_time
```

当前 `external/pylast/root/LactEventSource.cpp` 第一版已实现：

```text
config.schema_name/schema_version/profile/run_id 校验
camera_pixels -> dense pixel axis and CameraGeometry
telescopes/optics -> SubarrayDescription
corsika_events -> ArrayEvent.simulation.shower
observations.image_pe -> ArrayEvent.dl0 image and R1 single-sample fallback
observations.image_time_peak_ns -> ArrayEvent.dl0 peak_time
waveforms sparse COO -> ArrayEvent.r1 waveform matrix when present
triggered -> ArrayEvent.simulation.triggered_tels
```

暂未填充 `Pointing`、`EventMonitor` 和真实 `R0Event`。这是有意保留：LACT 输出目前是
p.e. 域结果，不是 ADC/R0 数据。

### image_pe 到 pylast 的映射

LACT integrated p.e. image 不是 raw ADC，不应映射成 `R0Event`。

建议映射：

```text
ArrayEvent.r1.tels[telescope_id].n_pixels  = camera pixel count
ArrayEvent.r1.tels[telescope_id].n_samples = 1
ArrayEvent.r1.tels[telescope_id].waveform  = Matrix<double>(n_pixels, 1)
ArrayEvent.r1.tels[telescope_id].waveform(pixel, 0) = image_pe[pixel]
ArrayEvent.r1.tels[telescope_id].gain_selection = VectorXi::Zero(n_pixels)
```

这样 pylast 后续如果已有 waveform extractor，可以把 `n_samples=1` 视为已经积分后的
single-sample p.e. image。

### timeseries_pe 到 pylast 的映射

建议映射：

```text
ArrayEvent.r1.tels[telescope_id].n_pixels  = camera pixel count
ArrayEvent.r1.tels[telescope_id].n_samples = n_time_bins
ArrayEvent.r1.tels[telescope_id].waveform  = Matrix<double>(n_pixels, n_time_bins)
waveform(pixel, time_bin) = pe
gain_selection = VectorXi::Zero(n_pixels)
```

如果下游算法只需要 image，可以对 `R1Camera.waveform` 沿 sample 轴积分得到
`image_pe`。reader 也可以提供一个 helper，但不应在文件中省略 `observations.image_pe`。

### monitor/calibration

LACT 输出已经是 p.e.，默认不需要 pedestal subtraction 或 ADC-to-p.e. calibration。

第一版建议：

```text
不写 R0Event
不强制写 EventMonitor
如 pylast 某些流程要求 monitor，则 LactEventSource 创建 identity monitor：
  pedestal_per_sample = 0
  dc_to_pe = 1
  calibration_mode = "already_pe"
```

## 11. 必填字段

第一版 `image_pe` profile 必填：

```text
config:
  schema_name, schema_version, profile, run_id, coordinate_convention, event_id_mode

telescopes:
  telescope_id, position, pointing, camera_id, optics_id

camera_pixels:
  camera_id, pixel_id, x_m, y_m, size_m, shape_code

optics:
  optics_id, equivalent_focal_length_m or effective_focal_length_m

corsika_events:
  event_id, shower_event_id, array_id, primary_type, energy_gev,
  theta_deg, phi_deg, core_x_north_m, core_y_west_m

observations:
  event_id, telescope_id, n_pixels_camera, pixel_id, image_pe, total_pe
```

第一版 `timeseries_pe` profile 额外必填：

```text
waveform_enabled
time_reference
time_bin_width_ns
time_window_start_ns
time_window_end_ns
n_time_bins
waveforms.event_id
waveforms.telescope_id
waveforms.pixel_id
waveforms.time_bin
waveforms.pe
observations.image_time_peak_ns
```

## 12. 实现状态

已完成：

```text
1. LACT_sim lact_event_root writer：支持 image_pe / timeseries_pe / debug_full profile。
2. LACT_sim writer：写 config、camera_pixels、optics、telescopes、corsika_events、
   observations、waveform_config、waveforms、trace_summary。
3. EventIO metadata：合并 EVTH 和 simtel MCShower，补 h_first_int/xmax/hmax/depth_start。
4. pylast LactEventSource：支持 image_pe -> DL0/R1 n_samples=1 fallback。
5. pylast LactEventSource：支持 timeseries_pe sparse COO -> R1 waveform matrix。
```

下一步验收标准：

```text
有 ROOT >= 6.24 的环境可以编译 LACT_sim ROOT writer 和 pylast lact_event_source
LACT_sim 产出 image_pe/timeseries_pe ROOT 文件
pylast LactEventSource 可以遍历 ArrayEvent 并支持 __getitem__
每个 ArrayEvent 有 simulation、subarray、dl0、r1
R1 integrated image 与 LACT_sim HDF5/CSV 中的 image_pe 数值一致
至少一个现有 pylast image/Hillas 流程可以直接消费该事件
```
