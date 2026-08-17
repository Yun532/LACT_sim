# 相机时间信息与 waveform proxy 输出说明

这部分用于在 CORSIKA/EventIO 相机输出中保留时间诊断信息。默认全部关闭，
不会增加正式大批量模拟的文件体积。

## 逐像素时间统计

配置：

```ini
output.write_pixel_time_stats=true
```

输出：

```text
/images/dense/time_mean_ns
/images/dense/time_rms_ns
```

这两个数据集的形状和 `/images/dense/pe` 一样，都是：

```text
[image_index, pixel]
```

时间均值和 RMS 按每个光子的最终 signal weight 加权。当前 signal weight
包含光学效率、光收集器接受情况、SiPM PDE 等已经接入的 integrated p.e.
权重。没有光子的像素写为 0。

## waveform proxy

当前程序还没有真实电子学 waveform。这里的 waveform 是一个 proxy，也就是把
到达相机或光收集器出口后的光子按时间 bin 累计成相机图像序列。

配置：

```ini
waveform.enabled=true
waveform.source=pe             # 可选 pe 或 photon_count
waveform.time_reference=image_first # 可选 absolute、image_mean 或 image_first
waveform.time_bin_width_ns=1
waveform.time_window_start_ns=-5
waveform.time_window_end_ns=20
```

输出：

```text
/waveforms/time_edges_ns
/waveforms/time_centers_ns
/waveforms/reference_time_ns
/waveforms/pixel_id_axis
/waveforms/pe                 # 当 waveform.source=pe
/waveforms/photon_count       # 当 waveform.source=photon_count
```

`/waveforms/pe` 和 `/waveforms/photon_count` 的形状是：

```text
[image_index, time_bin, pixel]
```

注意：waveform proxy 不做真实 SiPM/electronics 响应。如果同时开启 NSB 和
`waveform.source=pe`，当前 constant-rate NSB 会在每个 time bin 里独立做
Poisson 采样，并写入 `/waveforms/nsb_pe`；最终 `/images/dense/primary_nsb_pe`
等于这些 time bins 的积分。Cherenkov 光子如果落在配置的 waveform 时间窗
之外，不会进入 `/waveforms/cherenkov_pe`，因此需要精确检查 Cherenkov 时间
积分时，要把时间窗设得足够宽。

对 CORSIKA 相机 GIF，推荐：

```ini
waveform.time_reference=image_first
waveform.time_window_start_ns=-5
waveform.time_window_end_ns=20
```

这种模式在保存 HDF5 时先用每个 event/telescope 图像自己的第一个
Cherenkov 光子到达时间 `/images/index.time_first_ns` 作为 T0，再把光子填入
waveform bin。`/waveforms/reference_time_ns` 会逐 image 记录被减掉的 T0。
这样 HDF5 不用把绝对时间窗开得很宽，也不会因为某台望远镜的 CORSIKA 原始时间是负数而
把主峰裁掉。

如果设置：

```ini
waveform.source=electronics
```

程序会明确报错。这是给后续真实电子学 waveform 接口预留的位置。

## 光收集器逐光子 debug 输出

默认关闭，只建议小样本调试使用：

```ini
collector.debug_photon_output=true
collector.debug_photon_csv=run_logs/my_run/collector_debug_photons.csv
collector.debug_max_photons=20000
```

CSV 会保存通过相机响应模块后的逐光子诊断信息，包括：

```text
event_id,telescope_id,pixel_id,accepted,collector_reflections,
time_ns,collector_time_delay_ns,wavelength_nm,photon_weight,
relative_efficiency,signal_weight,exit_x_m,exit_y_m,exit_z_m,dir_u,dir_v,dir_w
```

目前光收集器内部没有额外传播时间模型，所以
`collector_time_delay_ns=0`。后续如果光收集器几何路径长度需要进入时间，
可以在这里扩展。

## 画图

逐像素时间图：

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/manual_checks/corsika_waveform_smoke/camera_waveform_dense.h5 \
  --event-id 46889802 \
  --telescope-id 3 \
  --quantity time_mean_ns \
  --output run_logs/manual_checks/corsika_waveform_smoke/event46889802_tel4_time_mean.png
```

proxy waveform 时间帧：

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_waveform_gif.py \
  run_logs/manual_checks/corsika_waveform_smoke/camera_waveform_dense.h5 \
  --event-id 46889802 \
  --telescope-id 3 \
  --quantity pe \
  --output-dir run_logs/manual_checks/corsika_waveform_smoke/waveform_frames \
  --stride 10
```

如果本地安装了 Pillow，可以加：

```bash
--gif run_logs/manual_checks/corsika_waveform_smoke/event46889802_tel4_pe.gif
```

## 快速测试配置

仓库提供了一个小样本配置：

```text
configs/experiments/corsika_waveform_smoke_camera.cfg
```

运行方式：

```bash
build_hessio_check/run_corsika_trace \
  configs/experiments/corsika_waveform_smoke_camera.cfg \
  /path/to/corsika.zst
```

正式服务器上通常使用：

```bash
./build/run_corsika_trace configs/experiments/corsika_waveform_smoke_camera.cfg /path/to/corsika.zst
```
