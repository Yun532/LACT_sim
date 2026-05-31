# 抛物面 CORSIKA 理想相机测试

这个测试使用 6 m 圆形抛物面镜，读取 CORSIKA/EventIO 光子，输出真实相机像素图像。它用于检查纯几何光追和相机像素收集；不考虑 PDE、镜面反射率曲线、滤光片透过率、NSB 或 trigger，但启用 MODTRAN 大气透过率，和 full-response 示例保持同一套大气处理。

## 配置文件

`configs/parabolic_corsika_examples/ideal_parabolic_camera.cfg`

关键配置：

- `mirror.config=../mirrors/mirror_6m_parabolic.cfg`：使用 6 m 圆形抛物面镜。
- `camera.config=../cameras/new_camera.cfg`：使用当前真实相机像素布局和光收集器。
- `sipm.config=../sipm/ideal_sipm.cfg`：SiPM/PDE 使用理想配置。
- `efficiency.config=../efficiency/ideal.cfg`：不引入波长效率曲线。
- `atmosphere.config=../atmosphere/modtran_4400_desert.cfg`：使用 MODTRAN 总光学厚度大气吸收。
- `trigger.config=../trigger/disabled.cfg`：不做 trigger 筛选。
- `waveform.source=photon_count`：时间序列保存相机收到的光子数，不使用 p.e. 权重。
- `waveform.time_reference=image_first`：每个 event/telescope 图像以第一个到达相机的 Cherenkov 光子作为 T0。
- `waveform.time_bin_width_ns=1`：保存 1 ns 一个时间 bin。
- `waveform.time_window_start_ns=-5` 和 `waveform.time_window_end_ns=20`：保存 `T0-5 ns` 到 `T0+20 ns`。
- `source.eventio_2d_input_plane_z_m=0`：若输入 CORSIKA/EventIO 是 2D photon bunch，则把光子 `x/y` 放在望远镜本地默认 `z=0 m` 输入平面。
- `source.eventio_2d_plane_mode=auto`：程序按输入平面位置自动选择追迹方式；这里不再使用之前调试用的 `z=-16 m` 平面。

默认会处理 CORSIKA 文件中的全部 shower event 和全部 array/core offset。若只想快速测试前几个 shower，可在 cfg 里临时加入：

```ini
source.max_shower_events=1
```

## 单独运行

```bash
make
BUILD_DIR=build tools/run_parabolic_corsika_examples.sh --corsika-file /path/to/input.zst --only ideal
```

主要输出：

- `run_logs/parabolic_corsika_examples/ideal/camera_ideal_dense.h5`
- `run_logs/parabolic_corsika_examples/ideal/plots/all_tel_photon_count/`
- `run_logs/parabolic_corsika_examples/ideal/plots/brightest_tel_photon_count.gif`
- `run_logs/parabolic_corsika_examples/ideal/plots/time_hist_photon_count.png`

一键脚本会对所有望远镜画静态图；GIF 默认只画该事件里信号最大的那台望远镜，避免 32 台望远镜生成几千张 1 ns frame。

## 手动画图

```bash
python3 python/plot_hdf5_camera.py run_logs/parabolic_corsika_examples/ideal/camera_ideal_dense.h5 \
  --event-id EVENT_ID --quantity photon_count \
  --output run_logs/parabolic_corsika_examples/ideal/plots/all_tel_photon_count

python3 python/plot_hdf5_waveform_gif.py run_logs/parabolic_corsika_examples/ideal/camera_ideal_dense.h5 \
  --event-id EVENT_ID --telescope-id TELESCOPE_ID --quantity photon_count \
  --output-dir run_logs/parabolic_corsika_examples/ideal/plots/waveform_photon_count_frames \
  --gif run_logs/parabolic_corsika_examples/ideal/plots/brightest_tel_photon_count.gif \
  --stride 1

python3 python/plot_hdf5_time_histogram.py run_logs/parabolic_corsika_examples/ideal/camera_ideal_dense.h5 \
  --event-id EVENT_ID --quantity photon_count \
  --output run_logs/parabolic_corsika_examples/ideal/plots/time_hist_photon_count.png
```
