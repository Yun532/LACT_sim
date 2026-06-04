# 抛物面 CORSIKA full-response + 遮挡测试

这个测试使用同一个 6 m 圆形抛物面镜，但打开更完整的响应链：遮挡模型、光收集器、镜面反射率、滤光片透过率、SiPM PDE、SkyCalc spectral NSB、小误差项和简单 trigger。图像和时间序列中的 `pe` 包含 Cherenkov p.e. 与 NSB p.e. 的总和；若需要分量检查，可看 HDF5 中的 `cherenkov_pe` 和 `nsb_pe`。

## 配置文件

`configs/parabolic_corsika_examples/full_response_obstruction_parabolic_camera.cfg`

关键配置：

- `mirror.config=../mirrors/mirror_6m_parabolic.cfg`：使用 6 m 圆形抛物面镜。
- `obstruction.config=../obstructions/raytrace_final_structure.cfg`：加入相机到镜面之间的 3D 遮挡结构。
- `sipm.config=../sipm/new_camera_sipm.cfg`：使用真实相机 SiPM 配置，PDE 在这里设置。
- `efficiency.config=../efficiency/curves_all.cfg`：使用镜面反射率和滤光片透过率曲线。
- `atmosphere.config=../atmosphere/modtran_4400_desert.cfg`：使用 MODTRAN 总光学厚度大气吸收。
- `nsb.config=../nsb/spectral_skycalc_dark_with_obstruction.cfg`：使用无月 SkyCalc LoNS 光谱 NSB，固定有效面积为 `22.606448 m^2`，程序算得约 `0.07438 p.e./ns/pixel`。
- `trigger.config=../trigger/example_simple_multiplicity.cfg`：计算简单 multiplicity trigger。
- `output.save_only_triggered=false`：示例保留所有望远镜图像，方便检查；如只想保存触发图像可改为 `true`。
- `waveform.source=pe`：时间序列保存经过效率/PDE 后的 p.e.。
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
BUILD_DIR=build tools/run_parabolic_corsika_examples.sh --corsika-file /path/to/input.zst --only full
```

主要输出：

- `run_logs/parabolic_corsika_examples/full_response_obstruction/camera_full_response_obstruction_dense.h5`
- `run_logs/parabolic_corsika_examples/full_response_obstruction/plots/all_tel_pe/`
- `run_logs/parabolic_corsika_examples/full_response_obstruction/plots/brightest_tel_pe.gif`
- `run_logs/parabolic_corsika_examples/full_response_obstruction/plots/time_hist_pe.png`

一键脚本会对所有望远镜画静态图；GIF 默认只画该事件里信号最大的那台望远镜，避免 32 台望远镜生成几千张 1 ns frame。

## 手动画图

```bash
python3 python/plot_hdf5_camera.py run_logs/parabolic_corsika_examples/full_response_obstruction/camera_full_response_obstruction_dense.h5 \
  --event-id EVENT_ID --quantity pe \
  --output run_logs/parabolic_corsika_examples/full_response_obstruction/plots/all_tel_pe

python3 python/plot_hdf5_waveform_gif.py run_logs/parabolic_corsika_examples/full_response_obstruction/camera_full_response_obstruction_dense.h5 \
  --event-id EVENT_ID --telescope-id TELESCOPE_ID --quantity pe \
  --output-dir run_logs/parabolic_corsika_examples/full_response_obstruction/plots/waveform_pe_frames \
  --gif run_logs/parabolic_corsika_examples/full_response_obstruction/plots/brightest_tel_pe.gif \
  --stride 1

python3 python/plot_hdf5_time_histogram.py run_logs/parabolic_corsika_examples/full_response_obstruction/camera_full_response_obstruction_dense.h5 \
  --event-id EVENT_ID --quantity pe \
  --output run_logs/parabolic_corsika_examples/full_response_obstruction/plots/time_hist_pe.png
```
