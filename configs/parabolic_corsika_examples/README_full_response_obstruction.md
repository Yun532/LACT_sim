# 抛物面 CORSIKA full-response + 遮挡测试

这个测试使用同一个 6 m 圆形抛物面镜，但打开更完整的响应链：遮挡模型、光收集器、镜面反射率、滤光片透过率、SiPM PDE、小误差项和简单 trigger。NSB 默认关闭，因此图像和时间分布只对应 Cherenkov 光经过光学系统后的 p.e.。

## 配置文件

`configs/parabolic_corsika_examples/full_response_obstruction_parabolic_camera.cfg`

关键配置：

- `mirror.config=../mirrors/mirror_6m_parabolic.cfg`：使用 6 m 圆形抛物面镜。
- `obstruction.config=../obstructions/raytrace_final_structure.cfg`：加入相机到镜面之间的 3D 遮挡结构。
- `sipm.config=../sipm/new_camera_sipm.cfg`：使用真实相机 SiPM 配置，PDE 在这里设置。
- `efficiency.config=../efficiency/curves_all.cfg`：使用镜面反射率和滤光片透过率曲线。
- `trigger.config=../trigger/example_simple_multiplicity.cfg`：计算简单 multiplicity trigger。
- `output.save_only_triggered=false`：示例保留所有望远镜图像，方便检查；如只想保存触发图像可改为 `true`。
- `waveform.source=pe`：时间序列保存经过效率/PDE 后的 p.e.。
- `waveform.time_reference=image_first`：每个 event/telescope 图像以第一个到达相机的 Cherenkov 光子作为 T0。
- `waveform.time_bin_width_ns=1`：保存 1 ns 一个时间 bin。
- `waveform.time_window_start_ns=-5` 和 `waveform.time_window_end_ns=20`：保存 `T0-5 ns` 到 `T0+20 ns`。
- `source.max_shower_events=1`：示例默认只跑第一个 shower event；批量运行时可以删掉或改大。

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
