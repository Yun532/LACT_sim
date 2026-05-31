# 抛物面 CORSIKA 示例测试

这里放两个独立的 6 m 圆形抛物面 CORSIKA 相机测试：

1. `ideal_parabolic_camera.cfg`：理想光学相机测试，启用 MODTRAN 大气透过率，画 `photon_count` 图像、1 ns/bin GIF 和光子到达时间 histogram。
2. `full_response_obstruction_parabolic_camera.cfg`：full-response + 遮挡测试，画 `pe` 图像、1 ns/bin GIF 和 p.e. 到达时间 histogram。

两个测试的 waveform 都使用 `waveform.time_reference=image_first`，也就是每个 event/telescope 以第一个到达相机的 Cherenkov 光子作为 T0，并保存 `T0-5 ns` 到 `T0+20 ns`。

两个示例都显式设置：

```ini
source.eventio_2d_input_plane_z_m=0
source.eventio_2d_plane_mode=auto
```

这表示如果输入 EventIO photon bunch 是 2D 格式，就把 `x/y` 放在望远镜本地默认 `z=0 m` 输入平面。`auto` 会按这个平面位置选择正常向前追迹，不再使用之前调试用的 `z=-16 m` 平面。

两个示例默认都会处理输入 CORSIKA 文件中的全部 shower event 和全部 array/core offset。若只想快速试跑，可在对应 cfg 中临时加入 `source.max_shower_events=1`。

推荐一键运行：

```bash
make
BUILD_DIR=build tools/run_parabolic_corsika_examples.sh --corsika-file /path/to/input.zst
```

只跑其中一个：

```bash
tools/run_parabolic_corsika_examples.sh --corsika-file /path/to/input.zst --only ideal
tools/run_parabolic_corsika_examples.sh --corsika-file /path/to/input.zst --only full
```

更详细说明见：

- `configs/parabolic_corsika_examples/README_ideal.md`
- `configs/parabolic_corsika_examples/README_full_response_obstruction.md`
