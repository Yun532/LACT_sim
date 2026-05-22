# 抛物面 CORSIKA 示例测试

这里放两个独立的 6 m 圆形抛物面 CORSIKA 相机测试：

1. `ideal_parabolic_camera.cfg`：理想光学相机测试，画 `photon_count` 图像、1 ns/bin GIF 和光子到达时间 histogram。
2. `full_response_obstruction_parabolic_camera.cfg`：full-response + 遮挡测试，画 `pe` 图像、1 ns/bin GIF 和 p.e. 到达时间 histogram。

两个测试的 waveform 都使用 `waveform.time_reference=image_first`，也就是每个 event/telescope 以第一个到达相机的 Cherenkov 光子作为 T0，并保存 `T0-5 ns` 到 `T0+20 ns`。

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
