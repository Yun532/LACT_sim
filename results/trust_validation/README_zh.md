# LACT_sim 可信链路检验结果

生成日期：2026-07-25
代码基线：GitHub `main` 的
`0048369cd06e87322d9ef3611a4a369a1152fe99`，加本轮工作区修改。

## 目录

### `muon_full_1000_final`

真实 CORSIKA 7.741 μ 子输入的全 1000-shower 结果：

- `corsika_trace.h5`：全量 shower、array event、稠密相机分量图像和 trigger；
- `lact_events.root`：同一运行的 ROOT 后端输出；
- `validation_summary.json`：可机器读取的形状、统计、sanity 和跨后端比较；
- `stdout.txt`：逐 shower/望远镜的光学计数；
- `stderr.txt`：输入读取阶段的状态信息，无运行错误。

全量 HDF5 为 1,124,582,303 byte，ROOT 为 209,898,581 byte。

### `muon_first_event_canonical_nsb`

同一输入第一 shower 的逐级诊断结果。除 HDF5、ROOT、日志和验证报告外，还包含：

- `trace_summary.csv`：每个 event-telescope 的光学链路汇总；
- `mirror_diagnostic.csv`：镜面求交和反射诊断；
- `collector_photons.csv`：collector 内部反射、光程、延迟和出口状态；
- `camera_pixels.csv`：相机像素表；
- `atmosphere_height_histogram.csv`：大气吸收高度统计。

该 HDF5 保存 31 个 2-ns bin 的稀疏 Cherenkov/NSB 波形，可独立复算 trigger。

## 验收结果

```text
全量 shower header                    1000
全量 array event                     13100
全量 event-telescope                 24706
ROOT/HDF5 trigger 比较               24706 / 24706 通过
第一 shower 保存波形 trigger 重算         12 / 12 通过
全量 sanity                          passed
```

完整解释见项目 `docs/trust_validation_muon_result_2026-07-25_zh.md`。

## SHA-256

```text
e64255e0882e2c2d0b2277c2a5f089caf8fc445f67e98e0ddf7ac1656a0aeeda  muon_full_1000_final/corsika_trace.h5
2f11d1f7fc37719374b3bcf0a8888c3d6b492a575d0be47ede3813f2e4475610  muon_full_1000_final/lact_events.root
0c3d22a8426f95ef094a3fcfec041074dd4445d5618e86bbf2f14e16bbeb5997  muon_full_1000_final/validation_summary.json
3d3b52b11d3ab429810518dd27bd2b773552a761f78ee88ddc399358794531d2  muon_first_event_canonical_nsb/corsika_trace.h5
b1e02a512df2e1db882f917d5f68fbc415bffc79fafaabaecdc4bf516a700ca6  muon_first_event_canonical_nsb/lact_events.root
0b19d6007d45f6c7179afde2774dafb0c598f84c774a514e280a5521c3fab42f  muon_first_event_canonical_nsb/validation_summary.json
```

## 复核

在具有 `h5py`、`numpy` 和 `uproot` 的 Python 环境中，从项目根目录运行：

```bash
python tools/summarize_trust_validation.py \
  results/trust_validation/muon_full_1000_final/corsika_trace.h5 \
  --root results/trust_validation/muon_full_1000_final/lact_events.root \
  --stdout-log results/trust_validation/muon_full_1000_final/stdout.txt \
  --pretty
```

全量配置为了控制文件大小，没有保存逐 bin 波形；需要检查波形时，对
`muon_first_event_canonical_nsb` 目录运行同一个命令。
