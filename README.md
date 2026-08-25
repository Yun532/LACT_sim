# LACT_sim

[**中文**](README.md) | [English](README_EN.md)

LACT_sim 是使用 LACT 望远镜默认镜片布局的 C++ 光学光线追迹程序，可用于：

- 平行光和有限距离点源光追；
- 自定义 Photon CSV 光子输入；
- CORSIKA/EventIO 光子输入；
- 镜面形变、随机误差和望远镜结构遮挡；
- 相机、light collector、SiPM、NSB 和 trigger 响应；
- 白板 CSV、相机 HDF5 及 LACT ROOT/pylast 输出。

本页只介绍用户最常用的安装、运行、cfg 和画图方法。完整流程见[用户使用手册](docs/user_guide_zh.md)，全部官方测试见[官方测试说明](docs/official_tests.md)。

## 安装和编译

```bash
git clone https://github.com/Yun532/LACT_sim.git
cd LACT_sim

# Ubuntu/Debian
sudo apt-get install cmake g++ make zlib1g-dev libhdf5-dev

# 基本绘图环境
python3 -m pip install numpy pandas matplotlib h5py

make
make test
```

`make` 会构建仓库内的 hessioxxx 和 LACT_sim。主要程序位于：

```text
build/run_optical_sim      平行光、点源和 Photon CSV
build/run_corsika_trace    CORSIKA/EventIO、相机、HDF5 和 ROOT
```

只使用人工光源、不需要 CORSIKA 时也可以执行：

```bash
make no-hessio
```

ROOT/pylast 环境和常见编译问题见[完整安装说明](docs/user_guide_zh.md#安装和编译)。

## 一键运行官方测试

不使用 CORSIKA：

```bash
tools/run_official_tests.sh --no-corsika
```

包含 CORSIKA/EventIO：

```bash
tools/run_official_tests.sh --corsika-file /path/to/input.zst
```

脚本会自动编译、运行 CTest、执行标准光追并画图。结果统一写入：

```text
run_logs/official_tests/
```

## 平行光

运行理想平行光白板测试：

```bash
mkdir -p run_logs/official_tests/perfect_parallel
build/run_optical_sim configs/official_tests/perfect_parallel_whiteboard.cfg
```

主要 cfg 内容：

```ini
telescope.config=telescope_1229_minimal.cfg
mirror.config=../mirrors/mirror_1229_imported.cfg
source.config=../sources/parallel_1M_on_axis.cfg
output.config=../outputs/whiteboard_f8.cfg
output.csv=run_logs/official_tests/perfect_parallel/hits.csv
```

画白板光斑：

```bash
MPLBACKEND=Agg python3 python/plot_spot_histogram.py \
  run_logs/official_tests/perfect_parallel/hits.csv \
  --config configs/official_tests/perfect_parallel_whiteboard.cfg \
  --output run_logs/official_tests/perfect_parallel/spot.png
```

有限距离点源使用：

```bash
build/run_optical_sim configs/official_tests/perfect_point_900m_whiteboard.cfg
```

结构遮挡、仰角形变和三维结构图见[人工光源官方测试](docs/user_guide_zh.md#人工光源官方测试)。

### 理想、支架形变与实测标定光学

原有基准均保留且互不覆盖：

- `configs/official_tests/perfect_parallel_whiteboard.cfg`：理想镜片、无误差；
- `configs/official_tests/deformation_parallel_whiteboard.cfg`：理想镜片，只加入原始随仰角支架形变；
- `configs/examples/lact2_measured_parallel.cfg`：20260622 实测标定光学白板入口；
- `configs/examples/lactroot_only.cfg`：实测标定光学生产默认入口，仅输出 pyLAST 所需 ROOT。

实测入口直接读取
`configs/calibrated/lact2_measured_20260622/mirror_elevation_series_20260622.csv`。该表在每个仰角、
每片镜子的同一行中保存支架形变、固定指向偏差和实测曲率半径。
全镜面统一粗糙度由配套的 `errors.cfg` 施加。
全局随波长实测反射率由
`configs/efficiency/mirror_reflectivity_dm0113_13point_mean.csv` 提供。完整字段和避免重复
施加误差的约定见
[`configs/calibrated/lact2_measured_20260622/README.md`](configs/calibrated/lact2_measured_20260622/README.md)。
逐镜 CSV、统一误差和效率乘数的完整优先级见
[`docs/facet_parameter_precedence_zh.md`](docs/facet_parameter_precedence_zh.md)。

### 强度干涉的完整单光子光学响应

长时间 HBT 模拟不需要把数小时内的每个光子反复追迹，但必须先用当前 `main`
校准单光子响应。下面的运行同时包含 LACT2 实测逐镜参数、70° 仰角插值、3D
结构遮挡、真实相机、集光器、反射率、滤光片和 PDE：

```bash
build/run_optical_sim configs/optics/lact2_measured_full_response_400nm.cfg
python3 tools/derive_full_optical_response.py \
  run_logs/sii_full_optics/lact2_measured_400nm_hits.csv \
  configs/optics/lact2_measured_full_response_400nm.csv \
  --input-photons 1000000 \
  --provenance-json configs/optics/lact2_measured_full_response_400nm.provenance.json \
  --source-config configs/optics/lact2_measured_full_response_400nm.cfg
```

生成后 `Instrument.from_repository()` 会自动优先使用该完整响应；没有生成时才回退到
旧的轴上理想误差时间核。光学响应决定探测概率、像素接收、PSF 和到达时间分布，
但 HBT 相关性仍必须在两台望远镜的入射热光统计中产生。线性被动光学不会凭空生成
两镜相关，也不会移动 UV 坐标；它改变 UV 点的误差和权重。

## Photon CSV

运行仓库示例：

```bash
build/run_optical_sim configs/examples/photon_csv_local_whiteboard.cfg
```

用户通常只需要修改：

```ini
source.mode=PhotonCsv
source.csv_path=/path/to/photons.csv
source.coordinate_frame=telescope_local
output.mode=hits
output.csv=run_logs/my_photons/hits.csv
```

CSV 最少包含：

```csv
x_m,y_m,z_m,dir_x,dir_y,dir_z
3.9014025878906251,0.16619310379028321,0,-0.33789920806884766,-0.014721360988914967,-0.94106716376520094
```

可直接参考
[`configs/sources/photon_csv_six_column_example.csv`](configs/sources/photon_csv_six_column_example.csv)。
用于作图的 event 1909、19 号望远镜输入是
[`configs/sources/event1909_tel19_minimal_photons.csv`](configs/sources/event1909_tel19_minimal_photons.csv)，
同样只有六列。
其中 `z_m=0` 保留 CORSIKA 2D bunch 的原始参考平面；两个示例 cfg 通过
`source.eventio_2d=true` 声明其来源，程序再按直接读取 EventIO 的相同规则，
自动使用望远镜本地 `-16 m` 输入面并进行有符号光线求交。

```bash
# 纯光学：输出白板命中和逐像素光子数，可画两张焦平面诊断图
build/run_optical_sim configs/examples/photon_csv_minimal_optics.cfg
python3 python/plot_minimal_photon_csv_outputs.py \
  --mode optics \
  --hits run_logs/examples/photon_csv_minimal/whiteboard_hits.csv \
  --photon-pixels run_logs/examples/photon_csv_minimal/camera_photon_counts.csv \
  --camera configs/cameras/new_camera_pixels_1664.csv \
  --output-dir run_logs/examples/photon_csv_minimal/plots

# 完整相机链：输出 ROOT，再由 pyLAST 画一张 p.e. 相机图
build/run_corsika_trace configs/examples/photon_csv_full_camera_root.cfg
python3 python/plot_photon_csv_root_pylast.py \
  run_logs/examples/photon_csv_full_camera/lact_events.root \
  --event-id 1909 \
  --output run_logs/examples/photon_csv_full_camera/camera_pe.png
```

两个 Photon CSV 示例默认使用 LACT2 20260622 实测镜面参数和补齐角部的 1664 像素相机；完整示例与 EventIO 输入共用 `run_corsika_trace` 的相机处理和 ROOT 输出。
CSV 未提供 `multiplicity` 时使用 cfg 的 `source.multiplicity`，本例设为 1，
即每行一个光子；示例还在 cfg 中统一设置 `400 nm`，并关闭 NSB、trigger
和 waveform。纯光学的两张图保持 LACT_sim 物理焦平面 `(u, v)`。
`LactEventSource` 在输入边界将它映射为 pyLAST 的 `camera_x=-v`、`camera_y=-u`。
ROOT 绘图脚本直接调用 pyLAST 原生 `EventVisualizer.plot_event()`，不再自行交换、
翻转或旋转相机坐标；白板和 HDF5 诊断图直接显示文件中保存的物理 `(u, v)`。
更完整的说明见[最简 Photon CSV 示例](docs/minimal_photon_csv.md)。

如需在输出中保存实际进入光追的位置和方向：

```ini
output.whiteboard_input_photon=true
```

完整列定义、其他输入坐标和 HTML 三维绘图见 [Photon CSV 格式](docs/photon_csv_format.md)及[用户使用手册](docs/user_guide_zh.md#photon-csv)。

## 电子学的两个最小测试

下面两个测试都使用仓库自带的 event 1909、19 号望远镜 Photon CSV，关闭 NSB，
因此不需要额外的 CORSIKA 文件。它们使用相同的光学、light collector、PDE 和随机种子，
区别只在于是否继续生成单 PE 波形，适合检查电子学开关没有改变上游结果。

### 1. 只考虑显式微单元饱和，不生成波形

```bash
build/run_corsika_trace \
  configs/examples/photon_csv_saturation_on_waveform_off_validation.cfg
```

输出位于：

```text
validation/measured_electronics/saturation_on_waveform_off/lact_events.root
validation/measured_electronics/saturation_on_waveform_off/lact_events.h5
```

标准结果是饱和前纯切伦科夫 `3198 PE`、饱和后 `3196 PE`，其中 `2 PE` 因同一微单元
已经触发而被拒绝；`waveform.enabled=false`，ROOT 中直接保存饱和后的积分 PE 图像。

### 2. 显式微单元饱和加实测单 PE 波形

```bash
build/run_corsika_trace \
  configs/examples/photon_csv_measured_spe_4ns_validation.cfg
```

输出位于：

```text
validation/measured_electronics/measured_waveform_no_nsb/lact_events.root
validation/measured_electronics/measured_waveform_no_nsb/lact_events.h5
validation/measured_electronics/measured_waveform_no_nsb/waveforms.csv
```

该测试保持相同的 `3198 -> 3196 PE` 微单元结果，然后为每次 fired PE 抽取实测相对
电荷并叠加实测单 PE 模板，最后生成 `4 ns`、单位为 `mV` 的 1664 像素波形。
单 PE 面积定标常数为 `84.03495572475859 mV ns / PE`。可进一步运行三种格式一致性检查：

```bash
python3 scripts/validate_measured_electronics.py \
  validation/measured_electronics \
  --json validation/measured_electronics/validation_report.json
```

波形 ROOT 可由 pyLAST 的 `FullWaveFormExtractor` 或带尾部补偿的
`LocalPeakExtractor` 读取；完整数据层级和 notebook 用法见
[实测单 PE 电子学说明](docs/measured_spe_electronics_zh.md)。

正式 CORSIKA 批量模拟可使用实测波形配置：

```bash
build/run_corsika_trace \
  configs/examples/corsika_lact_pylast_root_only_measured_waveform.cfg \
  /path/to/input.zst \
  -C output.lact_root_path=/path/to/lact_events.root
```

The production NSB configuration uses `nsb.window_ns=32`. This is the NSB
integration gate for the saved p.e. image, relative to the event reference
time: only NSB primary/fired p.e. in `[0, 32 ns)` enter `image_pe` and
`total_pe`. The measured-electronics waveform and voltage trigger still use
the complete configured waveform, including padded NSB outside this image
gate, so waveform edges remain stationary without inflating the saved image.

该配置开启暗夜 NSB、微单元饱和和实测单 PE 波形，只写 ROOT，并只保存通过相机触发的望远镜。触发条件是 20 ns 内至少 3 个像素的 4 ns 波形采样值达到 `19.267713 mV`；这个电压是实测模板经 4 ns 采样后的平均 8 PE 等效阈值。阵列触发保持关闭。

只保留微单元饱和、关闭 NSB 和波形的 CORSIKA 基线配置为：

```bash
build/run_corsika_trace \
  configs/examples/corsika_lact_pylast_root_only_microcell_no_nsb.cfg \
  /path/to/input.zst \
  -C output.lact_root_path=/path/to/lact_events.root
```

该配置使用 20 ns 内 `10 PE / 3 像素` 的计数触发，并保存所有最终积分图像非空的望远镜（包括低于触发阈值的望远镜；全零相机不保存）；ROOT 的 `observations.image_pe` 是无 NSB、经过显式微单元饱和后的稀疏积分图像，并由 pyLAST 自动补零后映射到 DL0。每台望远镜的触发结果保存在 `observations.triggered`。

正式配置默认使用紧凑 ROOT 输出（`output.lact_root_write_components=false`）。其中始终保留
`image_pe`、`image_time_peak_ns` 和纯切伦科夫真值
`image_primary_cherenkov_pe`；pyLAST 将后者映射为
`event.simulation.tels[tel_id].true_image`。开启波形时还会保留
`waveform_config` 与 `waveforms`。NSB/dark/fired 分量、gap/饱和损失、逐像素平均/均方根
时间和 `trace_summary` 都属于测试诊断量，只在
`output.lact_root_write_components=true` 时写入。逐 PE 和逐微单元树还需要同时开启各自的
`electronics.output.save_*` 开关；因此正式批量输出不会被测试记录撑大。

## CORSIKA 白板测试

白板测试适合在接入相机前检查 CORSIKA 光子是否正确读入：

```bash
mkdir -p run_logs/official_tests/corsika
build/run_corsika_trace \
  configs/official_tests/corsika_whiteboard.cfg \
  /path/to/input.zst
```

主要 cfg 内容：

```ini
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70
source.mode=EventIO
source.coordinate_frame=corsika_nwu_relative
output.format=csv
output.hits_csv=run_logs/official_tests/corsika/whiteboard_hits.csv
```

事件选择和白板绘图命令见 [CORSIKA 白板流程](docs/user_guide_zh.md#corsika-白板测试)。

## CORSIKA 相机和完整响应

理想相机 HDF5：

```bash
build/run_corsika_trace \
  configs/official_tests/corsika_new_camera.cfg \
  /path/to/input.zst
```

完整光学响应：

```bash
build/run_corsika_trace \
  configs/official_tests/corsika_full_response_camera.cfg \
  /path/to/input.zst
```

主要 cfg 内容：

```ini
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70
camera.config=../cameras/new_camera_1664.cfg
electronics.config=../electronics/explicit_microcell_saturation_only.cfg
waveform.enabled=false
source.mode=EventIO
source.coordinate_frame=corsika_nwu_relative
output.format=hdf5
output.hdf5_storage=dense
```

快速画第一幅相机图：

```bash
MPLBACKEND=Agg python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --image-index 0 \
  --quantity pe \
  --output run_logs/official_tests/corsika/camera_first.png
```

NSB、PE 计数 trigger、结构遮挡、阵列图和 ROOT/pylast 输出见[完整 CORSIKA 流程](docs/user_guide_zh.md#corsika-相机和完整响应)。
当前正式配置只做显式微单元饱和，不生成单 PE 或采样波形；波形配置仅保留在专用验证示例中。

## LACT ROOT / pylast

只保存 LACT ROOT：

```bash
build/run_corsika_trace \
  configs/examples/lactroot_only.cfg \
  /path/to/input.zst
```

默认输出：

```text
run_logs/lactroot_only/lact_events.root
```

配置默认写入饱和后的积分 p.e. 图像、CORSIKA truth 和望远镜/相机 metadata，
不写 waveform；默认处理输入文件中的全部 shower，并且只保存通过 PE 计数
trigger 的望远镜事件。快速画图：

```bash
python3 python/plot_lact_root_output.py \
  run_logs/lactroot_only/lact_events.root \
  --outdir run_logs/lactroot_only/root_quicklook
```

pylast 读取、ROOT tree 说明和 notebook 流程见
[ROOT 输出与 pylast 检查](docs/server_root_output_check_zh.md)。

完整的阵列/core、触发时延、相机图像、Hillas、方向/芯位重建和 3D SDP
流程见 [pyLAST Jupyter notebook](notebooks/lact_root_to_pylast_visualize.ipynb)。
同一完整结构的两种 NSB 示例同时保存在 LACT_sim 和 pyLAST：

- [无 NSB 单事件重建](notebooks/lact_event_reconstruction_no_nsb.ipynb)
- [加入 pyLAST 泊松 NSB 的单事件重建](notebooks/lact_event_reconstruction_with_nsb.ipynb)

两本 notebook 只在 Cell 1 的 NSB 参数上不同，其他读取、绘图、清理、
Hillas 和方向/芯位重建步骤保持一致。
LACT ROOT adapter 和这些可视化接口位于
[pyLAST `lact_sim` 分支](https://github.com/Yun532/pylast/tree/lact_sim)。

## 常用 cfg

| 用途 | cfg |
|---|---|
| 理想平行光 | `configs/official_tests/perfect_parallel_whiteboard.cfg` |
| 900 m 点源 | `configs/official_tests/perfect_point_900m_whiteboard.cfg` |
| 平行光加结构遮挡 | `configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg` |
| 平行光仅支架形变 | `configs/official_tests/deformation_parallel_whiteboard.cfg` |
| LACT2 实测光学平行光 | `configs/examples/lact2_measured_parallel.cfg` |
| Photon CSV | `configs/examples/photon_csv_local_whiteboard.cfg` |
| CORSIKA 白板 | `configs/official_tests/corsika_whiteboard.cfg` |
| CORSIKA 理想相机 | `configs/official_tests/corsika_new_camera.cfg` |
| CORSIKA + NSB + trigger | `configs/official_tests/corsika_nsb_trigger_camera.cfg` |
| CORSIKA 完整响应 | `configs/official_tests/corsika_full_response_camera.cfg` |
| LACT2 实测光学 ROOT-only | `configs/examples/lactroot_only.cfg` |

## 文档

- [文档总目录：按使用场景分类](docs/README.md)
- [用户使用手册：安装、运行、输出和画图](docs/user_guide_zh.md)
- [所有 official tests](docs/official_tests.md)
- [Photon CSV 格式](docs/photon_csv_format.md)
- [真实 LACT 三维坐标模型](docs/assets/lact-coordinate-system-3d.html)
- [交互式坐标系总图](docs/assets/coordinate-system-explorer.html)
- [坐标系中文说明](docs/coordinate_systems_zh.md)
- [坐标系英文参考](docs/coordinate_systems.md)
- [HDF5 输出格式](docs/hdf5_output_format.md)
- [电子学结果与 ROOT、HDF5、CSV 统一输出](docs/electronics_output_unification_zh.md)
- [相机时间和波形](docs/camera_timing_waveform_zh.md)
- [NSB 光谱模型](docs/nsb_spectral_model_zh.md)
- [ROOT/pylast 数据层级](docs/pylast_event_data_levels_zh.md)
