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
```

仓库中的 event 1909、19 号望远镜示例也只使用这六列：

```bash
# 纯光学：输出白板命中和逐像素光子数，可画两张焦平面诊断图
build/run_optical_sim configs/examples/photon_csv_minimal_optics.cfg
python3 python/plot_minimal_photon_csv_outputs.py \
  --mode optics \
  --hits run_logs/examples/photon_csv_minimal/whiteboard_hits.csv \
  --photon-pixels run_logs/examples/photon_csv_minimal/camera_photon_counts.csv \
  --camera configs/cameras/new_camera_pixels.csv \
  --output-dir run_logs/examples/photon_csv_minimal/plots

# 完整相机效率但不生成波形：只画一张期望 p.e. 相机图
build/run_optical_sim configs/examples/photon_csv_minimal_pe_camera.cfg
python3 python/plot_minimal_photon_csv_outputs.py \
  --mode camera \
  --pe-pixels run_logs/examples/photon_csv_minimal/camera_expected_pe.csv \
  --camera configs/cameras/new_camera_pixels.csv \
  --output-dir run_logs/examples/photon_csv_minimal/plots
```

纯光学的两张图保持 LACT_sim 焦平面方向；只有完整相机图采用 pyLAST
天空视角（`pix_x=-x_m`、`pix_y=-y_m`）。波长和时间都可从 CSV 省略：
纯光学关闭波长效率，完整示例则在 cfg 中统一设置 `400 nm`，且不生成波形。
更完整的说明见[最简 Photon CSV 示例](docs/minimal_photon_csv.md)。

如需在输出中保存实际进入光追的位置和方向：

```ini
output.whiteboard_input_photon=true
```

完整列定义、其他输入坐标和 HTML 三维绘图见 [Photon CSV 格式](docs/photon_csv_format.md)及[用户使用手册](docs/user_guide_zh.md#photon-csv)。

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
camera.config=../cameras/new_camera.cfg
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

NSB、trigger、结构遮挡、阵列图、波形和 ROOT/pylast 输出见[完整 CORSIKA 流程](docs/user_guide_zh.md#corsika-相机和完整响应)。

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

配置默认写入积分 p.e. 图像、CORSIKA truth、望远镜/相机 metadata 和稀疏
`timeseries_pe` waveform；默认处理输入文件中的全部 shower，并且只保存通过
trigger 的望远镜事件。快速画图：

```bash
python3 python/plot_lact_root_output.py \
  run_logs/lactroot_only/lact_events.root \
  --outdir run_logs/lactroot_only/root_quicklook
```

pylast 读取、ROOT tree 说明和 notebook 流程见
[ROOT 输出与 pylast 检查](docs/server_root_output_check_zh.md)。

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

- [用户使用手册：安装、运行、输出和画图](docs/user_guide_zh.md)
- [所有 official tests](docs/official_tests.md)
- [Photon CSV 格式](docs/photon_csv_format.md)
- [坐标系中文说明](docs/coordinate_systems_zh.md)
- [坐标系英文参考](docs/coordinate_systems.md)
- [HDF5 输出格式](docs/hdf5_output_format.md)
- [相机时间和波形](docs/camera_timing_waveform_zh.md)
- [NSB 光谱模型](docs/nsb_spectral_model_zh.md)
- [ROOT/pylast 数据层级](docs/pylast_event_data_levels_zh.md)
