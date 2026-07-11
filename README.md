# LACT_sim

[**中文**](README.md) | [English](README_EN.md)

LACT_sim 是面向 LACT 镜望远镜的 C++ 光学光线追迹程序，可用于：

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

## 常用 cfg

| 用途 | cfg |
|---|---|
| 理想平行光 | `configs/official_tests/perfect_parallel_whiteboard.cfg` |
| 900 m 点源 | `configs/official_tests/perfect_point_900m_whiteboard.cfg` |
| 平行光加结构遮挡 | `configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg` |
| Photon CSV | `configs/examples/photon_csv_local_whiteboard.cfg` |
| CORSIKA 白板 | `configs/official_tests/corsika_whiteboard.cfg` |
| CORSIKA 理想相机 | `configs/official_tests/corsika_new_camera.cfg` |
| CORSIKA + NSB + trigger | `configs/official_tests/corsika_nsb_trigger_camera.cfg` |
| CORSIKA 完整响应 | `configs/official_tests/corsika_full_response_camera.cfg` |
| pylast ROOT-only | `configs/examples/corsika_lact_pylast_root_only_full_response.cfg` |

## 文档

- [用户使用手册：安装、运行、输出和画图](docs/user_guide_zh.md)
- [所有 official tests](docs/official_tests.md)
- [Photon CSV 格式](docs/photon_csv_format.md)
- [坐标系说明](docs/coordinate_systems.md)
- [HDF5 输出格式](docs/hdf5_output_format.md)
- [相机时间和波形](docs/camera_timing_waveform_zh.md)
- [NSB 光谱模型](docs/nsb_spectral_model_zh.md)
- [ROOT/pylast 数据层级](docs/pylast_event_data_levels_zh.md)
