# LACT_sim

[**中文**](README.md) | [English](README_EN.md)

LACT_sim 是面向 LACT 1229 面镜望远镜的 C++ 光学光线追迹程序。目前支持：

- 理想光学条件下的平行光和有限距离点源；
- Photon CSV 自定义光子输入；
- CORSIKA/EventIO photon bunch 直接输入；
- 望远镜指向、阵列位置和多种输入坐标系；
- 随仰角变化的镜面结构形变、随机光学误差和结构遮挡；
- 反射率、滤光片、SiPM PDE、light collector、大气透过率和 NSB；
- 白板、相机像素、HDF5 以及面向 pylast 的 LACT ROOT 输出。

程序生成的构建目录、运行日志、CSV、HDF5、ROOT 和图片不纳入源码仓库，均可通过下述命令重新生成。

## 目录结构

```text
apps/        C++ 可执行程序和测试
configs/     望远镜、镜面、光源、相机、输出和官方测试配置
docs/        坐标系、输入输出和物理模型说明
include/     C++ 头文件
python/      绘图、检查和格式转换工具
src/         C++ 实现
tools/       构建及官方测试脚本
external/    随仓库提供的 hessioxxx 源码和元数据
```

## 构建

在仓库根目录执行：

```bash
make
make test
```

`make` 会先构建 `external/hessioxxx`，再使用 CMake 在 `build/` 中编译 LACT_sim。只运行平行光、点源和 Photon CSV，不需要 CORSIKA/EventIO 时，可以执行：

```bash
make no-hessio
ctest --test-dir build --output-on-failure
```

HDF5 相机输出需要系统提供 HDF5 C 开发库。例如：

```bash
# Ubuntu/Debian
sudo apt-get install cmake g++ make zlib1g-dev libhdf5-dev

# CentOS/RHEL 系列
sudo yum install cmake gcc-c++ make zlib-devel hdf5-devel
```

如果 CMake 能找到 ROOT 6.24 或更新版本，会自动启用 LACT ROOT 输出。ROOT 不包含在 `external/` 中，可使用服务器 module、系统安装或 conda 环境：

```bash
source /path/to/root/bin/thisroot.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DLACT_ENABLE_ROOT=ON \
  -DROOT_DIR="$(root-config --cmakedir)"
cmake --build build -j4
```

## 两条运行入口

- `build/run_optical_sim`：平行光、点源和 Photon CSV 光学调试。
- `build/run_corsika_trace`：CORSIKA/EventIO 输入、阵列事件、相机响应、HDF5 和 ROOT 输出。

不要用 `run_optical_sim` 直接读取 EventIO；它会拒绝 `source.mode=EventIO`。

## 平行光和有限距离点源

平行光示例：

```bash
build/run_optical_sim configs/official_tests/perfect_parallel_whiteboard.cfg
```

有限距离点源示例：

```bash
build/run_optical_sim configs/official_tests/perfect_point_900m_whiteboard.cfg
```

光源由 cfg 的 `source.mode` 决定：

```ini
source.mode=Parallel
source.coordinate_frame=telescope_local
```

或：

```ini
source.mode=Point
source.coordinate_frame=telescope_local
```

`beam_theta_deg` 和 `beam_phi_deg` 描述的是光子的传播方向，不是天空中光源方向。反射光学系统还会在焦平面形成倒像，因此判断离轴方向时应先阅读[坐标系说明](docs/coordinate_systems.md)。

## 望远镜位置和指向

```ini
telescope.position_m=0,0,0
telescope.pointing_az_deg=30
telescope.pointing_el_deg=50
```

其中：

```text
elevation = 90° - zenith
```

例如天顶角 40°、方位角 30°应设置：

```ini
telescope.pointing_az_deg=30
telescope.pointing_el_deg=50
```

如果输入位置使用望远镜本地或“相对望远镜”的坐标，不需要在 CSV 中重复加入望远镜全局位置。如果输入使用全局坐标，则必须正确配置 `telescope.position_m`，程序会在进入光追前完成平移和旋转。

## 输入坐标系

内部光追仍使用原有的望远镜本地光学坐标。`source.coordinate_frame` 只说明输入光子的位置和方向应该如何转换到这个内部坐标；它不会改变镜面、相机、输出平面或绘图坐标定义。

| `source.coordinate_frame` | 输入含义 | 位置处理 |
|---|---|---|
| `telescope_local` | 已经位于选定望远镜的本地光学坐标 | 不再做全局平移；推荐作为 Photon CSV 默认值 |
| `corsika_nwu_relative` | CORSIKA NWU 方向，位置已相对望远镜 | 使用 CORSIKA NWU 到望远镜本地的旋转 |
| `corsika_nwu_global` | CORSIKA NWU 阵列绝对坐标 | 先减一次望远镜位置，再旋转 |
| `lact_generic_global` | 原有 LACT 通用全局坐标，方位角从 `+x` 指向 `+y` | 减望远镜位置并使用原有通用框架转换 |

旧配置：

```ini
source.local_telescope_frame=true
```

仍可使用；`true` 等价于 `telescope_local`，`false` 等价于 `lact_generic_global`。新配置应显式使用 `source.coordinate_frame`。

## Photon CSV

最少需要六列：

```csv
x_m,y_m,z_m,dir_x,dir_y,dir_z
```

可选列：

```text
time_ns,wavelength_nm,weight,multiplicity,event_id,telescope_id,emission_altitude_km
```

方向向量不要求预先严格归一化，读取时会进行归一化。`weight * multiplicity` 是进入光学传播前的光子权重；`emission_altitude_km` 可供 MODTRAN 大气模型使用。

推荐从这个示例开始：

```bash
build/run_optical_sim configs/examples/photon_csv_local_whiteboard.cfg
```

对应的核心配置是：

```ini
source.mode=PhotonCsv
source.csv_path=configs/sources/photon_csv_center_test.csv
source.coordinate_frame=telescope_local
```

外部 CSV 的列名或单位不一致时，可使用：

```bash
python3 python/normalize_photon_csv.py external.csv photons.csv \
  --map x_m=x_cm,y_m=y_cm,z_m=z_cm,dir_x=ux,dir_y=uy,dir_z=uz \
  --scale-position 0.01 \
  --normalize-direction \
  --fail-on-nonfinite
```

完整格式见 [Photon CSV 说明](docs/photon_csv_format.md)。

## 保存进入光追时的光子

对白板 hits CSV 开启：

```ini
output.whiteboard_input_photon=true
```

输出会增加：

```text
input_local_x_m,input_local_y_m,input_local_z_m,
input_local_dir_x,input_local_dir_y,input_local_dir_z
```

这些字段直接保存光子完成输入坐标解释后、实际交给光追的望远镜本地位置和方向。输出阶段不会再做一次坐标转换。因此，无论原始输入选择哪一种 `source.coordinate_frame`，这些字段都可以直接用于比较、调试和三维绘图。

在望远镜三维结构图上叠加入射光子：

```bash
python3 python/plot_optical_layout_html.py \
  --config configs/examples/photon_csv_local_whiteboard.cfg \
  --input-photon-csv run_logs/examples/photon_csv_local_whiteboard/hits.csv \
  --output run_logs/examples/photon_csv_local_whiteboard/layout_with_photons.html
```

## CORSIKA/EventIO

CORSIKA 运行应使用：

```ini
source.mode=EventIO
source.coordinate_frame=corsika_nwu_relative
source.use_eventio_telescope_position=true
```

如果 CORSIKA 输入给出天顶角 `Z` 和方位角 `A`：

```ini
telescope.pointing_el_deg=90-Z
telescope.pointing_az_deg=A
```

运行：

```bash
build/run_corsika_trace my_corsika_run.cfg /path/to/input.zst
```

推荐从以下模板复制并修改输入、输出和指向：

```text
configs/templates/minimal_corsika_camera.cfg
configs/examples/corsika_new_user_full.cfg
```

大文件快速检查可设置：

```ini
source.max_shower_events=1
# 或
source.filter_shower_event_id=327666
```

## 大气吸收

CORSIKA 光子到达 EventIO 记录面之后，LACT_sim 可选择再施加望远镜附近或用户指定的额外大气透过率。不开启额外吸收：

```ini
atmosphere.transmission=none
```

简单常数或波长表：

```ini
atmosphere.transmission=0.92
atmosphere.transmission=configs/atmosphere/my_transmission.csv
```

随发射高度和波长变化的 MODTRAN 总光学深度：

```ini
atmosphere.config=configs/atmosphere/modtran_4400_desert.cfg
```

程序在望远镜光学传播前应用：

```text
T = exp(-tau_total)
```

EventIO photon bunch 优先使用自身的发射高度；Photon CSV 可通过 `emission_altitude_km` 提供。大气模型使用物理高度和光子传播信息，不会重新定义相机或绘图坐标。

## 输出

平行光、点源和 Photon CSV 常用：

```ini
output.mode=hits
output.csv=run_logs/my_run/hits.csv
```

相机像素汇总：

```ini
output.mode=pixel
output.pixel_csv=run_logs/my_run/camera_pixels.csv
```

CORSIKA 相机正式输出推荐 dense HDF5：

```ini
output.format=hdf5
output.hdf5_path=run_logs/my_corsika_run/corsika_trace.h5
output.hdf5_storage=dense
```

ROOT/pylast 输出示例：

```text
configs/examples/corsika_lact_pylast_root_only_full_response.cfg
configs/examples/corsika_lact_root_full_response.cfg
```

HDF5 结构见 [HDF5 输出格式](docs/hdf5_output_format.md)，ROOT/pylast 数据层级见 [pylast 数据层级说明](docs/pylast_event_data_levels_zh.md)。

## 官方测试

重新构建、运行 CTest、执行标准光学案例并生成结果：

```bash
# 不使用 CORSIKA
tools/run_official_tests.sh --no-corsika

# 使用指定 CORSIKA/EventIO 文件
tools/run_official_tests.sh --corsika-file /path/to/input.zst
```

输出位于：

```text
run_logs/official_tests/
```

每个 cfg 的独立命令、预期输出和检查方法见[官方测试说明](docs/official_tests.md)。

## 进一步阅读

- [程序整体流程](docs/program_overview_zh.md)
- [坐标系定义](docs/coordinate_systems.md)
- [Photon CSV 格式](docs/photon_csv_format.md)
- [CORSIKA/EventIO 适配](docs/corsika_eventio_adapter.md)
- [相机时间和波形](docs/camera_timing_waveform_zh.md)
- [NSB 光谱模型](docs/nsb_spectral_model_zh.md)
- [结构遮挡 primitives](docs/obstruction_primitives_csv_zh.md)
- [官方测试](docs/official_tests.md)

英文详细手册及全部逐项命令见 [README_EN.md](README_EN.md)。
