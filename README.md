# LACT_sim 用户版 v1.0

这是 LACT 光学模拟程序的简洁用户版本，保留人工光源和 CORSIKA/EventIO
两类常用流程，并提供四个直接运行的示例配置：

1. 平行光白板测试：`configs/examples/parallel_light.cfg`
2. Photon CSV 白板测试：`configs/examples/photon_csv_source.cfg`
3. CORSIKA/EventIO 完整相机响应：`configs/examples/full_response_corsika.cfg`
4. CORSIKA/EventIO -> LACT ROOT 输出：`configs/examples/lactroot_only.cfg`

用户版不包含开发测试、额外示例和详细开发文档。核心源码仍保留在 `src/`、`include/` 和 `apps/` 中，并保留了常用画图脚本用于检查输出结果。

## 目录说明

```text
apps/                         程序入口和小工具
  run_optical_sim.cpp          平行光、点源、CSV 光子等非 CORSIKA 光学模拟
  run_corsika_trace.cpp        CORSIKA/EventIO 输入的完整响应模拟
  compute_nsb_rate.cpp         从夜天光谱估算每像素 NSB rate

configs/examples/             用户直接运行和修改的示例 cfg
configs/mirror_1229_facets.csv 默认镜片布局文件
configs/cameras/              new_camera 像素布局
configs/efficiency/           SiPM PDE 曲线
configs/atmosphere/           MODTRAN 大气透过率表
configs/nsb/                  夜天光谱和 NSB 示例

external/hessioxxx/            CORSIKA/EventIO 读取所需的 hessioxxx
tools/build_hessio.sh          编译 hessioxxx 的辅助脚本
tools/ensure_hdf5.sh           检查 HDF5；缺失时尝试自动安装
python/                        平行光光斑和 CORSIKA 相机图像画图脚本
```

## 编译

在仓库根目录运行：

```bash
make
```

这会先编译 `external/hessioxxx/source`，然后在 `build/` 下生成：

```text
build/run_optical_sim
build/run_corsika_trace
build/compute_nsb_rate
```

如果只想编译平行光示例、不需要 CORSIKA/EventIO 支持，可以运行：

```bash
make no-hessio
```

这时只会生成 `build/run_optical_sim`，不会生成 `run_corsika_trace`。

默认 `make` 会检查 HDF5 C 库。若系统没有 HDF5，`tools/ensure_hdf5.sh` 会尝试用当前系统包管理器自动安装：

```text
macOS:          brew install hdf5
Ubuntu/Debian:  sudo apt-get install -y libhdf5-dev
Fedora/RHEL:    sudo dnf/yum install -y hdf5-devel
```

如果机器没有下载或安装权限，脚本会报错退出。正式完整响应建议使用 HDF5 输出。

如果暂时没有可用的 HDF5，或者 macOS 上 HDF5 与编译器架构不一致，可以先关闭 HDF5：

```bash
make LACT_ENABLE_HDF5=OFF
```

这种情况下 CORSIKA 示例的 `output.format` 需要改成 `csv`，不要使用 `hdf5`。

LACT ROOT 输出是可选功能。若系统能找到 ROOT 6.24 或更新版本，`make` 会自动为 `run_corsika_trace` 打开 ROOT 输出；如果找不到 ROOT，普通平行光和 HDF5/CSV 示例仍可编译运行，只是 `lactroot_only.cfg` 不能运行。也可以显式关闭 ROOT 检测：

```bash
make LACT_ENABLE_ROOT=OFF
```

在 macOS 上，Makefile 会默认使用当前机器架构设置 `CMAKE_OSX_ARCHITECTURES`。如果你明确需要指定架构，可以这样运行：

```bash
make CMAKE_OSX_ARCHITECTURES=arm64
```

或：

```bash
make CMAKE_OSX_ARCHITECTURES=x86_64
```

常见依赖：

```bash
# Ubuntu/Debian
sudo apt-get install cmake g++ make zlib1g-dev libhdf5-dev

# CentOS/RHEL-like
sudo yum install cmake gcc-c++ make zlib-devel hdf5-devel
```

画图脚本需要额外的 Python 包：

```bash
python3 -m pip install -r requirements.txt
```

macOS 或 Linux 如果运行时找不到 `libhessio`，可以手动设置库路径。

Linux:

```bash
export LD_LIBRARY_PATH="$PWD/external/hessioxxx/source/lib:${LD_LIBRARY_PATH:-}"
```

macOS:

```bash
export DYLD_LIBRARY_PATH="$PWD/external/hessioxxx/source/lib:${DYLD_LIBRARY_PATH:-}"
```

## 示例 1：平行光白板测试

运行：

```bash
build/run_optical_sim configs/examples/parallel_light.cfg
```

输出：

```text
run_logs/parallel_light/hits.csv
```

这个例子使用理想平行光照射默认镜片布局，然后把反射光打到焦平面白板上。它适合用来检查镜面布局、焦距、光斑位置和基本光追是否正常。

画出焦平面光斑：

```bash
python3 python/plot_spot_histogram.py \
  run_logs/parallel_light/hits.csv \
  --output run_logs/parallel_light/spot_histogram.png \
  --title "Parallel light spot"
```

如果提示缺少 Python 包，先运行 `python3 -m pip install -r requirements.txt`。

### parallel_light.cfg 怎么改

```ini
telescope.pointing_az_deg=0
telescope.pointing_el_deg=90
```

望远镜指向。`pointing_el_deg=90` 表示指向天顶。平行光源定义在望远镜本地坐标中，改指向后程序会整体变换坐标。

```ini
mirror.csv_path=configs/mirror_1229_facets.csv
```

镜片布局 CSV。一般不用改；如果要测试新的镜面排布，换成新的 CSV 路径。

```ini
source.n_bunches=1000000
source.beam_radius_m=4
source.beam_direction=0,0,-1
source.random_seed=1229
```

平行光设置。`n_bunches` 是光子束数量，越大统计越稳定但越慢。`beam_radius_m` 是入射光束半径。`beam_direction=0,0,-1` 表示沿望远镜光轴入射。

```ini
source.coordinate_frame=telescope_local
```

表示人工光源的位置和方向使用原有望远镜本地光学坐标。它只控制输入如何进入
光追，不改变镜面、白板或画图坐标。

如果想把这个例子改成轴上点源，可以参考 cfg 里的中文注释，把 `source.mode` 改成 `PointSource`，并设置：

```ini
source.source_position=0,0,50
source.aperture_z=0
source.aperture_radius_m=4
```

这表示点源在望远镜本地 `+z` 方向 50 m 处，光子采样目标平面是 `z=0`。

```ini
output.plane_point=0,0,-8
output.csv=run_logs/parallel_light/hits.csv
```

白板输出平面和输出 CSV。当前焦距是 8 m，所以白板放在 `z=-8`。如果修改焦距，要同步检查这里。

## 示例 2：Photon CSV 白板测试

这个示例从 CSV 读取用户给定的光子位置和方向，再使用与平行光示例相同的
望远镜、镜面和白板进行光追。

运行：

```bash
build/run_optical_sim configs/examples/photon_csv_source.cfg
```

输入和输出：

```text
configs/examples/photon_csv_source.csv
run_logs/photon_csv_source/hits.csv
```

画出白板光斑：

```bash
python3 python/plot_spot_histogram.py \
  run_logs/photon_csv_source/hits.csv \
  --output run_logs/photon_csv_source/spot_histogram.png \
  --title "Photon CSV spot"
```

CSV 至少需要六列：

```csv
x_m,y_m,z_m,dir_x,dir_y,dir_z
```

示例还给出了 `time_ns`、`wavelength_nm`、`weight`、`multiplicity`、
`event_id` 和 `telescope_id`。用户通常只需要复制 cfg 后修改：

```ini
source.csv_path=/path/to/my_photons.csv
source.coordinate_frame=telescope_local
output.csv=run_logs/my_photons/hits.csv
```

下面的开关会把完成输入坐标解释后、实际进入光追的望远镜本地光子位置和方向
直接追加到输出 CSV：

```ini
output.whiteboard_input_photon=true
```

对应列为：

```text
input_local_x_m,input_local_y_m,input_local_z_m,
input_local_dir_x,input_local_dir_y,input_local_dir_z
```

如需使用其他输入坐标，可将 `source.coordinate_frame` 设置为
`corsika_nwu_relative`、`corsika_nwu_global` 或 `lact_generic_global`；一般用户
自行生成的望远镜局部光子文件推荐保持 `telescope_local`。

## 示例 3：CORSIKA 完整相机响应

这个示例包含：

```text
CORSIKA/EventIO 输入
默认镜片布局光追
仰角相关结构形变
小幅随机光学误差
mirror reflectivity 曲线
filter transmission 曲线
new_camera 像素布局
Bezier light collector
SiPM PDE
simple multiplicity trigger
```

默认会计算 trigger，但不会只保存触发事件：

```ini
output.save_only_triggered=false
```

如果正式运行时只想保存触发事件，改成：

```ini
output.save_only_triggered=true
```

NSB 默认关闭。建议先确认无 NSB 版本跑通，再按 cfg 里的中文注释打开。

运行：

```bash
build/run_corsika_trace configs/examples/full_response_corsika.cfg /path/to/input.zst
```

第二个参数是 CORSIKA/EventIO 文件路径。也可以把路径直接写进 cfg 的 `source.eventio_path`，但推荐命令行传入，方便重复跑不同输入文件。

输出默认在：

```text
run_logs/full_response_corsika/corsika_trace.h5
run_logs/full_response_corsika/camera_pixel_image.csv
run_logs/full_response_corsika/corsika_trace_summary.csv
```

### full_response_corsika.cfg 怎么改

```ini
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70
```

望远镜指向。CORSIKA 常给出 `zenith` 和 `azimuth`，这里需要填：

```text
pointing_el_deg = 90 - zenith
pointing_az_deg = azimuth
```

例如 CORSIKA 输入是 `zenith=20 deg, azimuth=0 deg`，则写：

```ini
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70
```

```ini
mirror.csv_path=configs/mirror_1229_facets.csv
camera.csv_path=configs/cameras/new_camera_pixels.csv
```

望远镜结构文件。`mirror.csv_path` 是镜片布局，`camera.csv_path` 是相机像素布局。用户替换结构时主要改这两个路径。

```ini
error.structural_deformation_config=configs/mirrors/mirror_1229_elevation_series.cfg
error.facet_radial_position_sigma_m=0.002
error.facet_normal_sigma_deg=0.01
error.reflect_direction_sigma_deg=0.0179
```

结构形变和随机误差设置。`structural_deformation_config` 会根据当前望远镜仰角读取 `configs/mirror_1229_elevation_series.csv`。几个 `sigma` 参数是小幅随机扰动，默认值来自开发版 full-response 测试。

```ini
efficiency.mirror_reflectivity=configs/efficiency/mirror_reflectivity.csv
efficiency.filter_transmission=configs/efficiency/filter_transmission.csv
```

镜面反射率和滤光片透过率，按波长插值。

```ini
camera.collector=bezier
camera.collector_exit_size_m=0.0130
camera.collector_height_m=0.0297
```

光收集器设置。默认使用 new_camera 的 Bezier 型光收集器参数。

```ini
sipm.size_m=0.0130
sipm.pde=configs/efficiency/sipm_pde.csv
```

SiPM 设置。`sipm.pde` 是波长相关探测效率曲线。如果想先做理想光子计数，可以注释掉 `sipm.pde`。

```ini
atmosphere.transmission=none
```

大气额外透过率默认关闭。CORSIKA 光子通常已经是到达望远镜附近的光子；如果需要额外使用 MODTRAN 表，可以在 cfg 里改成：

```ini
atmosphere.config=configs/atmosphere/modtran_4400_desert.cfg
```

这个表会按光子波长和发射高度计算 `exp(-tau)`。如果输入 EventIO 里没有可用的发射高度，开启 MODTRAN 后程序会报错提醒。

```ini
nsb.enabled=false
nsb.model=constant_rate
```

NSB 默认关闭。简单常数模型可以直接设置：

```ini
nsb.enabled=true
nsb.model=constant_rate
nsb.rate_pe_per_ns_per_pixel=0.05
nsb.window_ns=16
```

如果希望从夜天光谱自动估算 rate，可以用 spectral 模型：

```ini
nsb.enabled=true
nsb.model=spectral_flux
nsb.spectrum_csv=configs/nsb/nsb_spectrum_skycalc_dark.csv
nsb.spectrum_unit=ph_s_nm_sr_m2
nsb.effective_area_m2=29.623570
nsb.pixel_solid_angle=auto
nsb.window_ns=16
```

可以先单独检查 spectral NSB 算出来的 rate：

```bash
build/compute_nsb_rate configs/examples/full_response_corsika.cfg
```

注意：`full_response_corsika.cfg` 默认 `nsb.enabled=false`。要检查 spectral NSB 时，先按 cfg 里的中文注释打开 `nsb.model=spectral_flux`，或者把 cfg 中的 NSB 部分替换为：

```ini
nsb.config=configs/nsb/spectral_skycalc_dark_no_obstruction.cfg
```

```ini
trigger.enabled=true
trigger.pixel_threshold_pe=10
trigger.camera_multiplicity=3
trigger.array_multiplicity=2
trigger.coincidence_window_ns=50
```

简单 multiplicity trigger。默认会计算触发结果，但 `output.save_only_triggered=false`，所以第一次检查时不会因为没触发而没有输出。

```ini
source.mode=EventIO
source.eventio_coordinate_frame=corsika_iact
source.use_eventio_telescope_position=true
source.missing_wavelength_model=cherenkov
source.eventio_2d_input_plane_z_m=0
source.eventio_2d_plane_mode=auto
```

CORSIKA/EventIO 输入设置。通常不用改。`corsika_iact` 表示按 hessio/CORSIKA IACT 光子束坐标约定读入，再根据望远镜指向旋转到 LACT 本地光学坐标。

`missing_wavelength_model=cherenkov` 用于处理 EventIO photon bunch 中 `lambda=0` 或缺失波长的情况。程序会按 `1/lambda^2` 随机补波长，波长范围优先从 CORSIKA input card 的 `CWAVLG` 行读取；如果输入文件没有 `CWAVLG`，可以在 cfg 里手动设置：

```ini
source.missing_wavelength_min_nm=260
source.missing_wavelength_max_nm=1000
source.missing_wavelength_seed=246813579
```

`eventio_2d_input_plane_z_m` 用于 2D EventIO photon bunch，表示记录平面在望远镜本地坐标的 z 位置。默认 `0`。`eventio_2d_plane_mode=auto` 会自动判断向前追迹还是回投到镜面；只有明确知道输入平面约定时才需要改成 `forward` 或 `backproject`。

```ini
source.max_shower_events=10
# source.filter_shower_event_id=327666
```

普通 CORSIKA 示例默认处理前 10 个 shower；需要时仍可筛选指定 shower。

```ini
output.format=hdf5
output.hdf5_path=run_logs/full_response_corsika/corsika_trace.h5
output.pixel_csv=run_logs/full_response_corsika/camera_pixel_image.csv
output.summary_csv=run_logs/full_response_corsika/corsika_trace_summary.csv
```

输出设置。正式结果建议用 HDF5。调试时可以把 `output.format` 改成 `csv` 或 `both`。

## 示例 4：LACT ROOT-only 输出

这个示例基于完整响应链，但只写 LACT ROOT 文件，适合后续接 pylast 或 ROOT quicklook：

```bash
build/run_corsika_trace configs/examples/lactroot_only.cfg /path/to/input.zst
```

默认输出：

```text
run_logs/lactroot_only/lact_events.root
```

这个 cfg 默认开启：

```text
MODTRAN atmosphere
spectral NSB
simple multiplicity trigger
timeseries_pe waveform
```

ROOT 输出相关参数：

```ini
output.format=root
output.lact_root_enabled=true
output.lact_root_path=run_logs/lactroot_only/lact_events.root
output.lact_profile=timeseries_pe
output.lact_root_write_components=false
```

`output.format=root` 表示不写 HDF5/CSV，只写 ROOT。`output.lact_profile=timeseries_pe` 会写适合后续 waveform/readout 检查的 p.e. 时间序列。

该 ROOT 输出示例默认设置 `source.max_shower_events=-1`，处理输入文件中的全部
shower。

默认只保存通过 trigger 的望远镜事件：

```ini
output.save_only_triggered=true
```

如果只是排查光追或 ROOT 写出链路，确实需要保留未触发事件，可以临时改成：

```ini
output.save_only_triggered=false
```

如果只想先检查 ROOT 输出链路、暂时不加 NSB，可以改：

```ini
nsb.enabled=false
```

如果没有 ROOT，运行这个 cfg 会提示 `output.lact_root_enabled=true` 但编译时没有启用 ROOT；这时需要安装 ROOT 后重新 `make`，或先运行前面几个非 ROOT 示例。

### 画 CORSIKA 相机图像

如果使用 CSV 输出，先确认 cfg 中包含：

```ini
output.format=csv
```

或：

```ini
output.format=both
```

然后画出某一个 event 的所有望远镜相机图像：

```bash
python3 python/plot_corsika_trace_output.py \
  run_logs/full_response_corsika/camera_pixel_image.csv \
  --event-number 1 \
  --camera-csv configs/cameras/new_camera_pixels.csv \
  --output-dir run_logs/full_response_corsika/plots/event_1
```

如果使用 HDF5 输出，可以直接从 `corsika_trace.h5` 画图：

```bash
python3 python/plot_hdf5_camera.py \
  run_logs/full_response_corsika/corsika_trace.h5 \
  --event-number 1 \
  --output run_logs/full_response_corsika/plots/hdf5_event_1
```

这两个脚本沿用开发版的画图风格，会按相机像素布局画相机图像；不指定 `--telescope-id` 时，会把所选 event 里的所有望远镜都画出来。
