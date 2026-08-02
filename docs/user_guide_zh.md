# LACT_sim 用户使用手册

[返回中文 README](../README.md) | [English guide](user_guide_en.md)

本文给出从安装、编译到运行和画图的完整用户流程。各物理模型和文件格式的内部细节由文末专题文档说明。

## 安装和编译

### 获取代码

```bash
git clone https://github.com/Yun532/LACT_sim.git
cd LACT_sim
```

### 系统依赖

Ubuntu/Debian：

```bash
sudo apt-get update
sudo apt-get install cmake g++ make zlib1g-dev libhdf5-dev python3 python3-pip
```

CentOS/RHEL 系列：

```bash
sudo yum install cmake gcc-c++ make zlib-devel hdf5-devel python3
```

基本绘图工具：

```bash
python3 -m pip install numpy pandas matplotlib h5py
```

ROOT 转换脚本还可能需要 `uproot`；LACT ROOT 原生输出需要服务器提供 ROOT 6.24 或更新版本。

### 标准编译

```bash
make
make test
```

`make` 会先编译 `external/hessioxxx`，再在 `build/` 中配置和编译 LACT_sim。成功后至少应有：

```text
build/run_optical_sim
build/run_corsika_trace
build/compute_nsb_rate
```

只运行平行光、点源或 Photon CSV：

```bash
make no-hessio
ctest --test-dir build --output-on-failure
```

禁用 ROOT 但保留 EventIO/HDF5：

```bash
make no-root
```

如果已安装 ROOT：

```bash
source /path/to/root/bin/thisroot.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DHESSIO_ROOT="$PWD/external/hessioxxx/source" \
  -DLACT_ENABLE_ROOT=ON \
  -DROOT_DIR="$(root-config --cmakedir)"
cmake --build build -j4
```

### 无图形界面的服务器

```bash
export MPLBACKEND=Agg
export MPLCONFIGDIR=/tmp
```

## 官方测试总入口

不需要 CORSIKA：

```bash
tools/run_official_tests.sh --no-corsika \
  2>&1 | tee run_logs/official_tests/run_no_corsika.log
```

包含 CORSIKA：

```bash
tools/run_official_tests.sh --corsika-file /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/run_with_corsika.log
```

该脚本会重新配置构建、运行 CTest、执行标准测试并生成图片。所有结果位于 `run_logs/official_tests/`。全部测试项目和检查标准见[官方测试说明](official_tests.md)。

## 人工光源官方测试

### 理想平行光

使用 `configs/official_tests/perfect_parallel_whiteboard.cfg`。它组合完整 54 片镜面、100 万轴上平行光和 8 m 白板：

```bash
mkdir -p run_logs/official_tests/perfect_parallel
build/run_optical_sim configs/official_tests/perfect_parallel_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/perfect_parallel/run.log

MPLBACKEND=Agg python3 python/plot_spot_histogram.py \
  run_logs/official_tests/perfect_parallel/hits.csv \
  --config configs/official_tests/perfect_parallel_whiteboard.cfg \
  --output run_logs/official_tests/perfect_parallel/spot.png \
  --max-bins 520 --dpi 350
```

主要结果：

```text
run_logs/official_tests/perfect_parallel/hits.csv
run_logs/official_tests/perfect_parallel/spot.png
run_logs/official_tests/perfect_parallel/run.log
```

### 900 m 点源

```bash
mkdir -p run_logs/official_tests/point_900m
build/run_optical_sim configs/official_tests/perfect_point_900m_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/point_900m/run.log

MPLBACKEND=Agg python3 python/plot_spot_histogram.py \
  run_logs/official_tests/point_900m/hits.csv \
  --config configs/official_tests/perfect_point_900m_whiteboard.cfg \
  --output run_logs/official_tests/point_900m/spot.png \
  --max-bins 520 --dpi 350
```

这个 cfg 将平行光源换成 900 m 轴上点源，并使用有限物距对应的输出平面。

### 平行光加三维结构遮挡

```bash
mkdir -p run_logs/official_tests/raytrace_structure_parallel
build/run_optical_sim \
  configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/raytrace_structure_parallel/run.log

MPLBACKEND=Agg python3 python/plot_spot_histogram.py \
  run_logs/official_tests/raytrace_structure_parallel/hits.csv \
  --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  --output run_logs/official_tests/raytrace_structure_parallel/spot.png

MPLBACKEND=Agg python3 python/plot_optical_layout_3d.py \
  --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  --show-obstruction \
  --output run_logs/official_tests/raytrace_structure_parallel/layout_3d.png

python3 python/plot_optical_layout_html.py \
  --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  --output run_logs/official_tests/raytrace_structure_parallel/layout_3d.html
```

遮挡由下面一行开启：

```ini
obstruction.config=../obstructions/raytrace_final_structure.cfg
```

### 仰角形变扫描

标准配置为：

```text
configs/official_tests/deformation_parallel_whiteboard.cfg
```

一键脚本会自动扫描多个仰角并生成光斑、镜面布局和汇总 CSV。独立扫描参数见[官方测试说明的仰角扫描章节](official_tests.md#5-支架形变仰角扫描)。

## Photon CSV

### 输入文件

最少需要：

```csv
x_m,y_m,z_m,dir_x,dir_y,dir_z
0.0,0.0,1.0,0.0,0.0,-1.0
```

推荐 cfg：

```text
configs/examples/photon_csv_local_whiteboard.cfg
```

运行：

```bash
build/run_optical_sim configs/examples/photon_csv_local_whiteboard.cfg
```

复制 cfg 后，用户主要修改：

```ini
source.csv_path=/path/to/photons.csv
source.coordinate_frame=telescope_local
output.csv=run_logs/my_photons/hits.csv
```

如果需要保存实际传入光追的位置和方向：

```ini
output.whiteboard_input_photon=true
```

然后可生成包含入射方向的三维 HTML：

```bash
python3 python/plot_optical_layout_html.py \
  --config configs/examples/photon_csv_local_whiteboard.cfg \
  --input-photon-csv run_logs/examples/photon_csv_local_whiteboard/hits.csv \
  --output run_logs/examples/photon_csv_local_whiteboard/layout_with_photons.html
```

可选列和外部文件归一化见 [Photon CSV 格式](photon_csv_format.md)；CORSIKA
NWU、东起始 ENU、望远镜本地坐标及画图方向见
[中文坐标系说明](coordinate_systems_zh.md)。

## CORSIKA 白板测试

这个测试不接相机，用来先检查 EventIO 读取、望远镜指向和白板光斑：

```bash
mkdir -p run_logs/official_tests/corsika/whiteboard_plots
build/run_corsika_trace \
  configs/official_tests/corsika_whiteboard.cfg \
  /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/whiteboard_run.log
```

主要输出：

```text
run_logs/official_tests/corsika/whiteboard_hits.csv
run_logs/official_tests/corsika/whiteboard_summary.csv
```

如果已知输出 event id：

```bash
MPLBACKEND=Agg python3 python/plot_corsika_trace_output.py \
  run_logs/official_tests/corsika/whiteboard_hits.csv \
  --summary-csv run_logs/official_tests/corsika/whiteboard_summary.csv \
  --event-id OUTPUT_EVENT_ID \
  --output-dir run_logs/official_tests/corsika/whiteboard_plots
```

一键脚本会自动选择一个共同事件并写入：

```text
run_logs/official_tests/corsika/plots/selected_event.env
```

## CORSIKA 相机和完整响应

### 理想相机 HDF5

```bash
build/run_corsika_trace \
  configs/official_tests/corsika_new_camera.cfg \
  /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_run.log
```

这个 cfg 使用标准镜面、真实相机像素和 light collector，但关闭额外大气、NSB、trigger 和 PDE 损失。输出：

```text
run_logs/official_tests/corsika/camera_dense.h5
```

快速画第一幅图：

```bash
MPLBACKEND=Agg python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --image-index 0 --quantity pe \
  --output run_logs/official_tests/corsika/camera_first.png
```

### NSB 和 trigger

```bash
build/run_corsika_trace \
  configs/official_tests/corsika_nsb_trigger_camera.cfg \
  /path/to/input.zst
```

该 cfg 主要增加：

```ini
sipm.config=../sipm/new_camera_sipm.cfg
electronics.config=../electronics/explicit_microcell_saturation_only.cfg
efficiency.config=../efficiency/curves_all.cfg
nsb.config=../nsb/spectral_skycalc_dark_no_obstruction.cfg
trigger.config=../trigger/camera_pe_count_array_off.cfg
output.hdf5_write_components=true
waveform.enabled=false
```

当前正式链路把 Cherenkov 和 NSB 的 primary p.e. 一起送入显式微单元饱和，
然后用饱和后的 PE 做相机触发；不生成单 PE 或采样波形。

### 完整响应

```bash
build/run_corsika_trace \
  configs/official_tests/corsika_full_response_camera.cfg \
  /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_full_response_run.log
```

这个 cfg 启用结构形变、随机误差、反射率、滤光片、MODTRAN 大气、SiPM PDE、结构遮挡和 trigger，输出：

```text
run_logs/official_tests/corsika/camera_full_response_dense.h5
```

相机图：

```bash
MPLBACKEND=Agg python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_full_response_dense.h5 \
  --image-index 0 --quantity pe \
  --output run_logs/official_tests/corsika/full_response_first.png
```

阵列图：

```bash
MPLBACKEND=Agg python3 python/plot_hdf5_array_layout.py \
  run_logs/official_tests/corsika/camera_full_response_dense.h5 \
  --shower-event-number 1 --array-id 0 --quantity pe --log-color \
  --output run_logs/official_tests/corsika/full_response_array.png
```

如果需要按事件画所有望远镜，应使用 `python/select_hdf5_event.py` 或一键脚本生成的 `selected_event.env`。完整参数见[官方测试说明](official_tests.md)。

### ROOT/pylast

ROOT-only 示例：

```bash
build/run_corsika_trace \
  configs/examples/lactroot_only.cfg \
  /path/to/input.zst
```

默认输出：

```text
run_logs/lactroot_only/lact_events.root
```

快速画图：

```bash
python3 python/plot_lact_root_output.py \
  run_logs/lactroot_only/lact_events.root \
  --outdir run_logs/lactroot_only/root_quicklook
```

该 cfg 默认使用 `timeseries_pe`，处理输入文件中的全部 shower，并且只保存通过
trigger 的望远镜事件。ROOT 文件结构、pylast 数据层级和 notebook 流程见：

- [ROOT 输出检查](server_root_output_check_zh.md)
- [pylast 数据层级](pylast_event_data_levels_zh.md)
- `notebooks/lact_root_to_pylast_visualize.ipynb`

## 用户修改 cfg 时常用的字段

```ini
# 望远镜指向
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70

# 输入
source.mode=EventIO
source.eventio_path=
source.max_shower_events=10

# 输出
output.format=hdf5
output.hdf5_path=run_logs/my_run/output.h5
output.save_only_triggered=false
```

EventIO 文件通常作为命令行第二个参数传入，因此 `source.eventio_path` 可以保持为空。人工光源和 Photon CSV 使用 `run_optical_sim`；EventIO 使用 `run_corsika_trace`。

## 常见问题

- 找不到 HDF5：安装 HDF5 C 开发包后删除 `build/` 并重新运行 `make`。
- 找不到 ROOT：先激活 ROOT 环境，再重新运行 CMake；不需要 ROOT 时使用 `make no-root`。
- 服务器不能画图：设置 `MPLBACKEND=Agg` 和 `MPLCONFIGDIR=/tmp`。
- CORSIKA 输出为空：先用 `corsika_whiteboard.cfg`，并暂时关闭 `output.save_only_triggered` 或减少筛选条件。
- 不确定 event id：运行一键脚本，或使用 `python/select_hdf5_event.py` 检查 HDF5。

## 专题文档

- [所有官方测试及检查标准](official_tests.md)
- [程序整体流程](program_overview_zh.md)
- [Photon CSV 格式](photon_csv_format.md)
- [坐标系中文说明](coordinate_systems_zh.md)
- [坐标系英文参考](coordinate_systems.md)
- [CORSIKA/EventIO 适配](corsika_eventio_adapter.md)
- [HDF5 输出格式](hdf5_output_format.md)
- [相机时间和波形](camera_timing_waveform_zh.md)
- [NSB 光谱模型](nsb_spectral_model_zh.md)
- [结构遮挡](obstruction_primitives_csv_zh.md)
