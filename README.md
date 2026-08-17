# LACT_sim 用户手册

## 编译

Ubuntu/Debian 常用依赖：

```bash
sudo apt-get install cmake g++ make zlib1g-dev libhdf5-dev
python3 -m pip install -r requirements.txt
make
```

`make` 会先构建仓库内的 hessioxxx，再生成：

```text
build/run_optical_sim
build/run_corsika_trace
build/compute_nsb_rate
```

只使用人工光源、不需要 EventIO 时：

```bash
make no-hessio
```

HDF5 可通过 `make LACT_ENABLE_HDF5=OFF` 关闭。LACT ROOT 是可选功能，编译器找到 ROOT 6.24+ 时自动启用；也可用 `make LACT_ENABLE_ROOT=OFF` 关闭。

两个 ROOT 示例的画图命令需要安装我们维护的 [pyLAST `lact_sim` 分支](https://github.com/Yun532/pylast/tree/lact_sim)。该分支包含 LACT ROOT adapter、触发时间和阵列可视化接口；安装后请确认下面的导入成功：

```bash
python3 -c "from pylast.io import LactEventSource"
```

## 1. 平行光白板

```bash
build/run_optical_sim configs/examples/parallel_light.cfg
python3 python/plot_spot_histogram.py \
  run_logs/parallel_light/hits.csv \
  --output run_logs/parallel_light/spot.png \
  --title "Parallel light"
```

常改参数是 `source.n_bunches`、`source.beam_radius_m`、`source.beam_direction` 和输出路径。
该配置显式使用 `camera.mode=whiteboard`；相机与电子学均不参与本测试。

## 2. Photon CSV 完整相机与 pyLAST

```bash
build/run_corsika_trace configs/examples/photon_csv_source.cfg

python3 python/plot_photon_csv_root_pylast.py \
  run_logs/photon_csv_source/lact_events.root \
  --event-id 1909 \
  --output run_logs/photon_csv_source/pylast_camera.png
```

该示例读取 event 1909、telescope 19 的 CORSIKA 二维 bunch 六列数据，通过 1664 像素相机链生成随机 p.e.，写入 `image_pe` LACT ROOT，再由 `pylast.io.LactEventSource` 和 `pylast.visualize` 绘图。波形关闭。
用户版默认保持线性 p.e. 响应（电子学关闭），但完整电子学模块仍可通过配置开启。

坐标约定是：LACT_sim 白板和探测器直接输出物理焦平面 `(u, v)`；`LactEventSource` 在 pyLAST 输入边界统一映射为 `pix_x=-v`、`pix_y=-u`。该脚本直接调用 pyLAST 原生 `EventVisualizer.plot_event()`，不再自行交换、翻转或旋转相机坐标。白板/HDF5 诊断图则直接显示文件中保存的物理 `(u, v)`。

从安装环境、拉取两个仓库到 CSV 字段说明和出图命令，见 [Photon CSV 从零运行指南](PHOTON_CSV_QUICKSTART.md)。

输入 CSV 的最小列为：

```csv
x_m,y_m,z_m,dir_x,dir_y,dir_z
```

示例输入保存在 `configs/sources/event1909_tel19_minimal_photons.csv`。最新 Photon CSV 读取器还支持每行的 `time_ns`、`wavelength_nm`、`weight`、`multiplicity`、`event_id`、`telescope_id`、`emission_altitude_km`、`eventio_2d` 和阵列位置字段。从 EventIO 二维 bunch 导出的 CSV 使用 `corsika_nwu_relative` 坐标、`eventio_2d=true` 和 `-16 m` 参考平面。

## 3. CORSIKA 完整相机响应

```bash
build/run_corsika_trace \
  configs/examples/full_response_corsika.cfg \
  /path/to/input.zst

python3 python/plot_hdf5_camera.py \
  run_logs/full_response_corsika/corsika_trace.h5 \
  --image-index 0 --quantity pe \
  --output run_logs/full_response_corsika/camera.png
```

必须令 `telescope.pointing_el_deg = 90 - CORSIKA zenith_deg`。示例使用 2026-06-22 实测光学和 1664 像素相机，波形关闭；处理前 10 个 shower、关闭 NSB，并保留所有 `total_pe > 0` 的望远镜（包括未触发望远镜），便于第一次核对输入；全零望远镜和全零事件不会写入图像表。该示例显式开启阵列 multiplicity 判决。需要生产筛选时可设 `source.max_shower_events=-1`、`output.save_only_triggered=true`；需要夜天光时参考下一个 ROOT 示例中的 `nsb.*`。
该用户示例默认关闭显式电子学，以保持原 `user_v2.0` 的输出定义。

## 4. LACT ROOT / pyLAST

```bash
build/run_corsika_trace \
  configs/examples/lactroot_only.cfg \
  /path/to/input.zst

python3 python/plot_lact_root_pylast.py \
  run_logs/lactroot_only/lact_events.root \
  --event-index 0 \
  --output run_logs/lactroot_only/pylast_cameras.png
```

这个配置默认处理全部 shower、加入光谱 NSB、只保留达到相机触发的望远镜，使用 2026-06-22 实测光学和 1664 像素相机，直接写 `image_pe`，不生成波形。为保持旧 `user_v2.0` ROOT 示例的实际保存行为，阵列级过滤显式关闭；`trigger.array_multiplicity` 不参与这个示例的输出筛选。绘图脚本使用 `LactEventSource → EventVisualizer.plot_event(image_level="dl0")`，显示随机 p.e. 积分图像。首次小规模检查可暂时把 `source.max_shower_events` 改为 `10`。
几何焦距仍为 8.0 m；写给 pyLAST 重建使用的有效焦距独立设为 8.1787 m。用户版默认关闭显式电子学饱和。

LACT ROOT 保留光学焦平面的原始 `u/v`。pyLAST `LactEventSource` 在输入边界统一映射为 `pix_x=-v`、`pix_y=-u`；用户绘图脚本和 notebook 不再增加 LACT 专用旋转。

完整的事件读取、相机图像、触发时延与 p.e.、Hillas 和 SDP 重建流程见 [完整 pyLAST Jupyter notebook](notebooks/lact_root_to_pylast_visualize.ipynb)。

所有相对路径都以仓库根目录为工作目录解析，请从仓库根目录运行命令。
