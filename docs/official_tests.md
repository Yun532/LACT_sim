# 官方测试说明

这份文档解释 `tools/run_official_tests.sh` 会运行哪些测试，以及每个
`configs/official_tests/*.cfg` 里的配置项具体代表什么。目标是让新用户只看
cfg 和这份文档，就能知道每个测试在测什么、怎么单独运行、应该检查哪些输出。

运行不依赖 CORSIKA/EventIO 的测试：

```bash
tools/run_official_tests.sh --no-corsika 2>&1 | tee run_logs/official_tests/run_no_corsika.log
```

运行完整测试，需要给一个 CORSIKA/EventIO 输入文件：

```bash
tools/run_official_tests.sh --corsika-file /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/run_all.log
```

默认编译目录是 `build/`。如果想用其它编译目录：

```bash
LACT_BUILD_DIR=build_release tools/run_official_tests.sh --no-corsika
```

## 共同配置规则

official cfg 基本都是“装配文件”：它们引用 telescope、mirror、source、
output、camera、sipm、NSB、trigger 等模块 cfg。每个 official cfg 只写本测试
真正需要改变或显式启用的配置。

常见配置项含义：

- `telescope.config`: 望远镜元信息和默认指向。官方测试一般用
  `configs/official_tests/telescope_1229_minimal.cfg`。
- `telescope.pointing_az_deg`: 望远镜指向方位角，单位 degree。
- `telescope.pointing_el_deg`: 望远镜指向仰角，单位 degree。
- `mirror.config`: 镜片布局。标准 1229 镜面是
  `configs/mirrors/mirror_1229_imported.cfg`。
- `source.config`: 非 CORSIKA 测试用的人工光源。
- `source.mode=EventIO`: CORSIKA/EventIO 输入模式。
- `source.eventio_path`: EventIO 文件路径。official cfg 里留空，因为运行时通过
  命令行传入，避免用户每次改 cfg。
- `source.event_id_mode=event_array100`: 输出事件编号采用
  `event_id = shower_event * 100 + array_id`。
- `source.eventio_coordinate_frame=corsika_iact`: 使用当前 EventIO/CORSIKA IACT
  坐标转换约定。
- `output.config`: 白板或相机焦平面位置。
- `output.format=csv`: 输出光子命中 CSV。
- `output.format=hdf5`: 输出 HDF5。
- `output.hdf5_storage=dense`: 每个图像保存完整 1616 像素，便于后续加 NSB/trigger。
- `output.hdf5_write_components=true`: 额外写出 `cherenkov_pe` 和 `nsb_pe`。
- `output.save_only_triggered=true`: HDF5 只保留通过 telescope trigger 的图像。
- `camera.config`: 相机像素和光收集器。`new_camera.cfg` 使用真实像素表、
  Bezier square-cone 光收集器和 `true_reflect` 材料。
- `sipm.config`: SiPM 有效面积和 PDE。理想测试用 `ideal_sipm.cfg`，PDE 不启用。
- `efficiency.config`: 镜面反射率、滤光片透过率等效率曲线。
- `atmosphere.config`: CORSIKA 到达望远镜之后的额外大气透过率。`ideal.cfg` 表示不额外衰减。
- `nsb.config`: 夜空背景光。`ideal.cfg` 表示关闭。
- `trigger.config`: 简单 multiplicity trigger。`disabled.cfg` 表示关闭。
- `obstruction.config`: 3D 遮挡模型。
- `error.config`: 镜片形变和随机误差项。
- `propagation.speed_of_light_m_per_ns`: 局部望远镜光路中的光速，默认真空光速。

注意：official cfg 现在不写 `electronics.config=../electronics/ideal_pe.cfg`。
目前 electronics 模块只是后续波形/电子学的 placeholder，没有实际响应；SiPM PDE
只通过 `sipm.pde` 设置，避免同一个 p.e. 转换逻辑出现两个入口。

## 一键脚本包含的测试

`tools/run_official_tests.sh --no-corsika` 会运行：

1. C++ 单元测试和小型回归测试。
2. 完美平行光白板。
3. 完美 900 m 点源白板。
4. 平行光 + 3D 遮挡白板。
5. 30 m 点源 + 3D 遮挡白板。
6. 支架形变 0 到 90 度仰角扫描。
7. 光收集器角响应测试。
8. 效率曲线验证测试。

如果没有 `--no-corsika`，脚本还会继续运行：

1. CORSIKA 白板 debug。
2. CORSIKA 完美相机 dense HDF5。
3. CORSIKA + NSB + trigger。
4. CORSIKA + 3D 遮挡 + NSB + trigger。
5. CORSIKA full-response smoke test。

CORSIKA 阶段的画图会先自动选择一个共同 event：

```bash
python3 python/select_hdf5_event.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5 \
  run_logs/official_tests/corsika/camera_obstruction_nsb_trigger_dense.h5 \
  run_logs/official_tests/corsika/camera_full_response_dense.h5 \
  --output-env run_logs/official_tests/corsika/plots/selected_event.env
```

选择规则：这个 `event_id` 必须在所有 camera HDF5 中都有可画图像，并且优先选择总
p.e. 较高的事件。后续 whiteboard、完美相机、NSB+trigger、遮挡+NSB+trigger、
full-response 和阵列芯位图都使用同一个 `event_id`。输出目录形如：

```text
run_logs/official_tests/corsika/plots/event_<event_id>/
```

## 0. 效率曲线验证测试

脚本会在非 CORSIKA 阶段运行：

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_efficiency_curves.py \
  --output-dir run_logs/official_tests/efficiency_curves \
  2>&1 | tee run_logs/official_tests/efficiency_curves/run.log
```

测试内容：分别检查镜面反射率、滤光片透过率、SiPM PDE 和大气透过率曲线，再检查它们的
逐波长乘积。默认大气项关闭，所以大气透过率图是一条 `1` 的理想曲线。

默认输入：

- `configs/efficiency/mirror_reflectivity.csv`: 镜面反射率。
- `configs/efficiency/filter_transmission.csv`: 相机滤光片透过率。
- `configs/efficiency/sipm_pde.csv`: SiPM PDE。
- `atmosphere_transmission=none`: 不额外加入大气衰减。

输出：

- `efficiency_summary.png/pdf`: 推荐检查图，所有单项效率和总效率画在同一张图上；同色虚线是
  输入表格理论线，同色实线是程序实际使用线。
- `mirror_reflectivity.png`: 原始镜面反射率点和程序插值曲线。
- `filter_transmission.png`: 原始滤光片点和程序插值曲线。
- `sipm_pde.png`: 原始 SiPM PDE 点和程序插值曲线。
- `atmosphere_transmission.png`: 大气项，未配置时为理想 `1`。
- `total_efficiency.png`: 所有已启用效率项的总乘积。
- `efficiency_curve_samples.csv`: 逐波长采样后的数值，方便人工核对。
- `efficiency_curve_report.txt`: 输入行数、重复波长、效率范围和总效率范围。

检查重点：灰点是输入表格点，橙色虚线是输入表格按波长连起来的理论/文件线，蓝线是程序
实际使用的曲线。同一波长有多个输入值时，程序会先取平均再插值；这个行为由 C++ 测试
`test_efficiency_curves` 锁定。更详细的说明见
`docs/efficiency_validation.md`。

## 1. 完美平行光白板测试

cfg:

```text
configs/official_tests/perfect_parallel_whiteboard.cfg
```

测试内容：理想 54 片镜子、无效率损失、无误差、无遮挡，用 100 万平行光打到 8 m
白板上，检查基础 PSF。

cfg 逐项解释：

- `telescope.config=telescope_1229_minimal.cfg`: 使用 1229 望远镜元信息。默认指向
  `az=0, el=90`，也就是竖直向上。
- `mirror.config=../mirrors/mirror_1229_imported.cfg`: 使用完整 54 片标准镜面。
- `source.config=../sources/parallel_1M_on_axis.cfg`: 平行光源，`n_bunches=1000000`，
  采样半径 `beam_radius_m=4`，方向 `beam_direction=0,0,-1`。
- `output.config=../outputs/whiteboard_f8.cfg`: 白板在本地 `z=-8 m`，法向量
  `0,0,-1`，图像坐标由 `plane_u_axis=1,0,0` 和 `plane_v_axis=0,1,0` 定义。
- `propagation.speed_of_light_m_per_ns=0.299792458`: 局部传播光速。
- `output.csv=run_logs/official_tests/perfect_parallel/hits.csv`: 输出白板命中光子。

单独运行：

```bash
mkdir -p run_logs/official_tests/perfect_parallel
build/run_optical_sim configs/official_tests/perfect_parallel_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/perfect_parallel/run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_spot_histogram.py \
  run_logs/official_tests/perfect_parallel/hits.csv \
  --output run_logs/official_tests/perfect_parallel/spot.png \
  --max-bins 520 --dpi 350
```

检查重点：`blocked_by_obstruction=0`，`hit_output_plane` 等于
`hit_output_before_obstruction`，spot 应该是理想光学 PSF。

## 2. 完美 900 m 点光源白板测试

cfg:

```text
configs/official_tests/perfect_point_900m_whiteboard.cfg
```

测试内容：理想镜面下，900 m 处轴上点光源在有限距离焦平面附近形成的光斑。

cfg 逐项解释：

- `telescope.config=telescope_1229_minimal.cfg`: 使用默认竖直望远镜。
- `mirror.config=../mirrors/mirror_1229_imported.cfg`: 完整 54 片标准镜面。
- `source.config=../sources/point_900m_on_axis.cfg`: 点源在本地 `z=900 m`，
  程序从 `aperture_z=0`、半径 4 m 的入瞳采样发射光线。
- `output.config=../outputs/whiteboard_point_900m_focus.cfg`: 900 m 物距对应的近似
  成像平面，`z=-8.0717488789 m`。
- `output.csv=run_logs/official_tests/point_900m/hits.csv`: 输出白板命中光子。

单独运行：

```bash
mkdir -p run_logs/official_tests/point_900m
build/run_optical_sim configs/official_tests/perfect_point_900m_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/point_900m/run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_spot_histogram.py \
  run_logs/official_tests/point_900m/hits.csv \
  --output run_logs/official_tests/point_900m/spot.png \
  --max-bins 520 --dpi 350
```

检查重点：有限距离点源光斑应明显不同于无限远平行光。

## 3. 平行光 + 3D 遮挡白板测试

cfg:

```text
configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg
```

测试内容：在理想平行光 PSF 基础上加入 3D 遮挡模型，检查支架、相机盒等结构对光子
数量和光斑形状的影响。

cfg 逐项解释：

- `telescope.config=telescope_1229_minimal.cfg`: 默认竖直望远镜。
- `mirror.config=../mirrors/mirror_1229_imported.cfg`: 完整 54 片标准镜面。
- `source.config=../sources/parallel_1M_on_axis.cfg`: 100 万平行光。
- `output.config=../outputs/whiteboard_f8.cfg`: 8 m 白板。
- `obstruction.config=../obstructions/raytrace_final_structure.cfg`: 启用由外部模型转换得到的
  3D primitives。当前会检查入射段和反射段遮挡。
- `output.csv=run_logs/official_tests/raytrace_structure_parallel/hits.csv`: 正常物理模式下，
  被遮挡光子不会写入这个 CSV。

单独运行：

```bash
mkdir -p run_logs/official_tests/raytrace_structure_parallel
build/run_optical_sim configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/raytrace_structure_parallel/run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_spot_histogram.py \
  run_logs/official_tests/raytrace_structure_parallel/hits.csv \
  --output run_logs/official_tests/raytrace_structure_parallel/spot.png \
  --max-bins 520 --dpi 350 \
  --title "Parallel beam with 3D obstruction"
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_optical_layout_3d.py \
  --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  --show-obstruction \
  --output run_logs/official_tests/raytrace_structure_parallel/layout_3d.png \
  --dpi 350
python3 python/plot_optical_layout_html.py \
  --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  --output run_logs/official_tests/raytrace_structure_parallel/layout_3d.html
```

检查重点：log 中应有 `blocked_by_obstruction`、
`output_transmission_after_obstruction` 和遮挡前后等效收集面积。

## 4. 30 m 点光源 + 3D 遮挡白板测试

cfg:

```text
configs/official_tests/point_30m_structure_whiteboard.cfg
```

测试内容：有限距离点源下遮挡结构的投影会变化，这个测试用来检查非平行光情况下的
遮挡光斑和镜面命中分布。

cfg 逐项解释：

- `source.config=../sources/point_30m_from_whiteboard_on_axis.cfg`: 点源距离标准白板
  30 m。因为白板在 `z=-8 m`，所以点源设为本地 `z=22 m`。
- `output.config=../outputs/whiteboard_f8.cfg`: 仍然使用标准 8 m 白板，不自动移动到点源
  最佳焦面。
- `obstruction.config=../obstructions/raytrace_final_structure.cfg`: 启用 3D 遮挡。
- 其它 telescope/mirror 配置与标准 1229 完美镜面一致。

单独运行：

```bash
mkdir -p run_logs/official_tests/raytrace_structure_point_30m
build/run_optical_sim configs/official_tests/point_30m_structure_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/raytrace_structure_point_30m/run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_spot_histogram.py \
  run_logs/official_tests/raytrace_structure_point_30m/hits.csv \
  --output run_logs/official_tests/raytrace_structure_point_30m/spot.png \
  --max-bins 520 --dpi 350 \
  --title "30 m point source with 3D obstruction"
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_mirror_hit_map.py \
  run_logs/official_tests/raytrace_structure_point_30m/hits.csv \
  --config configs/official_tests/point_30m_structure_whiteboard.cfg \
  --require-surface --overlay-facets \
  --output run_logs/official_tests/raytrace_structure_point_30m/mirror_hits_with_facet_outlines.png \
  --dpi 350 \
  --title "30 m point source: mirror hit points with facet outlines"
```

检查重点：spot 和镜面命中图都应体现有限距离和遮挡的共同影响。

## 5. 支架形变仰角扫描

cfg:

```text
configs/official_tests/deformation_parallel_whiteboard.cfg
```

测试内容：使用仰角相关镜面 series，扫描 0 到 90 度仰角下的平行光 PSF。

cfg 逐项解释：

- `telescope.config=telescope_1229_minimal.cfg`: 基础望远镜信息。
- `mirror.config=../mirrors/mirror_1229_imported.cfg`: 理想镜面作为基准。
- `error.config=../errors/structural_deformation_1229.cfg`: 启用结构形变 series。程序会根据
  当前 `telescope.pointing_el_deg` 插值/读取对应仰角下的镜片中心和法向。
- `source.config=../sources/parallel_1M_on_axis.cfg`: 平行光；扫描脚本会覆盖光子数。
- `output.csv=run_logs/official_tests/deformation_scan/hits.csv`: 基础输出路径，扫描时每个仰角会写到各自子目录。

单独运行：

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/run_elevation_parallel_scan.py \
  --config configs/official_tests/deformation_parallel_whiteboard.cfg \
  --run-binary build/run_optical_sim \
  --elevations 0,10,20,30,40,50,60,70,80,90 \
  --n-bunches 100000 \
  --output-dir run_logs/official_tests/deformation_scan
```

检查重点：每个仰角都有独立 spot 图和指标，RMS/质心应随仰角发生变化。

## 6. 光收集器角响应测试

这个测试不是 assembly cfg，而是专用 C++ 程序：

```text
apps/scan_light_collector_angular_response.cpp
```

测试内容：直接把光子打到单个光收集器入口，扫描 0 到 90 度入射角，检查几何接受率、
考虑材料权重后的接受率、多次反射次数。

固定模型：

- 单个 square pixel，入口边长 `2.44 cm`。
- Bezier square-cone 光收集器。
- 材料为 `true_reflect`，不是 100% 理想反射。
- SiPM 有效面边长 `1.30 cm`。
- 每个角度默认 `2000` 个光子，角度步长 `1 deg`。

单独运行：

```bash
mkdir -p run_logs/official_tests/collector_angular_response
build/scan_light_collector_angular_response \
  --photons-per-angle 2000 \
  --angle-step-deg 1 \
  --max-angle-deg 90 \
  --output run_logs/official_tests/collector_angular_response/collector_angular_response.csv \
  2>&1 | tee run_logs/official_tests/collector_angular_response/run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_collector_angular_response.py \
  run_logs/official_tests/collector_angular_response/collector_angular_response.csv \
  --output run_logs/official_tests/collector_angular_response/collector_angular_response.png \
  --dpi 350
```

输出 CSV 关键列：

- `geometric_acceptance`: 是否最终打到 SiPM 的几何比例。
- `weighted_acceptance`: 几何接受后再乘 collector 材料权重的比例。
- `mean_collector_weight`: 被接受光子的平均权重。
- `mean_reflections`: 被接受光子的平均反射次数。

## 7. CORSIKA 白板 debug 测试

cfg:

```text
configs/official_tests/corsika_whiteboard.cfg
```

测试内容：从 EventIO 读取 CORSIKA 光子，只做镜面光追到白板，不接相机。用于在相机
像素化前检查光子是否读入和坐标转换是否正常。

cfg 逐项解释：

- `telescope.pointing_az_deg=0`, `telescope.pointing_el_deg=70`: 望远镜指向，和 CORSIKA
  文件中的天顶角 20 度示例相匹配。
- `source.mode=EventIO`: 读取 EventIO 光子。
- `source.eventio_path=`: 留空，运行命令传入 zst 文件。
- `source.event_id_mode=event_array100`: 输出事件号为 `shower_event * 100 + array_id`。
- `source.eventio_coordinate_frame=corsika_iact`: 使用 CORSIKA IACT 坐标约定。
- `output.format=csv`: 输出白板光子 CSV。
- `output.hits_csv`: 白板命中光子。
- `output.summary_csv`: 每个 event/telescope 的简要统计。

单独运行：

```bash
build/run_corsika_trace configs/official_tests/corsika_whiteboard.cfg /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/whiteboard_run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_corsika_trace_output.py \
  run_logs/official_tests/corsika/whiteboard_hits.csv \
  --summary-csv run_logs/official_tests/corsika/whiteboard_summary.csv \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --output-dir "run_logs/official_tests/corsika/plots/event_${LACT_SELECTED_EVENT_ID}/whiteboard"
```

检查重点：每个 event/telescope 的 hit 数量、芯位、能量和方向应在 log 中可检查。

## 8. CORSIKA 完美相机 dense HDF5 测试

cfg:

```text
configs/official_tests/corsika_new_camera.cfg
```

测试内容：CORSIKA 光子经过理想镜面、真实 new_camera 像素和光收集器，输出 dense HDF5
相机图像。不加 NSB，不加 trigger，不加效率曲线和误差。

cfg 逐项解释：

- `telescope.pointing_az_deg=0`, `telescope.pointing_el_deg=70`: 望远镜指向。
- `mirror.config=../mirrors/mirror_1229_imported.cfg`: 标准 54 片理想镜面。
- `output.config=../outputs/focal_plane_f8.cfg`: 相机焦平面在 `z=-8 m`，法向指向镜面。
- `camera.config=../cameras/new_camera.cfg`: 使用真实相机像素表和 Bezier 光收集器。
- `sipm.config=../sipm/ideal_sipm.cfg`: SiPM 只提供有效面尺寸，PDE 关闭。
- `atmosphere.config=../atmosphere/ideal.cfg`: 不额外加大气透过率。
- `nsb.config=../nsb/ideal.cfg`: 背景光关闭。
- `trigger.config=../trigger/disabled.cfg`: trigger 关闭。
- `output.format=hdf5`: 输出 HDF5。
- `output.hdf5_path=run_logs/official_tests/corsika/camera_dense.h5`: 输出文件。
- `output.hdf5_storage=dense`: 每个 image 保存完整 1616 像素。
- `output.hdf5_write_components=false`: 只写最终 `pe/signal/photon_count`，不额外写分量。

单独运行：

```bash
build/run_corsika_trace configs/official_tests/corsika_new_camera.cfg /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" --quantity pe \
  --output "run_logs/official_tests/corsika/plots/event_${LACT_SELECTED_EVENT_ID}/camera/all_tel_pe"
```

检查重点：输出应包含所选 event 下可画的全部 telescope dense 相机图像，且没有
NSB/trigger 筛选。

## 9. CORSIKA + NSB + trigger 测试

cfg:

```text
configs/official_tests/corsika_nsb_trigger_camera.cfg
```

测试内容：在完美相机链路基础上加入常数 NSB 和简单 multiplicity trigger，检查 HDF5
中的 Cherenkov/NSB/final p.e. 分量和 trigger 表。

cfg 逐项解释：

- `nsb.config=../nsb/example_constant_rate.cfg`: 启用常数 NSB。
  当前参数为 `rate_pe_per_ns_per_pixel=0.05`、`window_ns=16`，所以平均
  `0.8 p.e./pixel`。
- `trigger.config=../trigger/example_simple_multiplicity.cfg`: 启用简单 trigger。
- `trigger.pixel_threshold_pe=10`: official smoke test 使用 10 p.e. 像素阈值。
- `output.hdf5_write_components=true`: 写出 `cherenkov_pe`、`nsb_pe` 和最终 `pe`。
- `output.save_only_triggered=true`: 只保存触发望远镜图像，减少输出体积。
- 其它镜面、相机、SiPM、EventIO 配置与完美相机测试相同。

单独运行：

```bash
build/run_corsika_trace configs/official_tests/corsika_nsb_trigger_camera.cfg /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_nsb_trigger_run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" --quantity pe \
  --output "run_logs/official_tests/corsika/plots/event_${LACT_SELECTED_EVENT_ID}/nsb_trigger/all_tel_final_pe"
```

检查重点：HDF5 里应有 `/images/dense/cherenkov_pe`、`/images/dense/nsb_pe` 和
`/trigger` 表。

## 10. CORSIKA + 3D 遮挡 + NSB + trigger 测试

cfg:

```text
configs/official_tests/corsika_obstruction_nsb_trigger_camera.cfg
```

测试内容：这是当前最接近完整链路的 smoke test：EventIO 光子、镜面光追、3D 遮挡、
真实相机、Bezier 光收集器、SiPM p.e.、NSB、trigger、dense HDF5。

cfg 逐项解释：

- `obstruction.config=../obstructions/raytrace_final_structure.cfg`: 启用 3D 遮挡模型。
  正常模式下，被遮挡光子在进入相机前被丢弃。
- `nsb.config=../nsb/example_constant_rate.cfg`: 加常数 NSB。
- `trigger.config=../trigger/example_simple_multiplicity.cfg`: 启用 trigger。
- `trigger.pixel_threshold_pe=10`: 像素阈值 10 p.e.。
- `output.hdf5_path=run_logs/official_tests/corsika/camera_obstruction_nsb_trigger_dense.h5`:
  独立输出文件，避免覆盖不带遮挡的 NSB+trigger 测试。
- `output.hdf5_write_components=true`: 写 Cherenkov/NSB/final 分量，方便比较遮挡前后。
- `output.save_only_triggered=true`: 只保存触发图像。

单独运行：

```bash
build/run_corsika_trace configs/official_tests/corsika_obstruction_nsb_trigger_camera.cfg /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_obstruction_nsb_trigger_run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_obstruction_nsb_trigger_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" --quantity pe \
  --output "run_logs/official_tests/corsika/plots/event_${LACT_SELECTED_EVENT_ID}/obstruction_nsb_trigger/all_tel_final_pe"
```

检查重点：log 中应输出遮挡统计；HDF5 里应有 final p.e. 图和 trigger 表。

## 11. CORSIKA full-response smoke test

cfg:

```text
configs/official_tests/corsika_full_response_camera.cfg
```

测试内容：加入目前已实现的“非理想响应”：支架形变、小随机误差、镜面反射率、滤光片
透过率、SiPM PDE 和 trigger。NSB 在这个测试中关闭，用于单独检查光学/效率响应。

cfg 逐项解释：

- `telescope.pointing_el_deg=70`: 望远镜仰角 70 度，结构形变会按这个仰角取值。
- `error.config=../errors/full_response_1229.cfg`: 启用仰角相关结构形变和小随机误差。
- `efficiency.config=../efficiency/curves_all.cfg`: 启用镜面反射率和滤光片透过率曲线。
- `sipm.config=../sipm/new_camera_sipm.cfg`: 启用 SiPM PDE 曲线。
- `nsb.config=../nsb/ideal.cfg`: NSB 关闭。
- `trigger.config=../trigger/example_simple_multiplicity.cfg`: trigger 开启。
- `trigger.pixel_threshold_pe=10`: official 阈值设为 10 p.e.。
- `output.hdf5_write_components=false`: 这个测试只关心最终 p.e.，不写 NSB 分量。
- `output.save_only_triggered=true`: 只保存触发图像。

单独运行：

```bash
build/run_corsika_trace configs/official_tests/corsika_full_response_camera.cfg /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_full_response_run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_full_response_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" --quantity pe \
  --output "run_logs/official_tests/corsika/plots/event_${LACT_SELECTED_EVENT_ID}/full_response/all_tel_pe"
```

检查重点：相对完美相机测试，图像 p.e. 会受到形变、误差和效率曲线影响。一键脚本还会
额外生成芯位/阵列分布图：

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_array_layout.py \
  run_logs/official_tests/corsika/camera_full_response_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" --quantity pe --log-color \
  --output "run_logs/official_tests/corsika/plots/event_${LACT_SELECTED_EVENT_ID}/layout/core_and_array_pe.png"
```

## 默认不在一键脚本里跑的诊断 cfg

这些 cfg 是调试用，不放进默认 official 脚本，避免测试时间和输出文件过大。

### 中间两圈镜片平行光

cfg:

```text
configs/official_tests/inner_two_rings_parallel_whiteboard.cfg
```

配置说明：

- `mirror.config=../mirrors/mirror_1229_inner_two_rings.cfg`: 只使用中心第一圈 6 片和第二圈
  12 片，共 18 片。
- 不启用遮挡、不启用误差。

### 最外圈镜片遮挡 mark-only 诊断

cfg:

```text
configs/official_tests/outer_ring_parallel_structure_mark_only_whiteboard.cfg
```

配置说明：

- `mirror.config=../mirrors/mirror_1229_outer_ring.cfg`: 只使用最外圈 18 片。
- `obstruction.mark_only=true`: 被遮挡光子只打标记，仍继续传播到白板。
- `mirror.mode` 和 `mirror.csv_path` 在这里显式写出，是为了让这个单独诊断 cfg 即使脱离
  mirror 组件文件也能清楚指向最外圈 CSV。

### 全镜面遮挡 mark-only 诊断

cfg:

```text
configs/official_tests/perfect_parallel_structure_mark_only_whiteboard.cfg
```

配置说明：

- 完整 54 片镜面。
- `obstruction.mark_only=true`: 用于定位具体遮挡来源，比如区分 incoming 和 reflected
  段遮挡。

mark-only 输出 CSV 会多出：

```text
obstruction_blocked
obstruction_blocked_incoming
obstruction_blocked_reflected
```

画被遮挡光子的镜面分布：

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_obstruction_marked_hits.py \
  run_logs/official_tests/raytrace_structure_parallel_mark_obstructed/hits.csv \
  --space mirror \
  --config configs/official_tests/perfect_parallel_structure_mark_only_whiteboard.cfg \
  --overlay-facets \
  --output run_logs/official_tests/raytrace_structure_parallel_mark_obstructed/mirror_hits_marked_obstructed.png
```
