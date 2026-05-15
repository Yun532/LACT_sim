# LACT_sim 程序中文总览

这份文档用于两个目的：

1. 帮助用户理解 LACT_sim 当前能做什么、怎么配置、输出代表什么。
2. 帮助新的开发对话快速接手程序结构，避免重新梳理坐标系、事件编号、HDF5、PDE 和 trigger 等约定。

本文档不是逐行 API 手册，而是模块级说明。若源码有更新，应同步更新这里的流程和关键假设。

## 1. 当前程序定位

LACT_sim 当前是一个面向 LACT 望远镜的光学和相机 integrated p.e. 模拟程序。主流程是：

```text
人工光源 / Photon CSV / CORSIKA EventIO
  -> 望远镜坐标变换
  -> 镜面命中和球面反射
  -> 可选支架形变、随机误差、3D 遮挡
  -> 白板或相机焦平面
  -> 可选 new_camera 像素映射和光收集器
  -> mirror/filter/atmosphere/SiPM PDE 权重
  -> 可选 NSB
  -> 可选 simple trigger
  -> CSV 或 HDF5 输出
```

当前电子学被简化为 integrated p.e. 图像：

```text
p.e. contribution =
  photon multiplicity
  * photon weight
  * mirror/filter/atmosphere efficiency
  * light collector acceptance/reflection weight
  * sipm.pde
```

也就是说，现在没有 waveform、采样、增益涨落、串扰、afterpulse、dark count 和真实 trigger board。它已经适合做光学验证、CORSIKA 到相机 p.e. 图像、NSB/trigger 流程打通，但还不是完整电子学模拟。

## 2. 两个主程序

### run_optical_sim

用途：人工光源、点光源、Photon CSV、白板测试、镜片布局和遮挡 benchmark。

典型命令：

```bash
build/run_optical_sim configs/official_tests/perfect_parallel_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/perfect_parallel/run.log
```

常见输出：

```text
hits.csv                 白板/输出平面逐光子命中
camera_pixel_image.csv   可选的像素聚合 CSV
run.log                  生效配置、统计和运行时间
```

### run_corsika_trace

用途：读取一个 CORSIKA/EventIO 文件，处理其中所有 shower event 和望远镜 stream，输出白板 CSV 或相机 HDF5。

典型命令：

```bash
build/run_corsika_trace configs/official_tests/corsika_new_camera.cfg input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_run.log
```

推荐正式输出是 dense HDF5：

```text
camera_dense.h5
camera_nsb_trigger_dense.h5
camera_full_response_dense.h5
```

白板调试仍推荐 CSV，因为它保存逐光子命中点，便于看 PSF 和 spot 结构：

```text
whiteboard_hits.csv
whiteboard_summary.csv
```

## 3. 配置系统

主 cfg 可以通过 `*.config=...` 包含子配置。常见目录是：

```text
configs/templates/        给新用户的入口模板
configs/official_tests/   官方 benchmark 配置
configs/mirrors/          镜片布局
configs/sources/          人工光源
configs/outputs/          白板/相机平面
configs/cameras/          new_camera 像素
configs/sipm/             SiPM 尺寸和 PDE
configs/efficiency/       mirror/filter/atmosphere 曲线
configs/errors/           误差和支架形变
configs/obstructions/     3D 遮挡 primitives
configs/nsb/              夜空背景
configs/trigger/          simple trigger
```

配置读取和展开由 `expandConfig()` 完成。每个子配置会被加上自己的前缀，例如 `mirror.config` 中的 `csv_path` 会展开为 `mirror.csv_path`。

推荐新用户从下面两个模板开始：

```text
configs/templates/perfect_whiteboard.cfg
configs/templates/minimal_corsika_camera.cfg
```

默认策略是：未写入的模块尽量保持理想或关闭。用户通常只需要修改输入文件、望远镜指向、镜片/相机结构和输出路径。

## 4. 坐标系总览

坐标系详细说明见 `docs/coordinate_systems.md`。这里保留最容易出错的结论。

### CORSIKA/EventIO 阵列坐标

CORSIKA IACT 水平坐标采用 NWU：

```text
+x : magnetic North
+y : West
+z : Up
```

HDF5 中显式保存为：

```text
array_x_north_m
array_y_west_m
array_z_up_m
core_x_north_m
core_y_west_m
```

绘制常规地图时才转换到 East-North 显示：

```text
plot_east_m  = -array_y_west_m
plot_north_m =  array_x_north_m
```

### 望远镜指向

CORSIKA 模式下：

```text
pointing_az_deg = 0    指向 magnetic North / +x
pointing_az_deg = 90   指向 East / -y
pointing_el_deg        地平高度角
zenith_deg             90 - pointing_el_deg
```

例如 CORSIKA 文件是 `zenith=20, azimuth=0`，则通常配置：

```ini
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70
```

### 望远镜本地坐标

光追在望远镜本地坐标中进行：

```text
local +z : 从镜面指向天空的光轴
local -z : 轴上入射光方向
local x/y: 相机和镜面横向坐标
```

`buildTelescopeFrame()` 根据望远镜位置和 az/el 生成全局到本地的旋转关系。`applyTelescopeFrame()` 用于把镜片、输出平面或光子变换到对应框架。

### 相机/白板平面

白板和相机都由 output plane 定义：

```ini
output.plane_point=0,0,-8
output.plane_normal=0,0,-1
output.plane_u_axis=1,0,0
output.plane_v_axis=0,1,0
```

`normal` 指向镜片侧。光线求交是双面的，所以 normal 正负不决定是否命中。图像坐标由 `u_axis` 和 `v_axis` 决定。

## 5. 事件编号

正式下游分析只使用一个编号：

```text
event_id
```

对于 CORSIKA/EventIO，默认：

```text
event_id = shower_event * 100 + array_id
```

例如：

```text
shower_event = 468898
array_id = 2
event_id = 46889802
```

这里的 `array_id` 是 CORSIKA 的 array/core offset stream，不是望远镜编号。望远镜编号是 `telescope_id`。

## 6. 核心数据结构和函数职责

下面列的是最重要的函数职责。它们不是所有函数，但足够帮助新开发者理解主线。

### 配置和组件构建

`readKeyValueConfig()`  
读取简单 key/value cfg。

`expandConfig()`  
展开 `telescope.config`、`mirror.config`、`source.config` 等子配置，并给子配置 key 加模块前缀。

`buildTelescopeConfig()`  
读取望远镜位置、az/el 指向、焦距和坐标系名称。

`buildTelescopeFrame()`  
根据望远镜配置生成本地坐标系。CORSIKA 模式下要特别注意 azimuth 是 North to East。

`buildFacetsFromConfig()`  
构造镜片列表。支持生成布局、导入 CSV、半边镜片、内外圈镜片和 elevation series。

`applyStructuralDeformation()`  
根据仰角相关的 series mirror 文件插值，得到支架形变后的镜片中心和朝向。

`applyFacetErrors()`  
添加随机误差，包括径向位置涨落、法向角误差、反射方向误差、曲率半径涨落和反射率 scale 涨落。

`buildOutputPlane()`  
读取白板/相机平面的点、normal、u/v 轴。

`buildCameraGeometry()`  
读取 new_camera 像素 CSV 或生成简单相机。

`buildLightCollector()`  
根据相机配置创建 Bezier square-cone 光收集器。入口边长来自 pixel size，出口尺寸来自 `sipm.size_m` 和 collector 配置。

`buildEfficiencyConfig()`  
读取 mirror reflectivity、filter transmission、atmosphere transmission、funnel acceptance 等效率配置。

`buildNsbConfig()` 和 `buildTriggerConfig()`  
读取 constant-rate NSB 和 simple multiplicity trigger。

### 光追核心

`OpticalTracer::traceToPlane()`  
单个光子的主光追函数。输入 photon、镜片 layout、输出平面和效率模型。输出 `OpticalSurfaceHit`。

`OpticalTracer::intersectMirror()`  
对每片镜子求交，选择最近的有效交点。

`intersectSphericalFacet()`  
计算光线和球面镜片的交点，并做孔径检查。

`insideAperture()`  
把交点转到单片镜子的局部孔径坐标中，判断是否落在六边形或矩形镜片内。六边形朝向由 `aperture_rotation_rad` 控制。

`OpticalTracer::intersectOutputPlane()`  
计算反射光线和白板/相机平面的交点，并给出 u/v 坐标。

`OpticalEfficiency::total()`  
按波长和入射角计算总效率因子。当前效率是乘法模型。

### 相机响应

`findContainingPixel()`  
根据焦平面 u/v 坐标找到包含该命中点的 camera pixel。

`traceLightCollector()`  
把命中 pixel 入口的光子送入 square-cone 光收集器，追踪多次反射。返回是否到达 SiPM、剩余 intensity 和反射次数。

`applyCameraResponse()`  
完成 pixel 查找、collector、SiPM PDE 和电子学 placeholder 处理，更新 `OpticalSurfaceHit` 中的 `pixel_id`、`relative_efficiency`、collector 信息和 accepted 状态。

`accumulatePixelHit()`  
把单个命中累计到 `(event_id, telescope_id, pixel_id)` 上，生成 photon_count、signal、pe 和时间统计。

### CORSIKA/EventIO

`readEventIOMetadata()`  
读取 EventIO 元数据，包括 input card、望远镜位置、shower header、array offsets。

`streamEventIOPhotonBunches()`  
顺序流式读取 photon bunch，不再一次性把全文件 photon preload 到内存。每个 bunch 会回调到 `run_corsika_trace` 的处理函数。

`transformEventIOBunchToTraceFrame()`  
把 EventIO photon bunch 转成光追使用的 photon。CORSIKA photon bunch 的 x/y/z 是相对望远镜的坐标，方向来自 EventIO direction cosines。

### 遮挡

`buildObstructionMask()`  
读取 2D mask 或 3D primitives。正式 3D 遮挡建议使用 primitives。

`incomingSegmentBlockedByObstruction()`  
检查光子从入射位置到镜面命中点之间是否被遮挡。

`segmentBlockedByObstruction()`  
检查镜面反射点到白板/相机平面之间是否被遮挡。

`segmentIntersectsCylinder()`、`segmentIntersectsAabb()`、`segmentIntersectsRegularPrism()`、`segmentIntersectsBoxWithHole()`  
分别处理支撑杆、盒子、正多边形柱和带孔连接块的几何相交。

### 输出

`writePixelCsv()`  
写像素聚合 CSV，主要用于调试。

`writeNativeTraceHdf5()`  
写正式 HDF5 输出，包括 static metadata、events、camera pixels、mirror facets、dense/sparse images、NSB components 和 trigger tables。

`writeSummaryCsv()`  
写 CORSIKA 白板或 CSV 调试的 event/telescope 摘要。

## 7. 输出文件含义

### whiteboard_hits.csv

逐光子白板命中。适合画 spot、PSF、镜面命中检查。文件大，但最直观。

### whiteboard_summary.csv

event/telescope 级摘要。绘图时用于按第几个 shower event 或 event_id 选择数据。

### camera_dense.h5

普通相机输出。`/images/dense/pe` 是每个 event/telescope 的 1616 像素 integrated p.e. 图像。

### camera_nsb_trigger_dense.h5

相机输出加 constant-rate NSB 和 simple trigger。若启用 `output.hdf5_write_components=true`，会额外保存：

```text
/images/dense/cherenkov_pe
/images/dense/nsb_pe
```

### camera_full_response_dense.h5

当前最接近完整响应的官方测试输出，包含支架形变、误差、效率曲线、SiPM PDE、光收集器和 trigger。它仍然不是 waveform 电子学。

## 8. HDF5 关键表

推荐分析 dense HDF5，而不是大量 CSV。

```text
/events/table
  event_index, event_id

/events/corsika
  event_id, shower_event_id, array_id, energy_gev, theta_deg, phi_deg,
  core_x_north_m, core_y_west_m, array_rotation_deg

/telescopes/table
  telescope_id, array_x_north_m, array_y_west_m, array_z_up_m,
  pointing_az_deg, pointing_el_deg, focal_length_m

/camera/pixels
  pixel_id, x_m, y_m, size_m, shape_code

/images/index
  image_index, event_id, telescope_id, total_photons, total_pe,
  total_signal, time_mean_ns, time_rms_ns

/images/dense/pe
  [n_images, n_pixels] final integrated p.e.

/images/dense/photon_count
  [n_images, n_pixels] collected photon/bunch count

/trigger/telescope
  event_id, telescope_id, triggered, n_pixels_above_threshold,
  total_pe, trigger_time_ns

/trigger/array
  event_id, array_triggered, n_triggered_telescopes
```

相机图像坐标直接来自 `/camera/pixels`，所以 HDF5 画图不需要额外 pixel CSV。

## 9. 当前物理边界

已经实现并可用于 benchmark：

```text
镜面 raytrace
标准/导入/变形镜片
人工光源和 CORSIKA EventIO
白板逐光子输出
new_camera 1616 像素相机
Bezier square-cone 光收集器
mirror/filter/atmosphere/SiPM PDE 权重
constant-rate NSB
simple multiplicity trigger
3D primitives 遮挡
HDF5 dense 相机输出
```

仍是简化或 placeholder：

```text
真实 waveform 电子学
SiPM 饱和、串扰、afterpulse、dark count
真实 trigger board 和邻接拓扑
真实 NSB sky model 或测量数据
大气传播重新模拟
大规模事件级并行
HDF5 chunk/compression 优化
```

## 10. 官方测试入口

完整官方测试脚本：

```bash
tools/run_official_tests.sh --corsika-file /path/to/input.zst
```

不跑 CORSIKA，仅跑光学和绘图：

```bash
tools/run_official_tests.sh --no-corsika
```

主要测试包括：

```text
完美平行光白板
900 m 点源白板
3D 遮挡平行光和 30 m 点源
支架形变 0 到 90 度扫描
半边/内外圈镜片测试
光收集器角响应
效率曲线验证
CORSIKA 白板
CORSIKA camera dense HDF5
CORSIKA + NSB + trigger
CORSIKA + obstruction + NSB + trigger
CORSIKA full response
```

每个测试的配置和命令见 `docs/official_tests.md`。

## 11. 常用绘图脚本

```text
python/plot_spot_histogram.py
  画白板光斑。

python/plot_mirror_hit_map.py
  画镜面命中点，可叠加镜片轮廓。

python/plot_optical_layout_3d.py
  画镜片、白板/相机、遮挡结构的静态 3D 图。

python/plot_optical_layout_html.py
  生成可交互 HTML 3D 检查图。

python/plot_hdf5_camera.py
  从 HDF5 画相机像素图。可画单台望远镜，也可省略 telescope_id 画所有望远镜。

python/plot_hdf5_array_layout.py
  画望远镜阵列分布、core 位置、arrival direction 和 telescope total p.e.。

python/select_hdf5_event.py
  从多个 HDF5 中选择一个共同可画的 event_id，官方测试用它保证不同输出画同一个事件。

python/plot_efficiency_curves.py
  验证 mirror/filter/atmosphere/SiPM PDE 曲线和总效率。
```

## 12. 新对话快速阅读顺序

如果以后开启新对话，让模型接手这个项目，建议按这个顺序读：

```text
1. docs/program_overview_zh.md
2. docs/coordinate_systems.md
3. docs/hdf5_output_format.md
4. docs/module_interfaces.md
5. docs/official_tests.md
6. include/app/OpticalSimCommon.hpp
7. apps/run_corsika_trace.cpp
8. apps/run_optical_sim.cpp
```

如果任务只涉及 CORSIKA 输出和画图，重点读：

```text
docs/corsika_eventio_adapter.md
docs/hdf5_output_format.md
python/plot_hdf5_camera.py
python/plot_hdf5_array_layout.py
```

如果任务涉及遮挡，重点读：

```text
docs/module_interfaces.md
configs/obstructions/raytrace_final_structure.cfg
configs/obstructions/raytrace_final_structure_primitives.csv
python/plot_optical_layout_html.py
```

## 13. 开发注意事项

1. 不要把 CORSIKA 的 `array_id` 当成 telescope id。
2. 不要把 CORSIKA NWU 坐标直接当 ENU 坐标用。
3. 相机图像方向由 output plane 的 u/v 轴决定，不由 normal 单独决定。
4. `sipm.pde` 是当前唯一推荐的 photon-to-p.e. 波长转换入口。旧的 `electronics.pe_conversion` 和 `efficiency.sipm_pde` 只作为兼容别名。
5. `output.save_only_triggered=true` 时，HDF5 `/images` 只保存触发望远镜，但 `/trigger` 表仍保存完整触发判断。
6. 白板 CSV 是逐光子调试文件，不适合大规模生产保存。
7. dense HDF5 是当前推荐的正式相机输出格式。
8. 物理验证优先用 official tests，不要只看单个事件图像下结论。

