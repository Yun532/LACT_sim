# LACT_sim 坐标系说明

本文说明 CORSIKA/EventIO、用户 Photon CSV、望远镜本地光追、阵列地图和
pyLAST 相机图分别使用什么坐标。英文参考版及更多输出字段说明见
[coordinate_systems.md](coordinate_systems.md)。

建议先打开[真实 LACT 三维坐标模型](assets/lact-coordinate-system-3d.html)：拖动可旋转、
滚轮可缩放；它使用正式 `mirror_1229` 镜面、遮挡结构和 1616 个真实相机像素，
同时标出输入 NWU、望远镜局部 `x/y/z`、天空/入射光方向、相机 `u/v` 及
pyLAST 符号关系。页面中的地面和底座仅用于辨认水平面与天空方向，不属于工程几何。

[交互式坐标系总图](assets/coordinate-system-explorer.html)则适合切换
`source.coordinate_frame`、望远镜方位/高度和天空偏移，同时对照输入 `x/y`、
局部光学坐标、相机 `u/v`、ROOT/HDF5 字段和 pyLAST 的最终画图轴。

三维 HTML 由仓库中的 Python 代码按实际 cfg 生成：

```bash
python python/plot_optical_layout_html.py \
  --config configs/official_tests/corsika_full_response_camera.cfg \
  --output docs/assets/lact-coordinate-system-3d.html \
  --coordinate-frame-mode source \
  --show-coordinate-axes \
  --show-camera-pixels \
  --show-ground \
  --trace-csv docs/assets/data/corsika-north-example-rays.csv \
  --trace-provenance docs/assets/data/corsika-north-example-provenance.json
```

## main 坐标审计与真实 CORSIKA 示例

这次按实际调用链检查了 EventIO 读取、输入适配、镜面/遮挡、本地光追、输出平面、
相机像素、CSV/HDF5/ROOT 字段、阵列画图和 pyLAST 边界。主数据流是：

```text
EventIO NWU photon bunch
  -> transformBunchToTelescopeLocal()
  -> telescope-local mirror + obstruction + output plane
  -> hit.u/v = dot(surface - plane_point, plane_u/v)
  -> camera_x/y = u/v
  -> ROOT: x_m/y_m = u/v
  -> pyLAST: pix_x=-v, pix_y=+u
```

审计中修正了三处会造成误读或画图不一致的地方：

1. `run_corsika_trace` 日志原来打印的是通用光学布局 frame，容易让人误以为它就是
   CORSIKA 输入适配器；现在直接打印真实 `source.coordinate_frame` 的本地轴在输入
   坐标中的方向，并明确光追几何始终在 telescope-local 中。
2. 白板 CSV 原来只保存输入位置，`dir_x/y/z` 实际是镜面反射后的方向；现在新增
   `input_dir_x/y/z`，可以无歧义地重建“输入锚点/方向 → 镜面命中 → 输出命中”。
3. 几个带 `--sky-up` 的 Python 图原来总使用通用布局 frame；现在统一通过
   `config_io.py` 选择真实 source adapter，CORSIKA 图与 HDF5 相机图使用同一套 NWU
   基底。静态坐标检查脚本也默认采用 source adapter。

三维页面中的浅黄色光线不是示意线，而是当前 C++ 程序实际运行结果。可复现实例使用：

```ini
config=configs/examples/corsika_coordinate_north_example.cfg
telescope.pointing_az_deg=0       # 指向磁北
telescope.pointing_el_deg=70
source.coordinate_frame=corsika_nwu_relative
```

输入文件为 `muon_E100_th0_run000001.zst`（SHA-256 和命令见
[`corsika-north-example-provenance.json`](assets/data/corsika-north-example-provenance.json)）。
实际 CORSIKA shower 1 的到达方位角为 `300.027133 deg`（北向东）、高度角为
`88.282787 deg`；这是事例真值，不会覆盖望远镜的 `az=0, el=70` 指向。

本次运行共保存 3385 个白板输出命中。页面从 `event_id=110, telescope_id=0` 的
1557 个输出命中中等距选取 64 条完整光路；该流有 520 个输入 bunch、
2474.320005 个加权输入光子、1885 个遮挡前镜面命中、328 个遮挡损失、1557 个
最终输出命中。逐光子原始字段见
[`corsika-north-example-rays.csv`](assets/data/corsika-north-example-rays.csv)，对应汇总见
[`corsika-north-example-summary.csv`](assets/data/corsika-north-example-summary.csv)。

这里把真实白板 `u/v` 命中叠加在真实 1616 像素几何上；白板运行没有启用相机响应，
所以它展示的是“光学输出落在相机平面的哪里”，不是 PE、触发或电子学输出。

## 先区分五层坐标

| 层次 | `x` 方向 | `y` 方向 | 说明 |
|---|---|---|---|
| CORSIKA/EventIO 原始阵列坐标 | 磁北 North | 西 West | NWU；CORSIKA 默认输入 |
| `enu_east_*` Photon CSV | 东 East | 北 North | ENU；用户显式选择 |
| `telescope_local` Photon CSV | 望远镜本地横轴 | 望远镜本地纵轴 | 纯光学输入默认值 |
| 阵列地图显示 | 东 East | 北 North | 只在画图时由 NWU 转换 |
| 焦平面/相机原始坐标 | output plane 的 `u` | output plane 的 `v` | 通常等于 telescope-local `x/y` |

最容易混淆的是：阵列地图横轴虽然是 East，但 CORSIKA 文件里存储的
`x` 仍然是 North。转换关系为：

```text
plot_x_East  = -CORSIKA_y_West
plot_y_North =  CORSIKA_x_North
```

`telescope.coordinate_system=array` 当前主要是望远镜配置中的说明性元数据。
输入行采用哪套坐标，由 `source.coordinate_frame` 决定。

## 默认选择规则

### CORSIKA/EventIO

`run_corsika_trace` 的正常 CORSIKA 输入默认等价于：

```ini
source.mode=EventIO
source.coordinate_frame=corsika_nwu_relative
```

旧配置名 `source.eventio_coordinate_frame=corsika_iact` 也会归一化为
`corsika_nwu_relative`。

### 用户 Photon CSV 和纯光学测试

Photon CSV 没有写 `source.coordinate_frame` 时，默认是：

```ini
source.coordinate_frame=telescope_local
```

手写平行光、单色光和简单光学测试推荐使用该默认值。从 CORSIKA photon
bunch 裁剪出来的 CSV 应显式使用：

```ini
source.coordinate_frame=corsika_nwu_relative
```

不要先把同一行转到本地坐标，又让程序按 `corsika_nwu_relative` 转一次。

## CORSIKA/EventIO 默认流程

### 1. 原始水平坐标

CORSIKA IACT 使用磁北-西-上坐标，即 NWU：

```text
+x : magnetic North
+y : West
+z : Up
```

普通 EventIO 2D photon bunch 给出：

```text
x, y   : 相对当前望远镜中心的光子位置 [cm]
cx, cy : 光子传播方向余弦
```

程序把位置从厘米换成米，并补出向下传播的分量：

```text
dir_z = -sqrt(1 - dir_x^2 - dir_y^2)
```

因此 `dir_x/y/z` 表示光子正在传播的方向，不是光源来向。下降光子通常有
`dir_z < 0`。2D bunch 的 `z=0` 是输入参考平面，不是发射高度或海拔。

EventIO bunch 保存相对各台望远镜的 `x、y`，但没有逐 bunch 的 `z`。把这个锚点
旋转到望远镜局部坐标后，程序用原有的单值配置将它平移到光学模型的局部原点：

```ini
source.eventio_reference_z_m=-16
```

当前默认值为 `-16 m`。这是当前导入 LACT 光学模型的坐标映射：镜面顶点约在
`z=-16 m`，相机约在 `z=-8 m`；而 sim_telarray 对应光学坐标中镜面顶点约为
`z=0`、相机约为 `z=+8 m`。它不表示光子在镜面下方 16 m 产生，也不改变
CORSIKA 的全局观测高度。

仍可显式设为 `0` 做坐标诊断，但对于当前导入模型它不是等价的原点。只修改这个
值会保持 EventIO 的 `x、y`、方向和时间不变，从而在局部光学模型中定义另一条
平行光线，所以相机结果本来就会不同。

默认的 `source.eventio_2d_plane_mode=auto` 对所有二维 EventIO bunch 都按完整
入射直线求交，交点参数 `t` 可以为正或负；随后仍用光子传播方向检查是否从镜片
正面入射。只有当输入位置确实是物理上游起点而非记录面锚点时，才应显式使用
`source.eventio_2d_plane_mode=forward`。

CORSIKA 中各望远镜的位置（包括不同安装高度）仍是阵列/全局元数据，用于簇射
波前和时间。由于 bunch 的 `x、y` 已经相对各自望远镜给出，使用同一套光学模型
的望远镜都采用同一个局部 `-16 m` 平移；不要再把望远镜安装高度加到这个值上。

从 CORSIKA 裁剪成六列 Photon CSV 时保留原始 `z_m=0`，并在 cfg 设置
`source.eventio_2d=true`，不要机械地把 NWU 的 z 列改为 `-16`。

### 2. 望远镜指向

CORSIKA 模式的方位角从磁北开始向东增加：

```text
az =   0 deg : North，NWU +x
az =  90 deg : East， NWU -y
az = 180 deg : South，NWU -x
az = 270 deg : West， NWU +y

el            : 地平高度角
zenith        : 90 deg - el
```

例如 `zenith=20 deg, azimuth=0 deg` 对应：

```ini
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70
```

光追使用配置文件中的 `telescope.pointing_az_deg/el_deg`。EventIO 事件方向
会写入输出元数据，但不会隐式覆盖配置中的望远镜指向。

### 3. NWU 转望远镜本地坐标

令 `A=azimuth`、`E=elevation`。程序在 NWU 中构造：

```text
e_x = ( sin(A),         cos(A),        0     )
e_y = (-sin(E) cos(A),  sin(E) sin(A), cos(E))
e_z = ( cos(E) cos(A), -cos(E) sin(A), sin(E))
```

三条轴的含义是：

```text
local +x : 水平横向轴；A=0 时指向西，与方位角增加方向相反
local +y : 高度角增加方向，朝天顶（sky-up）
local +z : 望远镜光轴，从镜面指向天空
local -z : 正入射光子的传播方向
```

每个 photon bunch 通过点积旋转：

```text
p_local = (dot(p_NWU, e_x),
           dot(p_NWU, e_y),
           dot(p_NWU, e_z))

d_local = normalize((
           dot(d_NWU, e_x),
           dot(d_NWU, e_y),
           dot(d_NWU, e_z)))
```

`corsika_nwu_relative` 中的位置已经相对望远镜，不再减望远镜阵列位置。
只有显式选择：

```ini
source.coordinate_frame=corsika_nwu_global
```

才使用：

```text
p_local = R^T (p_NWU - telescope_position_NWU)
```

方向向量只旋转，不平移。完成这一步后，CORSIKA 光子、镜面和相机都在
望远镜本地光学坐标中追迹，不再套用第二次通用全局旋转。

### 4. 焦平面和相机

标准输出平面通常设置为：

```ini
output.plane_point=0,0,-8
output.plane_normal=0,0,-1
output.plane_u_axis=1,0,0
output.plane_v_axis=0,1,0
```

因此：

```text
camera_x = hit.u = local x
camera_y = hit.v = local y
```

这些是原始焦平面坐标。镜面反射、像素查找、光学效率和 PE 累积均使用
该坐标，显示方向不会改变像素编号或模拟结果。

## 纯光学默认的 telescope-local 坐标

`telescope_local` 输入直接按望远镜本地坐标解释：

```text
local +z : 从镜面指向天空的光轴
local -z : 从天空射向镜面的正入射方向
local +x : 焦平面/相机第一条横向轴
local +y : 与 local x/z 组成右手系的第二条横向轴
```

最简单的轴上光子可以写成：

```csv
x_m,y_m,z_m,dir_x,dir_y,dir_z
0,0,10,0,0,-1
```

在该模式下，CSV 坐标不旋转，程序只归一化方向。`run_optical_sim` 如需把
整套本地几何放进通用全局展示坐标，会使用：

```text
e_z = ( cos(E) cos(A),  cos(E) sin(A), sin(E))
e_x = (-sin(A),         cos(A),        0     )
e_y = e_z cross e_x
```

这里方位角只是从通用 `+X` 朝 `+Y` 增加；不要把这两个通用轴自动解释成
CORSIKA 的 North/West。

## 显式 ENU East-start 输入

需要“东为 `x`”的用户 Photon CSV 可以选择：

```ini
source.coordinate_frame=enu_east_relative
# 或
source.coordinate_frame=enu_east_global
```

定义为：

```text
+x : East
+y : North
+z : Up

az =   0 deg : East
az =  90 deg : North
az = 180 deg : West
az = 270 deg : South
```

`_relative` 表示位置已经相对望远镜；`_global` 表示绝对阵列位置，程序会
减一次 `telescope.position_m`。

CORSIKA NWU 转为相同物理方向的 ENU：

```text
x_ENU = -y_NWU
y_ENU =  x_NWU
z_ENU =  z_NWU

dir_x_ENU = -dir_y_NWU
dir_y_ENU =  dir_x_NWU
dir_z_ENU =  dir_z_NWU

az_ENU_east_start = 90 deg - az_NWU_north_to_east
```

角度可归一化到 `[0, 360)`。该变换不包含磁北到真北的磁偏角修正。

## 画图层的方向

### 阵列地图

阵列位置和 shower core 在文件中继续保存为 NWU；地图显示时才转换为：

```text
horizontal = East  = -West
vertical   = North
```

### 光学白板和光子数相机图

纯光学白板、原始 photon-count 相机图显示 output plane 的原始 `u/v`，
不做 pyLAST 天空视图反转。

### 完整相机和 pyLAST

LACT_sim ROOT 输出保存原始焦平面坐标。`LactEventSource` 在进入 pyLAST
相机几何时转换成 source-offset/天空视图坐标：

```text
pyLAST pix_x = -LACT_sim camera_y
pyLAST pix_y = +LACT_sim camera_x
```

这只是完整相机图的显示约定，不是 CORSIKA NWU 输入转换，也不改变
`pixel_id` 和每个像素的 PE。

当前 `pylast.visualize.plot_camera_image()` 还采用“第二个坐标画在水平轴、
第一个坐标画在垂直轴”的显示方式：

```text
Matplotlib horizontal = pyLAST pix_y = +LACT_sim u
Matplotlib vertical   = pyLAST pix_x = -LACT_sim v
```

若需要用同一个 pyLAST event 检查 LACT_sim 原始焦平面 `u/v`，运行：

```bash
python3 python/plot_photon_csv_root_pylast.py \
  run_logs/examples/photon_csv_full_camera/lact_events.root \
  --event-id 1909 --telescope-id 19 \
  --coordinate-view lact-uv \
  --output run_logs/examples/photon_csv_full_camera/camera_uv.png
```

`lact-uv` 只改变绘图坐标与轴标签，不改变 pixel 顺序或 PE 数组。

## Event 1909 示例

仓库中的 Event 1909 最简 CSV 显式配置为：

```ini
source.coordinate_frame=corsika_nwu_relative
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70
```

所以 CSV 中：

```text
x_m     : 相对望远镜的磁北方向位置
y_m     : 相对望远镜的西方向位置
z_m=0   : CORSIKA 2D 输入参考平面
dir_x   : 北向传播分量
dir_y   : 西向传播分量
dir_z   : 向上传播分量，下降光子为负
```

程序按正常 CORSIKA NWU 适配器将其转到 telescope-local；辅助脚本没有
提前转换一次，因此不会发生重复旋转。
