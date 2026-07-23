# LACT_sim 坐标系说明

本文说明 CORSIKA/EventIO、用户 Photon CSV、望远镜本地光追、阵列地图和
pyLAST 相机图分别使用什么坐标。英文参考版及更多输出字段说明见
[coordinate_systems.md](coordinate_systems.md)。

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
e_x = (-sin(E) cos(A),  sin(E) sin(A), cos(E))
e_y = (-sin(A),        -cos(A),        0     )
e_z = ( cos(E) cos(A), -cos(E) sin(A), sin(E))
```

三条轴的含义是：

```text
local +x : 高度角增加方向，朝天顶
local +y : 方位角从北向东增加方向
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
pyLAST pix_x = -LACT_sim camera_x
pyLAST pix_y = -LACT_sim camera_y
```

这只是完整相机图的显示约定，不是 CORSIKA NWU 输入转换，也不改变
`pixel_id` 和每个像素的 PE。

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
