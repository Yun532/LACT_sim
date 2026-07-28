# LACT_sim 坐标系说明

完整公式、转换顺序、数值验证和代码出处见[坐标与转换审计笔记](coordinate_transform_audit_zh.md)。交互核查页见[坐标系与真实事例核查](assets/lact-coordinate-system-3d.html)。

## 四个页面

1. **全局坐标定义**：只解释阵列、望远镜本地、镜片、遮挡、输出面和相机坐标，不绑定事例。
2. **平行光**：望远镜固定 `az=0° / el=70°`，分别运行 `el=71°`、`el=69°`、`az=-1°`、`az=+1°` 四组真实输入；选择相机图时同步切换三维光路、镜面反射点和遮挡诊断点。
3. **不同天顶角平行光**：按真实形变 CSV 切换镜片位置，用颜色表示相对理想镜片中心的位移，并显示同一次程序运行产生的光路与完整输出面 `u/v`。
4. **CORSIKA 事例**：同一 event 的芯位、到达方向、阵列、望远镜指向、二维 photon bunch、`zem` 反演发射点和相机图一起显示。

每个含地面的页面都固定采用北、东、南、西与天空方向。俯视时北在上、东在右；拖动只改变观察相机，不改变程序坐标或方位角定义。

## 必须牢记的坐标边界

```text
EventIO NWU photon bunch
  -> transformBunchToTelescopeLocal()
  -> telescope-local mirror + obstruction + output plane
  -> hit.u/v = dot(surface - plane_point, plane_u/v)
  -> camera_x/y = u/v
  -> ROOT: x_m/y_m = u/v
  -> pyLAST: pix_x=-v, pix_y=+u
```

| 数据层 | x | y | z / 说明 |
|---|---|---|---|
| CORSIKA / EventIO 原始输入 | North | West | Up；NWU |
| 地图显示 East/North | East = `-West` | North = `x` | 只用于地图显示 |
| 望远镜本地光学 | 水平横向 `local x` | 天空向上 `local y` | `local +z` 为光轴，指向天空 |
| 镜片与遮挡 CSV | local x | local y | local z；直接进入本地光学几何 |
| 输出面 / LACT 相机 | `u = local x` | `v = local y` | 当前输出面约为 local `z=-8 m` |
| pyLAST `LactEventSource` | `pix_x=-v` | `pix_y=+u` | 与最新版 `main` 一致 |
| pyLAST 当前画布 | 横轴 `pix_y=+u` | 纵轴 `pix_x=-v` | `plot_camera_image()` 的参数顺序 |

`buildTelescopeFrame()` 用于通用全局显示；`buildCorsikaNwuTelescopeFrame()` 用于 NWU 输入适配。两者的全局轴名称不同，但进入镜片、遮挡、输出面和相机的规范光学本地语义必须一致：`+x` 为横向、`+y` 为天空向上、`+z` 为光轴。网页不得再用旧的“CORSIKA local x=仰角、local y=方位”基底展开望远镜结构。

## theta / phi 的两种定义

程序中有两套不能混用的 `theta/phi`：

1. 合成平行光的 `source.beam_theta_deg / source.beam_phi_deg` 是望远镜本地角。实际代码计算 `d_local=(sinθ cosφ, sinθ sinφ, -cosθ)`，所以 `θ=0°` 沿 local `-z` 入射；`φ=0°/90°/180°/270°` 分别偏向 local `+x/+y/-x/-y`。
2. CORSIKA event header 的原始 `theta_deg / phi_deg` 是 shower 头字段。`theta` 从天顶量，因此 `altitude=90°-theta`；原始 `phi` 不能直接画成北起地图方位。程序使用 `A=(array_rotation-phi+180°) mod 360°` 得到 `azimuth_north_to_east_deg`，其中 `A=0°/90°/180°/270°` 分别为北/东/南/西。

四页左下角共用同一个固定地图：始终北上、南下、西左、东右，橙色箭头显示天空来向，方位 `A` 从 North 向 East 顺时针增加；地图下方同时给出 `A` 和从天顶量起的 `θz=90°-altitude`。第 2–4 页右上角的精简角度卡再显示当前事例的原始值、局部光束公式或 CORSIKA 换算公式。CORSIKA 页固定地图中的橙色箭头画的是程序换算后的 `A`，不是原始 `phi`。

固定地图在 CORSIKA 页严格使用程序 CORSIKA NWU 坐标：`x=North、y=West、z=Up`，为了“北上东右”，画布右侧 East 对应 `-y`。第 1–3 页的固定地图只是通用全局坐标的方位显示参考，不是 CORSIKA 输入 `x/y`；网页已在地图标题下明确标出这一区别。

固定地图旁边保留一个紧凑的“屏幕方向”框，它才随观察视角旋转；重合于视线的轴使用 `⊗屏幕内 / ⊙屏幕外`。三维投影恢复使用完整画布高度。望远镜角度、事例和观察视角控制合并在底部“视角 / 参数”面板中并默认收起，避免长操作提示和仰角滑块遮挡镜片、支架或遮挡点。

## CORSIKA 二维 bunch 的三维反演边界

EventIO 二维 bunch 提供记录面上的 `x/y`、传播方向、时间和 `zem`，其中 `z=0` 是输入参考面，不是发射高度。页面沿原始传播直线，用 `zem` 高度反推三维发射点：

```text
h = zem - observation_altitude - telescope_z
s = (h - anchor_z) / (-dir_z)
emission = anchor - s * direction
```

反演点是派生量，不是 EventIO 显式存储的三维坐标。生成器必须验证方向模长、反演高度和直线共线残差，并保留原始单位和来源信息。

## 输出面与相机

默认页面直接显示 LACT_sim 原始 `u/v`：横轴 `u`、纵轴 `v`，不平移、不旋转、不取反。选择 pyLAST 显示时只改变坐标表达，复用同一批 output 点、像素 id 和信号值。

物理相机范围严格按相机 CSV 的全部像素显示；光斑放大只改变取景范围。两种视图都不能移动物理零点或改变数值。

## 支架形变

理想镜片来自镜片配置，仰角系列来自支架形变 CSV。页面使用：

```text
delta_center_local = deformed_center_local - ideal_center_local
```

镜片颜色表示位移模长，箭头表示原始 local `dx/dy/dz`。箭头倍率只用于观察，不改变镜片位置和统计值。

## 生成与验证

```powershell
python python/build_coordinate_diagnostics_html.py --output docs/assets/lact-coordinate-system-3d.html
node python/test_coordinate_diagnostics_runtime.js docs/assets/lact-coordinate-system-3d.html
```

生成器必须核验镜片 id、形变差值、CORSIKA raw input、surface 与 `u/v`、输出行数、相机像素 id、event/array/telescope 元数据和来源版本。真实输出重新生成后，审计笔记会记录 Git commit、源码树哈希、配置、输入文件哈希、运行命令和结果哈希。

平行光三维图可省略大部分光路线，但必须保留全部真实镜面反射点和全部遮挡记录的诊断端点。当前 CSV 没有保存与遮挡 primitive 的精确相交位置，因此这些点只能称为“遮挡诊断端点”，不能解释为精确遮挡交点。
