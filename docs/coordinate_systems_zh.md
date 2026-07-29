# LACT_sim 坐标系说明

完整公式、转换顺序、数值验证和代码出处见[坐标与转换审计笔记](coordinate_transform_audit_zh.md)。交互核查页见[坐标系与真实事例核查](assets/lact-coordinate-system-3d.html)。

## 四个页面

1. **全局坐标定义**：只解释阵列、望远镜本地、镜片、遮挡、输出面和相机坐标，不绑定事例。
2. **平行光**：望远镜固定 `az=0° / el=70°`，四组光源都与光轴精确相差 `4°`：上/下沿本地 `±v`，东/西沿本地 `±u`。东、西不是把天空方位角简单加减 `4°`；考虑 70° 仰角后，它们等效为 `az=±11.5550088° / el=69.6199950°`。四组均由最新程序独立运行；选择相机图时同步切换三维光路、镜面反射点和遮挡诊断点。
3. **不同天顶角平行光**：按真实形变 CSV 切换镜片位置，用颜色表示相对理想镜片中心的位移，并显示同一次程序运行产生的光路与完整输出面 `u/v`。
4. **CORSIKA 事例**：进入页面首先显示完整簇射总览，默认南朝屏幕外；同一 event 的芯位、到达方向、32 台阵列/指向、二维 photon bunch 与 `zem` 反演发射点一起显示。右上角可在 tel19/16/21/20 四幅 ROOT 相机图和 pyLAST 阵列/芯位图之间切换。点击相机或芯位图中的望远镜只把完整场景自动缩放到该台并把旋转中心设为该台，不删除其他 31 台模型或其他 photon 路径；“返回完整簇射”只恢复初始尺度和视角。

每个含地面的页面都固定采用北、东、南、西与天空方向。俯视时北在上、东在右；拖动只改变观察相机，不改变程序坐标或方位角定义。

## 必须牢记的坐标边界

```text
EventIO NWU photon bunch
  -> transformBunchToTelescopeLocal()
  -> telescope-local mirror + obstruction + output plane
  -> hit.u/v = dot(surface - plane_point, plane_u/v)
  -> camera_x/y = u/v
  -> ROOT: x_m/y_m = u/v
  -> pyLAST: pix_x=-v, pix_y=-u
```

| 数据层 | x | y | z / 说明 |
|---|---|---|---|
| CORSIKA / EventIO 原始输入 | North | West | Up；NWU |
| 地图显示 East/North | East = `-West` | North = `x` | 只用于地图显示 |
| 望远镜本地光学 | 水平横向 `local x`，朝方位角增加方向；`A=0°` 时为 East | 天空向上 `local y` | `local +z` 为光轴，指向天空 |
| 镜片与遮挡 CSV | local x | local y | local z；直接进入本地光学几何 |
| 输出面 / LACT 相机 | `u = local x` | `v = local y` | 当前输出面约为 local `z=-8 m` |
| pyLAST `LactEventSource` | `pix_x=-v` | `pix_y=-u` | 与 LACT 坐标修复 `4fb44f8` 配套 |
| pyLAST 当前画布 | 横轴 `pix_y=-u` | 纵轴 `pix_x=-v` | `plot_camera_image()` 的参数顺序 |
| pyLAST 阵列/芯位图 | East=`-array_y_west` | North=`array_x_north` | `plot_event_cores()` 的地面轴 |

`buildTelescopeFrame()` 用于通用全局显示；`buildCorsikaNwuTelescopeFrame()` 用于 NWU 输入适配。两者的全局轴名称不同，但进入镜片、遮挡、输出面和相机的规范光学本地语义必须一致：`+x/+u` 朝方位角增加方向，`+y/+v` 为天空向上，`+z` 为光轴。网页不得再用旧的“CORSIKA local x=West”基底展开望远镜结构。

CORSIKA 方位角从 North 向 East 增加。令 `A=azimuth`、`E=elevation`，修正后的程序基底是：

```text
e_x = (-sin(A),        -cos(A),        0     )   # local +x / +u
e_y = (-sin(E) cos(A),  sin(E) sin(A), cos(E))   # local +y / +v
e_z = ( cos(E) cos(A), -cos(E) sin(A), sin(E))   # 光轴
```

NWU 的物理轴顺序 `North-West-Up` 是左手排列，所以这套物理右手光学基在 NWU 数值分量中满足 `e_x × e_y = -e_z`。不能为了让数值叉积变成 `+e_z` 而翻转 `e_x`；旧实现正是因此让同一个向东偏移的光源在通用入口与 CORSIKA 入口产生相反的相机 `u`。程序回归测试现已直接比较两种入口的“向东偏 1°”结果。

## theta / phi 的两种定义

程序中有两套不能混用的 `theta/phi`：

1. 合成平行光的 `source.beam_theta_deg / source.beam_phi_deg` 是望远镜本地角。实际代码计算 `d_local=(sinθ cosφ, sinθ sinφ, -cosθ)`，所以 `θ=0°` 沿 local `-z` 入射；`φ=0°/90°/180°/270°` 分别偏向 local `+x/+y/-x/-y`。
2. CORSIKA event header 的原始 `theta_deg / phi_deg` 是 shower 头字段。`theta` 从天顶量，因此 `altitude=90°-theta`；原始 `phi` 不能直接画成北起地图方位。程序使用 `A=(array_rotation-phi+180°) mod 360°` 得到 `azimuth_north_to_east_deg`，其中 `A=0°/90°/180°/270°` 分别为北/东/南/西。

四页左下角共用同一个固定地图：始终北上、南下、西左、东右，橙色箭头显示天空来向，方位 `A` 从 North 向 East 顺时针增加；地图下方同时给出 `A` 和从天顶量起的 `θz=90°-altitude`。第 2–4 页右上角的精简角度卡再显示当前事例的原始值、局部光束公式或 CORSIKA 换算公式。CORSIKA 页固定地图中的橙色箭头画的是程序换算后的 `A`，不是原始 `phi`。

固定地图在 CORSIKA 页严格使用程序 CORSIKA NWU 坐标：`x=North、y=West、z=Up`，为了“北上东右”，画布右侧 East 对应 `-y`。第 1–3 页的固定地图只是通用全局坐标的方位显示参考，不是 CORSIKA 输入 `x/y`；网页已在地图标题下明确标出这一区别。

固定地图旁边保留一个紧凑的“屏幕方向”框，它才随观察视角旋转；重合于视线的轴使用 `⊗屏幕内 / ⊙屏幕外`。三维投影恢复使用完整画布高度。望远镜角度、事例和观察视角控制合并在底部“视角 / 参数”面板中并默认收起，避免长操作提示和仰角滑块遮挡镜片、支架或遮挡点。

## CORSIKA 望远镜位置与 pyLAST 芯位图

三维模型始终用 ROOT `telescopes` 表中的原始 NWU 位置 `(North, West, Up)`。地面视图在每台模型旁显示 `tel_id` 和 North/West 坐标；其他尺度至少在所选模型旁显示三位小数坐标，避免完整簇射尺度下 32 组长标签重叠。

右上角的“pyLAST 芯位图”不是嵌入 PNG，而是逐条模仿 `EventVisualizer.plot_event_cores()` 重画同一 event：望远镜横轴坐标为 `East=-West`，纵轴为 `North`，红星为 ROOT 真芯位，红箭头为 event header 的真到达方向。pyLAST 原函数把望远镜标成 `T{tel_id+1}`，所以图中的 `T20` 对应网页和 ROOT 的 `tel19`；标题同时显示两种编号及所选台精确坐标。为提高可读性，网页采用纯黑背景和高对比线性色带，但不改变 PE 归一化、坐标或芯位数值。

该 ROOT 的 32 个 `triggered` 标志均为 `false`，但 32 台都有真实 `signal_pe`。芯位图明确采用 `plot_event_cores(..., include_non_triggered=True)` 的显示语义，使颜色表示原始 PE 而不是把全部望远镜误画成零；没有补造触发状态。点击图中的任一望远镜会选择对应 ROOT `telescope_id`，并把三维视图聚焦到同一坐标。

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

光学程序里一直同时出现过两层字段名，但它们不是两套坐标。光线与输出面的交点先按输出面基底投影成 `hit.u_m/hit.v_m`；进入探测器像素查找时，代码直接令 `camera_x_m=hit.u_m`、`camera_y_m=hit.v_m`。相机配置 CSV 和 ROOT 相机几何再把像素中心写成 `x_m/y_m`。因此当前配置下完整对应关系是：

```text
输出面       hit.u_m      hit.v_m
探测器接口   camera_x_m   camera_y_m
相机几何     x_m          y_m
LACT 图      横轴 u       纵轴 v
              ↑数值相等    ↑数值相等
```

这里的相机 `x_m/y_m` 是相机面字段，不能与 CORSIKA 阵列 NWU 的 North/West `x/y` 混为一谈。

CORSIKA 页的 LACT u/v 也不是由网页根据 photon bunch 的方向反算。相机几何直接来自同一个 ROOT 的 `camera_pixels.x_m/y_m/size_m`，依据上面的程序赋值链有 `x_m≡u`、`y_m≡v`；所选望远镜的信号来自 `observations.image_cherenkov_pe`，按 `pixel_id` 与这 1616 个几何像素连接。生成时会把 ROOT 几何与相机配置逐像素比较。

该 ROOT 中 32 台望远镜的 `triggered` 标志全部为 `false`。因此右上角原来的四幅图只能按原值准确称为“主要成像望远镜 / `image_cherenkov_pe` 前四”，不能改写成 ROOT 已触发；它们按 PE 总和排序正好是 tel19、tel16、tel21、tel20。完整簇射总览固定显示这四台；若下拉框选择其他望远镜，聚焦视图的右下角会替换为所选台的真实 ROOT 图像，前三幅仍保留 PE 前三，返回总览后恢复 PE 前四。所有 32 台都有各自的 `image_cherenkov_pe`，这里没有合成光斑。

四台主图使用固定红/蓝/黄/绿，其他望远镜也按 `telescope_id` 使用固定颜色；相机、三维模型、pointing 与 photon 路径复用同一个颜色。总览中的光子线是从每台完整 EventIO bunch 中均匀抽出的真实记录，保留完整条数与抽样出处；聚焦时仍保留全部 32 组当前已载入路径。ROOT 没有保存逐光子的镜面反射折线，所以网页不伪造反射段。

物理相机范围严格按相机 CSV 的全部像素显示；光斑放大只改变取景范围。两种视图都不能移动物理零点或改变数值。

LACT 和 pyLAST 像素图现在统一使用同一套暗色背景、像素形状、强度色标与物理/光斑缩放方式。统一的是网页画法，不是坐标数值：LACT 仍直接画 `(u,v)`；pyLAST 仍严格画 reader 和 renderer 得到的 `(plot_x,plot_y)=(+u,-v)`。

`pix_x=-v, pix_y=-u` 的作用是把 LACT 的物理焦平面基底适配到 pyLAST 的源偏移基底。LACT 中 `u` 是方位型、`v` 是仰角型；pyLAST reconstruction 中 `pix_x` 是仰角型、`pix_y` 是方位型，所以必须交换轴。坐标修复 `4fb44f8` 后，LACT 的 `+u` 随 North→East 方位增加，而 pyLAST TelescopeFrame 的正方位基底朝 West，因此 `pix_y=-u`；天空向上在焦平面 `v` 上经过镜面成像反向，因此 `pix_x=-v`。边界矩阵为 `[[0,-1],[-1,0]]`。这是两个程序坐标基底之间的适配，不是网页为了视觉一致而翻图。

真实 4° 平行光给出同样的符号：上方光源的 LACT 质心 `v=-0.571319 m`，转成 pyLAST 后 `pix_x=+0.571319 m`；东侧光源的 LACT 质心 `u=-0.571243 m`，转成 pyLAST 后 `pix_y=+0.571243 m`。这些值由页面直接读取程序输出，不是示意数值。

event 1909 的端到端复算进一步锁定了第二个负号：在同一份固定版 LACT ROOT、同一批 10 台望远镜和同一重建流程下，`(-v,-u)` 得到方向误差 `0.05450°`、芯位误差 `13.04 m`；保留旧的 `(-v,+u)` 则变成 `1.57150°` 和 `165.06 m`。因此 LACT 的 `4fb44f8` 修复与 pyLAST reader 的 `pix_y=-u` 必须作为一组变更。

## 支架形变

理想镜片来自镜片配置，仰角系列来自支架形变 CSV。页面使用：

```text
delta_center_local = deformed_center_local - ideal_center_local
```

第 3 页可在“指向偏转”和“中心位移”之间切换，默认显示指向偏转。指向模式的颜色及橙色/蓝色曲线分别表示每片镜子的法向夹角、54 片中的最大夹角和夹角 RMS，单位为 `mdeg`；箭头以实际形变法向为 1×，100×/500× 只放大相对理想法向的角差。位置模式的颜色及曲线表示中心位移模长，箭头表示原始 local `dx/dy/dz`。两种模式都让 0°–90° 的全部镜片共用从零到各自全局最大值的固定物理标尺。选择器只改变诊断显示，光线追迹始终同时采用形变 CSV 中真实的 center 和 normal。

第 2、3 页的平行光与第 4 页 event 1909 都已在 `main@48ea631` 加坐标修复 `4fb44f8` 的同一份源码上重新运行；event 1909 重新读取了原始 4 GB EventIO 文件。各类数据的来源提交、二进制、输入和结果哈希均在网页数据与审计笔记中记录，不把旧二进制结果冒充为新运行。

## 生成与验证

```powershell
python python/build_coordinate_diagnostics_html.py --output docs/assets/lact-coordinate-system-3d.html
node python/test_coordinate_diagnostics_runtime.js docs/assets/lact-coordinate-system-3d.html
```

生成器必须核验镜片 id、形变差值、CORSIKA raw input、surface 与 `u/v`、输出行数、相机像素 id、event/array/telescope 元数据和来源版本。真实输出重新生成后，审计笔记会记录 Git commit、源码树哈希、配置、输入文件哈希、运行命令和结果哈希。

平行光三维图可省略大部分光路线，但必须保留全部真实镜面反射点和全部遮挡记录的诊断端点。当前 CSV 没有保存与遮挡 primitive 的精确相交位置，因此这些点只能称为“遮挡诊断端点”，不能解释为精确遮挡交点。
