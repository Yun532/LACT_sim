# LACT_sim 坐标与转换审计笔记

这份笔记是 `docs/assets/lact-coordinate-system-3d.html` 的唯一坐标定义来源。页面只能展示本笔记中已经找到代码出处的定义；没有程序原值的内容必须标记为“派生”或“仅用于坐标演示”。

## 1. 三套不能混用的三维基底

### 1.1 CORSIKA / EventIO 原始坐标：NWU

| 轴 | 程序含义 | 单位 |
|---|---|---|
| `+x` | magnetic North | m（EventIO 输入先由 cm 乘 `0.01`） |
| `+y` | West | m |
| `+z` | Up / 天空 | m |

二维 EventIO photon bunch 的程序读取规则：

```text
position_raw_nwu = (b.x*0.01, b.y*0.01, 0)
direction_raw_nwu = normalize(b.cx, b.cy, -sqrt(1-b.cx²-b.cy²))
time_ns = b.ctime
emission_altitude_km = b.zem*1e-5
```

代码出处：

- `src/io/EventIOPhotonSource.cpp:48`：`downwardDirZ()` 强制二维下降光子的 `dir_z<=0`。
- `src/io/EventIOPhotonSource.cpp:553`：`makeBunch()`；位置、方向、时间和 `zem` 的实际赋值在约 `570-590` 行。
- `src/io/EventIOPhotonSource.cpp:619`：三维 EventIO bunch 使用显式 `x/y/z/cx/cy/cz`，不能套用二维反演规则。

### 1.2 CORSIKA NWU 输入适配到望远镜本地

设 CORSIKA 方位角 `A` 从磁北向东增加，仰角 `E` 从地平线向天空增加。程序构造：

```text
e_x = (-sin(E)cos(A),  sin(E)sin(A), cos(E))
e_y = (-sin(A),       -cos(A),       0)
e_z = ( cos(E)cos(A), -cos(E)sin(A), sin(E))
```

含义：

```text
local +x = 仰角增加方向
local +y = 方位角 North→East 增加方向
local +z = boresight，镜面指向相机/天空
```

这是右手系：`e_x × e_y = e_z`。

转换规则：

```text
relative NWU position: p_local = (dot(p_nwu,e_x), dot(p_nwu,e_y), dot(p_nwu,e_z))
global NWU position:   p_local = (dot(p_nwu-origin,e_x), ...)
direction:             d_local = normalize(dot(d_nwu,e_x), dot(d_nwu,e_y), dot(d_nwu,e_z))
```

代码出处：

- `src/app/OpticalSimCommon.cpp:214`：`rotateVectorToLocal()` 使用三个点积。
- `src/app/OpticalSimCommon.cpp:218`：`pointToLocal()` 先减 `origin`。
- `src/app/OpticalSimCommon.cpp:584`：`buildCorsikaNwuTelescopeFrame()` 定义上面的三个基向量。
- `src/app/OpticalSimCommon.cpp:690` 附近：`transformBunchToTelescopeLocal()`；`corsika_nwu_relative` 只旋转，`corsika_nwu_global` 先减望远镜位置再旋转。
- `apps/test_coordinate_frames.cpp:31-70`：对 az=0、el=70° 的轴方向和 relative/global 往返做了数值断言。

### 1.3 通用/合成光线追迹全局基底

`run_optical_sim` 的平行光先在 `telescope_local` 中生成，然后使用另一套历史基底 `buildTelescopeFrame()` 把光子、镜片和输出面整体转到全局：

```text
generic e_z = (cos(E)cos(A), cos(E)sin(A), sin(E))
generic e_x = (-sin(A), cos(A), 0)
generic e_y = e_z × e_x
```

代码出处：

- `src/app/OpticalSimCommon.cpp:555`：`buildTelescopeFrame()`。
- `src/app/OpticalSimCommon.cpp:205`：`rotateVector()` 把本地向量展开到全局。
- `src/app/OpticalSimCommon.cpp:209`：`pointToGlobal()` 加上望远镜原点。
- `apps/run_optical_sim.cpp:26-45`：读取 telescope、先应用结构形变，再将镜片和输出面送入通用 frame。
- `apps/run_optical_sim.cpp:410-430`：输入先转入本地，再由 `applyTelescopeFrame(photon, telescope_frame)` 转入通用全局追迹坐标。

重要审计结论：`buildTelescopeFrame()` 和 `buildCorsikaNwuTelescopeFrame()` 不是同一横轴定义。页面第一部分必须同时画出并明确命名，不能把合成平行光的 generic frame 冒充 CORSIKA NWU input frame。

## 2. EventIO 二维参考面的 z 平移

EventIO 二维 bunch 经 NWU→local 旋转后，程序再执行：

```text
p_local.z += source.eventio_reference_z_m
```

当前 CORSIKA 示例配置为 `-16 m`。它只改变本地 `z`，不改变物理 `x/y`，也不是光子发射高度。

代码出处：

- `src/app/OpticalSimCommon.cpp:749`：`applyEventIOReferenceZOffset()`。
- `apps/test_coordinate_frames.cpp:192-205`：二维和三维输入都只平移 `z` 的断言。
- `configs/examples/corsika_coordinate_north_example.cfg`：显式配置 `source.eventio_reference_z_m=-16`。

## 3. 镜片、遮挡和形变坐标

镜片基础布局、形变布局和 obstruction primitives 都使用同一套望远镜本地 `x/y/z`。光线追迹前，程序再根据运行模式把它们整体放入追迹 frame。

数据出处：

- 理想镜片：`configs/mirror_1229_facets.csv`。
- 仰角形变：`configs/mirror_1229_elevation_series.csv`。
- 遮挡与支架：`configs/obstructions/raytrace_final_structure_primitives.csv`。

结构形变规则：

```text
当前仰角 = telescope.pointing_el_deg
center = 邻近仰角锚点的线性插值
normal = 邻近仰角锚点的 SLERP
实际镜片 center/normal = 形变 series 覆盖理想镜片的 center/normal
```

代码出处：

- `src/app/OpticalSimCommon.cpp:1641`：`slerpUnitVectors()`。
- `src/app/OpticalSimCommon.cpp:1655`：`buildElevationSeriesFacets()`；约 `1750-1790` 行执行中心插值和法向 SLERP。
- `src/app/OpticalSimCommon.cpp:3473`：`applyStructuralDeformation()`；约 `3500` 行覆盖 center/normal。
- `configs/errors/structural_deformation_1229.cfg`：把 error 配置连接到 elevation series。

页面中的形变量必须直接计算：

```text
delta_local = deformed_center_local - ideal_center_local
```

镜片颜色使用 `|delta_local|`；只有箭头端点允许乘显示倍率。

## 4. 输出面、相机 u/v 和法向

当前 8 m 输出面配置：

```text
plane_point  = (0,0,-8)
plane_normal = (0,0,-1)
plane_u_axis = (1,0,0)
plane_v_axis = (0,1,0)
```

输出坐标：

```text
rel = surface_point - plane_point
u_m = dot(rel, plane_u_axis)
v_m = dot(rel, plane_v_axis)
```

因此在当前配置中：

```text
u = local x
v = local y
camera CSV x_m = u
camera CSV y_m = v
camera plane normal = -local z（从相机朝镜面）
boresight / 天空方向 = +local z
```

代码出处：

- `configs/outputs/whiteboard_f8.cfg:2-5`：输出面点、法向和 u/v 轴。
- `src/optics/OpticalTracer.cpp:183-184`：实际 `u_m/v_m` 点积。
- `src/app/OpticalSimCommon.cpp:1150` 附近：相机 CSV 的 `x_m/y_m` 读入 `CameraPixel.center`。
- `include/io/SurfaceHitCsvWriter.hpp:28-58`：白板 CSV 同时保存 surface、u/v 和 camera 字段。

## 5. pyLAST 边界

LACT_sim 到 pyLAST 的已实现约定：

```text
pyLAST pix_x = -LACT_sim u
pyLAST pix_y = -LACT_sim v
```

pyLAST 当前绘图函数又把第二个坐标放在水平轴、第一个坐标放在垂直轴。诊断页面默认画 LACT_sim 原始 `u/v`；只有用户显式选择“pyLAST 相机显示”时才应用该变换，不得静默翻转。

之所以视觉上反常，是因为这里连续发生了两件事：ROOT reader 为得到 pyLAST 的源偏移定义而对 LACT 焦面传播坐标双轴取负；visualizer 再交换字段在 Matplotlib 画布上的横纵顺序。最终画布是 `horizontal=-v、vertical=-u`，但三维望远镜和物理光路没有因此再旋转一次。

诊断页面现在提供显式坐标选择器。`LACT 原始 u/v` 保持横轴 `u`、纵轴 `v`；`pyLAST 相机显示` 严格使用本机核对的 pyLAST 0.0.4 `root/LactEventSource.cpp:245-246` 与 `src/pylast/visualize/visualize.py:82-83`，即横轴 `pix_y=-v`、纵轴 `pix_x=-u`。选择器只改变坐标表达，必须复用同一 output 点、像素 id 和 `image_cherenkov_pe`，不得重新模拟或改变事件内容。

代码出处：

- `python/compare_minimal_csv_to_corsika_pylast.py:193-206`：用 `actual_x+config_x`、`actual_y+config_y` 验证双轴取负。
- `python/plot_photon_csv_root_pylast.py:31-55`：记录 pyLAST 的绘图轴顺序并提供 LACT `u/v` 视图。

## 6. CORSIKA 芯位、到达轴和 photon bunch 反演

`MC_TELOFF` 保存“阵列相对 shower core 的偏移”，所以程序输出芯位为：

```text
core_x_north_m = -MC_TELOFF.xoff
core_y_west_m  = -MC_TELOFF.yoff
```

代码出处：

- `src/io/EventIOPhotonSource.cpp:530-550`：读取并由 cm 转 m。
- `apps/run_corsika_trace.cpp:230-273`：`outputEventMetadata()` 按 array id 取负得到芯位。
- `apps/run_corsika_trace.cpp:2792-2795`：终端显示 West 时转换为 `core_E=-core_y_west`。

二维 bunch 没有显式三维发射点。诊断页面只做直线反演：

```text
height_relative_m = zem_km*1000 - observation_altitude_m - telescope_z_m
s = (height_relative_m-anchor_z)/(-direction_z)
emission_derived = anchor - s*direction
path_length_m = |emission_derived-anchor|
emission_time_derived_ns = record_time_ns-path_length_m/c
```

其中 `anchor/direction/zem/time` 是程序读取原值，`emission_derived` 和 `emission_time_derived` 必须标为派生。time 不能单独把二维位置唯一恢复成三维位置。

芯位与 photon anchor 的物理定义不同：芯位是 shower header 主轴在阵列 `z=0` 的交点；每个 anchor 是 EventIO telescope 表位置加该 photon bunch 的相对落点。诊断页因此不要求二者重合，而是做两项一致性检查：光子传播方向与 header 到达轴的夹角，以及反演发射点到同高度主轴点的横向距离。当前北向事例方向夹角中位数为 0.4698°、95 分位为 0.6867°，空间残差中位数为 6.23 m；对准事例分别为 0.3990°、0.6874° 和 6.36 m。这说明“看起来不重合”来自单光子横向分布，不是把 core 与 photon 坐标混用了。

实现出处：

- `python/prepare_event1909_coordinate_case.py:42-98`：event 1909 的原始 bunch 读取与 `zem` 直线反演。
- `src/io/EventIOPhotonSource.cpp:570-590`：C++ 中二维 bunch 的位置、方向、time 与 `zem` 赋值。

### 6.1 prod1 event 1909 的反演边界

页面第 4 部分直接读取本地 4,002,090,371 字节 prod1 EventIO 文件中的 `shower 19 / reuse 10`（程序输出编号 `event 1909 / array_id 9`），并与该输入生成的 ROOT 输出逐台对照。程序日志记录：

```text
input_photon_format = 2d
input_bunches_2d = 86550
input_bunches_3d = 0
read_emitter_info = false
```

这里的 `read_emitter_info=false` 只表示没有读取可选的 STORE-EMITTER 附加记录，不表示二维 photon bunch 的 `zem` 字段不存在。用 Python `eventio 2.1.1` 直接检查该事例，bunch dtype 包含 `x, y, cx, cy, time, zem, photons, wavelength`。例如 tel19 第一条 bunch 的原值为：

```text
x=390.14026 cm, y=16.61931 cm
cx=-0.3378992, cy=-0.014721361
time=89.006516 ns, zem=1600121 cm, photons=4.921779
```

页面分别从 tel19 / tel16 / tel21 / tel20 均匀抽样 240 条。每条的 `x/y/cx/cy/time/zem/photons` 保留 EventIO 原值；仅做程序已有的 cm→m 和二维向下方向补全，再用第 6 节公式派生三维发射点。发射点是 `anchor + direction + zem` 的直线反演结果，不是 EventIO 中直接储存的三维坐标，也不是任意长度的显示线段。

观测高度不再由网页或生成脚本硬编码。脚本直接读取所选 shower header 的 `n_observation_levels=1` 和 `observation_height[0]=440000 cm`，换算为 4400 m。对页面实际嵌入的全部 960 个抽样 bunch，生成期逐点验证方向为单位向量、反推点的 array-z 精确复现 `zem-observation_height`、以及 `emission + path_length*direction` 回到原始接收面 anchor。只有这些校验全部通过才写出网页数据。

所有 CORSIKA 三维坐标最终都以米存储：EventIO `x/y/zem` 乘 0.01，ROOT 阵列位置分支本身就是米。网页 `fit()` 对 x/y/z、全部望远镜、接收点、发射点只计算一个共同 `span`，所以没有逐轴或逐望远镜的独立缩放。散点半径是为了可见性而固定的屏幕像素，不代表光团的物理半径。

四台的 EventIO bunch 总数必须逐台等于同一 ROOT `trace_summary.input_bunches`；tel19 为 27,159 条。ROOT 中的阵列位置、芯位、shower header、逐台追迹统计和 `image_cherenkov_pe` 都绑定同一个 event 1909。

数据出处：

- 本地 4 GB `.zst`：用 `eventio.IACTFile` 直接选择 shower 19 / reuse 10 并读取四台 photon bunch。
- `python/prepare_event1909_coordinate_case.py:42-98`：均匀抽样、二维方向补全和基于原始 `zem` 的三维发射点反演。
- ROOT `corsika_events`：event、能量、到达角、芯位。
- ROOT `telescopes`：32 台 `array_x_north_m / array_y_west_m / array_z_up_m`。
- ROOT `trace_summary / observations`：逐台光学统计和相机 `image_cherenkov_pe`。

## 7. 四个页面的数据边界

1. **全局坐标定义**：只画一台完整理想望远镜，并在同一个原点叠加 NWU/CORSIKA 输入基底与 generic trace frame；不读取 event。完整镜片、遮挡、相机和 u/v 必须使用实际 `run_optical_sim` 的 `buildTelescopeFrame()` 展开，CORSIKA NWU 基底只画输入轴，不能驱动结构。程序角度默认 `az=0° / el=70°`，运行时测试要求其 54 块理想镜片的全部顶点与第 2 页同角度结构逐点一致。默认显示相机取 `view_el=0° / view_az=180°`，视线位于地平面内，因此地面严格投影为一条水平地平线；该显示视角不改变望远镜 pointing。
2. **平行光**：望远镜固定 `pointing_az=0° / pointing_el=70°`，读取四个独立的真实 `run_optical_sim` 运行。天区上/下光源为 `az=0° / el=71°、69°`；天区左/右光源为 `az=-1°、+1° / el=70°`。全局天空方向先用 `buildTelescopeFrame()` 的点积转换为程序所需本地传播向量，再原值写入 `source.beam_direction`。右上角四张白板 `u/v` 图分别按各自原始包围盒自动放大，保留绝对 u/v 刻度、正方向和米单位，并作为统一 case 选择器；选择后，三维区只显示该 run 的光路、镜面反射点和遮挡段。热图逐行统计完整 CSV 的 output `u/v`，完整 hit-camera `u/v` 保留用于强校验但不重复覆盖在光斑上。
3. **不同天顶角平行光**：每个角度有独立 `run_optical_sim` 输出，并与同角度的形变镜片绑定。默认显示相机同样取 `view_el=0° / view_az=180°`，所以地面只是一条地平线；显示相机角与望远镜仰角/天顶角是两个独立量。
4. **CORSIKA 事例**：固定绑定 4 GB prod1 的 event 1909；画 32 台真实阵列位置，以及 tel19 / tel16 / tel21 / tel20 的原始 photon anchor、direction 和 `zem` 派生发射点；四台各从完整 bunch 列表均匀抽取 240 条并用不同颜色的大散点显示。同一 ROOT 中信号最大的四台相机同步显示。逐光子镜面反射折线未写入 ROOT，页面不得伪造。

每个包含地面的页面都同时画地面三维轴和两套方向指示。固定地图始终“北上、东右”，不接受三维视角旋转；相邻的“当前屏幕投影”才随观察视角变化。CORSIKA 页面按 NWU 标为 `North=+x、West=+y、Sky=+z`，因此地图 East=`-y`；通用平行光页面按程序通用全局轴把 `North=+x、East=+y` 明确标成“显示参考”，避免把它误称为 CORSIKA 地理坐标。某方向与观察视线重合时，投影指示用 `⊗ 屏幕内 / ⊙ 屏幕外` 表示，因此即使默认地面侧视成一条地平线也不会丢失南北方向。页首说明、左侧原值与校验、相机图说明及形变图说明均默认折叠；平行光白板采用与页面一致的深色背景，这些仅是界面层变化，不修改数据或坐标。

三维观察相机也必须直接使用程序方位基底。通用页的屏幕右、屏幕上、深度轴分别等于 `buildTelescopeFrame()` 的 `x_axis/y_axis/z_axis`。CORSIKA 适配器的本地横轴定义不同：local `+x` 增加仰角、local `+y` 增加方位角，因此 CORSIKA 页的屏幕右、屏幕上、深度轴分别等于 `buildCorsikaNwuTelescopeFrame()` 的 `y_axis/x_axis/z_axis`。这样 `az=0°` 朝 North，`az=90°` 在通用显示参考中朝 `+y`、在 CORSIKA NWU 中朝 East=`-y`；从天空向地面看都表现为顺时针向东。鼠标拖动只连续修改观察 az/el，不允许更换叉乘符号或恢复另一套屏幕右轴。

页面不允许相机与三维区各自保留不同事例。平行光右上角图像可以被点击，但它修改的是唯一的 `selectedCase`，相机、三维光路、镜片统计和左侧说明必须同步切换。

平行光数据来自当前工作树源码在隔离目录中的真实 C++ 编译与运行，源归档 SHA-256 为 `afd630aedc69825f55eea3961a942ae7904c2b0bb6061f016775f3f10bc49afe`，基线/仰角结果归档 SHA-256 为 `47493ae480d4e5f00d38c656b39926f956c47da35415f7836c4802e6300cc77f`，新四方向结果与配置的组合 SHA-256 为 `8e584fdc9036fd980b5370f1176e298a19199f990d7d4076c0300f5ffb5af955`。每个 case 都保存自己的 basis、天空光源角、程序本地传播向量、raw input、镜面反射点、遮挡段标志、global surface、原始 `u/v`、相机像素和 C++ summary。

四方向结果所用服务器源码与当前工作区逐文件 SHA-256 再核对一致：`apps/run_optical_sim.cpp=13bb4068…`、`src/app/OpticalSimCommon.cpp=9984bfd5…`、`include/io/SurfaceHitCsvWriter.hpp=de7885de…`；运行二进制 SHA-256 为 `5c486b2e…`。当前分支后来只修改了 CORSIKA 输出代码，没有修改这三个平行光执行路径文件。

右上角相机核对不再使用 400 点抽样：四组分别嵌入全部 `14676 / 14619 / 14632 / 14674` 个物理 output-plane `u/v`，并严格仿照程序 `plot_whiteboard()`：每组按自身包围盒中心取最大边长构成正方形范围，四周留 8% 边距，以 `140×140` 分箱绘制二维 count/bin 热图；红、蓝、黄、绿主色与对应三维光路一致，白叉为全量质心，白圈为 R68。LACT 模式显示未平移的绝对 u/v 原值和正方向；pyLAST 模式对同一批点显示横轴 `-v`、纵轴 `-u`。取景中心不强制为 `(0,0)`，零点只有位于当前范围内才画零线。全部 `11243 / 11212 / 9100 / 9133` 个 `hit_camera` 的 `camera_x/y` 仍保留在事例数据和强校验中，但不再作为第二层散点遮住光斑密度。完整 hits CSV 按 `pixel_id` 重计数后，必须逐像素等于完整 camera CSV 的 `photon_count`。三维区也直接嵌入并绘制全部 `hit_mirror` 镜面反射点（包括后来在反射段被遮挡的行），以及全部 incoming/reflected 遮挡记录对应的诊断端点；只有光路线段为避免遮挡而分层抽取 48 条。

相机画布不从三维图重新投影落点。二维热图逐行读取每个 run 的完整 hits CSV `u_m/v_m` 并统计 count/bin；完整 `camera_x_m/camera_y_m` 与 camera CSV 像素计数继续参与生成期校验，但不作为第二层散点覆盖光斑。生成期逐行验证所有 camera hit 的 `camera_x_m==u_m、camera_y_m==v_m`，同时用该 run 保存的 `buildTelescopeFrame` 基底把 global surface 投影回局部面并复现同一 `u/v`。LACT 模式不做变换；pyLAST 模式只在画布入口按源码双轴取负并交换绘图轴，不回写或修改原始数组。

为保留遮挡诊断行，四方向配置显式设置 `obstruction.mark_only=true`。CSV 的 `obstruction_blocked_incoming/reflected` 是 C++ 原始标志；`mirror_point` 和 `surface_point` 也是 tracer 输出原值。程序没有输出光线与 obstruction primitive 的精确交点，因此页面将全部 incoming-blocked 行画成其理论镜面端点红叉，将全部 reflected-blocked 行画成其理论输出端点紫叉。它们是完整的遮挡记录诊断点集合，但不能冒充 primitive 的精确相交位置。

## 8. 生成期强制校验

- 镜片与形变表 id 完整一致。
- `delta = deformed-ideal` 逐坐标成立。
- `surface_x==u`、`surface_y==v` 对当前输出配置成立。
- camera pixel id 存在于 1616 像素几何。
- 每个 case 的 ray/output/camera 标识一致。
- CORSIKA 的 event/array/telescope/pointing 一致。
- 每个真实平行光 case 的 global surface 用该 run 基底投影后必须复现原始 `u/v`。
- camera `x/y` 必须等于原始 `u/v`，像素光子总数必须等于 C++ `hit_camera`。
- 平行光必须同时保留当前源码归档与结果归档的 SHA-256。
- 四方向平行光必须精确包含天空光源 `(az,el)=(0,71),(0,69),(-1,70),(+1,70)`，并验证其转换后的本地传播向量。相机质心必须体现反射倒像：天区上/下落在 `-v/+v`，天区左/右落在 `+u/-u`。
- `mark_only_rows` 必须等于 `hit_output_before_obstruction`，无遮挡 CSV 行数必须等于 `hit_output_plane`；像素和仍必须等于 `hit_camera`。
- 四方向相机图嵌入的完整 output `u/v` 行数必须等于 `hit_output_plane`，完整 camera-hit `u/v` 行数必须等于 `hit_camera`；按 `pixel_id` 重计数必须逐像素等于 camera CSV。
- 四方向三维图嵌入的完整镜面反射点数必须等于 `hit_mirror`（包括后来在反射段被遮挡的光子），incoming/reflected 遮挡诊断端点数必须分别等于 C++ `blocked_incoming/blocked_reflected`；只有光路线段允许抽样。
- event 1909 必须包含 32 台望远镜；四台 EventIO raw bunch 数必须逐台等于 ROOT `input_bunches`，其中 tel19 为 27,159 条；相机图只读同一 ROOT 的 `image_cherenkov_pe`。
- event 1909 的观测高度必须直接来自所选 EventIO event header；每个显示发射点必须同时复现原始 `zem` 高度并与原始 anchor/direction 共线。
- CORSIKA 的 EventIO cm→m、ROOT m 和网页单一等比例 span 必须显式记录；仅标记半径允许使用屏幕像素。
- 两套 telescope frame 都保持单位长度、正交、右手性；页面标签不得互换。
