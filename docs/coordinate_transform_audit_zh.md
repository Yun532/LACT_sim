# LACT_sim 坐标与转换审计笔记

这份笔记是 `docs/assets/lact-coordinate-system-3d.html` 的坐标审计记录。页面只能展示已经由实际代码、数值测试或真实输出证明的定义；没有程序原值的内容必须标记为“派生”或“仅用于坐标演示”。注释只能帮助定位，不能单独作为结论。

本轮审计基线是远端 `main` 提交 `da51c09e9f93df44a3dc956a6a40eb028efb25f5`（`Update base0 elevation-dependent mirror geometry`）。判定顺序是：读取实际函数实现 → 独立计算已知角度的数值基向量 → 运行 C++ 回归测试 → 用最新版二进制重新生成输出 → 由网页运行时测试逐顶点比较结构。任何一项不一致都不得标记为通过。

## 1. 三层必须明确边界的三维基底

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
e_x = (-sin(A),        -cos(A),        0)
e_y = (-sin(E)cos(A),   sin(E)sin(A), cos(E))
e_z = ( cos(E)cos(A), -cos(E)sin(A), sin(E))
```

含义：

```text
local +x = 水平横向；A=0 时指向 East，与方位角增加方向一致；相机 u=local x
local +y = 仰角增加方向 / sky-up
local +z = boresight，镜面指向相机/天空
```

NWU 的物理轴顺序 `North-West-Up` 是左手排列。光学本地坐标本身仍按物理右手语义定义，但写成 NWU 数值分量时满足 `e_x × e_y = -e_z`。不能为了让数值叉积看起来为 `+e_z` 而翻转 `e_x`，否则同一个向东偏移光源会在通用入口和 CORSIKA 入口得到相反的相机 `u`。

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
- `apps/test_coordinate_frames.cpp`：对 az=0、el=70° 的三轴数值、NWU 物理手性、relative/global 往返、71° 光源映射，以及同一向东偏 1° 光源跨通用/CORSIKA 入口的 `u/v` 符号一致性做断言。

### 1.3 通用/合成光线追迹全局基底

`run_optical_sim` 的平行光先在 `telescope_local` 中生成，然后使用通用历史全局基底 `buildTelescopeFrame()` 把光子、镜片和输出面整体转到全局：

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

重要审计结论：两者的全局坐标语义仍不同——通用全局只定义 `az: +x→+y`，不能自动把其 `+x/+y` 宣称为 CORSIKA 的 North/East；但进入镜片、遮挡、输出面和相机后的规范光学本地语义必须一致：`local +x/+u` 朝方位角增加方向，`local +y/+v` 为 sky-up，`local +z` 为光轴。早期网页曾把 CORSIKA `local x/y` 交换；随后网页虽跟随了程序，却也如实暴露出程序把水平轴取成 West 的缺陷。两者现在都由跨入口向东偏 1° 的数值回归锁定。

对 `A=0°、E=70°`，最新版函数必须实际返回：

```text
CORSIKA NWU e_x = (0, -1, 0)  # 物理 East
CORSIKA NWU e_y = (-sin70°, 0, cos70°)
CORSIKA NWU e_z = ( cos70°, 0, sin70°)
```

网页运行时测试不读取注释，而是调用页面公式计算上述数值；比较通用 NEU 与 CORSIKA NWU 几何时先应用唯一的物理轴边界 `East=-West`，再把 54 块镜片逐顶点减去各自镜面顶点比较，最大误差必须小于 `1e-9 m`。测试还直接构造同一物理向东偏 1° 的光源，要求两种入口得到完全相同的 local `u/v` 符号。

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

第 3 页默认以镜片法向夹角 `acos(normal_ideal·normal_deformed)` 着色，并绘制 54 片镜子的最大夹角/RMS 曲线；用户可切换到 `|delta_local|` 中心位移着色及其最大值/RMS 曲线。指向模式的箭头在 1× 时使用真实形变法向，100×/500× 只放大相对理想法向的角差；位置模式的箭头端点按同一倍率放大 `delta_local`。颜色归一范围不得按当前仰角自动变化：两种模式分别让 0°–90° 的全部镜片共用从零到各自全局最大值的固定范围。选择器不得修改实际镜片或光路数据。

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
pyLAST pix_x = -LACT_sim v
pyLAST pix_y = -LACT_sim u
```

pyLAST 当前绘图函数把第二个坐标放在水平轴、第一个坐标放在垂直轴。诊断页面默认画 LACT_sim 原始 `u/v`；只有用户显式选择“pyLAST 相机显示”时才应用该变换，不得静默翻转。

最终画布是 `horizontal=pix_y=-u、vertical=pix_x=-v`。这是 LACT `u/v` 到 pyLAST 字段的边界映射加上绘图参数顺序，不会让三维望远镜或物理光路再旋转一次。

诊断页面提供显式坐标选择器。`LACT 原始 u/v` 保持横轴 `u`、纵轴 `v`；`pyLAST 相机显示` 使用与 LACT 坐标修复 `4fb44f8` 配套并经 event 1909 重建核对的 `pix_x=-v、pix_y=-u`，即横轴 `-u`、纵轴 `-v`。选择器只改变坐标表达，必须复用同一 output 点、像素 id 和 `image_cherenkov_pe`，不得重新模拟或改变事件内容。

代码出处：

- `python/compare_minimal_csv_to_corsika_pylast.py`：数值核对 LACT `u/v` 与 pyLAST `pix_x/pix_y`。
- `python/plot_photon_csv_root_pylast.py:31-58`：明确 reader 与绘图轴顺序；当前配套映射为 `pix_x=-v、pix_y=-u`。

### 5.1 `theta/phi` 角度方向审计

代码中同名角度存在两个作用域，网页必须明确标注，不能互相替代：

```text
合成平行光（望远镜本地）：
d_local = (sin(theta) cos(phi), sin(theta) sin(phi), -cos(theta))
theta=0° -> local -z
phi=0°/90°/180°/270° -> 偏向 local +x/+y/-x/-y
```

代码出处为 `src/app/OpticalSimCommon.cpp:1936-1947`。`SyntheticPhotonSource` 最终直接使用归一化后的 `beam_direction`，见 `src/io/SyntheticPhotonSource.cpp:92-103`。因此这里的 `phi` 是望远镜本地 x/y 平面的角，不是北起地图方位。

CORSIKA shower header 则按以下顺序保留原值并生成地图方位：

```text
theta_deg = MC_EVTH[10]                    # 从天顶量
phi_deg = MC_EVTH[11]                      # CORSIKA 原始 phi
altitude_deg = 90° - theta_deg
A_north_to_east = (array_rotation - phi + 180°) mod 360°
```

代码出处为 `src/io/EventIOPhotonSource.cpp:300-312`；ROOT 同时保存 `theta_deg`、`phi_deg`、`azimuth_north_to_east_deg` 和 `altitude_deg`，见 `src/io/LactEventRootWriter.cpp:863-867`。固定角度参考的橙色箭头只画最终 `A_north_to_east`：`0°=North、90°=East、180°=South、270°=West`。原始 `phi` 只作为原值和换算输入显示，不能单独贴到北上东右地图上。

event 1909 的实际原值为：`theta=19.3797875057°`、`phi=-179.1167112851°`、`array_rotation=-1.4299999943°`；代入公式得到 `A=357.6867155597°`，与 ROOT 的 `azimuth_north_to_east_deg` 相同，`altitude=70.6202124943°`。

固定角度参考的视觉方向在四页完全相同：俯视始终北上东右，`A/φ_map` 从北向东顺时针；侧视始终天顶在上、地平在右，`θ_sky` 从天顶向地平量。坐标边界仍分页面：CORSIKA 页 `x=North、y=West、z=Up` 是代码输入/ROOT 输出的 NWU 原值，地图 East 因而是 `-y`；全局、平行光和仰角页的 `North=+x、East=+y` 仅为通用全局坐标的显示参考，不是 CORSIKA 输入 `x/y`。网页标题下必须显示当前属于哪一类。

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

1. **全局坐标定义**：只画一台完整理想望远镜，不读取 event。完整镜片、遮挡、相机和 `u/v` 使用实际 `run_optical_sim` 的 `buildTelescopeFrame()` 展开；CORSIKA NWU 轴只用于解释另一种输入坐标，不能把两个全局 `x/y` 的名字直接等同。程序角度默认 `az=0° / el=70°`，运行时测试要求其 54 块理想镜片的全部顶点与第 2 页同角度结构逐点一致。默认观察为地平侧视，因此地面投影为一条水平线；观察角不改变望远镜 pointing。
2. **平行光**：望远镜固定 `pointing_az=0° / pointing_el=70°`，读取四个独立的真实 `run_optical_sim` 运行。天区上/下光源为 `az=0° / el=71°、69°`；天区左/右光源为 `az=-1°、+1° / el=70°`。全局天空方向先用 `buildTelescopeFrame()` 的点积转换为程序所需本地传播向量，再原值写入 `source.beam_direction`。右上角四张白板 `u/v` 图分别按各自原始包围盒自动放大，保留绝对 u/v 刻度、正方向和米单位，并作为统一 case 选择器；选择后，三维区只显示该 run 的光路、镜面反射点和遮挡段。热图逐行统计完整 CSV 的 output `u/v`，完整 hit-camera `u/v` 保留用于强校验但不重复覆盖在光斑上。
3. **不同天顶角平行光**：每个角度有独立 `run_optical_sim` 输出，并与同角度的形变镜片绑定。默认显示相机同样取 `view_el=0° / view_az=180°`，所以地面只是一条地平线；显示相机角与望远镜仰角/天顶角是两个独立量。
4. **CORSIKA 事例**：固定绑定 4 GB prod1 的 event 1909；画 32 台真实阵列位置，以及 tel19 / tel16 / tel21 / tel20 的原始 photon anchor、direction 和 `zem` 派生发射点；四台各从完整 bunch 列表均匀抽取 240 条并用不同颜色的大散点显示。同一 ROOT 中信号最大的四台相机同步显示。逐光子镜面反射折线未写入 ROOT，页面不得伪造。

每个包含地面的页面都同时画地面三维轴和两套方向指示。左下固定地图始终“北上、南下、西左、东右”，不接受三维视角旋转；相邻的紧凑“屏幕方向”才随观察视角变化。CORSIKA 页面按 NWU 标为 `North=+x、West=+y、Sky=+z`，因此地图 East=`-y`；通用平行光页面按程序通用全局轴把 `North=+x、East=+y` 明确标成“显示参考”，避免把它误称为 CORSIKA 地理坐标。某方向与观察视线重合时，投影指示用 `⊗ 屏幕内 / ⊙ 屏幕外` 表示，因此即使默认地面侧视成一条地平线也不会丢失南北方向。`fit()` 仍以全部场景点的包围球半径确定统一 span，但三维投影恢复使用完整画布高度。望远镜角度、事例和观察视角控件合并到默认收起的“视角 / 参数”面板；展开后才显示简短的 `az/el`、`el/θz` 和拖动操作，避免原来的长提示与滑块覆盖镜片、支架和遮挡点。页首说明、左侧原值与校验、相机图说明及形变图说明也默认折叠；这些仅是界面层变化，不修改数据或坐标。

三维观察相机必须服从所显示的全局坐标。通用页按 `buildTelescopeFrame()` 的抽象全局基底观察；CORSIKA 页按 NWU 观察，其中 `North=+x、East=-y、Up=+z`。最新版 CORSIKA 光学轴为 `local x=水平横向、local y=sky-up、local z=boresight`；俯视屏幕为了保持北上东右，屏幕右方向在 NWU 中是 East=`-y`，因此它是 `-local x` 而不是把 `local x/y` 交换。鼠标拖动只修改观察相机，不修改程序基向量。

页面不允许相机与三维区各自保留不同事例。平行光右上角图像可以被点击，但它修改的是唯一的 `selectedCase`，相机、三维光路、镜片统计和左侧说明必须同步切换。

本轮平行光和 CORSIKA 数据都来自坐标修复提交 `4fb44f876f37a7c920f0705abc3a8552533588e9`（基于当时最新 `main@48ea631`）的同一份隔离源码编译。源码归档 SHA-256 为 `e528d7f24e1410c2daf0822e9e920da82cfd53518b0a1a9bcd57fb9fb60917e7`；包含基线、10 个仰角、四方向、event 1909 ROOT/CSV 和运行日志的结果归档 SHA-256 为 `db2e5e5a7b02fe3f511f6c7564713ea9841a8ce9ced0d602594846b54bc57b0d`；四方向输入/输出组合 SHA-256 为 `e860355e4cf1a52f71a65a9e4fd0bd324e6d1625fedca914a5cfcbe2405eea12`。结果包下载后再次计算得到相同 SHA-256。

实际执行文件 SHA-256：`test_coordinate_frames=fd5e28bca7dc3d43deaf0883cfd0ed10255066d9e5f18cfac3086b59e73e8c36`、`run_optical_sim=4058e7a9c0cba8e2f3d662f646043f0fb4c30ed75a009070ee573af09800cddb`、`run_corsika_trace=d53f84de2db69de785c5f9e5b7e84b1bba8407694a90c26b588c41615bec12ff`。`test_coordinate_frames` 在该构建中返回 0。服务器实际执行路径文件的 SHA-256 为：`apps/run_optical_sim.cpp=ebf65379…`、`apps/run_corsika_trace.cpp=6682d808…`、`src/app/OpticalSimCommon.cpp=9685f610…`、`include/io/SurfaceHitCsvWriter.hpp=df7a1de1…`、`apps/test_coordinate_frames.cpp=91963c3c…`。

4 GB EventIO 原文件大小为 `4,002,090,371` 字节，SHA-256 为 `3feee5b7f3a001858201eea2cf75ba3f5f0277283e29900b5f259bd2c9bc4220`；本轮 event 1909 ROOT 的 SHA-256 为 `36876c1e586c3c6dcc1ff5d2f7a49cc66430a0cfce3db019454dde3684bf29bc`，运行日志 SHA-256 为 `6133e36b173121087cc1b5b066db128284b782f239a9cd50fa6dfdf832f7887f`。网页生成器会拒绝缺少这些哈希、来源提交不一致或平行光/CORSIKA 源归档不一致的数据。

### 7.1 `da51c09` 镜片仰角序列更新的实测影响

从上一轮基线 `f2b6617` 到本轮 `da51c09`，程序仓库中唯一变化的文件是 `configs/mirror_1229_elevation_series.csv`。逐镜片比较显示：不同仰角下中心位置最大变化为 `2.732–3.956 mm`，中心位置 RMS 变化为 `2.088–3.195 mm`，法向最大夹角变化为 `41.674–155.472 arcsec`。因此第 3 页必须重跑，不能只替换网页里的 CSV。

下表是同一批平行光输入、同一模拟流程的旧版 → 新版结果；质心和 RMS 都直接由程序 output-plane 原始 `u/v` 统计：

| 仰角 | `hit_output_plane` 旧→新 | 质心 `(u,v)` mm 旧→新 | 质心位移 | RMS 旧→新 |
|---:|---:|---:|---:|---:|
| 0° | 8,556→8,496 | `(-0.010,+0.065)`→`(-0.008,+0.133)` | 0.067 mm | 4.039→4.010 mm |
| 10° | 8,462→8,687 | `(-0.010,-0.073)`→`(-0.006,+0.273)` | 0.346 mm | 4.015→4.021 mm |
| 20° | 8,710→8,837 | `(-0.023,-0.248)`→`(+0.009,+0.450)` | 0.699 mm | 4.027→4.069 mm |
| 30° | 8,984→9,051 | `(-0.049,-0.471)`→`(+0.036,+0.674)` | 1.149 mm | 4.077→4.151 mm |
| 40° | 9,294→9,327 | `(-0.086,-0.804)`→`(+0.075,+1.009)` | 1.820 mm | 4.162→4.263 mm |
| 50° | 9,673→9,845 | `(-0.133,-1.158)`→`(+0.123,+1.365)` | 2.537 mm | 4.280→4.402 mm |
| 60° | 10,260→10,358 | `(-0.189,-1.702)`→`(+0.178,+1.908)` | 3.629 mm | 4.421→4.558 mm |
| 70° | 10,685→10,676 | `(-0.251,-2.089)`→`(+0.241,+2.297)` | 4.413 mm | 4.581→4.730 mm |
| 80° | 10,789→10,898 | `(-0.318,-2.596)`→`(+0.309,+2.805)` | 5.438 mm | 4.757→4.914 mm |
| 90° | 10,919→11,063 | `(-0.389,-3.061)`→`(+0.380,+3.272)` | 6.380 mm | 4.932→5.087 mm |

独立重跑后，第 2、3 页的 baseline、四方向与 10 个仰角结果在排除来源元数据后均与上一轮逐值一致，说明本次修复没有改变通用平行光入口。第 4 页的 event、pointing、array、原始 bunch、`zem` 反演点和空间尺度不变，但相机水平方向按修复后的 CORSIKA 入口重新追迹；不能再声称其 camera views 与旧版逐值一致。

四方向真实输出摘要如下；质心是 LACT 原始 `(u,v)`，不是网页重投影值：

| 光源 | `hit_mirror` | incoming / reflected 遮挡 | `hit_output_plane` | `hit_camera` | 原始质心 `(u,v)` m |
|---|---:|---:|---:|---:|---:|
| 上 `az=0°, el=71°` | 14,803 | 2,973 / 127 | 14,676 | 11,243 | `(0.000023, -0.143107)` |
| 下 `az=0°, el=69°` | 14,754 | 2,988 / 135 | 14,619 | 11,212 | `(0.000015, +0.143108)` |
| 左 `az=-1°, el=70°` | 14,748 | 2,952 / 116 | 14,632 | 9,100 | `(+0.048908, -0.000399)` |
| 右 `az=+1°, el=70°` | 14,795 | 2,904 / 121 | 14,674 | 9,133 | `(-0.048943, -0.000386)` |

网页运行时测试直接比较三页结构，而不是比较标签文字：全局/平行光理想镜片最大顶点误差为 `4.681018003915913e-10 m`；通用 NEU 与 CORSIKA NWU 先按 `East=-West` 还原成同一物理轴后，镜面中心归一的最大顶点误差为 `9.769163457519441e-11 m`；平行光/仰角页基底误差为 0；CORSIKA `az=0°/el=70°` 数值基底与最新版 C++ 公式误差为 0；同一向东偏 1° 光源跨入口的 local `u/v` 误差为 0；三维 x/y/z 等比例缩放测试通过。

右上角相机核对不使用抽样：四组分别嵌入全部物理 output-plane `u/v`，并按程序 `plot_whiteboard()` 的方式以原始包围盒取景、留 8% 边距、使用 `140×140` count/bin。LACT 模式显示未平移的绝对 `u/v`；pyLAST 模式对同一批点显示横轴 `+u`、纵轴 `-v`。具体行数必须从本轮最新版真实输出重新写入来源记录，不能沿用旧二进制结果。完整 camera-hit `x/y` 继续参与逐行和逐像素强校验。三维区绘制全部 `hit_mirror` 反射点以及全部 incoming/reflected 遮挡记录诊断端点；只有光路线允许抽样。

相机画布不从三维图重新投影落点。二维热图逐行读取每个 run 的完整 hits CSV `u_m/v_m` 并统计 count/bin；完整 `camera_x_m/camera_y_m` 与 camera CSV 像素计数继续参与生成期校验。生成期逐行验证 `camera_x_m==u_m、camera_y_m==v_m`，并用该 run 保存的 `buildTelescopeFrame` 基底把 global surface 投影回局部面复现同一 `u/v`。LACT 模式不做变换；pyLAST 模式只在画布入口应用 `horizontal=+u、vertical=-v`，不回写原始数组。

为保留遮挡诊断行，四方向配置显式设置 `obstruction.mark_only=true`。CSV 的 `obstruction_blocked_incoming/reflected` 是 C++ 原始标志；`mirror_point` 和 `surface_point` 也是 tracer 输出原值。程序没有输出光线与 obstruction primitive 的精确交点，因此页面将全部 incoming-blocked 行画成其理论镜面端点红叉，将全部 reflected-blocked 行画成其理论输出端点紫叉。它们是完整的遮挡记录诊断点集合，但不能冒充 primitive 的精确相交位置。

### 7.2 CORSIKA `u` 方向问题的版本时间线与实测影响

这不是最近一次网页配色或相机画法单独引入的问题，时间线如下：

| 时间 / 提交 | 程序或网页状态 | 结论 |
|---|---|---|
| 2026-07-27 19:20，程序 `f2b6617` | CORSIKA 本地轴由“`x=sky-up、y=East`”改成“`x=West、y=sky-up`” | 修正了 `x/y` 语义次序，但水平轴符号选反；相机 `u=local x` 从此与通用入口相反 |
| 2026-07-28 12:17，网页 `ae94c96` | 网页仍画旧的 `x=sky-up、y=East` | 网页与当时程序轴次序不一致 |
| 2026-07-28 13:25，网页 `693aaa9` | 网页改为 `x=West、y=sky-up` | 网页开始忠实复现程序，但也把程序的 `u` 反号完整显示出来 |
| 2026-07-28 23:00，程序 `4fb44f8` | CORSIKA 改为 `x=East、y=sky-up`，新增向东偏 1° 跨入口回归 | `u/v` 的定义与第 1–3 页统一；网页与真实输出同步重建 |

因此，“CORSIKA 与 1–3 页的 `u` 不一致”是 `f2b6617` 的最近一次坐标重构引入的程序问题；网页后来的变化只是先后经历了“轴次序画错”和“正确展示程序反号”两种状态，并不是网页在最后几次界面优化中把已有正确结果翻坏。

使用同一个 4 GB 输入、同一个 event 1909、同一相机像素几何比较修复前后 ROOT，相机 PE 加权质心如下。`v` 基本不变，`u` 精确呈镜像翻转，正是本次水平轴修复应有的特征；PE 总数的极小差异来自方形像素边界上的重新归属。

| 望远镜 | 旧版 `(u,v)` m | 修正版 `(u,v)` m | PE 总数 旧→新 |
|---|---:|---:|---:|
| tel19 | `(-0.141253, +0.080178)` | `(+0.141268, +0.080214)` | 51,438 → 51,417 |
| tel16 | `(+0.172656, +0.064227)` | `(-0.172621, +0.064210)` | 15,161 → 15,166 |
| tel21 | `(-0.376956, +0.055549)` | `(+0.377051, +0.055642)` | 9,444 → 9,445 |
| tel20 | `(-0.348826, -0.225677)` | `(+0.348708, -0.225621)` | 9,324 → 9,334 |

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
- 两套 telescope frame 都保持单位长度和正交；通用数值基满足 `x×y=z`，CORSIKA NWU 数值基因 `North-West-Up` 的左手物理轴顺序满足 `x×y=-z`。CORSIKA 的实际数值必须与最新版 C++ 公式一致，页面不得翻转或交换 `local x/y`。
- 平行光与 CORSIKA 在 `az=0°/el=70°` 下的 54 块镜片必须在各自镜面顶点归一后逐顶点一致，最大误差 `<1e-9 m`。
- pyLAST 选择器必须对同一批点实现 `horizontal=+u、vertical=-v`，且切换前后三维几何逐值不变。
