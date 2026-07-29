# LACT 坐标系代码审计

本笔记只依据可执行代码、真实运行输出和数值回归建立，不把项目说明文档或代码注释当作坐标定义。注释仅在表达式核对完成后用于交叉检查。

## 1. 审计基线

- LACT_sim 远端 `main`：`48ea63115ffff4dd840e792149f79a2549f3abb0`
- 坐标统一修复：`4fb44f876f37a7c920f0705abc3a8552533588e9`
- pyLAST 远端 `lact_sim`：`953d07392a542a3362a0ce35f14bcdbe4c347c9a`
- 2026-07-28 重新查询 GitHub 后，以上三个提交分别仍是对应远端分支的最新提交。
- 网页中的平行光和 CORSIKA 数据来自同一个 LACT_sim 源码归档；生成器会比较源码归档、二进制和 ROOT 的 SHA-256，不一致时拒绝构建。

代码位置均写成 `仓库@提交:文件:行`。本网站分支只保存网页、绘图数据、生成脚本和本笔记，不复制 LACT_sim / pyLAST 主程序源码。

## 2. 网页的固定物理方向

网页把物理方向统一显示为：北、东、南、西、天空、地面。俯视时固定为上北、下南、左西、右东。

四页默认观察相机统一取 `view_az=0°`：North 的观察深度为正（背离屏幕、屏幕内），South 的观察深度为负（朝向屏幕、屏幕外）。第 1、3 页同时取 `view_el=0°`，因此地面严格投影成一条水平地平线；第 2、4 页只增加负的观察仰角以保留结构深度，不改变南/北的前后关系。这里的 `view_az/view_el` 只控制网页观察相机，不是望远镜 pointing。

CORSIKA 阵列使用 NWU 分量：

- `+X = North`
- `+Y = West`
- `+Z = Up / Sky`
- 因此物理 East 是 `−Y`。

通用光学全局框架的水平分量按 North/East 展示：

- `+X = North`
- `+Y = East`
- `+Z = Up / Sky`

网页不会假装这两组数组分量相同。它们只共享同一组物理方向标签；从 CORSIKA NWU 到屏幕物理 East 时，必须显示 `East = −West`。

## 3. 通用望远镜框架

代码在 `LACT_sim@4fb44f8:src/app/OpticalSimCommon.cpp:555-581` 构造：

```text
z = (cos(el) cos(az), cos(el) sin(az), sin(el))
x = (-sin(az), cos(az), 0)
y = z × x
```

局部向量与点的正反变换由 `OpticalSimCommon.cpp:206-223` 实现：

```text
global_vector = x_axis*local.x + y_axis*local.y + z_axis*local.z
global_point  = origin + global_vector
local.x/y/z  = dot(global-origin, x_axis/y_axis/z_axis)
```

数值含义：当 `az=0°` 时，光轴 `local +z` 朝北，`local +x` 朝东，`local +y` 是仰角增加方向。网页第 1、2、3 页的望远镜结构、镜片、遮挡和相机面均使用这一套基底展开；没有为某一页另写镜片朝向。

## 4. CORSIKA NWU 到望远镜局部框架

代码在 `LACT_sim@4fb44f8:src/app/OpticalSimCommon.cpp:584-609` 构造：

```text
x_local_in_NWU = (-sin(A), -cos(A), 0)
y_local_in_NWU = (-sin(E) cos(A), sin(E) sin(A), cos(E))
z_local_in_NWU = ( cos(E) cos(A),-cos(E) sin(A), sin(E))
```

当 `A=0°` 时：

- `local +z` 指北并具有仰角 `E`；
- `local +x=(0,-1,0)`，在 NWU 中就是物理 East；
- `local +y` 是天空/仰角增加方向。

因此通用入口与 CORSIKA 入口的物理光学局部轴一致。`apps/test_coordinate_frames.cpp:40-85` 对“向东偏移产生相同 local u/v 符号”做了数值回归。

输入适配由 `OpticalSimCommon.cpp:691-750` 决定：

- `telescope_local`：位置不变，方向只归一化；
- `corsika_nwu_relative`：把相对望远镜的 NWU 位置、方向投影到上述局部基底；
- `corsika_nwu_global`：先减望远镜原点，再投影；
- ENU 两种入口先做 `East→−West` 和方位角换算，再进入同一 NWU 物理框架；
- `lact_generic_global`：使用第 3 节的通用框架做逆变换。

`OpticalSimCommon.cpp:752-760` 的 EventIO reference offset 只加局部 `z`，不修改局部 `x/y`。

## 5. 两条完整光学执行链

### 5.1 `run_optical_sim`

入口在 `LACT_sim@4fb44f8:apps/run_optical_sim.cpp`：

1. `:33` 构造通用望远镜框架；
2. `:37` 应用支架形变；
3. `:39` 构造输出面；
4. `:47` 把镜片和输出面从望远镜局部展开到通用全局；
5. `:414-424` 把光子先适配为望远镜局部，再展开到同一通用全局；
6. 遮挡检查把全局线段用同一个框架投回局部，实际表达式见 `OpticalSimCommon.cpp:3376,3394-3395,3412,3470-3471`；
7. `run_optical_sim.cpp:484` 把追迹命中的原始 `u/v` 交给相机响应。

所以网页第 2、3 页不能套用 CORSIKA NWU 基底去摆放光学结构；它们必须使用通用框架。

第 2 页四组 3° 天区输入均由同一程序公式换到望远镜局部传播方向，随后运行最新 `run_optical_sim`，不是在网页中移动旧光斑：

| 名称 | 天区输入 | `source.beam_direction` 原值 | 局部 `θ / φ` |
|---|---|---|---|
| 上 | `az=0°, el=73°` | `(0, -0.0523359562, -0.9986295348)` | `3°, 270°` |
| 下 | `az=0°, el=67°` | `(0, +0.0523359562, -0.9986295348)` | `3°, 90°` |
| 西 | `az=-3°, el=70°` | `(+0.0178999513, -0.0004404590, -0.9998396860)` | `1.0259569321°, 358.5904234021°` |
| 东 | `az=+3°, el=70°` | `(-0.0178999513, -0.0004404590, -0.9998396860)` | `1.0259569321°, 181.4095765979°` |

西/东的局部 `θ` 小于 3° 是球面几何的正常结果：在仰角 70° 处，方位角差投影到切平面会乘以约 `cos(70°)`；天区输入方位角本身仍是完整的 `±3°`。

### 5.2 `run_corsika_trace`

入口在 `LACT_sim@4fb44f8:apps/run_corsika_trace.cpp`：

1. `:168-170` 把 EventIO 光子转换为望远镜局部并应用 reference-z；
2. `:3358-3366` 应用形变、构造输出面；CORSIKA 模式构造 NWU 望远镜框架只用于输入/物理放置解释；
3. 镜片、遮挡和输出面本身在望远镜局部追迹；不会先套通用全局框架；
4. `:3946` 使用追迹产生的原始 `u/v` 做相机响应。

网页第 4 页把局部望远镜结构放到 NWU 阵列地图时，使用第 4 节的程序基底；这只是三维位置表达，不修改相机输出数组。

## 6. 镜片与遮挡的 x/y

理想镜片 CSV 的 `center_x/y/z` 和 `normal_x/y/z` 由 `OpticalSimCommon.cpp:1870-1910` 读取为望远镜局部量。镜片孔径轴由 `include/geometry/MirrorFacet.hpp:61-71` 确定：若 CSV 明确给轴就直接使用；否则由法向量和参考轴构造正交 `u/v`。

遮挡 primitive 同样以望远镜局部坐标参与相交测试。网页把镜片、支架、相机和遮挡作为同一个局部整体变换，禁止分别交换 x/y。

## 7. 输出面 u/v 与相机 x/y

输出面默认基底由 `LACT_sim@4fb44f8:include/optics/OutputPlane.hpp:21-30` 构造：

```text
u_axis = normalize(reference × normal)
v_axis = normal × u_axis
```

若配置显式给 `output.plane_u_axis/v_axis`，`OpticalSimCommon.cpp:2187-2205` 检查它们与法向量及彼此正交后直接采用。

真正写入命中的值在 `src/optics/OpticalTracer.cpp:179-184`：

```text
relative = surface_point - plane.point
u_m = dot(relative, plane.u_axis)
v_m = dot(relative, plane.v_axis)
```

相机边界在 `OpticalSimCommon.cpp:2724-2726`：

```text
camera_x_m = hit.u_m
camera_y_m = hit.v_m
pixel = findContainingPixel(camera, camera_x_m, camera_y_m)
```

结论：LACT 相机白板横轴直接画输出 `u_m`，纵轴直接画输出 `v_m`。网页允许按原始包围盒选择显示窗口、做等比例缩放和 count/bin 分箱，但不允许对数据换轴、翻符号、重设零点或旋转。

## 8. ROOT 到 pyLAST 的真实边界

LACT_sim ROOT writer 在 `src/io/LactEventRootWriter.cpp:775-788` 把相机几何写成：

```text
camera_pixels.x_m = pixel.center.x = u
camera_pixels.y_m = pixel.center.y = v
```

`observations.image_cherenkov_pe` 的填充与写入见 `LactEventRootWriter.cpp:301-496,893,1028`。

最新 pyLAST reader 在 `pyLAST@953d073:root/LactEventSource.cpp:245-265` 执行：

```text
pix_x = -root_camera_pixels.y_m = -v
pix_y = +root_camera_pixels.x_m = +u
```

`root/LactEventSource.cpp:696-697` 对 `image_cherenkov_pe` 做最近整数舍入后存入 `true_image`。`src/pylast/visualize/event_visualizer.py:393-401,2216-2218` 再把 `pix_y` 放到画布横轴、`pix_x` 放到画布纵轴。

网页现在按这两个可执行边界在 Canvas 中逐点重画，而不嵌入原始 PNG：

```text
reader:    (u, v) -> (pix_x, pix_y) = (-v, +u)
renderer:  (pix_x, pix_y) -> (plot_x, plot_y) = (pix_y, pix_x)
combined:  (plot_x, plot_y) = (+u, -v)
```

这里的 `(+u,-v)` 不是网页为“看起来一致”增加的翻转，而是严格依次执行 reader 与 `EventVisualizer._camera_to_plot_xy()` 的结果。实现中保留了 `pylastReaderCoordinates()` 和 `pylastPlotCoordinates()` 两个独立函数，运行时测试分别核对，避免把两层边界压成无出处的显示技巧。

像素颜色也模仿最新版 `EventVisualizer._image_norm()`：零值像素为白色，非零像素使用 plasma 色带；输入 PE 先按 `LactEventSource.cpp:696-697` 的 `round()` 规则形成 `true_image`。网页只重画几何和着色，不更改用于绘制的原始 LACT/ROOT 输出值。

本次数值审计结果：

- ROOT SHA-256：`36876c1e586c3c6dcc1ff5d2f7a49cc66430a0cfce3db019454dde3684bf29bc`
- pyLAST 坐标最大误差：`0.0 m`
- 全部 32 台望远镜的 `true_image` 都直接读取同一 ROOT 的 `image_cherenkov_pe`；网页没有按信号排序后丢弃其余望远镜。该事例的值本来就是整数，所以 pyLAST 最近整数规则的最大舍入差是 `0 PE`。
- 第 2 页四个平行光 pyLAST 面板读取各自 `run_optical_sim` 的真实像素 CSV；第 4 页读取同一 ROOT event 1909 的 `image_cherenkov_pe`。
- 切换 LACT/pyLAST 只改变相机 Canvas 的坐标表达和配色，不修改三维望远镜、光路或 CORSIKA 发射点。

由于网页自行重画，标题直接使用 ROOT/输出中的真实 `telescope_id`，不会继承 `event_visualizer.py:1950-1958` 的 `tel_id + 1` 显示偏移。

## 9. EventIO 光子、芯位和发射点

二维 photon bunch 在 `LACT_sim@4fb44f8:src/io/EventIOPhotonSource.cpp:570-590` 解码：

```text
position_m = (x*0.01, y*0.01, 0)
direction  = (cx, cy, negative_z_from_direction_cosines)
zem_km     = zem*1e-5
```

三维 bunch 在 `EventIOPhotonSource.cpp:648-667` 读取三维位置、方向和 emission altitude。阵列方位角的可执行换算在 `EventIOPhotonSource.cpp:309`；ROOT 芯位和阵列位置分别由 `LactEventRootWriter.cpp:83-93,818-843` 写成 North/West/Up 原值。

第 4 页的彩色“发射点”不是 EventIO 直接保存的三维出射点，而是诊断派生量。对二维 bunch，网页数据准备程序使用原始接收平面 anchor、原始方向和原始 `zem` 解直线与该高度面的交点：

```text
h = zem_asl - observation_altitude - telescope_z
s = (h - anchor_z) / (-direction_z)
emission = anchor - s*direction
```

该派生只用于显示簇射空间范围；相机图不使用这个反推点。网页校验方向模长、回代共线误差和回算高度，且 x/y/z 使用同一空间 scale。

## 10. 结构形变的“理想值”

形变序列在 `OpticalSimCommon.cpp:1530-1790` 读取并按仰角插值：中心线性插值，法向量用球面插值。`applyStructuralDeformation()` 在 `OpticalSimCommon.cpp:3476-3504` 用当前仰角得到的 `center` 和 `normal` 直接替换理想镜片值；CSV 保存的是该仰角的绝对镜片状态，不是 delta。

网页的两种颜色是诊断派生：

- 指向偏转：`acos(normal_ideal · normal_deformed)`；
- 中心位移：`center_deformed − center_ideal` 的模。

这里的 `ideal` 明确指同一次构建读取的 `configs/mirror_1229_facets.csv`，不是零度形变文件，也不是相机 R68。所有仰角共用一个颜色 scale。

## 11. 四页数据来源和允许操作

| 页面 | 三维结构 | 光路/数据 | 相机图 |
|---|---|---|---|
| 1 全局定义 | 从镜片、遮挡、相机 CSV 原值按程序通用基底放置 | 不放事例 | 不放事例；只标 LACT 与 pyLAST 字段方向 |
| 2 平行光 | 通用基底，az=0° / el=70° | 四组最新 C++ 输出：上 `el=73°`、下 `el=67°`、西 `az=-3°`、东 `az=+3°`；每次显示所选组的全量反射点、遮挡端点和抽样光路 | LACT：完整 `output u/v` 原值分箱；pyLAST：同组真实像素输出按 reader/renderer 规则重画 |
| 3 天顶角 | 每个仰角的绝对形变镜片状态 | 对应仰角真实 C++ 平行光输出 | 对应仰角完整 `output u/v` 原值分箱 |
| 4 CORSIKA | ROOT 阵列 NWU 原值；32 台均按各自 ROOT pointing 绘制镜片、支架/遮挡和相机外框；所选台另绘完整 1616 像素 | 32 台 event 1909 原始 anchor/direction/zem；发射点为注明公式的派生诊断 | 右上角选择 32 台中的任一台；LACT 直接用 ROOT `image_cherenkov_pe`，pyLAST 按同一数据与最新版代码规则重画 |

允许的网页显示操作：观察相机旋转/平移、等比例缩放、投影、抽样光路线、按原值窗口分箱。

禁止的数据操作：为“看起来一致”任意交换 `u/v`、改变正负号、把光斑质心当成 `(0,0)`、让切换相机坐标改变三维望远镜结构。pyLAST 模式唯一允许的换轴/符号变化是第 8 节逐行对应最新版 reader 与 renderer 的确定性映射。

## 12. 构建时自动拒绝的不一致

`python/build_coordinate_diagnostics_html.py` 在生成网页前检查：

- 平行光和 CORSIKA 是否来自同一 LACT_sim 基线与源码归档；
- 每个平行光事例是否包含全部输出面 `u/v` 和全部相机命中 `u/v`；
- `surface_point` 投影是否逐条复现输出 `u/v`；
- `camera_x/y` 是否逐条等于 `u/v`，像素计数之和是否等于 C++ `hit_camera`；
- 四方向平行光是否严格为 `el=73/67°` 和 `az=-3/+3°`，且配置方向向量与代码公式一致；
- pyLAST reader/renderer 的两步映射是否对给定 `u/v` 产生精确的 `pix_x=-v, pix_y=+u, plot=(+u,-v)`。

`python/test_coordinate_diagnostics_runtime.js` 还会数值比较第 1/2/3/4 页镜片顶点方向、固定地图、四页默认南朝屏幕外/北朝屏幕内、地平线投影、通用/CORSIKA 东向偏移符号、CORSIKA 角度换算、三维等比例 scale，并确认平行光和 CORSIKA 切换到 pyLAST 重画时不会修改三维结构。

## 13. cfg 原值与 CORSIKA 望远镜选择

数据准备脚本现在解析 cfg 中未注释的 `key=value`，将原始字符串作为 `provenance.run_config_values` 写入 JSON。网页“说明 / 原值 / 校验”直接读取这些字符串，不根据网页几何反推配置。

平行光四方向的共同配置是：

```text
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70
source.n_bunches=30000
obstruction.mark_only=true
output.mode=both
output.whiteboard_input_photon=true
```

东西两例的关键 cfg 原值为：

```text
East sky az=+3°, el=70°:
source.beam_direction=-0.017899951255297582,-0.00044045903963302324,-0.9998396860201599

West sky az=-3°, el=70°:
source.beam_direction=+0.017899951255297582,-0.00044045903963302324,-0.9998396860201599
```

这里的 East/West 名称来自 cfg 注释记录的天空方位和代码审核得到的物理方向；网页三维光线直接读取运行输出中的方向列，不重新用名称生成向量。

CORSIKA 数据包现在包含 32 个 ROOT observation 相机和 32 个 EventIO photon group，每台均匀抽样 120 条，仅用于三维线段显示，共 3840 条。完整 photon bunch 数仍保留在逐台 `raw_bunch_count` 中；相机像素不抽样。完整簇射/地面视图绘制全部 32 台结构。右上角选择望远镜后：

1. 三维场景只保留所选台的详细结构、完整相机像素和该台抽样入射线；
2. `fit()` 只接收以所选 ROOT 位置为中心的局部点集，因此旋转中心随望远镜切换；
3. 局部入射线是原始 `anchor_array_nwu_m` 和 `direction_nwu` 所定义直线在望远镜附近的截取，不是新增光学追迹；
4. 右侧相机同时切到同一 `telescope_id` 的 ROOT `image_cherenkov_pe`。

运行时测试显式检查：阵列场景含 `32 × 54 = 1728` 个镜片面、选择器含 32 项、选择 tel0 后场景只含 tel0 的 54 个详细镜片，并且三维拟合中心的 North/West 坐标等于 tel0 的 ROOT 位置。
