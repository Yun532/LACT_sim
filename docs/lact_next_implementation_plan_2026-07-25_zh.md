# LACT_sim 下一步源码修改说明

> 状态说明：本文是基于修改前基线 `0048369` 写出的实施设计。当前实际落地状态、
> 配置语义和验证进度见
> `lact_trustworthy_optics_implementation_2026-07-25_zh.md`。

> 基线：GitHub 远端 `main`
> `0048369cd06e87322d9ef3611a4a369a1152fe99`
> 原则：只修复 LACT 自身的正确性、内部自洽性、数值安全和可验证性；不要求实现方式与 sim_telarray 一致。

## 1. 当前程序真实执行流程

### 1.1 配置和静态光学对象

`apps/run_corsika_trace.cpp::main()` 首先：

1. 读取主配置并展开 telescope、mirror、camera、SiPM、efficiency、atmosphere、
   error、obstruction、trigger 等组件配置；
2. 创建望远镜配置；
3. 创建 nominal mirror facets；
4. 应用结构形变和逐镜效率；
5. 创建输出平面；
6. 创建相机像素和 collector；
7. 创建带随机误差的逐望远镜镜面缓存；
8. 创建效率、大气、光速、PE response 和输出对象。

对应代码：

```text
run_corsika_trace.cpp:3289-3368
```

重要调用关系：

```text
readKeyValueConfig
  -> expandConfig
  -> buildTelescopeConfig
  -> buildFacetsFromConfig
  -> applyStructuralDeformation
  -> applyFacetEfficiencyScales
  -> buildOutputPlane
  -> buildCameraGeometry
  -> buildLightCollector
  -> TelescopeOpticsCache
  -> buildEfficiencyConfig
  -> buildAtmosphereTransmissionConfig
  -> OpticalTracer
```

### 1.2 EventIO 解码

二维 bunch 在 `EventIOPhotonSource.cpp::makeBunch()` 中转换：

```text
x,y: cm -> m
z:   暂设 0
dir: (cx,cy,-sqrt(1-cx²-cy²))
time: ctime ns
weight: source default
multiplicity: b.photons × default multiplicity
```

三维 bunch 在 `makeBunch3d()` 中转换：

```text
x,y,z: cm -> m
dir: (cx,cy,cz)
time/multiplicity/wavelength: 与二维相同语义
```

二维和三维此时都仍在输入数据的坐标语义中，还没有完成 LACT 光学原点转换。

对应代码：

```text
EventIOPhotonSource.cpp:544-571
EventIOPhotonSource.cpp:598-624
```

### 1.3 输入坐标到 telescope-local

每个 bunch 进入：

```text
transformEventIOBunchToTraceFrame
  -> transformBunchToTelescopeLocal
```

`transformBunchToTelescopeLocal()` 根据 `source.coordinate_frame` 区分：

- 已经是 telescope-local；
- CORSIKA NWU 相对坐标；
- CORSIKA NWU 全局坐标；
- ENU 相对/全局坐标；
- 旧 LACT generic global。

相对坐标只旋转向量；全局坐标先减望远镜位置再旋转。方向只做旋转，不做平移。

对应代码：

```text
OpticalSimCommon.cpp:643-758
run_corsika_trace.cpp:142-174
```

### 1.4 `z=0 -> z=-16`

输入 frame 转换完成后，代码执行：

```cpp
if (bunch.eventio_2d) {
    bunch.photon.pos.z += eventio_2d_input_plane_z_m;
}
```

当前默认值为 `-16 m`。

它的物理意义不是把光子传播 16 m，而是把 EventIO fiducial/reference 原点表达为
LACT 导入几何的光学原点。因为只是换坐标表示：

- x/y 不需要沿光线传播；
- `ctime` 不应改变；
- 镜面约在 `z=-16 m`；
- 相机约在 `z=-8 m`。

二维使用 full-line/backprojection 与这个原点平移是两个独立概念：

- 原点平移回答“这个物理点在 LACT 坐标里是多少”；
- full-line 回答“二维记录点只是光线锚点时，怎样稳定找到镜面”。

### 1.5 PE 候选和预几何效率

`PhotonResponseSampler` 有两个模式：

#### expectation

一个 bunch 只追一条代表光线：

```text
photon.weight *= multiplicity
```

大气和光学效率继续作为权重相乘。

#### stochastic_pe

根据 multiplicity 展开 represented photon，使用：

```text
shower_event_id
array_id
telescope_id
source_bunch_index
represented_index
```

生成稳定的 `stream_id`，先随机决定波长、大气和预几何效率是否存活。

对应代码：

```text
PhotonResponseSampler.cpp:90-113
PhotonResponseSampler.cpp:121-215
```

### 1.6 主镜追迹

`OpticalTracer::traceToPlaneImpl()`：

1. 遍历镜片；
2. 按镜片类型做平面、球面或抛物面解析求交；
3. 做真实口径裁剪；
4. 把法线定向到镜片前表面；
5. 拒绝背面入射；
6. 选择光线参数顺序中的第一块有效镜片；
7. 计算镜面反射方向；
8. 可选加入反射方向散射；
9. 与输出平面求交；
10. 计算 `u/v`、到达时间和相对效率。

对应代码：

```text
OpticalTracer.cpp:111-207
OpticalTracer.cpp:210-445
```

命中时间当前为：

```text
photon.time + (mirror_intersection_t + mirror_to_plane_t) / c
```

对于二维 full-line，第一段是相对 EventIO 时间锚点的有符号时间修正；第二段是镜面到焦面的正向传播时间。

### 1.7 遮挡

光追先得到不含结构遮挡的镜面和焦面交点，随后主程序分别检查：

```text
物理上游入射射线 -> 镜面
镜面 -> 输出平面有限线段
```

对应代码：

```text
run_corsika_trace.cpp:3798-3857
OpticalSimCommon.cpp:3060-3168
```

### 1.8 相机和 collector

到达输出平面后：

1. 用 `u/v` 查找包含该点的 pixel；
2. 没有 collector 时，pixel containment 就是相机接受；
3. 有 collector 时，把位置和方向转换到单个 square cone 局部坐标；
4. 在四个壁面和出入口之间迭代反射；
5. 检查出口点是否落到 SiPM 有效面积；
6. 把 collector intensity 乘进相对效率。

对应代码：

```text
OpticalSimCommon.cpp:2306-2452
LightCollectorSquareCone.hpp:213-302
LightCollectorSquareCone.hpp:521-648
```

### 1.9 后几何接受、累计和输出

stochastic 模式在几何完成后，再根据剩余的 collector/逐镜等效率随机接受。

随后：

```text
pixel PE 累积
  -> waveform proxy 时间分箱
  -> 可选 NSB
  -> 简化 camera/array trigger
  -> CSV/HDF5/ROOT
```

当前电子学 metadata 明确写的是：

```text
integrated_pe_placeholder
```

这不是完整模拟波形电子学。

## 2. 第一项必须修改：重构 Bezier collector 求交

### 2.1 当前错误位置

Bezier 壁面求交先得到曲线参数 root 和交点：

```text
LightCollectorSquareCone.hpp:225-267
```

但是返回值只有：

```cpp
std::pair<Position, bool>
```

曲线参数 root 被丢掉。

之后为了求法线，`get_normal_vector(pos)` 又分别从交点的 x 和 z 反解二次方程，并尝试在四个根中找到两个近似相等的值：

```cpp
auto it = std::find_if(roots.begin(), roots.end(), ...);
auto t = *it;
```

位置：

```text
LightCollectorSquareCone.hpp:270-299
```

当没有找到近似相等根时，`it==roots.end()`，第 291 行仍解引用，ASan 已经确认发生 heap-buffer-overflow。

### 2.2 另一个算法错误

四个 cone 壁面选择最近交点时使用：

```cpp
return d1 < d2 || fabs(d1 - d2) < 1e-6;
```

位置：

```text
LightCollectorSquareCone.hpp:549-556
```

当两个距离很接近时：

```text
comp(a,b) == true
comp(b,a) == true
```

这违反 `std::min_element` comparator 的严格弱序要求。

### 2.3 不建议只做的补丁

只增加：

```cpp
if (it == roots.end()) return ...;
```

虽然可以避免越界，但没有解决“求交时已经知道参数，求法线时却再次不稳定反解”的根本设计问题。

### 2.4 建议修改结构

让每个 Surface 一次返回完整交点：

```cpp
struct SurfaceHit {
    bool hit = false;
    double ray_t = 0.0;
    double surface_u = 0.0;  // Bezier parameter
    Position point;
    DirectionVecter normal;
    std::size_t surface_index = 0;
};
```

Bezier 求交得到 `u` 后立即计算：

```text
x(u), z(u)
dx/du, dz/du
```

Bezier cylindrical surface 可以参数化为：

```text
S(u,y) = (x(u), y, z(u))
```

两个切向量为：

```text
dS/du = (dx/du, 0, dz/du)
dS/dy = (0, 1, 0)
```

法线由它们的叉积得到。这样交点和法线共享同一个 `u`，不再反解。

最近交点使用：

```cpp
if (candidate.ray_t < best.ray_t) ...
```

距离近似相等时使用明确、确定性的 seam/corner 规则，例如按壁面 index 打破平局；不能把 epsilon 写进 `<` comparator。

### 2.5 collector 同时应补的时间

当前 collector 返回：

```text
intensity, exits, exit_position, exit_direction, reflections
```

没有返回 collector 内部累计光程。因此 `hit.time_ns` 停留在 collector 入口/输出平面的到达时间，反射后更长的 collector 光程没有计入 PE 时间。

建议 `SurfaceHit` 保存每段 `ray_t`，`ray_trace_impl()` 累加：

```text
collector_path_length_mm
```

然后：

```text
hit.time_ns += collector_path_length_m / propagation_speed
```

同时把该长度写入调试输出。29.7 mm 直达路径约为 0.1 ns，多次反射会更长；如果做亚纳秒定时，这不能长期省略。

### 2.6 collector 其他边界

一并明确：

- `number <= 500` 当前最多可执行到 501 次，改成显式 `max_reflections` 和 `<`；
- 超过最大反射次数返回明确状态，不能与普通“未命中出口”混为一类；
- 无交点、退化法线、NaN、近平行和切线都返回失败状态；
- corner 上的法线没有唯一物理定义，当前平均相邻壁面法线可以作为模型，但必须写成明确策略并测试；
- 如果 transform 不保证纯旋转，法线必须用 inverse-transpose 变换。

### 2.7 验收标准

至少满足：

1. `test_light_collector_response` 在 ASan/UBSan 下通过；
2. 入口中心直达；
3. 单次和多次反射；
4. 四个壁面旋转对称；
5. 对角 seam/corner；
6. 切线和近平行；
7. 无实根和退化根；
8. 最大反射次数；
9. 随机位置/角度 fuzz 不崩溃、不产生 NaN；
10. collector 光程和时间随反射路径增长。

## 3. 第二项必须修改：一般化 EventIO 光学原点

### 3.1 当前问题

当前字段是：

```cpp
double eventio_2d_input_plane_z_m = -16.0;
```

当前转换只在：

```cpp
if (bunch.eventio_2d)
```

成立时应用。

但 EventIO 3D 的 `(x,y,z)` 也相对同一个 telescope/fiducial center。3D 读取后只是 cm→m，没有把该原点映射到 LACT 的 `z=-16 m`。

所以目前：

```text
2D EventIO fiducial origin -> LACT z=-16
3D EventIO fiducial origin -> LACT z=0
```

这不是两种合理模型，而是同一个输入参考原点被解释了两次。

### 3.2 建议的新配置

引入 telescope-local 三维平移：

```cpp
Vec3 eventio_reference_to_optical_offset_m{0.0, 0.0, -16.0};
```

转换顺序：

```text
输入坐标 -> telescope-local 旋转/去望远镜位置
         -> 加 eventio_reference_to_optical_offset_m
```

平移要在变到 telescope-local 后应用，因为该 offset 是 LACT 光学本地坐标中的量。

二维和三维都应用同一平移。

旧配置：

```text
source.eventio_2d_input_plane_z_m
```

可以暂时作为兼容别名映射到新 offset 的 z，并输出 deprecated warning。

不建议从镜面 `z_min` 自动猜 offset。原点关系属于输入/几何契约，应在配置中明确，不能依赖某个镜面恰好放在哪里。

### 3.3 不要混淆 3D 和 full-line

修复 3D 原点不意味着让 3D 自动使用二维的 full-line 模式。

- 二维没有真实 z，使用 line anchor/backprojection；
- 三维有完整到达位置，应继续按它声明的传播方向追 forward ray。

要改的是共同的坐标原点，不是共同的求交方向。

### 3.4 验收标准

新增测试：

1. 2D `(0,0)` 仍映射到 `z=-16`，保持现有结果；
2. 3D `(0,0,0)` 映射到 `z=-16`；
3. 3D `(x,y,z)` 映射到 `(x,y,z-16)`；
4. 相对和全局 CORSIKA frame 都正确；
5. 斜入射 3D 光线命中解析已知镜片；
6. 只整体平移 source、mirror、plane 时局部交点不变；
7. 现有全 2D EventIO 样例统计不变。

## 4. 第三项必须修改：统一稳定的随机流身份

### 4.1 当前散射问题

`OpticalTracer.cpp:168-176` 把以下量混进散射 seed：

```text
facet ID
photon.pos.x/y/z
photon.dir.x/y/z
photon.random_stream_id
```

因此只改变坐标原点，`photon.pos` 就变化，粗糙度散射抽样也变化。

此外 `hashDouble()` 使用 `std::hash<double>`，它不是严格的跨标准库可复现协议。

### 4.2 当前已有的好基础

`PhotonResponseSampler::photonResponseStreamId()` 已经使用：

```text
shower + array + telescope + source_bunch + represented_index
```

这比使用浮点坐标作为身份更合理。

但现在只有 stochastic 分支在：

```cpp
out.photon.random_stream_id = out.stream_id;
```

expectation 分支在这之前就返回，因此 `random_stream_id` 保持默认 0。

### 4.3 建议修改

建立一个公共、固定算法的随机流派生工具，例如基于明确实现的 SplitMix64：

```text
deriveSeed(base_seed, physical_stream_id, stage, optional_object_key)
```

stage 使用固定枚举：

```text
MissingWavelength
PreGeometryAcceptance
MirrorRoughness
PostGeometryAcceptance
FacetPositionError
FacetNormalError
FacetRadiusError
FacetReflectivityError
Nsb
```

在 `PhotonResponseSampler::candidate()` 中，无论 expectation 还是 stochastic，都先生成：

```text
stream_id
photon.random_stream_id
```

镜面散射 seed 改成：

```text
error/random base seed
+ photon.random_stream_id
+ persistent facet identity
+ MirrorRoughness stage
```

移除绝对位置、方向和 `std::hash<double>`。

如果希望显示 ID 可以重编号而随机 realization 仍完全不变，需要额外的持久
`facet.random_key`；否则应明确 facet ID 就是物理镜片的稳定身份。不同程序之间不需要共享这个 key。

### 4.4 镜片误差的顺序依赖

`applyFacetErrors()` 当前用一个 RNG 顺序遍历 `facets`：

```text
OpticalSimCommon.cpp:3257-3294
```

同一镜片集合只重排 CSV 行，就会把随机误差分给不同镜片。

应为每个：

```text
telescope ID + persistent facet identity + error stage
```

派生独立 seed。不同误差项使用不同 stage，避免增加一个新误差项后改变所有旧误差抽样。

### 4.5 验收标准

1. scatter 开启时整体平移不改变局部结果；
2. 整体刚体旋转后结果只做相同旋转；
3. 镜片 CSV 行重排不改变固定 seed realization；
4. expectation 中不同 source bunch 有不同但稳定的散射流；
5. stochastic 中 represented photon 的随机流互不相同；
6. 单线程、并行和不同容器遍历顺序结果一致；
7. 明确是否要求跨编译器 bitwise 一致，并用固定 seed golden vector 测试。

## 5. 第四项：补齐相机和 collector 配置契约

### 5.1 当前 loader 的缺口

镜片有独立的 `validateMirrorFacets()`，会检查重复 ID、有限值、非零法线、尺寸和曲率。

相机 `readCameraCsv()` 只是：

```text
stoi/stod
camera.addPixel(pixel)
```

没有检查：

- pixel ID 是否非负、唯一；
- x/y/size 是否有限；
- size 是否大于 0；
- 像素是否重叠；
- collector 是否与 pixel shape 兼容。

负 ID 与 `-1` 未命中 sentinel 冲突。重复 ID 会使 ROOT/HDF5 的
`pixel_id -> column` map 覆盖已有列。

### 5.2 单一 collector 尺寸假设

`buildLightCollector()` 只创建一个 cone。

如果没有显式 `collector_entrance_size_m`，入口尺寸来自：

```cpp
camera.pixels().front().size
```

这个 cone 随后用于所有像素。

当前正式相机的 1,616 个像素都是 `0.0244 m`，所以当前数据没有问题；但通用代码隐含假设所有像素尺寸相同。

处理方式二选一：

1. 当前阶段明确要求 collector-enabled camera 所有像素 shape、size 相同，否则报错；
2. 未来支持不同像素时，按 collector geometry 建 cache，每个像素选择对应 cone。

第一种更适合现在，简单且不会偷偷算错。

### 5.3 analytic funnel 与真实 collector

`OpticalEfficiency` 还可以启用 `funnel_acceptance=cos(theta)`，而相机又可以启用真实 collector。

当前正式配置把 analytic funnel 关闭，所以没有双算；但配置系统应拒绝或明确警告：

```text
analytic funnel acceptance + explicit ray-traced collector
```

同时启用，否则可能把 collector 接受算两次。

### 5.4 验收标准

- 当前正式相机 CSV 通过；
- 重复/负 ID 拒绝；
- NaN/Inf/非正尺寸拒绝；
- collector 开启时不同 pixel size 默认拒绝；
- pixel CSV 行重排不改变 ID 命中；
- 重叠像素要么拒绝，要么有文档化且确定性的归属规则；
- analytic funnel 与显式 collector 的组合有明确行为。

## 6. 第五项：效率和输入概率必须 fail-fast

### 6.1 当前效率曲线没有范围校验

`EfficiencyCurve::loadCsv()` 接受任意数值，没有检查：

- wavelength 是否有限且为正；
- efficiency 是否有限；
- 物理概率是否在 `[0,1]`。

`parseEfficiencyFactor()` 的 constant 也没有相同校验。

stochastic 路径的 `probability()` 会静默 clamp 到 `[0,1]`；expectation 路径则可能保留负数或大于 1 的效率。

这会让同一个非法配置在两种 response mode 下产生不同含义。

### 6.2 建议区分概率和标度

以下字段应是物理概率，要求 `[0,1]`：

```text
mirror reflectivity curve
filter transmission
atmosphere transmission
SiPM PDE
collector reflectivity
```

以下字段如果允许大于 1，应明确叫 calibration/weight scale，而不能伪装成概率：

```text
constant_scale
per-facet relative scale
```

即使单项 scale 可大于 1，进入 stochastic 接受前的最终概率仍必须在 `[0,1]`。
超过范围应报错，不能静默 clamp。

### 6.3 Photon CSV/EventIO 输入校验

同时补：

- 位置、方向、时间、波长和权重有限；
- 方向非零；
- weight、multiplicity 非负；
- 2D `cx²+cy²` 只允许浮点容差范围内略大于 1；
- 3D `cx,cy,cz` 有限且非零；
- expectation 与 stochastic 使用相同的基础合法性检查。

当前 `downwardDirZ()` 对任何 `cx²+cy²>1` 都夹成 `cz=0`，严重坏数据可能被静默变成水平光。

## 7. 第六项：消除大气双算歧义

当前同时存在：

1. `OpticalEfficiencyConfig::atmosphere_transmission`：只依赖波长的效率因子；
2. `AtmosphereTransmission`：依赖波长、发射高度和方向的 MODTRAN/tau 模型。

主流程把它们相乘。这个数学行为没有隐藏，但物理语义容易误配成同一大气衰减两次。

建议：

- 给第一项改成不含糊的名称，例如 `additional_spectral_transmission`；
- 保留旧 key 作为兼容别名并警告；
- 两者同时开启时要求配置显式写：

```text
atmosphere.allow_additional_spectral_factor=true
```

- 输出 metadata 分别记录两项，而不是都叫 atmosphere；
- 保存每个曲线文件的 hash。

## 8. 第七项：补齐光程和运行 provenance

### 8.1 光程诊断

`OpticalSurfaceHit` 当前只保存最终时间，没有保存：

- 输入锚点到镜面的有符号距离；
- 镜面到输出平面的距离；
- collector 内部距离；
- 总几何/光学路径。

建议增加这些字段。这样 `z=-16`、full-line 时间、collector 延时和空气/真空光速都能独立检查，而不是只能看最终时间。

### 8.2 输出可复现性

ROOT 当前：

```text
producer_version = unknown
source_sha256 = empty
```

虽然保存了展开配置文本，但还应保存：

```text
git commit + dirty flag
compiler/build flags
dependency versions
source EventIO hash
mirror/camera/efficiency/PDE/atmosphere/obstruction/error 文件 hash
随机 seed 和流派生版本
坐标与单位约定
启用/关闭的模型清单
```

自定义 ROOT/HDF5 格式本身没有问题；缺的是复现一次运行所需的证据。

## 9. 目前不建议修改的部分

### 9.1 不改成 sim_telarray 的坐标和编号

LACT 当前的 telescope-local、镜片 ID、pixel ID 和输出列映射可以继续保留。

应增加内部不变量测试，不应做无物理收益的大规模坐标/编号重写。

### 9.2 不把真空光速直接改成 sim_telarray 的空气光速

当前 `propagation.speed_of_light_m_per_ns` 是明确可配置模型。先保存光程并量化真空/空气差异，再根据定时精度目标决定是否增加折射率/群速度模型。

### 9.3 不优先实现完整电子学

当前核心阻塞是 collector、3D 原点和随机流。SPE、模拟脉冲、ADC、SiPM
saturation/crosstalk/afterpulse/dark count 和硬件 trigger 可以在光学核心可信后作为独立模块实现。

配置名中的 `full_response` 应先改成更准确的：

```text
optics_pe_nsb_simple_trigger
```

避免能力声明超过实际实现。

## 10. 推荐实施顺序

### PR 1：Collector correctness

- 重构 `SurfaceHit`；
- Bezier 交点和法线共用参数；
- 修 strict comparator；
- 增加 collector 光程；
- ASan/UBSan、角度扫描和 fuzz。

这是第一优先级，因为当前正式配置会进入已证实的越界路径。

### PR 2：EventIO reference-to-optical transform

- 新增通用三维 offset；
- 2D/3D 共用；
- 旧配置兼容；
- 增加 3D 和坐标平移测试。

### PR 3：Deterministic random streams

- 公共固定 seed 派生；
- expectation/stochastic 都设置 photon stream；
- scatter 去掉绝对坐标 hash；
- 逐镜误差按 telescope/facet/stage 派生；
- 平移和行重排不变量测试。

### PR 4：Validation contracts

- camera geometry validator；
- collector uniform-geometry 检查；
- Photon/EventIO 数值合法性；
- efficiency/probability 范围；
- 禁止或确认 collector/funnel、大气双算。

### PR 5：Path/provenance/scope

- 保存分段光程；
- 完整 hash/build/random metadata；
- 更新配置命名和能力说明；
- 端到端可复现测试。

## 11. 完成判据

完成以上修改后，才能按下列顺序重新判断：

1. 解析单镜和坐标不变量证明几何正确；
2. ASan/UBSan 和 collector fuzz 证明数值安全；
3. 2D/3D EventIO 证明输入语义一致；
4. 固定随机身份和行重排证明可复现；
5. 概率范围和大气/collector 组合证明效率没有双算；
6. 输出 hash、光程和配置证明一次运行可独立复现；
7. 最后再用 sim_telarray 或其他程序做公共物理模型的独立旁证。

这套顺序不会把 LACT 改成 sim_telarray，而是把 LACT 自己的物理契约逐项闭合。
