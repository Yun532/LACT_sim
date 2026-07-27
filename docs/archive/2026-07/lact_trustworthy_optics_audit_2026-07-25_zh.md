# LACT_sim 可信光学追踪独立审计

> 状态说明：本文记录的是修改前基线 `0048369` 的审计结论。文中“当前缺陷”
> 不应被理解为修改后仍然存在；对应实施状态和新测试见
> `lact_trustworthy_optics_implementation_2026-07-25_zh.md`。

> sim_telarray 只作为成熟程序的旁证和遗漏提示，不是实现规范。
> 审计日期：2026-07-25
> LACT_sim 基线：GitHub 远端 `main`，`0048369cd06e87322d9ef3611a4a369a1152fe99`
> sim_telarray 参考：MPIK 官方 2025-02-20 发布包，SHA-256
> `98889ecc104da160894b2982c534261c0beebabf33b5fae99e5a72e86aaa9723`

## 1. 审计原则

本审计只用以下标准判断 LACT_sim：

1. 输入字段、单位、时间和坐标语义是否明确；
2. 同一物理系统在合法的坐标平移、刚体旋转、数据重排和 ID 重编号后，物理结果是否保持应有的不变性；
3. 镜面求交、口径裁剪、法线、反射、遮挡、焦面和 collector 是否符合 LACT 自己声明的物理模型；
4. 效率、随机误差和光电子抽样是否有明确概率语义并可复现；
5. 数值算法是否安全，异常输入是否被拒绝而不是静默产生错误结果；
6. 输出是否足以复现本次运行；
7. 测试是否能证明上述性质，而不只是证明某组样例恰好通过。

因此，以下差异本身都不是问题：

- LACT 与 sim_telarray 使用不同的内部坐标轴方向；
- 镜片、像素或数组使用不同的 ID；
- 使用不同的数据结构、求交算法、随机数发生器和输出格式；
- LACT 只实现自己需要的单反射面光学，而不实现 sim_telarray 的全部望远镜类型；
- 光学追踪器不包含完整电子学。

只有当不同约定在 LACT 内部被混用，或者改变了 LACT 声明要模拟的物理结果时，才记为问题。

## 2. 总结论

### 2.1 “我们现在是否正确”

不能把当前整个 `main` 判定为已经完全正确，但可以更精确地说：

- **理想单反射面主链已经有较强正确性证据**：二维 EventIO 输入、坐标旋转、当前 `z=-16 m` 原点转换、镜片求交、前表面选择、反射、焦平面求交、像素映射以及 expectation/stochastic PE 主链的结构基本合理；
- **当前仍有一个会实际越界的 P0 collector 缺陷**；
- **3D EventIO 的光学原点转换与 2D 不一致**；
- **启用反射方向散射后，随机结果依赖绝对坐标原点**；
- 若使用自定义相机 CSV，当前缺少足够的几何和 ID 校验；
- 配置允许把两种大气透过率同时相乘，但没有明确防止误用。

所以目前最准确的状态是：

> 理想主镜到焦面的核心光路接近可信；
> 启用当前 Bezier collector 的正式相机链还不能判定可信；
> 3D EventIO 支持还不能判定正确。

### 2.2 “我们现在是否完整”

若“完整”指可信的 LACT 光学追踪核心，则应至少覆盖：

```text
输入光线
  -> 单位与坐标转换
  -> 镜片几何与误差
  -> 入射遮挡
  -> 镜面求交、口径、法线与反射
  -> 反射后遮挡
  -> 焦面
  -> 相机像素与 collector
  -> 光程时间
  -> 波长相关效率
  -> 可复现输出和验证
```

这些模块在代码中基本都有，但 collector 安全性、3D 原点、随机不变量、输入校验和运行 provenance 仍有缺口，因此光学核心尚未达到“完整且可信”。

完整电子学不是光学追踪器成立的前提。当前程序中的 NSB、PE 时间分箱和 multiplicity trigger 可以作为简化的下游模型使用，但不应称为完整探测器电子学。

## 3. 从输入到输出的当前实现

`run_corsika_trace` 的实际主链是：

```text
EventIO
  -> PhotonBunch
  -> CORSIKA/阵列坐标到 telescope-local
  -> EventIO 光学原点转换
  -> 波长解析与大气透过
  -> expectation 或 stochastic PE 候选
  -> 入射遮挡
  -> 镜片求交、反射和光学散射
  -> 反射后遮挡
  -> 输出平面
  -> camera pixel / square collector / SiPM 面
  -> PE 累积
  -> 可选 NSB、时间分箱和简化 trigger
  -> CSV / HDF5 / lact_event ROOT
```

模块级判断如下：

| 模块 | 当前判断 | 主要依据或缺口 |
|---|---|---|
| EventIO 2D 字段和单位 | 基本正确 | cm→m、ns、multiplicity、波长分支明确 |
| EventIO 2D 光学原点 | 当前几何下正确 | `z=0` 是源格式占位；转换后 `z=-16` |
| EventIO 3D 光学原点 | 不正确/未完成 | 3D 仍保留输入原点，未映射到 LACT 光学原点 |
| 坐标旋转 | 较强证据 | 有 NWU/global/local 路径和坐标测试 |
| 镜片 CSV | 较强证据 | 有限值、法线、尺寸、曲率、重复 ID 校验 |
| 平面/球面/抛物面求交 | 较强证据 | 解析求交、物理球冠筛选、口径裁剪和单元测试 |
| 前表面与反射 | 较强证据 | 法线定向、前表面检查、反射公式清楚 |
| 遮挡 | 基本完整 | 分入射路径和反射路径处理 |
| 输出平面 | 较强证据 | 正交基显式校验，使用局部 `u/v` |
| 光程时间 | 模型明确 | 默认真空光速；空气折射率是可选精化，不是当前错误 |
| 相机像素 | 当前数据可用，通用校验不足 | 当前 1616 像素 ID 唯一且尺寸为正；loader 未强制保证 |
| Bezier collector | 不可信 | ASan 证实越界；交点比较器也违反严格弱序 |
| 效率/PE | 结构基本正确 | expectation/stochastic 分支和预应用效率标志明确 |
| 随机误差 | 统计模型可用，不变量不足 | 散射依赖绝对坐标；镜片误差依赖 CSV 顺序 |
| 输出 | 数据结构可用，provenance 不足 | ROOT 仍写 `producer_version=unknown`、空 `source_sha256` |
| 电子学 | 明确是占位/简化模型 | 不是完整波形与硬件读出模拟 |

## 4. 输入、坐标和 `z=0 -> z=-16`

### 4.1 二维 EventIO 的字段语义

vendored hessio `mc_tel.h::struct bunch` 定义：

- `photons`：bunch 中代表的光子数；
- `x,y`：相对望远镜的到达位置，单位 cm；
- `cx,cy`：光子传播方向余弦；
- `ctime`：到达时间，单位 ns；
- `zem`：发射高度，单位 cm；
- `lambda`：nm；零值表示未指定波长。

LACT 的读取逻辑：

- `x,y` 乘 `0.01` 变成 m；
- `z` 暂时设为 0，因为二维记录没有 z 字段；
- `cz=-sqrt(1-cx²-cy²)`，随后归一化；
- `ctime` 保持 ns；
- `photons` 进入 multiplicity；
- 正、零、负波长分别进入实际波长、缺失波长模型和已预应用光学效率语义。

这部分的主要语义是自洽的。

### 4.2 为什么当前应当转换到 `z=-16`

这里必须区分两个动作：

1. `makeBunch()` 中的 `z=0`：只是在 EventIO 源坐标中补齐缺失的第三维；
2. `transformEventIOBunchToTraceFrame()` 中加 `-16 m`：把 EventIO 的望远镜/fiducial 参考原点表达成 LACT 导入几何所用的光学原点。

当前 LACT 镜面 CSV 的镜片中心约为：

```text
z_min = -15.95807029 m
z_max = -15.44059938 m
```

相机平面为 `z=-8 m`。因此在当前坐标约定下，EventIO 的望远镜参考原点应映射到 LACT 的 `z=-16 m` 附近。

这不是“因为 sim_telarray 也这么做”，而是因为：

- 输入数据的原点有自己的定义；
- LACT 几何文件的原点有自己的定义；
- 两者之间存在已知的 16 m 平移。

sim_telarray 只是帮助暴露了这个原点差异。

所以结论是：

> 对当前 EventIO 和当前镜面/相机几何，`z=-16 m` 是正确的。
> 它不应被理解为所有输入和所有镜面模型的通用常数。

更稳妥的长期设计是把它明确命名成“EventIO 原点到 LACT 光学原点的变换”，并让 2D、3D 共用，而不是把它命名成只属于二维记录面的 z 参数。

### 4.3 3D EventIO 的真正问题

`mc_tel.h::struct bunch3d` 明确说明：

```text
x,y,z = relative to telescope,
fiducial sphere center at (0,0,0)
```

LACT 读取 3D bunch 时把 cm 转为 m，但之后只对 `eventio_2d=true` 的 bunch 加 `-16 m`。于是：

- 2D：fiducial 原点被映射到 LACT 的 `z=-16 m`；
- 3D：fiducial 原点仍停留在 LACT 的 `z=0`；
- 镜面和相机却仍在 `-16/-8 m`。

这是 LACT 内部对同一个 EventIO 参考原点使用了两种解释，与 sim_telarray 的实现无关。

现有服务器样例的首个 shower 共 1,360 个 bunch，全部是 2D，因此现有样例不受影响；但代码声称支持 3D EventIO，所以必须补齐一般化的原点变换和 3D 回归测试。

### 4.4 输入异常值还需要更严格处理

目前仍应增加以下 fail-fast 校验：

- Photon CSV 的位置、方向、时间、波长、权重和 multiplicity 必须有限；
- 方向不能为零向量；
- 权重和 multiplicity 不能为负；
- 2D EventIO 若 `cx²+cy²` 明显大于 1，应报错；只能对浮点量化造成的微小超界做容差夹紧；
- 3D EventIO 的方向必须有限且非零。

现在 `downwardDirZ()` 会把任何负的 `1-cx²-cy²` 直接夹到 0，严重损坏的输入可能被静默变成水平光线；expectation 模式也没有像 stochastic 模式一样完整拒绝非法 multiplicity。

## 5. 坐标轴、数组编号和 ID

### 5.1 不需要与 sim_telarray 对齐

LACT 可以自由选择：

- `+x/+y/+z` 的方向；
- 方位角零点和旋转方向；
- 镜片 ID、像素 ID 和望远镜 ID；
- 容器中的行顺序。

只要所有 LACT 输入、几何、效率表、遮挡、相机映射和输出都使用同一套约定，物理结果就没有问题。

因此，之前把“LACT 局部 x/y 与 sim_telarray 不同”列为待修问题是错误的，现已撤销。

### 5.2 编号什么时候才会真的出问题

编号本身不会造成物理错误。只有代码隐含假设以下等式时才会出错：

```text
physical_id == vector_index == CSV_row_number == output_column
```

当前 LACT 在若干关键位置已经避免了这种假设：

- 镜片校验拒绝重复 ID；
- 逐镜反射率文件按 ID 匹配，并检查未知、重复和缺失 ID；
- ROOT/HDF5 相机输出建立 `pixel_id -> column` 映射，而不是直接把 ID 当数组下标。

这说明不同 ID 体系本身不是当前问题。

当前正式相机 CSV 的 1,616 个像素也没有重复 ID，尺寸均为正。

不过 `readCameraCsv()` 本身没有强制检查：

- 重复或非法 pixel ID；
- NaN/Inf 坐标；
- 非正尺寸；
- 像素几何重叠造成的多义归属。

重复 pixel ID 会使输出列映射覆盖已有像素，因此通用相机 loader 仍需补齐校验。这里要修的是 LACT 的输入契约，不是改成 sim_telarray 的编号。

### 5.3 应验证的内部不变量

可信性测试应直接验证：

- 整个光子、镜面、遮挡和焦面同时平移后，局部命中位置与光程差不变；
- 整个系统刚体旋转后，结果只做同样旋转；
- 镜片 CSV 行重排不改变结果；
- 像素 CSV 行重排不改变 `pixel_id` 命中；
- 只对 ID 做一一重编号，并同步所有 ID 引用后，光学结果不变；
- 使用不同内部坐标轴但表达同一物理场景时，结果通过已知变换对应。

这些测试比“内部编号与 sim_telarray 一致”更能证明程序可信。

## 6. 镜面、反射、遮挡和焦面

### 6.1 镜面几何

当前支持：

- generated dish；
- 直接镜片 CSV；
- 随 elevation 插值的镜片序列；
- 平面、球面和抛物面求交；
- 圆形、六边形和方形口径；
- 逐镜反射率、粗糙度和失调参数。

镜片输入的校验较完整：

- ID 非负且唯一；
- 中心和法线有限，法线非零；
- 尺寸和曲率合法；
- 不支持的 Polynomial 明确拒绝；
- 效率和误差参数有限且非负。

这是当前代码中做得较可靠的一部分。

### 6.2 求交和反射

当前核心实现的正确点包括：

- 平面镜用解析直线/平面交点；
- 球面镜解二次方程，并额外裁掉同一完整球面的远端球冠；
- 抛物面使用局部正交基中的解析方程；
- 求交后按真实口径裁剪；
- 法线统一到镜片前表面；
- 拒绝从背面入射的交点；
- 反射方向按几何法线计算；
- 只接受沿反射方向前方的输出平面交点。

二维 EventIO 使用完整有向直线查找镜面，是因为记录点只是该物理光线在输入参考系中的锚点。镜片选择用参数 `t` 的顺序，而不是用 `|t|`，因此只沿同一光线移动锚点不会改变镜片选择。相关 near/far cap、前后表面和锚点回归已有测试。

### 6.3 遮挡

程序区分：

- 入射光到镜面之前的上游遮挡；
- 镜面到焦面之间的反射后遮挡；
- 标记型和真正阻断型结构；
- 结构在 telescope-local 与 trace frame 之间的变换。

这种分段方式符合单反射面光学追踪的物理边界。后续完整性重点应是用解析结构和非对称结构验证遮挡率，而不是复制其他程序的结构数据格式。

### 6.4 输出平面和时间

输出平面具有显式 point、normal、u、v。自定义 `u/v` 必须同时给出，并校验互相正交且垂直于 normal。

命中时间计算为：

```text
t_hit = t_input + (incoming_path + reflected_path) / c_model
```

默认 `c_model=0.299792458 m/ns`，即真空光速。这是一个明确、可配置的模型选择，并不因为 sim_telarray 使用空气折射率修正就自动成为错误。

如果 LACT 的目标精度进入亚纳秒绝对定时，应新增随高度/介质变化的群速度模型并量化误差；在此之前应在输出元数据中明确记录采用的光速，而不是强行换成某个参考程序的值。

## 7. Collector：当前最严重的问题

正式相机配置使用：

```text
camera.collector=bezier
camera.collector_material=true_reflect
```

ASan 在 `test_light_collector_response` 中确认：

```text
heap-buffer-overflow
SquareCone::ray_trace_impl
  -> get_reflect_direction
  -> Bezier2CylindricalSurface::get_normal_vector
```

直接原因位于 `LightCollectorSquareCone.hpp`：

```cpp
auto it = std::find_if(roots.begin(), roots.end(), ...);
auto t = *it;
```

代码假设四个候选根中一定能找到一对近似相等的根，但数值误差或无对应根时 `it==roots.end()`，随后发生越界读取。

同一文件中选择最近相交面的 comparator：

```cpp
return d1 < d2 || fabs(d1 - d2) < 1e-6;
```

对近似相等的 `d1,d2`，`comp(a,b)` 和 `comp(b,a)` 都可能为 true，违反 `std::min_element` 要求的严格弱序，会产生未定义或不稳定选择。

这两个问题完全由 LACT 自身的数值安全标准判定，不需要 sim_telarray 证明。

修复时不应只加一个 `if (it==end)` 就结束，还应：

1. 直接从交点求解时得到的 Bezier 参数 `t` 传给法线计算，避免再次分别从 x/z 反解并寻找“相等根”；
2. comparator 只使用严格的 `d1 < d2`；
3. 对边界、角点、切线、近平行光、最大反射次数和无根情况做性质测试；
4. 在 ASan/UBSan 下对 collector 做密集角度扫描和随机 fuzz；
5. 明确无有效法线时是拒绝该光线，不能继续使用未初始化结果。

在这些修复和验证完成前，关闭 collector 的白板/焦面结果可以使用；启用 Bezier collector 的正式相机结果不能作为已经可信的依据。

## 8. 光学效率、波长和随机模型

### 8.1 效率链

未预应用效率的光子使用：

```text
constant
× mirror reflectivity(λ)
× filter transmission(λ)
× atmosphere factor(λ)
× SiPM PDE(λ)
× facet reflectivity scale
× collector/funnel acceptance
```

负波长 `CEFFIC` 语义通过 `optical_efficiency_preapplied` 避免重新施加波长相关预筛选，只保留 LACT 本地仍应施加的逐镜和 collector 因子。expectation 和 stochastic PE 路径也分别处理权重与离散接受。

公共理想配置下，18,940 个匹配事件的独立交叉检查得到：

- stochastic 总 true PE 比：`0.9974932`；
- stochastic 事件相关系数：`0.9998604`；
- expectation 总量比：`0.9974188`；
- expectation 事件相关系数：`0.9998689`。

这说明公共理想模型的统计实现没有明显偏差，但它只是一项旁证，不是 LACT 设计必须复制的数值。

### 8.2 两种“大气”可能被重复相乘

代码同时允许：

- `efficiency.atmosphere_transmission` 或旧别名 `atmosphere.transmission`：仅波长相关的效率因子；
- `AtmosphereTransmission`：使用发射高度、波长和天顶角的 MODTRAN 光学深度模型。

如果用户同时启用两者，代码会把两项都乘上。这可能是有意表达两个独立损耗，也可能是把同一个大气衰减算了两次。

可信程序不应依赖用户猜测，应当：

- 为两项使用不混淆的名称；
- 同时启用时要求显式确认其物理含义，或默认报错；
- 在输出中分别记录两项及其输入文件 hash；
- 增加“只开 A、只开 B、同时开 A+B”的乘法测试。

### 8.3 散射随机数依赖绝对坐标

反射方向散射的 seed 当前混入：

```text
global seed
+ facet ID
+ photon.pos.x/y/z
+ photon.dir.x/y/z
+ photon.random_stream_id
```

只把整个物理系统平移一个常量，几何交点本应只做同样平移，但 `photon.pos` 的 hash 会改变，于是抽到不同的粗糙度散射方向。

另外，`std::hash<double>` 不提供跨标准库实现一致的序列保证，因此即使用户
seed 相同，当前做法也不能作为严格的跨平台可复现随机流规范。

这违反坐标原点不变量。它不是统计分布错误，但会导致：

- 同一物理场景换一个合法原点后，逐光子结果不再对应；
- `z=0/-16` 这样的纯坐标表达变化也会改变散射随机实现；
- 回归结果依赖任意坐标表示。

建议使用稳定的物理/数据身份生成随机流：

```text
run/event/telescope/source_bunch/represented_photon
+ facet ID
+ random stage
+ user seed
```

而不要使用绝对位置作为随机身份。随后增加 scatter 开启状态下的整体平移和旋转不变量测试。

### 8.4 镜片误差依赖 CSV 行顺序

`applyFacetErrors()` 用一个顺序 RNG 遍历 `facets`。同一组镜片只改变 CSV 行顺序，就会把不同随机误差分配给不同镜片。

统计分布不因此改变，但固定 seed 的物理 realization 和复现性会改变。更稳妥的做法是按：

```text
user seed + telescope ID + facet ID + error stage
```

派生独立随机流。这样 ID 可以与 sim_telarray 不同，但 LACT 自己的镜片身份和随机 realization 是稳定的。

## 9. 相机、SiPM、NSB、电子学和 trigger

### 9.1 光学追踪核心内应负责的部分

以下属于当前光学/光电转换链，应当正确：

- 焦面局部坐标；
- pixel containment；
- collector 的入口、壁面反射和 SiPM 出口；
- 镜面、滤光片、collector 和 PDE；
- 光子到 PE 的 expectation 或 stochastic 语义；
- PE 到达时间。

其中 collector 仍是当前阻塞项。

### 9.2 当前电子学实际范围

代码自己已经明确标记：

- `ElectronicsConfig` 为空；
- metadata 写入 `integrated_pe_placeholder`；
- `waveform.source=electronics` 会直接报错；
- 当前 waveform 是 photon-count 或 PE 时间 histogram proxy；
- 未包含真实 SPE 脉冲、增益涨落、transit-time jitter、模拟带宽、ADC/FADC、pedestal/noise、饱和、crosstalk、afterpulse、dark count 和硬件 trigger board；
- 当前 trigger 是简化的 PE 窗口 multiplicity 模型。

这不妨碍 LACT 成为可信光学追踪器，但配置名中的 `full_response` 容易让使用者误以为已经模拟完整探测器响应。建议改成更准确的名称，例如：

```text
optics_pe_nsb_simple_trigger
```

如果未来项目目标扩展为端到端探测器模拟，再单独定义电子学接口、脉冲模型和硬件验证标准，不需要把 sim_telarray 的电子学实现原样搬过来。

## 10. 输出与可复现性

CSV、HDF5 和自定义 ROOT 都可以是可信输出格式；不需要改成 sim_telarray EventIO。

当前输出的优点：

- 保存 event/telescope/pixel 身份；
- pixel ID 与列号分离；
- 可保存命中、像素积分、稀疏/稠密 waveform proxy；
- HDF5 中记录多项模型配置。

当前 provenance 缺口：

- ROOT writer 的 `producer_version` 仍为 `unknown`；
- `source_sha256` 为空；
- 没有系统记录所有 include 后的最终配置文本；
- 镜面、相机、效率、PDE、大气、遮挡和误差 CSV 的 hash 不完整；
- 缺少编译器、依赖版本、随机 seed 派生规则和 dirty-tree 状态。

要让一次光追可以被独立复现，至少应输出：

```text
git commit + dirty flag
build/compiler/dependency versions
fully resolved config
all external input paths + SHA-256
random seeds and stream policy
coordinate-frame declaration
units
enabled/disabled model inventory
```

## 11. 测试证据与还缺的验证

### 11.1 已执行证据

在与远端 `main` 相同代码树的服务器构建中：

- 含 EventIO/ROOT 的普通 CTest：`27/27` 通过；
- 完整基础 UBSan 构建：`26/26` 通过；
- ASan：`25/26`，唯一失败是 `test_light_collector_response`，并给出确定的 heap-buffer-overflow；
- 另一个所谓“full UBSan”目录有 3 个测试可执行文件未构建，属于构建不完整，不能算成运行通过或运行失败。

普通测试全通过不能覆盖已经被 ASan 发现的内存错误，因此当前 release/normal 测试绿色不等于 collector 正确。

### 11.2 为“可信光学追踪器”必须补的测试

优先级最高：

1. Bezier collector 的 ASan/UBSan 边界、角点、切线和随机 fuzz；
2. 3D EventIO 原点转换的解析场景；
3. scatter 开启时的整体平移、刚体旋转不变量；
4. 镜片和像素 CSV 行重排不变量；
5. 非法输入 fail-fast 测试；
6. 同时启用两套大气损耗时的显式行为测试。

随后补充：

- 单镜解析基准：平面、球面、抛物面；
- DC 多镜中心光线、轴外 spot、时间 spread；
- 非对称镜片效率、非对称遮挡和非对称相机；
- 大角度、掠入射、近平行和数值容差扫描；
- 总权重/接受概率范围检查；
- 固定输入与固定 seed 的跨编译器可复现策略；
- 输出 round-trip 和 provenance 完整性。

## 12. sim_telarray 在本审计中的正确用法

本次参考它只做了三件事：

1. 提醒检查 EventIO fiducial 原点，由此发现并解释了当前 `z=-16`；
2. 在双方都关闭额外模型的公共理想场景下提供独立统计交叉检查；
3. 提醒逐项盘点波长、大气、遮挡、collector、光程、PE、NSB、trigger 和输出 provenance。

它没有被用作以下要求：

- 不要求 LACT 改坐标轴；
- 不要求 LACT 改镜片/像素/望远镜 ID；
- 不要求使用相同数组布局；
- 不要求使用相同随机数；
- 不要求输出 sim_telarray EventIO；
- 不要求复制它的电子学；
- 不要求支持它的所有望远镜光学类型。

## 13. 修复顺序

### P0：阻止当前正式 collector 产生未定义行为

1. 重构 Bezier 交点/法线参数传递；
2. 修复 `min_element` comparator；
3. ASan/UBSan + 角度扫描 + fuzz 全通过；
4. 在完成前，不把 Bezier collector 输出标为可信正式结果。

### P1：修复 LACT 自身的物理不自洽

1. 把 EventIO→LACT 光学原点变换一般化到 2D 和 3D；
2. 用真实或合成 3D bunch 验证位置、交点和时间；
3. 将散射随机流从绝对坐标中解耦；
4. 补全相机 CSV 的 ID、有限值、尺寸和歧义校验。

### P2：完成可信性工程

1. 让逐镜误差对 CSV 行顺序稳定；
2. 消除两套大气配置的误用空间；
3. 完整记录 provenance；
4. 增加坐标/顺序/ID 不变量测试；
5. 把 `full_response` 改成符合当前真实范围的名称。

## 14. 最终判断

不需要、也不应该把 LACT 改成 sim_telarray。

当前最值得保留的独立设计包括：

- LACT 自己清楚的 telescope-local 光学坐标；
- 显式输出平面和 `u/v`；
- 镜片 CSV 的严格验证；
- ID 到输出列的显式映射；
- expectation/stochastic 两种响应模式；
- 自定义 HDF5/ROOT 输出；
- 可配置的光速和模块化效率。

当前真正阻止“可完全信赖”的不是双方内部实现不同，而是 LACT 自身仍存在：

- 一个已证实的 collector 越界；
- 一个 2D/3D EventIO 原点不一致；
- 一个散射随机数依赖绝对坐标的问题；
- 若干输入契约和 provenance 缺口。

先修复这些 LACT 自己的正确性和可验证性问题，再用解析解、内部不变量和独立程序交叉检查共同验证，才是把它做成可信光学追踪程序的正确路线。
