# LACT / pyLAST 坐标系代码审计

本文只依据执行代码、真实输出和数值回归；旧注释与旧网页不作为坐标定义证据。

审计基线：

- LACT_sim：远端 `main@48ea63115ffff4dd840e792149f79a2549f3abb0`，在该 main 工作区直接应用本次坐标改动；
- pyLAST：远端 `main@7a65cf1f12e37e0e437968c37da34b9a7fc4dd08`，线性纳入 LACT ROOT reader 至 `953d07392a542a3362a0ce35f14bcdbe4c347c9a`，再直接应用 reader 改动；
- 没有以旧 `fix` 分支作为代码基线。

## 1. CORSIKA / 阵列全局坐标

CORSIKA/EventIO 和 ROOT 阵列真值使用 NWU：

```text
+X = North
+Y = West
+Z = Up / Sky
East = -Y
```

望远镜方位角 `A=0°` 指北，正角度顺时针转向东；仰角 `E` 从地平面向天空增加。

代码出处：

- `src/io/EventIOPhotonSource.cpp`：EventIO `x/y` 从 cm 转为 m，下降光子的 `dir_z<0`；
- `src/app/ArrayTimingCorrection.cpp::corsikaNwuViewingDirection()`：天空方向 `(cosE cosA, -cosE sinA, sinE)`；
- `src/io/LactEventRootWriter.cpp`：真芯位字段 `core_x_north_m/core_y_west_m`；
- pyLAST `include/Coordinates.hh`：同样计算 `(cos az cos alt, -sin az cos alt, sin alt)`。

## 2. 镜片、遮挡与望远镜结构坐标

镜片和遮挡 CSV 的局部坐标按“从相机看向镜面”的结构视角解释。望远镜指北时：

```text
mirror local +x = West
mirror local +y = 仰角增加 / sky-up
mirror local +z = 镜面指向相机与天空的光轴
```

在 NWU 中的实际代码公式是：

```text
m_x = ( sinA,          cosA,         0    )
m_y = (-sinE cosA,     sinE sinA,    cosE )
m_z = ( cosE cosA,    -cosE sinA,    sinE )
```

代码出处：`src/app/OpticalSimCommon.cpp::buildTelescopeFrame()` 与
`buildCorsikaNwuTelescopeFrame()`。本次让两个入口逐分量使用同一公式，并由
`apps/test_coordinate_frames.cpp` 数值比较，防止平行光与 CORSIKA 结构再次翻轴。

单块镜片还有自己的 aperture `u_f/v_f`。若 CSV 未显式给出，代码在
`include/geometry/MirrorFacet.hpp::apertureFrame()` 中由镜片法向构造；它不等于整台望远镜的 mirror `x/y`。

## 3. 白板、焦平面和相机 LACT `u/v`

白板/相机是“从镜片看向相机”的观察方向，与结构坐标的水平观察方向相反：

```text
+u = -mirror local x   # 指北时向 East
+v = +mirror local y   # sky-up
normal = -mirror local z
```

标准配置原值：

```ini
output.plane_point=0,0,-8
output.plane_normal=0,0,-1
output.plane_u_axis=-1,0,0
output.plane_v_axis=0,1,0
```

代码链：

```text
OpticalTracer.cpp: u/v = dot(surface-plane_point, plane_u/v)
OpticalSimCommon.cpp: camera_x_m=u, camera_y_m=v
LactEventRootWriter.cpp: ROOT camera_pixels x_m/y_m = u/v
```

`buildOutputPlane()` 现在强制 `u × v = normal`。旧组合
`normal=-z, u=+x, v=+y` 是左手系，会在程序启动时被拒绝，而不是只改网页标签。

## 4. pyLAST 边界

pyLAST 的相机坐标表示天空/source-offset；LACT 的 `u/v` 表示硬件焦平面反射光斑。ROOT reader 只转换一次：

```text
pix_x = -LACT v
pix_y = -LACT u
```

pyLAST 原生画图再使用：

```text
screen horizontal = pix_y
screen vertical   = pix_x
```

代码出处：`root/LactEventSource.cpp`、`root/test_lact_event_source.cpp`、
`test/test_coordinates.cpp` 与 `src/pylast/visualize/visualize.py`。

## 5. 四方向端到端符号

望远镜 `az=0°/el=70°`，天源偏离光轴 4°：

| 天源方向 | 光子传播横向 | LACT 真实输出 | pyLAST 天空偏移 |
|---|---|---|---|
| 东 | `+mirror x`（向西传播） | `u=-0.571054 m` | `pix_y=+0.571054 m` |
| 西 | `-mirror x` | `u=+0.571243 m` | `pix_y=-0.571243 m` |
| 上 | `-mirror y` | `v=-0.571319 m` | `pix_x=+0.571319 m` |
| 下 | `+mirror y` | `v=+0.571366 m` | `pix_x=-0.571366 m` |

数值是四次 `run_optical_sim` 对无遮挡物理输出按程序 `weight×relative_efficiency` 计算的质心；完整 `hit_surface` 行、全量反射点和全部遮挡标志也一并保存在网页数据中。网页不交换、不翻号、不重设光斑中心。

## 6. CORSIKA event 1909

网页第 4 页重新运行了 4 GB EventIO 文件中的 event 1909，并保存新的 ROOT：

- 32 台望远镜；
- `86550` photon bunch；
- `197576` 个无遮挡输出面命中；
- `116995` 个相机命中 / PE；
- ROOT SHA-256：`546bdafd4d8ce31218cbc717c3c8b7463c4a48c80aebd044b6286a11f7ca174b`；
- EventIO SHA-256：`3feee5b7f3a001858201eea2cf75ba3f5f0277283e29900b5f259bd2c9bc4220`。

二维 bunch 的三维“发射点”不是 EventIO 显式字段，只是诊断派生量：以原始接收平面 anchor、原始 direction 和原始 `zem` 求直线与高度面的交点。相机图始终读取 ROOT `camera_pixels` 与 `image_cherenkov_pe`，不使用派生发射点反算相机坐标。

## 7. 编译与回归证据

服务器隔离目录实际编译并运行：

- `test_coordinate_frames` 通过；
- `test_off_axis_orientation` 真实 1229-facet 光追通过：传播 `+mirror x` 得 `u=-0.142653 m`，传播 `-mirror x` 得 `u=+0.142654 m`；
- pyLAST `test_coordinates`：9 cases / 30 assertions 通过；
- pyLAST `test_lact_event_source` 通过；
- `run_optical_sim` SHA-256：`20737903cdfeca7a7851ff32d31a07dbb664d758ad8339e1ce815f2d7025a914`；
- `run_corsika_trace` SHA-256：`db3946f33a8bc4d198d63191fe09f43b9d850aaeed5f22772698af3161afb3c9`；
- 同源审计归档 SHA-256：`e7335bba7e65863fbbeff6c17e06d753a954e213302a69557f698dfe43cfb2b3`；
- 全部最新结果归档 SHA-256：`9778624e091498d5ec021e132e47b0f8cb88f354632cd9bed0bc19aa5c43f67d`。

网页允许的变换只有观察相机旋转、缩放、平移以及 pyLAST reader/renderer 的明确边界映射；LACT 原始 `u/v`、ROOT PE、CORSIKA NWU 和镜片/遮挡原坐标不得在绘图脚本中另行改写。
