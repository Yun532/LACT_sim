# LACT 坐标系代码审计（main）

本文只记录从执行代码和数值测试得到的坐标约定；项目旧注释和旧说明不作为判定依据。

审计基线：`LACT_sim main` 提交 `48ea63115ffff4dd840e792149f79a2549f3abb0`。

## 1. 全局阵列与 CORSIKA 输入

CORSIKA/EventIO 数组坐标是 NWU：

- `+X`：磁北 North；
- `+Y`：西 West；
- `+Z`：天顶 Up；
- 方位角 `A=0` 指北，`A` 增大时顺时针转向东；
- 高度角 `E` 从地平面向天空增加。

代码出处：

- `src/io/EventIOPhotonSource.cpp`：EventIO 的 `x/y` 由厘米乘 `0.01` 写入米，方向 `z` 取下降光子的负根；
- `src/app/ArrayTimingCorrection.cpp::corsikaNwuViewingDirection()`：天空方向为
  `(cosE cosA, -cosE sinA, sinE)`；
- `src/io/LactEventRootWriter.cpp`：ROOT 真值写成 `core_x_north_m`、`core_y_west_m`。

## 2. 镜片、遮挡和望远镜局部轴

镜片 CSV 与遮挡 CSV 的局部坐标按“相机看向镜面”的视角解释：

- `+m_x`：望远镜指北时朝西；
- `+m_y`：俯仰增加方向，即 sky-up；
- `+m_z`：从镜面指向天空的光轴。

在 NWU 中：

```text
m_x = ( sinA,          cosA,         0    )
m_y = (-sinE cosA,     sinE sinA,    cosE )
m_z = ( cosE cosA,    -cosE sinA,    sinE )
```

且 `m_x × m_y = m_z`。代码出处是
`src/app/OpticalSimCommon.cpp::buildTelescopeFrame()` 与
`buildCorsikaNwuTelescopeFrame()`。二者必须逐分量相同，分别供普通光学入口和 CORSIKA 入口使用。

`transformBunchToTelescopeLocal()` 用点积把 NWU 光子位置和传播方向投影到这三根轴；相对坐标不减望远镜位置，全局坐标只减一次。

## 3. 白板、焦平面和相机原始坐标

焦平面从镜面看向相机。由于这个观察方向与镜片 CSV 的观察方向相反，水平轴定义为：

```text
+u = -m_x   # 指北时朝东
+v = +m_y   # sky-up
normal = -m_z
```

因此标准 8 m 焦平面配置是：

```ini
output.plane_point=0,0,-8
output.plane_normal=0,0,-1
output.plane_u_axis=-1,0,0
output.plane_v_axis=0,1,0
```

代码出处：

- `src/optics/OpticalTracer.cpp`：`u=dot(hit-plane_point, u_axis)`，`v=dot(hit-plane_point, v_axis)`；
- `src/app/OpticalSimCommon.cpp::applyCameraResponse()`：`camera_x_m=u`、`camera_y_m=v`；
- `src/io/LactEventRootWriter.cpp`：相机像素几何的 `x_m/y_m` 原样写入 ROOT。

`buildOutputPlane()` 现在强制 `u × v = normal`。它能直接拒绝旧的
`normal=-z, u=+x, v=+y` 左手组合，避免只改图中文字却不改模拟。

## 4. 东、西、上、下的物理符号检查

望远镜指北时：

| 天空光源 | 光子传播横向 | LACT 原始焦平面 | pyLAST 天空偏移 |
|---|---:|---:|---:|
| 东侧 | `+m_x`（向西传播） | `u < 0` | `pix_y > 0` |
| 西侧 | `-m_x`（向东传播） | `u > 0` | `pix_y < 0` |
| 指向上方 | `-m_y` | `v < 0` | `pix_x > 0` |
| 指向下方 | `+m_y` | `v > 0` | `pix_x < 0` |

这里“焦平面正负”由真实镜面求交和反射结果得到，不是绘图层人为翻转。
回归代码在 `apps/test_off_axis_orientation.cpp`。

## 5. pyLAST 边界

pyLAST 的 `AltAzFrame` 同样使用 NWU，且天空方位角从北向东为正。它的相机坐标是天空/source-offset，而 LACT ROOT 的 `u/v` 是硬件焦平面命中坐标，因此只在 ROOT 读取边界转换一次：

```text
pyLAST pix_x = -LACT v
pyLAST pix_y = -LACT u
```

随后像素编号、PE 和波形数组不重排。pyLAST 自己的绘图函数把 `pix_y` 放在屏幕横轴、`pix_x` 放在屏幕纵轴。

对应代码在 pyLAST main 的：

- `include/Coordinates.hh`；
- `include/CoordFrames.hh` 和 `src/CoordFrames.cpp`；
- `root/LactEventSource.cpp::load_telescopes()`；
- `src/pylast/visualize/visualize.py::plot_camera_image()`。

## 6. 必须通过的回归

- 普通平行光 frame 与 CORSIKA frame 三根轴逐项一致；
- 指北时 `m_x=West`、`m_y=sky-up`；
- 输出面满足 `u=-m_x`、`v=m_y`、`u×v=normal`；
- 东侧光源得到 `u<0`；
- pyLAST 东侧天空方向得到 `pix_y>0`；
- LACT ROOT 像素 `(u,v)` 读成 `(-v,-u)`；
- 平行光和 CORSIKA 通过同一相机像素几何，不能在画图脚本再做 LACT 专用翻转。
