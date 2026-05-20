# 遮挡模型 CSV 说明

本文档说明 `configs/obstructions/raytrace_final_structure_primitives.csv`
的格式和调参方法。这个 CSV 是当前 3D 遮挡模型的正式输入，供
`run_optical_sim` 和 `run_corsika_trace` 读取。

## 坐标和单位

- 单位：米。
- 坐标系：望远镜本地坐标系。
- 标准 1229 模型中，镜面附近约在 `z=-16 m`，白板/相机附近约在
  `z=-8 m`。
- CSV 文件头部的注释记录了原始模型到 LACT 坐标的转换，例如：

```text
coordinate_transform=swap source x/y, center on lens_reference, z(lens_reference)=-16 m
```

这表示当前版本已经完成了原始模型的 `x/y` 对调、单位 mm 到 m 的转换，以及
以 `lens_reference` 对齐镜面参考位置。

## 行类型

`type` 决定一行代表哪类几何体。

### cylinder

有限长圆柱，用于支撑杆。

必需字段：

```text
type=cylinder
name
role
x0_m,y0_m,z0_m
x1_m,y1_m,z1_m
radius_m
```

含义：

- `(x0_m,y0_m,z0_m)`：圆柱一端中心。
- `(x1_m,y1_m,z1_m)`：圆柱另一端中心。
- `radius_m`：圆柱半径。
- `role=support_strut`：支撑杆。反射段遮挡目前只考虑这类结构。

常见微调：

- 平移支撑杆：同时给 `x0/x1`、`y0/y1`、`z0/z1` 加同样偏移。
- 改粗细：改 `radius_m`。
- 改长度或方向：只改某一端的 `x/y/z`。

### box

轴对齐长方体，可选带正多边形孔。当前用于相机支撑连接块
`camera_support_gap_box`。

必需字段：

```text
type=box
name
role
center_x_m,center_y_m,center_z_m
half_x_m,half_y_m,half_z_m
```

可选孔字段：

```text
hole_radius_m
hole_rotation_rad
hole_sides
```

含义：

- `center_*`：box 中心。
- `half_x/y/z_m`：box 在三个方向的半宽、半高、半厚。
- `hole_radius_m`：中间孔的外接圆半径。
- `hole_sides=8`：八边形孔。
- `hole_rotation_rad`：孔在 `x-y` 平面内的旋转角。

光追判定时，光线先命中 box，再检查命中段是否穿过孔；孔内不算实心遮挡。

### polygon_prism

正多边形柱体。当前用于相机主体八边柱 `camera_body`。

必需字段：

```text
type=polygon_prism
name
role
center_x_m,center_y_m,center_z_m
radius_m
height_m
rotation_rad
sides
```

含义：

- `center_*`：棱柱中心。
- `radius_m`：正多边形外接圆半径。
- `height_m`：沿 `z` 方向的高度。
- `sides=8`：八边柱。
- `rotation_rad`：横截面在 `x-y` 平面内的旋转角。

## role 和方向规则

遮挡检查分两段：

1. 入射段：光源到镜片。
2. 反射段：镜片到白板/相机。

当前方向规则：

- 局部光线方向 `dir.z < 0`：认为光线从相机/天空侧打向镜面，检查所有实心结构。
- 局部光线方向 `dir.z > 0`：认为光线从镜面反射回白板/相机，只检查
  `role=support_strut` 的支撑杆。

因此：

- `support_strut` 会参与入射段和反射段遮挡。
- `camera_body` 和 `camera_support_gap_box` 只参与入射方向相关的遮挡，不参与
  镜面反射到白板/相机的遮挡。

## cfg 入口

遮挡 cfg 文件示例：

```ini
enabled=true
mode=primitives
primitives_csv=raytrace_final_structure_primitives.csv
check_incoming=true
check_reflected=true
```

字段：

- `enabled=true`：启用遮挡。
- `mode=primitives`：使用 3D primitives。
- `primitives_csv`：遮挡 CSV 路径。
- `check_incoming`：检查入射段遮挡。
- `check_reflected`：检查反射段遮挡。
- `mark_only=true`：调试模式。被遮挡光子只打标记，不丢弃。

## 推荐调参流程

不要直接覆盖标准文件，建议复制一份：

```bash
cp configs/obstructions/raytrace_final_structure_primitives.csv \
   configs/obstructions/raytrace_final_structure_primitives_tuned.csv
cp configs/obstructions/raytrace_final_structure.cfg \
   configs/obstructions/raytrace_final_structure_tuned.cfg
```

然后把 tuned cfg 改为：

```ini
primitives_csv=raytrace_final_structure_primitives_tuned.csv
```

单独画遮挡模型：

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_obstruction_primitives_3d.py \
  configs/obstructions/raytrace_final_structure_primitives_tuned.csv \
  --output run_logs/manual_checks/obstruction_tuned/layout.png
```

带镜片和白板一起检查：

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_optical_layout_3d.py \
  --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  --obstruction-primitives configs/obstructions/raytrace_final_structure_primitives_tuned.csv \
  --show-obstruction \
  --output run_logs/manual_checks/obstruction_tuned/layout_with_mirror.png
```

运行光追检查遮挡影响：

```bash
build/run_optical_sim configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  2>&1 | tee run_logs/manual_checks/obstruction_tuned/run.log
```

如果只想看哪些光子被遮挡，而不是直接丢掉它们，可以在 cfg 中加：

```ini
obstruction.mark_only=true
```

然后用 `python/plot_obstruction_marked_hits.py` 画被遮挡光子的分布。
