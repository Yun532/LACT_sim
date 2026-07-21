# 单台望远镜逐镜参数与统一参数的优先级

这套规则用于避免同一个物理量同时在镜面 CSV、`errors.cfg` 和
`efficiency.cfg` 中重复施加。固定偏转、随机误差和效率乘数分开处理。

## 推荐职责

- `mirror`：镜片中心、法向、曲率半径、口径等确定性几何；
- `error`：未逐镜测量时使用的统一随机误差；
- `efficiency`：随波长反射率曲线以及统一或逐镜反射率乘数。

固定逐镜指向偏转应直接合入 `normal_x/y/z`。`misalign_sigma_rad` 表示随机法向
误差的标准差，不是固定偏转量。

## 字段优先级

### 曲率半径

`radius_of_curvature` 可以出现在静态镜面 CSV 或仰角序列 CSV 中：

1. CSV 值大于 `0`：使用逐镜值；
2. CSV 缺列、留空或值为 `0`：使用 `mirror.radius_of_curvature`；
3. cfg 也未设置：使用程序默认值 `16 m`；
4. 负值或非有限值：直接报错。

因此理想光学仍可只写 `radius_of_curvature=16`；实测配置在 CSV 已完整给出
曲率半径时，不需要再写一遍统一值。

### 随机反射方向与随机法向误差

逐镜列分别为：

- `roughness_sigma_rad`：反射方向散射；
- `misalign_sigma_rad`：镜片法向随机偏差。

对应的统一参数分别为：

- `error.reflect_direction_sigma_deg`；
- `error.facet_normal_sigma_deg`。

每一列独立按整列判断：

1. 列不存在、全部留空或所有数值均为 `0`：使用 `errors.cfg` 的统一值；
2. 只要任意镜片值大于 `0`：启用逐镜模式，不再叠加统一值；
3. 逐镜模式中的 `0` 是该镜片明确的零误差；
4. 负值或非有限值直接报错。

这种设计兼容旧 CSV 中整列填零的占位写法，同时避免统一 sigma 和逐镜 sigma
重复计算。

### 反射率乘数

绝对反射率仍是
`efficiency.mirror_reflectivity` 指定的全局随波长曲线。无量纲逐镜乘数按以下
顺序选择：

1. `efficiency.mirror_reflectivity_scale_csv`：推荐的逐镜表，格式为
   `id,reflectivity_scale`，必须精确覆盖当前全部镜片；
2. 旧镜面 CSV 中显式提供的 `reflectivity_scale`：保留兼容；
3. `efficiency.mirror_reflectivity_scale`：统一回退值；
4. 均未设置时默认 `1`。

逐镜效率表存在时，它整体覆盖旧镜面 CSV 和统一回退值。`0` 表示该镜片完全
失效，是有效物理值，绝不触发回退；负值或非有限值直接报错。最终单镜反射率为

```text
全局随波长反射率 × 该镜片 reflectivity_scale
```

## 最小示例

理想镜面可以保持简洁：

```ini
mode=elevation_series
radius_of_curvature=16
```

实测镜面 CSV 只需保留真实存在的逐镜字段：

```text
elevation_deg,id,center_x,center_y,center_z,normal_x,normal_y,normal_z,radius_of_curvature
```

统一粗糙度写在 `errors.cfg`：

```ini
reflect_direction_sigma_deg=0.00725416907642
```

统一反射率乘数写在 `efficiency.cfg`：

```ini
mirror_reflectivity_scale=1
```

不要为了“字段齐全”向实测 CSV 增加无意义的全零列；程序仍能读取旧表中的这些
可选列，但新配置只写实际要逐镜控制的参数。
