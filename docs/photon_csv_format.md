# LACT_sim PhotonCsv Format

`PhotonCsv` is the common photon-table format used before the optical ray
tracer. External sources such as CORSIKA/EventIO adapters should convert their
output into this table first.

## Required Columns

```text
x_m,y_m,z_m,dir_x,dir_y,dir_z
```

- `x_m,y_m,z_m`: photon position in meters.
- `dir_x,dir_y,dir_z`: photon direction vector. It does not have to be exactly
  normalized; the C++ loader normalizes it before tracing.

## Optional Columns

```text
time_ns,wavelength_nm,weight,multiplicity,event_id,telescope_id,emission_altitude_km
```

If an optional column is absent, the source config supplies the default value:

```ini
time_ns=0
wavelength_nm=400
photon_weight=1
multiplicity=1
event_id=0
telescope_id=0
```

`weight * multiplicity` is applied to the photon before optical propagation.
`emission_altitude_km` is used by MODTRAN tau atmosphere tables; omit it when
atmosphere modeling is disabled or when the atmosphere config supplies an
explicit default emission altitude.

`eventio_2d` is also optional. For a six-column file cut from CORSIKA 2D
bunches, prefer keeping the CSV minimal and setting `source.eventio_2d=true`
in the cfg. LACT_sim then retains raw `z_m=0`, applies the normal local
`-16 m` EventIO input-plane default after coordinate rotation, and uses signed
ray intersection. Leave this setting off for literal hand-written positions.

## Coordinate Convention

推荐的默认输入坐标是 PhotonCsv 光学调试原来使用的望远镜本地坐标：

```ini
source.coordinate_frame=telescope_local
```

这保持原有 PhotonCsv 行为：CSV 行按配置的望远镜本地光学坐标解释，再通过
`buildTelescopeFrame()` 放入既有光追坐标。镜面、相机、输出和画图轴定义
均不改变。

还支持以下显式输入坐标：

```ini
source.coordinate_frame=corsika_nwu_relative
source.coordinate_frame=corsika_nwu_global
source.coordinate_frame=enu_east_relative
source.coordinate_frame=enu_east_global
source.coordinate_frame=lact_generic_global
```

- `corsika_nwu_relative`：CORSIKA NWU 向量，位置已经相对选定望远镜；
  与正常 EventIO photon bunch 相同。
- `corsika_nwu_global`：CORSIKA NWU 阵列绝对坐标；程序只减一次
  `telescope.position_m`。
- `enu_east_relative`：ENU 坐标，`+x=East,+y=North,+z=Up`，位置已经
  相对选定望远镜；方位角从 East 开始向 North 增加。
- `enu_east_global`：相同的 ENU 轴定义，但位置是阵列绝对坐标；程序只减
  一次 `telescope.position_m`。
- `lact_generic_global`：LACT 全局 NWU；`+x` 向北、`+y` 向西、`+z` 向上，方位角从北向东增加。

这些输入坐标、CORSIKA 默认流程和画图层方向的区别见
[中文坐标系说明](coordinate_systems_zh.md)。

```ini
source.local_telescope_frame=true
```

保留为兼容写法：`true` 等价于 `telescope_local`，`false` 等价于
`lact_generic_global`。新 cfg 应使用 `source.coordinate_frame`。

Raw CORSIKA/EventIO photon bunches should be read with `run_corsika_trace`.
Do not feed raw EventIO-derived rows into `run_optical_sim` as generic global
PhotonCsv unless they have first been explicitly converted into this LACT frame.

## Array Distribution

For multi-telescope input, include `telescope_id` in the CSV. During an array
run, `python/run_array_sim.py` writes:

```ini
source.filter_telescope_id=<current telescope_id>
```

so each telescope traces only its own rows.

For event-level slicing, `run_array_sim.py` can also write:

```ini
source.filter_event_id=<current event_id>
```

Run multiple events with:

```bash
python3 python/run_array_sim.py ... --event-ids 1,2,3
```

Each event is written into its own `eventNNN/` subdirectory, and the top-level
`array_run_summary.csv` combines all event/telescope rows.

## Normalizing External Files

Use `python/normalize_photon_csv.py` when an external file has different column
names or units:

```bash
python3 python/normalize_photon_csv.py external.csv photons.csv \
  --map x_m=x_cm,y_m=y_cm,z_m=z_cm,dir_x=ux,dir_y=uy,dir_z=uz,telescope_id=tel_id,event_id=event \
  --defaults wavelength_nm=400,weight=1,multiplicity=1,time_ns=0 \
  --scale-position 0.01 \
  --normalize-direction \
  --fail-on-nonfinite
```

The output can then be used with:

```ini
mode=PhotonCsv
csv_path=...
coordinate_frame=telescope_local
```

## ENU East-Start Input

PhotonCsv also accepts:

```ini
source.coordinate_frame=enu_east_relative
source.coordinate_frame=enu_east_global
```

Both use `+x` East, `+y` North, and `+z` Up. Telescope pointing azimuth starts
at East and increases toward North, so `pointing_az_deg=0` points East and
`pointing_az_deg=90` points North. The relative mode expects positions already
relative to the telescope. The global mode expects absolute ENU array
positions and subtracts `telescope.position_m` exactly once.
