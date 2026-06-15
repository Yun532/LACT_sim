# pylast event 数据层级与可视化取数说明

## 版本结论

当前仓库内用于验证和修改的是 `external/pylast`，其 `pyproject.toml`
声明版本为 `0.0.4`。`/Users/yun/Downloads/pylast-main` 声明版本为
`0.0.1`，且缺少当前 notebook 使用的 `pylast.visualize` 模块，因此两者
不是同一版本。

## ArrayEvent 结构

`pylast` 的单个事件是 `ArrayEvent`。主要字段如下：

| 字段 | 含义 |
| --- | --- |
| `simulation` | MC 真值和模拟相机信息，包括 shower 真值、每台望远镜的 true image、impact parameter 等。 |
| `r0` | 原始读出波形层，按望远镜保存双增益 waveform 和 waveform sum。 |
| `r1` | 校准后波形层，保存单个 waveform 矩阵和 gain selection。 |
| `dl0` | 从波形抽取出的相机图像层，保存每个像素的积分 image 和 peak time。 |
| `dl1` | 图像清理和参数化层，保存清理 mask、masked image 对应的 Hillas 参数等。 |
| `dl2` | 阵列级重建层，保存方向、芯位等重建结果。 |
| `pointing` | 望远镜阵列指向。 |

## 各层做了什么

### simulation

simtelarray 输入中，`event.simulation.shower` 保存真值：
`energy`、`alt`、`az`、`core_x`、`core_y`、`x_max`、`h_first_int` 等。
`event.simulation.tels[tel_id].true_image` 是模拟真值相机图像。

原始 `visual_new_diy.py` 的 `read_event_data()` 读取的就是这一层：
`event.simulation.tels[tel_id].true_image` 和 `true_image_sum`。

### R0

`event.r0.tels[tel_id]` 保存原始读出：
`waveform[0]`、`waveform[1]` 为两个增益通道的二维矩阵，形状通常是
`n_pixels x n_samples`；`waveform_sum[0/1]` 为每个像素的波形和。

### R1

`event.r1.tels[tel_id].waveform` 是校准后的波形矩阵，
`gain_selection` 记录每个像素选用的增益通道。simtel 示例中，触发望远镜
可从 `r1.tels.keys()` 得到。

### DL0

`Calibrator(source.subarray)(event)` 从 R0/R1 波形抽取 raw integrated image，
写入 `event.dl0.tels[tel_id].image` 和 `peak_time`。notebook 的
`plot_raw_images()`、`plot_event_cores()`、`plot_event_sdp_planes()` 默认读取
这一层，即 `image_level="dl0"`。

### DL1

`ImageProcessor(source.subarray)(event)` 对 DL0 图像做 tail-cuts clean，并计算
Hillas 参数，写入 `event.dl1.tels[tel_id]`：
`image`、`peak_time`、`mask`、`image_parameters`。

原始 `visual_new_diy.py` 的 `read_dl1_data()` 读取的是：
`event.dl1.tels[tel_id].image * event.dl1.tels[tel_id].mask`。
当前 `pylast.visualize.plot_clean_images()` 也是同样取数逻辑。

### DL2

`ShowerProcessor(source.subarray, config)(event)` 使用 DL1 Hillas 参数做方向和
芯位重建，结果在 `event.dl2.geometry["HillasReconstructor"]`。如果参与重建的
望远镜数量不足或质量筛选太严，`is_valid` 会是 `False`。

## 当前 pylast.visualize 的取数约定

`pylast.visualize.event_visualizer.read_event_data()` 通过 `image_level` 选择数据层：

| `image_level` | 读取数据 |
| --- | --- |
| `"simulation"` | `event.simulation.tels[tel_id].true_image` |
| `"dl0"` | `event.dl0.tels[tel_id].image` |
| `"dl1"` | `event.dl1.tels[tel_id].image * mask` |

notebook 里的默认绘图对应关系：

| 图 | 函数 | 数据层 |
| --- | --- | --- |
| 阵列 core 图 | `plot_event_cores()` | `dl0` |
| SDP 平面图 | `plot_event_sdp_planes()` | `dl0` |
| raw camera 图 | `plot_raw_images()` | `dl0` |
| clean camera + Hillas 图 | `plot_clean_images()` | `dl1` |
| 3D SDP | `plot_event_sdp_planes_3d()` / interactive | `dl0` + `dl2` |

## LACT ROOT adapter 测试记录

使用 `/Users/yun/Downloads/lact_events.root`、`run_logs/lact_root_only_full_response/lact_events.root`
和 `run_logs/lact_root_full_response/lact_events.root` 测试：

- `LactEventSource` 能正常读出 `simulation`、`r1`、`dl0`、`pointing`。
- ROOT 文件的 `observations.triggered` 分支存在，例如 event 103 的触发望远镜为
  telescope 1。
- 当前 `LactEventSource` Python wrapper 会统一提供 `get_triggered_tels(event)`：
  新 binding 优先使用 `simulation.triggered_tels`，旧本地 build 缺少该属性时才在
  source 层从 ROOT `observations.triggered` 兜底读取。notebook 和 visualizer 不再
  自己扫描 ROOT tree。
- LACT ROOT adapter 可以把 `image_pe` 复制到 `simulation.tels[*].true_image`
  作为兼容字段。`ImageProcessor` 默认保留 source/adapter 已经给出的
  `triggered_tels`；只有显式设置 `poisson_noise > 0` 时，才从 `true_image`
  加 Poisson noise 并用 fake trigger 重算触发望远镜。
- 为了和 simtelarray 的默认结构一致，LACT adapter 应默认只把
  `observations.triggered == true` 的 telescope 放进 `r1` 和 `dl0`；
  非触发 observation 可以留在 `simulation.tels` 中作为兼容/诊断信息。
