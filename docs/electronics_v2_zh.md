# 电子学 v2 链路

本实现从 `main` 的 `c75cb18` 开始。该基线已经包含坐标修复、实测约束光学、真实光收集器追迹和反射率边界修复；默认相机切换为已有的 1664 像素布局。

## 数据流

```text
CORSIKA/EventIO
  -> 镜面、遮挡、光收集器
  -> SiPM PDE（只应用一次）
  -> PrimaryPeHit
  -> 可选 NSB 合并
  -> 显式微单元（无恢复，每个微单元最多触发一次）
  -> FiredCellHit
  -> 可选单 p.e. 模板叠加
  -> 8 通道直接求和
  -> 4 ns 采样
  -> p.e. 计数或电压相机触发
  -> 阵列触发接口（默认关闭）
```

`uniform_interleaved` 是当前可替换的逻辑通道布线：一个 13 mm × 13 mm SiPM 使用 528 × 512 网格，共 270336 个微单元；`global_cell % 8` 给出通道号。获得真实器件布线后，只需替换映射，不改变后续接口。

当前明确不模拟恢复时间、串扰、后脉冲和暗计数。NSB 已可作为逐个 primary p.e. 在饱和前合并。

## 输出开关

```ini
electronics.pipeline.enabled=true
electronics.microcell.enabled=true
electronics.single_pe.enabled=true
waveform.enabled=true
waveform.source=electronics

electronics.output.save_primary_sequence=false
electronics.output.save_fired_sequence=false
electronics.output.save_microcell_decisions=false
```

- 关闭 `microcell.enabled`：最终积分图等于最初 primary p.e. 图。
- 开启微单元、关闭 `single_pe.enabled`：内部和可选输出都是饱和后的 fired p.e. 序列。
- 再关闭序列保存：ROOT 只保留积分图，不保存逐 hit 树。
- 开启 `single_pe.enabled`：生成采样波形；`unit=mV` 时使用电压触发。

详细序列默认不保存，以避免文件大小剧增。验证配置
`configs/examples/corsika_electronics_v2_validation.cfg` 才会显式开启全部诊断树。

ROOT 数据级别：

- `observations.image_cherenkov_pe`：饱和前、纯切伦科夫真值；
- `observations.image_fired_cherenkov_pe`：饱和后、纯切伦科夫真值；
- `observations.image_fired_pe` / `image_pe`：饱和后的最终探测器图；
- `primary_pe_hits`、`microcell_decisions`、`fired_pe_hits`：可选逐 hit 诊断；
- `waveforms.sample_value`：4 ns 采样值，单位见 `waveform_config.sample_unit`。

电子学波形的规范序列化目前使用 `lact_event` ROOT。若请求
`waveform.source=electronics` 与 HDF5 波形同时输出，程序会明确报错，不会写出空波形。

## 运行

独立电子学：

```bash
./build/run_camera_electronics \
  configs/examples/electronics_v2_voltage_trigger.cfg \
  configs/examples/electronics_v2_primary_hits.csv \
  validation_outputs/voltage_trigger
```

纯 NSB：

```bash
./build/run_camera_electronics \
  configs/examples/electronics_v2_nsb_only.cfg \
  none \
  validation_outputs/nsb_only
```

真实 CORSIKA gamma：

```bash
./build/run_corsika_trace \
  configs/examples/corsika_electronics_v2_validation.cfg \
  /path/to/gamma_eventio.zst
```

验证 notebook：

`validation/electronics_chain_v2/electronics_chain_v2_validation.ipynb`

第一格集中放路径、事件、望远镜和像素等可改参数；后续格只读取上述程序输出。
