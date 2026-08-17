# 实测单PE电子学链路

## 配置入口

- `configs/electronics/explicit_microcell_saturation_only.cfg`：正式的仅饱和基线。
- `configs/electronics/measured_spe_4ns.cfg`：微单元饱和、实测单PE模板、经验电荷涨落和4 ns mV采样。
- `configs/electronics/devices/s17351_tiled_2x4.cfg`：固定器件参数。通道数、总微单元数、逻辑网格和封装尺寸由程序自动推导。
- `configs/examples/electronics_measured_spe_4ns.cfg`：独立电子学验证入口。
- `configs/examples/photon_csv_measured_spe_4ns_validation.cfg`：1909号事例、19号望远镜的可重复光学到电子学回归入口。
- `configs/examples/photon_csv_saturation_off_waveform_off_validation.cfg`：原始PDE后PE基线。
- `configs/examples/photon_csv_saturation_on_waveform_off_validation.cfg`：只加入显式微单元饱和。
- `configs/examples/photon_csv_measured_spe_nsb_only_validation.cfg`：切伦科夫效率置零的NSB-only链路。
- `configs/examples/photon_csv_measured_spe_cherenkov_nsb_validation.cfg`：切伦科夫与NSB共同进入同一电子学链路。

所有正式文件使用`measured`命名，不在运行接口中暴露分析Notebook的历史版本号。

器件派生结果是8通道、33792微单元/通道、270336微单元/pixel、
`528 x 512`逻辑网格和`13.4 mm x 13.4 mm`封装范围。

## 当前包含和不包含的物理过程

包含：PDE后primary PE、显式微单元映射、无恢复硬饱和、实测单PE公共模板、逐雪崩经验电荷涨落、8通道直接求和、4 ns平均电压采样、PE-count或mV相机触发。

不包含：微单元恢复、基线涨落、电子学噪声、SiPM本征时间抖动、串扰、后脉冲、暗计数。NSB仍是独立的primary-PE来源，可单独开关；纯切伦科夫验证默认关闭。

## 输出量

逐pixel积分真值始终区分：

- `primary_*_pe`：PDE后、微单元饱和前；
- `fired_*_pe`：微单元实际触发次数；
- `gap_lost_pe`与`saturation_lost_pe`：几何gap和重复占用损失。

开启波形时，额外输出4 ns平均电压样本。ROOT `waveform_config`和HDF5 `/waveforms`元数据记录：

- `sample_unit=mV`；
- `single_pe_area_mv_ns=84.0349557248`；
- `template_time_reference=peak`；
- `charge_fluctuation_enabled=true`；
- `time_jitter_enabled=false`。

ROOT `waveform_config` 与 HDF5 `/waveforms` 还各自保存一份
`reference_pulse_time_ns` 和 `reference_pulse_amplitude`。它们是文件级标定元数据，
供 pyLAST 的局部积分窗口计算脉冲包含比例，不随事件、pixel 或望远镜重复保存。

验证模式可保存逐primary PE、逐fired PE、`charge_factor`、`time_jitter_ns`和8通道波形；正式模式默认不保存这些大体量中间量。

## pylast波形积分

LACT_sim的mV样本是每个4 ns区间的平均电压，所以：

```python
charge_mv_ns = waveform_mv.sum(axis=-1) * sample_width_ns
reconstructed_fired_pe = charge_mv_ns / single_pe_area_mv_ns
```

当前基线固定为0。实测电荷因子均值归一为1，因此波形重建PE在统计平均上等于`fired_pe`，但单个pixel会有物理涨落。这个标定不能把`fired_pe`自动反演成饱和前`primary_pe`。

NSB primary PE的生成范围会根据单PE模板的前后时间支撑自动向读出窗两侧扩展。
当前峰值对齐模板为`[-40, 180] ns`、保存窗为`[-40, 220) ns`，因此NSB在
`[-220, 260) ns`撒点，但ROOT/HDF5仍只保存原来的65个4 ns样本。窗外padding
hit参与微单元占用和波形卷积，却不计入窗口内的NSB积分图像，避免抬高PE真值。
靠近读出窗边界的单PE脉冲仍可能只在保存窗内留下部分电荷，因此NSB模式的
“所有fired hit电荷”不应被用作有限波形全积分的闭合真值；
无NSB的1909/19标准事例把完整脉冲包含在窗口内，用它执行严格积分闭合。

## 独立运行

```bash
build/run_camera_electronics \
  configs/examples/electronics_measured_spe_4ns.cfg \
  configs/examples/electronics_v2_primary_hits.csv \
  run_logs/electronics_measured_spe_4ns
```

## 命令行批量覆盖

生产程序支持与 sim_telarray 相同习惯的可重复 `-C key=value`。命令行值在主cfg和
所有组件cfg展开后再次应用，因此优先级最高；同一个键重复出现时最后一个值生效。

```bash
build/run_corsika_trace configs/examples/photon_csv_measured_spe_4ns_validation.cfg \
  -C electronics.microcell.saturation_enabled=true \
  -C electronics.single_pe.enabled=false \
  -C waveform.enabled=false \
  -C nsb.enabled=false \
  -C response.seed=1909 \
  -C output.lact_root_path=run_logs/batch/event1909_tel19.root
```

`-Ckey=value`、`--set key=value`和`--set=key=value`也是等价写法。组件入口也能覆盖；
组件路径与cfg内路径一样相对主cfg解析，例如此处主cfg位于`configs/examples`时使用
`-C electronics.config=../electronics/measured_spe_4ns.cfg`。这样批量脚本只需
保留一个基准cfg，不需要为每个随机种子、物理开关和输出路径复制文件。

显式S17351微单元几何依赖光收集器出口坐标。如果同时关闭`camera.collector`，程序现在
会在启动阶段报错，而不是把默认坐标`(0,0)`静默映射到中央channel gap。

## Notebook验证

LACT_sim三种格式一致性验证：

```bash
python scripts/validate_measured_electronics.py \
  validation/measured_electronics \
  --json validation/measured_electronics/validation_report.json
```

pyLAST功能分支中的`notebooks/lact_measured_electronics_validation.ipynb`读取这些实际
ROOT/HDF5/CSV输出，检查FullWaveFormExtractor、LocalPeakExtractor、五种模式、
NSB独立波形，以及1909号事例/19号望远镜的三层相机图。
