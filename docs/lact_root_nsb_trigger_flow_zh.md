# LACT ROOT 输出中的 NSB 与触发流程

本文档描述 `run_corsika_trace` 当前用于 ROOT/pylast 输出的时间响应、NSB、
触发和写出顺序。

## 模型边界

当前 waveform 是 **pre-electronics p.e. proxy**，表示经过大气、镜面、光学
效率、light collector 和 PDE 后的离散 p.e. 到达时间序列。它不是 ADC
waveform，尚未包含 SPE 脉冲、增益涨落、transit-time jitter、afterpulse、
crosstalk、基线噪声、饱和及高低增益通道。

这个 p.e. 序列是未来真实电子学的输入边界，而不是需要被删除的临时格式。
真实电子学接入后应形成：

```text
Cherenkov/NSB p.e. arrivals
-> SiPM + front-end electronics response
-> analog/digital waveform
-> discriminator/trigger primitives
-> FADC/ADC output
```

`waveform.source=pe` 应继续保留，作为光学验证和电子学前后对照模式；真实
电子学使用独立的 `waveform.source=electronics` 实现。

## 时间触发与输出 profile 分离

代码分别判断：

```text
waveform_pe_available = waveform.enabled && waveform.source == pe
write_time_series     = profile 是 timeseries_pe 或 debug_full
evaluate_time_series  = waveform_pe_available &&
                        (write_time_series || trigger.enabled)
```

因此 `output.lact_profile` 只控制 ROOT 序列化内容，不改变物理触发：

- `image_pe`：按完整 p.e. 时间序列触发，但不写 `waveforms` 树。
- `timeseries_pe`：使用相同触发，并写稀疏 p.e. 时间序列。
- `debug_full`：使用相同触发，并写完整调试内容。

如果没有可用的 p.e. 时间序列，积分图像才会被视为单个时间 bin。

## Event 结束后的处理顺序

每个 event 的光学追迹先积累：

```text
summaries
pixels
waveforms / raw_waveform_hits
```

随后 ROOT writer 对每个候选望远镜执行：

1. 根据 `waveform.time_reference` 将 Cherenkov p.e. 放入时间 bin。
2. 对配置的全部像素和全部时间 bin 生成同一份确定性 NSB realization。
3. 将 Cherenkov 与 NSB 合并为最终 p.e. 时间序列和积分图像。
4. 使用 `trigger.camera_coincidence_window_ns` 做滑动相机触发。
5. 使用望远镜触发时间和 `trigger.array_coincidence_window_ns` 做阵列符合。
6. 最后应用 `output.save_only_triggered`。
7. 仅在 profile 要求时把稀疏 waveform 写入 ROOT。

NSB cell 由 `event_id + telescope_id + pixel + time_bin + seed` 决定，同一配置
下重复调用得到相同结果。ROOT 和 HDF5 因而可以使用同一物理 realization，
不会在触发和写出阶段重新抽取另一份 NSB。

## 为什么不再做 Cherenkov 峰附近的局部 NSB 预判

旧实现为了加速 `save_only_triggered=true`，只在 Cherenkov 峰附近抽取一段
NSB；局部预判失败时便不生成完整 NSB。这会漏掉其他时间位置的纯 NSB
accidental trigger，并让输出 profile 改变触发结果。

当前实现始终先评价完整时间范围，再筛选输出。代价是 NSB 开启时增加
`候选望远镜数 × 像素数 × 时间 bin 数` 量级的触发计算，但 `image_pe` 不会
因此写出 waveform，ROOT 文件尺寸仍保持轻量。

如果未来性能基准证明完整评价成本不可接受，可以增加明确标记的近似模式；
近似模式不得伪装成与严格 ROOT/HDF5 或 sim_telarray 触发等价。

## 当前仍未包含的内容

- 真实 SiPM/电子学脉冲响应与 ADC waveform。
- 电子学噪声、增益涨落、饱和、高低增益和硬件 discriminator。
- afterpulse、crosstalk 和 dark count。
- 逐像素、星场和离轴相关的 NSB rate。

这些内容应在统一电子学响应层实现，不应再次分别塞入 ROOT 和 HDF5 writer。
