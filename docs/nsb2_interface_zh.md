# nsb2 与 LACT_sim 接口

本接口保持一条 `nsb2 -> PhotonCsv -> LACT_sim -> ROOT -> pyLAST` 数据链，
通过两种模式分别解决长曝光期望图和短时真实波形。

## 1. 长曝光期望图

```bash
python tools/nsb2_photoncsv.py \
  --mode expectation \
  --output-dir run_logs/nsb2_crab \
  --exposure-ns 1000000000 \
  --rays 1000000
```

生成的 PhotonCsv 不含 `time_ns`。固定数量的加权光线代表完整曝光时间，配合：

```ini
response.mode=expectation
electronics.enabled=false
waveform.enabled=false
trigger.enabled=false
```

运行 `configs/examples/crab_nsb2_expectation_root.cfg`。ROOT
`observations.image_pe` 是连续的期望 p.e.；除以 `config.integration_time_ns`
即可得到 p.e./ns，乘 `1e9` 得到 p.e./s。

## 2. 时间波形

```bash
python tools/nsb2_photoncsv.py \
  --mode timed \
  --output-dir run_logs/nsb2_crab \
  --integration-start-ns 0 \
  --integration-end-ns 2000 \
  --guard-pre-ns 250 \
  --guard-post-ns 250
```

此模式对入瞳光子数做 Poisson 抽样，每行都是独立物理光子，固定
`weight=1,multiplicity=1,origin=nsb`。运行
`configs/examples/crab_nsb2_waveform_root.cfg`。

保存窗口为 `[0,2000) ns`，但生成窗口为 `[-250,2250) ns`。guard 内光子参与
SPE 下降沿、SiPM 饱和与恢复，却不进入积分图像。稳定 NSB 必须使用
`waveform.time_reference=absolute`。

nsb2 已包含大气消光和散射，LACT_sim 应使用 `atmosphere/ideal.cfg`，避免重复
施加大气响应。

## 3. 分量选择

默认保存 nsb2 的所有分量。可使用：

```bash
--exclude-stars
--include-component TOKEN
--exclude-component TOKEN
```

所有实际采用的分量都会写入 manifest 与 rate cube。`--rate-cube` 可以复用
已有预测结果，只重新生成 PhotonCsv。

## 4. pyLAST

触发关闭时需显式读取保存的未触发数据：

```python
from pylast.io import LactEventSource

source = LactEventSource("lact_nsb2_waveform.root", read_untriggered=True)
event = source[0]
raw = source.get_raw_waveform(event, telescope_id=0)
```

原簇射工作流的默认 `read_untriggered=False` 不变。
