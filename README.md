# LACT 单像素恒星强度干涉模拟

使用LACT光学与电子学响应，生成双镜光电子波形、全阵列平方可见度测量，并进行GLS估计与图像重建。

- [物理原理与适用范围](docs/SII_PHYSICS_ZH.md)
- [代码实现、参数接口与复现](docs/SII_IMPLEMENTATION_ZH.md)
- [完整可执行Notebook](notebooks/sii_complete_waveform_report.ipynb)
- [验证结果](validation/sii_science/summary.json)

## 安装与运行

在仓库根目录执行：

```bash
python -m pip install -r requirements.txt
python tools/execute_notebook.py notebooks/sii_complete_waveform_report.ipynb --cwd . --timeout 3600
python -m pytest -q tests
python tools/check_sii_science_artifacts.py
```

需要与已验证环境保持相同直接依赖版本时，使用`requirements-sii-validated.txt`。

## Python接口

```python
from pathlib import Path
import sys
import pandas as pd

root = Path.cwd()
sys.path.insert(0, str(root / "python"))
import sii_unified as sii

# 读取本地仓库的仪器配置及实测响应。
instrument = sii.Instrument.from_repository(root)
source = sii.BinarySource()
observation = sii.Observation()
layout = pd.read_csv(root / "configs/arrays/layout_0803_reco32_coordinates.csv")

# 仪器响应、星等或背景改变后，需要重新标定。
calibration, _ = sii.simulate_waveform_gls_calibration(
    instrument, source_ab_magnitude=source.ab_magnitude, seed=1)
result = sii.run_sii_pipeline(
    layout, source, observation, instrument, seed=1,
    estimator="waveform_gls", waveform_calibration=calibration,
    do_reconstruction=False)

# 重建只使用测量、误差和采样信息。
image = sii.reconstruct_uv(
    result.measurements, grid_size=32, fov_mas=0.70,
    support_radius_mas=0.32, starts=3, max_iter=8000)
```

主要模块：`python/sii_unified.py`负责模拟和GLS标定，`python/sii_reconstruction.py`负责重建，`python/sii_validation.py`提供独立验证，`python/sii_observation.py`提供多镜共享波形、时延补偿及分块相关。

主Notebook第4.1节实际执行三镜观测链、同区间插值收敛与24/96 μs曝光对照；也可运行`python tools/validate_sii_observation.py`。结果保存在`validation/sii_observation/`，输入格式和适用范围见[实现说明](docs/SII_IMPLEMENTATION_ZH.md#41-多镜观测入口)。

## 参数与输出

默认仪器入口为`configs/examples/corsika_lact_pylast_root_only_measured_waveform.cfg`，由组件配置读取光谱响应、SPE、电荷分布和微单元几何。参数覆盖及数据格式见[实现说明](docs/SII_IMPLEMENTATION_ZH.md)。

缺失电子学参数保留接口，默认`0`且关闭；SPE长尾不作为微单元恢复时间。尚无实现或未标定的效应设为非零会明确报错。2 nm通带和固定源透过率属于场景假设。

Notebook输出保存在`validation/sii_science/`，图片位于`docs/sii_science_figures/`。短波形实际生成光电子与采样电压；长曝光使用经标定的统计模型，不生成整夜ADC。结果用于给定假设下的性能预测，重建的收敛状态和形态稳定性需分别检查。

## C++光学与电子学入口

在Linux环境构建并运行：

```bash
make
build/run_optical_sim configs/optics/lact2_measured_single_pixel_400nm.cfg
build/run_camera_electronics configs/sii/single_pixel_electronics.cfg INPUT_PE.csv OUTPUT_DIR
build/run_corsika_trace configs/examples/corsika_lact_pylast_root_only_measured_waveform.cfg INPUT.zst
```

首次`make`会解压并构建仓库内的`external/hessioxxx/hessioxxx.tar.gz`。
