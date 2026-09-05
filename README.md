# LACT 单像素恒星强度干涉模拟

这个分支只保留两条真正使用的链路：

1. `main`：C++ 光学追迹、SiPM/单 p.e. 波形、NSB、ADC 与输出；
2. `SII`：32 台望远镜的 UVW、带噪 `|V|²`、双镜波形和无相位图像重建。

强度干涉只使用每台望远镜的一个像素，不模拟完整相机图像。完整夜间观测也不会
逐采样点保存数小时波形：短波形用来标定光学和电子学响应，长曝光在平稳噪声等假设下
模拟相应的相关统计量。这是经过有限闭合检查的统计近似，不是数小时完整ADC波形验证。

## 文件在哪里

```text
apps/                         C++ 可执行程序入口
include/ + src/               main 的光学、电子学和输入输出主体
python/sii_unified.py         SII 源模型、UVW、噪声测量和一键流程
python/sii_reconstruction.py  独立的 |V|² 无相位重建
python/sii_validation.py      独立热光场、解析响应与参数区间验证
python/config_io.py           读取 main 的 cfg/CSV
configs/                      当前流程实际引用的配置和实测数据
notebooks/sii_complete_waveform_report.ipynb
                              当前主报告：波形GLS、独立验证、参数区间与图像稳定性
notebooks/lact_sii_paper_simulation.ipynb
                              保留的早期模拟报告，不代表当前主结果
docs/SII_PHYSICS_ZH.md         物理因果链、公式、输入来源与结果边界
docs/SII_IMPLEMENTATION_ZH.md  函数、数据接口与复现说明
SII_COMPLETE_WORKFLOW_ZH.md    当前说明及计算证据入口
tests/test_sii_*.py            SII自动检查
```

SII核心Python文件包含中文模块说明、关键物理步骤注释和函数docstring。
两份说明分别讲[物理原理](docs/SII_PHYSICS_ZH.md)和[代码实现](docs/SII_IMPLEMENTATION_ZH.md)。
Notebook保存实际执行的参数、独立波形验证、1000次参数区间覆盖率检验和图像重复结果。

## 物理流程

```text
双星参数 + 32 镜坐标
        ↓
地球自转投影：每条基线生成 (u,v,w)
        ↓
源模型得到理论 |V(u,v)|²
        ↓
main 参数：镜面/遮挡/滤光片/PDE/NSB/光学时间核/SPE/SiPM
        ↓
短波形标定 + 独立统计噪声 + 全体UV共享的标定增益误差
        ↓
模拟每个 UV 点的 |V|² 与不确定度
        ↓
正值、有限支撑、平滑正则、多起点无相位重建
```

两台望远镜实际记录的是单像素电流或 ADC 波形 `I1(t)`、`I2(t)`。用

```text
tau_g = dot(B12, source_direction) / c
```

校正几何时延后计算

```text
g2(0)-1 = <delta I1(t) delta I2(t)> / (<I1><I2>)
```

当前主报告使用多个延迟点的相关曲线，通过标定模板与噪声协方差做GLS估计，得到
该基线的 `P_hat` 与 `sigma_P`，不只是读取零延迟峰值。重建默认最小化绝对误差对应的
高斯负对数似然，并剖面拟合共同标定增益，预测值保留与观测相同的时间/光谱/UV平均；平滑强度由测量留出集选择。
上式给出通常的强度相关定义；程序实际相关向量以两条波形标准差的乘积归一化，
其到平方可见度的转换包含在标定模板中，不能把原始相关系数直接当成`g2-1`。
强度干涉不直接测傅里叶相位，非负、归一化和有限支撑等约束仍然必要，结果不保证唯一。
旧经验权重方法仅通过 `likelihood="legacy"` 用于历史对照。
当前主重建直接优化非负像素并归一化，默认使用`parameterization="flux"`和`1e-12`停止容差；
`parameterization="softmax"`保留旧参数化对照。数值停止、局部驻点和形态正确性分别检查。

## SII 公共接口

```python
from pathlib import Path
import sys

root = Path.cwd()
sys.path.insert(0, str(root / "python"))

from sii_unified import (
    Instrument,
    generate_uvw,
    simulate_waveform_gls_calibration,
    simulate_uv_observation,
    reconstruct_uv,
    run_sii_pipeline,
)

# 直接读取 main 当前使用的实测光学和电子学配置。
instrument = Instrument.from_repository(root)

# 实测SPE、光学核、星等或背景改变后，需要重新标定，不能复用旧误差。
calibration, _ = simulate_waveform_gls_calibration(
    instrument, source_ab_magnitude=source.ab_magnitude, seed=1
)
uvw = generate_uvw(layout, observation, instrument)
measurements, metadata = simulate_uv_observation(
    uvw, source, observation, instrument, seed=1,
    estimator="waveform_gls", waveform_calibration=calibration
)
image = reconstruct_uv(measurements)

# 或使用等价的一键流程：
result = run_sii_pipeline(
    layout, source, observation, instrument,
    estimator="waveform_gls", waveform_calibration=calibration
)
```

`generate_uvw`、`simulate_uv_observation` 和 `reconstruct_uv` 相互独立。只改变星等、
NSB或电子学参数时可以复用UVW；改变阵列、源赤纬、时角、时间分段或波长时重新生成。
相同轨迹的独立多夜重复可以合并统计量，但不能把改变曝光时长等同于增加新的UV覆盖。
上例假设已定义阵列、源和观测参数；可直接运行的完整设置见主notebook。
为兼容早期调用，省略`estimator`时仍走旧解析噪声分支，不等于当前报告的波形GLS流程。

## main 参数如何传到 SII

`Instrument.from_repository()` 默认读取
`configs/examples/corsika_lact_pylast_root_only_measured_waveform.cfg`，再沿组件配置读取：

| 内容 | 文件 |
|---|---|
| 32 台坐标 | `configs/arrays/layout_0803_reco32_coordinates.csv` |
| LACT2 实测镜面 | `configs/calibrated/lact2_measured_20260622/` |
| 单像素光学时间核 | `configs/optics/lact2_measured_single_pixel_400nm.csv` |
| 镜面、滤光片、PDE | `configs/efficiency/` |
| 暗夜 NSB 光谱 | `configs/nsb/` |
| 实测单 p.e. 波形和电荷涨落 | `configs/electronics/parameters/` |
| S17351 微单元几何 | `configs/electronics/devices/s17351_tiled_2x4.cfg` |

这些 cfg/CSV 改变后，下一次创建 `Instrument` 会重新读取，不需要同步修改 SII 代码。
当前没有实测来源的电子学量保留为接口且设为`0`，包括微单元恢复、本征抖动、串扰、
后脉冲、暗计数、额外电子噪声和仪器增益/时钟/透明度漂移。实测SPE长尾不被当成恢复时间。
获得实测值后通过`Instrument.from_repository(..., parameter=value)`覆盖并重新标定；
尚无时间/电荷模型的串扰、后脉冲及未标定GLS漂移若设为非零会明确报错。
2 nm通带和固定源透过率仍是场景假设。有限Monte Carlo标定误差来自程序估计，另行传播。

## 运行

Python 环境：

```bash
python -m pip install -r requirements.txt
python tools/execute_notebook.py notebooks/sii_complete_waveform_report.ipynb --cwd . --timeout 3600
python -m pytest -q tests
python tools/check_sii_science_artifacts.py
```

C++ main（Linux）：

```bash
make

# 单像素 400 nm 平行光，生成真实光学响应
build/run_optical_sim configs/optics/lact2_measured_single_pixel_400nm.cfg

# 两台望远镜的逐 p.e. CSV 可交给该电子学入口
build/run_camera_electronics \
  configs/sii/single_pixel_electronics.cfg \
  INPUT_PE.csv OUTPUT_DIR

# 完整 CORSIKA/EventIO 主流程
build/run_corsika_trace \
  configs/examples/corsika_lact_pylast_root_only_measured_waveform.cfg INPUT.zst
```

`make` 首次运行会从仓库内 `external/hessioxxx/hessioxxx.tar.gz` 解压并构建依赖；供应商
源码不再逐文件放在 GitHub 中。

## 当前结论的边界

NSB率随滤光片等设置变化，具体数值以主notebook开头的仪器参数表为准。
最长投影基线和波长给出条纹角尺度，但可可靠恢复的细节还取决于UV覆盖、SNR及先验；
检出极限受星等、积分时间、NSB、带宽和系统标定误差共同控制。
notebook 给出的是可复现的研究预测，不是望远镜已经验收的性能。滤光片角度响应、实际
过压/温度下暗计数、串扰、后脉冲和跨镜时钟稳定度有实测数据后，应通过现有接口补入。
