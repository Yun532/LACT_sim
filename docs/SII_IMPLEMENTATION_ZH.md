# SII程序实现、数据接口与复现

当前主入口是[`sii_complete_waveform_report.ipynb`](../notebooks/sii_complete_waveform_report.ipynb)。它从仪器配置开始，实际生成标定波形，做独立验证，再生成整阵列测量、参数区间和图像重建。本文说明程序如何计算；物理推导见[原理说明](SII_PHYSICS_ZH.md)。不要用历史Notebook的图片或早期标定数值替代当前结果。

## 1. 文件职责与调用关系

| 文件 | 负责什么 | 不应混用的层次 |
|---|---|---|
| [`python/sii_unified.py`](../python/sii_unified.py) | 读配置、源模型、几何、光子流、波形、GLS标定、长曝光数据和数据整理 | 波形GLS与旧解析SNR分支分别命名 |
| [`python/sii_reconstruction.py`](../python/sii_reconstruction.py) | 只根据测量和采样重建；处理共同增益先验、正则化与收敛诊断 | 真值只能在拟合完成后用于评价 |
| [`python/sii_validation.py`](../python/sii_validation.py) | 独立热光场、连续响应解析对照、留出波形、参数剖面和覆盖率工具 | 解析参考不调用共享光电子对生成器 |
| [`python/sii_layout.py`](../python/sii_layout.py) | 从TELESCOPE input读取原始NWU厘米坐标，转换ENU米并保留原值及行号 | 镜号按input顺序，标签保留；半径400 cm表示4 m半径 |
| [`python/sii_observation.py`](../python/sii_observation.py) | 天体逐波长相干、多镜共享事件及波形、时延追踪、分块相关、独立多镜热光模 | 主Notebook第4.1至4.3节实际执行；第5节整阵列长曝光仍是独立统计分支 |
| [`python/sii_performance.py`](../python/sii_performance.py) | 整数时延对齐、相位模板、段内误差积分、圆瞳对照和角直径似然 | 第9节最终性能入口；没有替换前面的探索性图像估计器 |
| [`tools/build_sii_science_notebook.py`](../tools/build_sii_science_notebook.py) | 用nbformat生成主Notebook的代码及中文解释 | 重新生成会清空该Notebook的执行输出，随后必须从头运行 |
| [`tools/execute_notebook.py`](../tools/execute_notebook.py) | 用nbclient在指定目录完整执行并保存Notebook | 单元报错时保留已执行输出，进程仍返回失败 |
| [`tests`](../tests)中的SII测试文件 | 物理恒等式、统计目标、数据独立性及数值边界检查 | 自动测试通过不等于所有科学主张通过 |
| [`validation/sii_science`](../validation/sii_science) | 本次计算结果表、标定、图像数组及汇总 | 图像重复、区间覆盖和波形检验各有不同样本含义 |

```text
Instrument.from_repository
    ├─ 模板/电荷样本/时间核/同通带光谱积分
    └─ simulate_waveform_gls_calibration
         ├─ simulate_hbt_primary_pe
         ├─ render_pe_waveform
         ├─ waveform_cross_correlation
         └─ 分组估计模板、协方差、增益和误差

run_sii_pipeline
    ├─ generate_uvw
    ├─ simulate_uv_observation(estimator="waveform_gls")
    │    ├─ segment_uv_samples + segment_sampling_weights
    │    ├─ source_visibility → 时间/通带平均的P
    │    └─ sample_waveform_gls_visibility2 → 随机测量
    └─ prepare_reconstruction_uv → UvData
         └─ reconstruct_uv_data
              ├─ power_sampling_kernel → 与数据相同的平均算子
              ├─ statistical_loss → 剖面共同标定增益
              ├─ 测量留出选择平滑项 + 多起点优化
              └─ 驻点、预测残差、峰数和拟合后评价
```

## 2. 从main读取的是什么

`Instrument.from_repository(ROOT)`读取**当前本地检出的配置文件**，不会在每次调用时联网取GitHub。此次已核对main提交`2926031b14c8aa2164a4d9233f5c2d0a75324127`中的10个关键参数文件，其内容与当前工作树一致。清单存于[`main_parameter_manifest.json`](../configs/sii/main_parameter_manifest.json)，使用统一换行后的SHA-256，避免Windows与Linux换行差异造成假变更。

另同步保留了main的[`spe_model_measured.yaml`](../configs/electronics/parameters/spe_model_measured.yaml)。运行时的数值来自cfg和CSV；YAML是SPE参数来源补充，不用其中的长尾解释替代微单元恢复测量。

默认响应入口为`configs/examples/corsika_lact_pylast_root_only_measured_waveform.cfg`，由既有`config_io`展开组件。数据流为：

1. 读取镜面、滤光片和PDE曲线，乘用户设置的SII通带，积分星光率、NSB率及相干面积。
2. 从电子学组件读取4 ns采样、SPE路径、电荷样本和微单元几何。电荷样本强制归一化为均值1，并由样本计算二阶矩。
3. 优先采用`lact2_measured_single_pixel_400nm.csv`和对应光追溯源JSON。有效探测面积已含光追损失，不重复乘PDE、镜面等因素。
4. 若完整光追核不可用，库保留旧核回退行为。因此正式报告另记录实际使用的输入哈希，不把回退核说成完整实测光学响应。
5. 最后应用调用者明确给出的覆盖参数。常数`source_transmission_scale=0.7836336`始终标为场景假设。

主Notebook记录源代码、布局和光追缓存哈希、随机种子及Python/NumPy/SciPy版本。光追缓存的原始逐光子文件位于历史溯源记录中的路径；本次没有重新运行百万光子的完整C++光追，不能将输入参数一致性核对说成重新验证了整条光学链。

### 2.1 阵列input的转换

当前阵列输入为[`lact36_20260906.input`](../configs/arrays/lact36_20260906.input)，来源是用户提供的36行TELESCOPE卡片，独立于main仪器参数清单。`read_corsika_layout`只读取TELESCOPE卡片，检查有限坐标、正半径和唯一标签；按出现顺序生成1起始镜号和0起始索引，保留原始NWU厘米、转换ENU米、半径与input行号。`generate_uvw`和`run_sii_pipeline`可直接接收`.input`路径，也兼容已有CSV及DataFrame。DataFrame中的非有限ENU位置会明确报错。

转换遵循[`EventIOPhotonSource.cpp`](../src/io/EventIOPhotonSource.cpp)的厘米到米，以及[`CorsikaTraceEventIOInput.cpp`](../src/io/CorsikaTraceEventIOInput.cpp)中的`(-y, x, z)`。例如Tel.1得到ENU=(-434.92063, 417.46032, 0) m，Tel.36得到(266.66667, 136.50794, 0) m；所有镜的半径为4 m。旧32台CSV是历史布局，当前入口不再默认使用它。需要重新导出可运行：

```powershell
python python/sii_layout.py configs/arrays/lact36_20260906.input configs/arrays/lact36_20260906_coordinates.csv
```

测试和结果检查脚本同时核对raw input与派生CSV，检查630条唯一基线和11340条分段测量，避免只更新坐标文件却继续使用旧阵列结果。

## 3. 缺失参数的零值约定

`Instrument`是不可变dataclass，修改使用`dataclasses.replace`或`from_repository(..., **overrides)`。字段名保留英文，物理用途和关键算法用中文注释。

| 接口 | 当前值 | 零值行为 | 将来非零时的处理 |
|---|---:|---|---|
| `electronic_noise_rms_mv` | 0 | 不产生加性噪声 | 当前支持独立白噪声；有色谱需扩展模型 |
| `dark_count_rate_hz` | 0 | 不增加暗计数 | 当前支持独立泊松暗计数并重新标定 |
| `microcell_recovery_time_ns` | 0 | 所有恢复比例严格为1 | 已有指数恢复模型；仅有独立实测恢复时间后启用 |
| `intrinsic_time_jitter_ns` | 0 | 不扰动光电子时刻 | 当前支持独立高斯时间抖动，单位ns |
| `adc_bits`、`adc_full_scale_mv` | 0、0 | 不量化、不剪裁 | 必须同时给有效位数及满量程，重新做波形标定 |
| `sipm_crosstalk_probability` | 0 | 无串扰 | 非零显式报未实现，等待实测时间/电荷模型 |
| `sipm_afterpulse_probability` | 0 | 无后脉冲 | 非零显式报未实现，等待实测延迟/幅度分布 |
| 每镜/每夜增益、基线零点、时延残差、透明度及NSB漂移 | 0 | 不额外扰动长曝光测量 | GLS分支对未标定非零项报错，避免套用平稳标定 |

这些0是“关闭效应”，不是伪造零误差测量。实测SPE、电荷涨落和光追展宽仍然存在。采样率、通带宽度、总光子率等必要物理输入不能用0占位。

最小使用方式：

```python
from pathlib import Path
from dataclasses import replace
import sys

ROOT = Path.cwd()  # 在仓库根目录运行
sys.path.insert(0, str(ROOT / "python"))
import sii_unified as sii

instrument = sii.Instrument.from_repository(ROOT)
# 有本征时间抖动实测后才填写；当前保持0。
instrument = replace(instrument, intrinsic_time_jitter_ns=0.0)
calibration, diagnostics = sii.simulate_waveform_gls_calibration(
    instrument, null_records=2048, signal_records=512, seed=20260905)
layout = ROOT / "configs/arrays/lact36_20260906.input"
result = sii.run_sii_pipeline(
    layout, sii.BinarySource(), sii.Observation(), instrument,
    estimator="waveform_gls", waveform_calibration=calibration,
    do_reconstruction=False, seed=20261305)
# 重建只传测量表；生成函数保留真值列供后验评价，整理函数不读取它。
image = sii.reconstruct_uv(result.measurements, grid_size=32,
    fov_mas=0.70, support_radius_mas=0.32, starts=3, max_iter=8000)
```

波形路径中的`electronics_bandwidth_hz`须与采样Nyquist频率一致，实际前端形状由SPE体现。`with_matched_effective_bandwidth`是旧解析SNR分支的辅助接口，不能先应用它再进入波形模拟，以免重复计算带宽损失。

## 4. 短波形与GLS标定的实现

`simulate_hbt_primary_pe`生成两镜的共享对与独立光电子，保持单镜平均率不变。`hbt_pair_rate_scale`只用于已知注入倍率，输出同时记录物理对率及注入对率。与main逐光电子接口兼容的`origin="cherenkov"`是原接口的信号标签，并不表示这里模拟的是切伦科夫光。

`render_pe_waveform`在每个SPE支持范围内计算连续平移模板，在ADC采样中心插值求和，用分块数组限制内存。它不先把光子时刻取整到采样格。微单元恢复关闭时不抽随机单元；电噪声关闭时不抽噪声数组，保持真正无效应的零值语义。

`waveform_cross_correlation`用FFT计算两条去均值波形的互相关，按重叠样本数修正滞后边缘，并按波形方差归一化。有限样本比值偏差通过解析对照和独立波形检查评估。

`WaveformGLSCalibration`包含以下数值，而不只是一个SNR常数：

| 字段 | 内容 |
|---|---|
| `lags_ns` | 完整相关峰的滞后坐标 |
| `null_mean` | 采用的零信号期望向量 |
| `peak_per_visibility2` | 单位平方可见度的峰模板 |
| `covariance_per_block`、`block_duration_s` | 单个标定块的协方差与时长 |
| `star_rate_hz`、`background_rate_hz` | 标定所对应的源与背景率 |
| `instrument_signature` | 仪器数值参数及响应内容签名 |
| `response_relative_uncertainty` | 固定选定权重下的有限响应增益不确定度 |
| `sigma_relative_uncertainty` | 噪声尺度估计的有限样本相对误差 |

`simulate_waveform_gls_calibration`的样本划分为：零信号一半训练、四分之一选收缩、四分之一定尺度；信号一半训练、四分之一校正响应、最后四分之一检验。候选收缩量为`1e-4, 1e-3, 1e-2, 5e-2`。模型假定平稳，将协方差按滞后差做Toeplitz平均，并给特征值设置相对`1e-10`下限。

`waveform_gls_weights`返回固定峰权重及条件统计误差。`sample_waveform_gls_visibility2`只生成GLS投影后的标量高斯统计量，并为整批数据抽一次共同增益。它不会生成6小时的250 MS/s数组。参数或光子率改变后，`simulate_uv_observation`检查签名和率，不匹配就要求重新标定。

保存的`waveform_calibration.npz`可用`np.load(..., allow_pickle=False)`读取，再将零维数组转为Python标量构造`WaveformGLSCalibration`。加载后仍需执行同样的仪器/率一致性检查；不能仅凭文件名断言有效。

### 4.1 多镜观测入口

`sii_observation.py`提供独立短记录入口，沿用main仪器响应。`simulate_array_photon_times`返回逐镜时间数组及对率、边缘率、填充信息；`simulate_array_waveforms`每镜只渲染一次，返回`adc_mv[镜,样本]`。它们接受复相干矩阵，先检查Hermitian、单位对角和半正定条件。弱对模型按每镜参与的全部对率扣除独立星光率，超额注入导致负率时明确报错。背景率显式传入，调用者将NSB与暗计数相加一次。

`geometric_arrival_delays_ns`返回相对参考镜的**到达时差**，即负的位置投影除以光速；`align_waveforms`执行`x_i(t+d_i)`。它返回共同有效样本中心、实际时长、首个输入索引和丢弃样本数。整数偏移精确索引，分数偏移默认16点半宽；改变半宽、采样率或时差后需匹配处理后的标定。没有共同数据时抛出异常。

分数插值使用有限线性相关实现FIR求和，交给SciPy选择直接计算或补零FFT，避免Python逐抽头循环。`mode='valid'`与显式公共裁剪确保输出只使用实际输入支持区间；这不是对原始记录作周期延拓。

`correlate_blocks`对所有镜使用同一组完整块，返回`block_correlations[块,基线,滞后]`、等曝光平均、基线顺序、块数和实际曝光；不足一块的尾数据计数后丢弃，恒定通道拒绝归一化。它不把相邻块自动视为统计独立，块相关须由重复数据检查。下面是可运行的最小例子，其中单位矩阵表示零跨镜相干的校验场景：

```python
import sys
from pathlib import Path
import numpy as np

root = Path.cwd()
sys.path.insert(0, str(root / "python"))
from sii_unified import Instrument, detected_star_rate_hz
from sii_layout import read_corsika_layout
from sii_observation import (simulate_array_waveforms, align_waveforms,
                             correlate_blocks, tracking_geometry)

instrument = Instrument.from_repository(root)
layout = read_corsika_layout(root / 'configs/arrays/lact36_20260906.input')
positions = layout.iloc[[0, 1, 5]][['east_m', 'north_m', 'up_m']].to_numpy()
state = tracking_geometry(positions, 0., .3, np.deg2rad(29.3586), elapsed_s=3600.)
delay_ns = state['arrival_delays_ns']
rate_ns_per_s = state['arrival_delay_rates_ns_per_s']
raw = simulate_array_waveforms(
    np.random.default_rng(42), 24_000., detected_star_rate_hz(2., instrument),
    instrument.detected_nsb_rate_hz + instrument.dark_count_rate_hz,
    np.eye(3), instrument, arrival_delays_ns=delay_ns,
    arrival_delay_rates_ns_per_s=rate_ns_per_s)
aligned = align_waveforms(raw["adc_mv"], delay_ns, instrument.sample_width_ns,
                          arrival_delay_rates_ns_per_s=rate_ns_per_s)
result = correlate_blocks(aligned["adc_mv"], instrument.sample_width_ns,
                          block_samples=aligned["adc_mv"].shape[1] // 4)
print(result["mean_correlation"].shape, result["effective_duration_s"])
```

`joint_thermal_mode_counts`是另外一条复高斯联合场→光强→泊松计数路径，用于验证二阶、三阶及边缘方差，不调用共享对生成器。它返回离散模平均的计数与光强，不输出真实纳秒波形。缺失时变增益、透明度等参数保持0；在当前联合波形入口设为非零会明确报尚未实现，避免被静默忽略。

### 4.2 随时角更新几何及块内时延率

`tracking_geometry`接受参考时角、赤纬、纬度、经过的SI秒和参考镜索引，返回该时刻的`arrival_delays_ns`、`arrival_delay_rates_ns_per_s`及`curvature_bound_ns_per_s2`。每条短记录的局部时间从0开始，避免用整夜绝对纳秒数参与插值而损失精度。曲率上界乘记录时长平方的一半，是一阶时延展开的误差界。

光子入口与补偿入口都接受可选的`arrival_delay_rates_ns_per_s`。默认省略或全0时沿用固定时延算法；非零时，光子参考时间先作线性伸缩和延迟，然后才抽光学时间响应。元数据中的`star_rate_hz`按参考时间定义，`received_star_rate_hz`明确包含时间映射的Jacobian；背景按接收时间生成。非有限时延率或非递增的时间映射会报错。

动态补偿按ADC半样本中心计算逐点读取位置，用同样的Lanczos-sinc权重，每批最多4096个输出点。公共裁剪确保每个通道的所有抽头都落在原始数据中；它不会把已计算的相关峰平移来冒充波形追踪。非零时延率下，部分通道即使恰好整数延迟，也按共同支持区间保守裁剪。

[`validate_sii_tracking.py`](../tools/validate_sii_tracking.py)在过中天前后3小时的七个时刻，分别生成96条训练、128条注入留出和128条零信号记录。每条24 μs，总计2464条三镜记录，原始共同时间累计0.059136 s；这不是6小时连续ADC，三个通道也不把曝光乘3。每个时刻的两种处理使用同一批ADC和共同参考时间区间。追踪响应只从该时刻训练集估计；初始时延对照除以追踪增益，不用接近零的自身响应放大噪声。

结果保存在[`validation/sii_tracking`](../validation/sii_tracking)：`geometry.csv`为36镜每分钟几何，`epochs.csv`为七个波形时刻的有效单条曝光及展开误差界，`records.csv`保留两种处理的逐条投影，`response.csv`记录训练增益、留出结果、零信号及误差。固定相干矩阵是处理检验场景，未随UV变化；该表不能直接作为恒星六小时光变或36镜重建测量。

本次共同裁剪后单条曝光为20.984至22.400 μs，2464条训练/留出/零信号记录合计有效时间0.054198144 s。统计总时长时按`epoch, case, record`去重，不能再次累加同一记录的三条基线。

### 4.3 从天体亮度构造通带内的联合相干

`source_coherence_spectrum`从实际ENU镜心位置形成`r_j-r_i`基线，调用已有`uvw_from_enu`和`source_visibility`，返回`coherence[波长,镜,镜]`、`wavelength_nm`、`spectral_weights`及`pair_visibility2[镜,镜]`。它支持已有双星、圆盘、椭圆和静态遮挡盘，沿用仪器的HBT谱积分节点；单色场景使用一个节点。每个波长的矩阵都检查Hermitian、单位对角和半正定。

`band_averaged_coherence_power`计算各节点模平方的加权和，不平均复振幅，也不把功率开方后伪装成相干矩阵。权重必须有限、非负且和为1。光子与波形入口新增可选`spectral_weights`，既有二维矩阵调用保持兼容；元数据保留`band_averaged_visibility2`。弱对生成器使用该功率决定对率，但仍不模拟多色热光的完整高阶时间过程。`joint_thermal_mode_counts`继续只接受单个合法矩阵。

下面续接第4.1节已经定义的`positions`、`state`、`instrument`及相关变量，改为实际圆盘源模型：

```python
from sii_unified import BinarySource
from sii_observation import source_coherence_spectrum

spectrum = source_coherence_spectrum(
    positions, state['hour_angle_rad'], .3, np.deg2rad(29.3586),
    BinarySource(primary_diameter_mas=.16), instrument, 'single_disk')
raw = simulate_array_waveforms(
    np.random.default_rng(43), 24_000., detected_star_rate_hz(2., instrument),
    instrument.detected_nsb_rate_hz + instrument.dark_count_rate_hz,
    spectrum['coherence'], instrument, arrival_delays_ns=delay_ns,
    arrival_delay_rates_ns_per_s=rate_ns_per_s,
    spectral_weights=spectrum['spectral_weights'])
aligned = align_waveforms(raw['adc_mv'], delay_ns, instrument.sample_width_ns,
                          arrival_delay_rates_ns_per_s=rate_ns_per_s)
```

主Notebook第4.3节实际运行`python tools/validate_sii_tracking.py --source-case single_disk --seed 20260908`，输出到[`validation/sii_source_tracking`](../validation/sii_source_tracking)。每个时刻用独立未分辨点源训练响应，`calibration_truth=1`；待测圆盘的`truth`只用于生成和拟合后检验，不参与求训练增益。点源为同星等、同通带、同几何状态的理想参考源，没有假称获得了实测标定星。

该目录沿用第4.2节的逐记录、响应和曝光表，并增加`source_predictions.csv`：七个时刻×630条基线共4410个36镜通带功率预测。`epochs.csv`中的`source_power_change_over_record`记录该24 μs窗口两端的最大功率差；信号强度在短记录内冻结，时延保留一阶变化。源矩阵、光谱权重与半正定检查不同于完整36镜波形或完整联合误差；这两项尚未由该表提供。

## 5. 测量表和重建器之间的边界

`generate_uvw`采用ENU坐标和西向为正的时角。时间推进使用`SIDEREAL_DAY_S`；曝光和协方差缩放使用秒。每段9个时间中点乘5个谱节点，组成45个实际子采样。

| 测量列 | 含义和使用范围 |
|---|---|
| `u_lambda`、`v_lambda` | 段中心、中心波长的代表坐标，用于索引和显示 |
| `uv_samples_u`、`uv_samples_v`、`uv_samples_weight` | 实际时间/光谱采样，用于正向平均与重建 |
| `visibility2_measured` | 随机测量，不做0至1裁剪 |
| `sigma_visibility2` | 独立统计标准差，用于逆方差合并 |
| `calibration_relative_sigma`、`calibration_id` | 共同增益先验及其标定身份 |
| `visibility2_true`、`visibility2_center` | 仅模拟及拟合后评价；不属于重建输入 |

`sigma_visibility2_calibration`是便于显示的逐行估算，**不能**代替共同协方差，也不能平方相加后假装所有基线独立。旧字段`segment_time_smearing_delta`目前包括时间与光谱平均相对中心值的合并差异。

`prepare_reconstruction_uv`以逆方差合并同一UV格的测量，统计误差按合并权重计算；把每一行内的真实采样权重乘该行在合并格中的权重，存入`UvData.sampling`。共同增益先验保持原值。它拒绝混合多个独立标定ID，因为那需要多组共同增益的似然，不能静默合成一个先验。

主Notebook导出的`binary_measurements.csv`为易读测量表，省去了变长子采样列。**不要直接用这张简表调用通用CSV重建入口来复现精确平均的主结果。** 主结果应从Notebook中的配置、观测及随机种子重新生成完整DataFrame。通用`read_uv_measurements`针对单点UV CSV，不会凭空恢复被省去的积分轨迹或共享标定先验。

## 6. 重建目标和独立性保护

`power_sampling_kernel`利用图像自相关表示平方可见度，提前平均每个像素位移处的余弦核。`power_from_image`用FFT得到图像自相关后乘该核；`_power_gradient`给出同一算子的伴随导数。算子本身只依赖几何和平均权重，不含源真值。

`statistical_loss`采用绝对误差，保留完整逆方差尺度。共同增益由`profile_calibration_gain`解析剖面消去。若直接提供完整`UvData.covariance`，代码检查对称性、正定性及其对角线与`sigma²`一致。当前自动交叉验证只支持独立统计误差加单一共同增益。

验证评分使用训练集的增益后验。若验证预测向量为$m$，统计对角阵为$D$，训练后的增益标准差为$s$，则验证预测协方差是$D+s^2mm^T$。评分包含二次型及其秩一行列式项，不在验证集重新拟合增益。

图像参数默认是有界非负流量再归一化，避免softmax让暗像素难以增长。旧`legacy`模式只用于历史算法对照，保留明确名称，不作为当前科学报告结果。多起点包括集中、弥散和随机平滑形态。

图像结果同时保存优化器状态、迭代次数、各起点目标函数、统计卡方、标定先验惩罚和单纯形驻点间隙。当前默认迭代上限为8000，达到上限仍报告失败，不把正常返回对象当成收敛。驻点阈值是明确数值诊断，不是全局最优证明。

拟合后图像配准只允许平移与180°中心反演，不使用周期卷绕。配准保留通量另作记录。`_peak_diagnostic`寻找局部极大值；单峰不输出第二颗“星”。真值只用于拟合后比较，相关自动测试还会污染真值列，确认测量整理不依赖它。

## 7. 独立验证为何不是程序自己证明自己

| 检查 | 独立输入或方法 | 能排除的问题 |
|---|---|---|
| `thermal_mode_counts` | 复高斯场→光强→条件泊松，独立于共享对生成器 | 热光二阶相关和超泊松方差公式/实现错误 |
| `analytic_waveform_calibration` | 连续SPE自相关、光学传递、Bartlett协方差 | 两个波形函数共享同一错误却互相闭合 |
| `waveform_records` | 固定新种子、完全独立记录，不训练权重 | 训练数据上的噪声低估和响应偏差 |
| 不同注入倍率、块长 | 重新生成波形而非缩放既有结果 | 注入非线性及块长归一化错误 |
| 时延率与变化插值 | 自转几何有限差分、独立解析信号的逆时间映射、共享光子逆映射 | 时延导数符号、ns/s单位、半样本位置或边缘率错误 |
| 天体联合相干 | 双点源场Gram矩阵、圆盘正亮度像素积分、单色与非相干谱叠加对照 | 基线相位符号、矩阵合法性及先平均振幅造成的功率错误 |
| 梯度有限差分 | 对标量目标直接作扰动 | 共同增益剖面导数或平均功率伴随错误 |
| `profile_model_grid` | 逐模型解析剖面，与显式似然优化对照 | 批量参数似然和区间实现错误 |
| 1000次直径区间 | 独立增益及噪声实现，报告二项区间 | 在指定场景中的区间覆盖不足 |
| 网格/合并/图像重复 | 固定测量改数值设置，或固定设置改噪声 | 数值误差与随机误差混淆、形态不稳定 |

独立解析对照只适用于线性探测和弱信号。ADC剪裁、微单元恢复、串扰和后脉冲开启时，它显式拒绝计算。白噪声和高斯本征抖动可在该解析对照中传播。解析网格0.05 ns与0.025 ns另作收敛检查。

圆盘剖面区间采用0.12至0.20 mas、801个网格点，并报告碰边和不连通情况。新36台布局的1000次检验为959次覆盖，比例0.959，95%二项区间约[0.9448, 0.9704]，与名义95%相容。噪声尺度改变约±6.13%后，覆盖率为0.942和0.971；有限标定误差仍会影响推断。

## 8. 如何完整复现

在仓库根目录运行，Python环境需包含[`requirements.txt`](../requirements.txt)中的依赖。本次使用Python 3.13.13、NumPy 2.4.4、SciPy 1.17.1；跨平台优化器停止位置可能有小差异，应比较误差及目标函数而不要求所有像素逐位相等。

[`requirements-sii-validated.txt`](../requirements-sii-validated.txt)固定本次直接依赖版本，可用于建立同版本环境；它不锁定操作系统、BLAS或全部传递依赖。`ipykernel`已列为显式依赖，避免干净环境只有nbclient却没有可执行的Python内核。

```powershell
python -m pip install -r requirements.txt
python -m pytest -q tests
python tools/build_sii_science_notebook.py
python tools/execute_notebook.py notebooks/sii_complete_waveform_report.ipynb --cwd . --timeout 3600
```

运行时间主要由几千条短波形和多次图像优化决定。Notebook固定随机种子，并将BLAS线程数设为1，避免小矩阵过度并行。正式运行不要同时修改导入模块，也不要手工编辑中间CSV后只执行末尾绘图单元。

主要输出为：

- `parameters.csv`、`thermal_modes.csv`：实际输入数值及独立热光检验。
- `response.csv`、`covariance.csv`、`waveform_calibration.npz`：本次GLS标定。
- `independent_null.csv`、`injection_scale.csv`、`block_duration.csv`：独立波形结果。
- `binary_measurements.csv`：供检查和画图的长曝光测量简表。
- `array_baselines.csv`：由36台input转换位置生成的630条唯一基线及长度。
- `diameter_profile.csv`、`diameter_coverage.csv`、`time_quadrature.csv`：参数区间、覆盖率和时间积分检验。
- `reconstruction.csv`、`convergence.csv`、`image_alignment.csv`、`image_repeats.csv`和各源`*_image.npy`：重建与稳定性。
- `summary.json`：实际运行版本、输入哈希、种子、所有核心结果及适用范围。
- `docs/sii_science_figures`：Notebook直接生成的科学图，不是手工示意图。

`tools/validate_sii_gls.py`是较小的独立标定/重建检查入口，可用于修改后的快速定位，不能替代主Notebook的完整证据。

主Notebook第4.1节通过当前Python解释器实际执行三镜观测验证，显示结果表与两幅验证图。该节的子进程失败会使Notebook失败；结果摘要的统一换行哈希随主摘要保存，避免旧图表与新计算混用。也可单独运行：

```powershell
python tools/validate_sii_observation.py
python tools/check_sii_science_artifacts.py
```

它默认生成256条训练记录、各384条注入留出/物理倍率/零信号记录，共1408条24 μs三镜原始波形；另生成384条96 μs物理倍率波形，并做30万条独立热光模计算。结果位于[`validation/sii_observation`](../validation/sii_observation)：`records.csv`保留固定解析权重的逐记录投影，`response.csv`保存各处理方法及基线的训练增益和标准误，`baseline_covariance.csv`保存基线联合协方差和近似相关区间，`thermal_moments.csv`为独立物理矩检验，`summary.json`记录输入、种子、代码哈希和适用范围。图由本次运行直接绘制。

`kernel_records.csv`保留16、32、64、128、256、512点半宽在相同19.232 μs区间内的逐记录投影；`kernel_convergence.csv`报告相对512点的配对响应差、均值标准误和噪声标准差比。比较时不对每个核单独除以新训练增益，否则会掩盖数值响应差。512点只是有限宽度参考，不能把参考曲线的零差、零误差棒解释为模拟真值精确已知。

`long_records.csv`保存独立生成的96 μs记录；`exposure_scaling.csv`将其噪声与原24 μs物理倍率样本比较。计算使用裁剪后的95.024和23.024 μs，标准差比乘曝光比平方根；表中95%区间来自近似正态块统计量的F分布。自动异常检查采用预先固定的双侧99.8%区间，图中仍展示较窄的95%区间，二者不能混称。通过这一检查不证明微秒到整夜的全部外推成立。

此入口采用固定解析GLS投影，再用独立训练数据校准补偿后的响应；尚未重新优化插值后的完整滞后权重。噪声尺度区间以训练增益固定为条件，增益误差在`response.csv`另报；注入恢复误差棒包含训练与留出两部分。物理倍率记录不与放大注入混合估计噪声，更不能用放大注入的基线相关外推整夜。处理链检验与主Notebook的长曝光灵敏度、重建结果有各自的代码哈希，保持清楚的证据范围。

主Notebook包含30个单元，其中16个代码单元；自动测试共98项。第4.2、4.3节分别执行固定相干与圆盘源的七时刻验证，第9节执行36镜性能验证，主摘要绑定这些结果摘要的哈希。完整执行状态及结果版本由下述检查脚本核对。
运行`python tools/check_sii_science_artifacts.py`可再检查说明链接、公式分隔符、Notebook执行状态和结果中的源代码/输入哈希。
新36台布局的36次图像重复均由优化器正常结束，并通过驻点阈值。逐次误差及配准保留通量仍公开在`image_repeats.csv`；通过数值诊断不等于图像唯一、全局最优或遮挡结构已检出。程序仍保留失败记录，不按成功状态筛选稳定性汇总。

## 9. 收到下一批真实参数时怎样更新

先保存现有分支或标签，将真实参数文件同步到工作树，记录main提交及文件哈希。`verify_main_parameters`会拒绝与旧清单不一致的文件，提醒来源已经改变；更新清单前应实际核对新文件，不是为使测试通过而覆盖哈希。

SPE形状、电荷分布、采样、通带、光子率或光学核变化后重新生成标定和全部相关结果。新增非线性/相关电子学数据时，先实现对应分布与时间结构，在短波形级验证后再允许长曝光使用。若是每镜各不相同的响应，需要按镜对或响应分组标定；当前示例是同型镜共享一套标定。

当前实现的适用条件及尚未验证的效应见[原理说明中的适用范围](SII_PHYSICS_ZH.md)。

## 10. 36镜角直径性能入口及SiPM规格书

本版固定交付范围是新36镜布局下的SII灵敏度与均匀圆盘角直径性能。完整Notebook实际执行这一验证集；只重算它时运行：

```powershell
python tools/validate_sii_performance.py
python tools/check_sii_science_artifacts.py
```

第二个命令检查整本Notebook及全部证据。修改数值代码或仪器输入后，需要重跑完整Notebook，不能仅更新摘要哈希。原始备份标签仍为`backup/sii-gls-before-fixes-20260905`。

### 10.1 从相位匹配的共享波形到长曝光

`integer_align`每镜仅平移最接近的整数样本、裁剪公共支持，返回镜对剩余分数时延。它不插值ADC、不周期卷绕。`analytic_waveform_calibration(..., residual_delay_ns=...)`将剩余时延放入连续相关峰模板，以相同零信号协方差求GLS权重。`phase_template_bank`使用65个相位节点，覆盖−4至4 ns。

正式验证用全部36镜、630条基线，在三个时角各生成96条点源训练、192条独立零信号和96条圆盘留出记录。共1152条24 μs原始记录，原始时间0.027648 s，公共有效时间0.026164224 s；它们不是连续6小时ADC。点源和圆盘注入倍率300，只用于响应检验；物理零信号数据用于噪声学习，放大注入协方差不进入最终噪声模型。圆盘真值只用于生成和评价，不用于校准。

630条基线共享同一记录的36通道，训练误差因此按每条记录的阵列均值计算，不能当630次独立训练。统一增益为1.00480，共享相对不确定度1.0804%；零信号标准差相对解析预测为1.00628，学习不确定度0.2621%。跨三个时刻的离散也保守计入学习误差。24个相位组的圆盘留出及零信号最大标准化偏差分别为2.341、2.267。预定5个标准误阈值是异常检查，不是逐项95%区间。

`tracked_segment_precision`沿实际基线、恒星日时角和段内相位积累方差。每段等时长短块平均与现有等时间可见度模型对应，因此使用平均块方差；不能改为逆方差平均却仍拟合等时间模型。1200个相位积分节点加倍到2400后，标准差最大变化0.00187%。连续观测假设缓冲区覆盖整数时延，不把独立测试记录的首尾损失重复扣在实际每个ADC块上。

`compressed_profiles`只用候选圆盘模型建立QR子空间，抽样其中的高斯充分统计量；检查每个候选的表示误差。正交补对固定噪声尺度的所有候选只贡献共同卡方常数，因而可以消去。真值不参与基底选择。共同增益逐次抽样并剖面；噪声学习误差作为一次观测共享的正值尺度抽样。拟合误差保守增加1.96倍噪声尺度标准误，并乘弱光Bartlett协方差界对应的标准差修正。

候选直径固定为0至0.64 mas、间隔0.0005 mas；平滑似然细化到0.00005 mas。33个预定非网格中点用直接模型检查，整体模型差最大约为统计误差的0.00000111。27场景各500次，共13500次观测实现；`profile_gain`和`fixed_gain`在每个场景使用相同模拟数据。边界及不连通诊断保留，不筛除困难结果。

| 输出，位于`validation/sii_performance` | 内容 |
|---|---|
| `raw_projections.npz` | 每个时刻训练、零信号、留出的逐记录630基线投影，及相位、真值和解析误差 |
| `raw_epochs.csv`、`raw_baselines.csv`、`raw_phase_holdout.csv` | 标定与误差学习、1890条基线结果、相位分组留出 |
| `phase_templates.csv` | 剩余相位与24 μs条件误差 |
| `performance.csv` | 27场景×2方法的偏差、RMSE、覆盖率和边界诊断 |
| `covariance_budget.csv`、`approximation_budget.csv` | 共镜修正、热光计数四阶对照、圆瞳及数值积分对照 |
| `datasheet_scenarios.csv` | 独立的25℃暗计数情景；不改变默认仪器 |
| `summary.json`、`performance.png` | 参数、种子、输入/代码/产物哈希及性能图 |

产物哈希按实际字节保存，本次Windows生成的性能CSV通过`.gitattributes`固定为CRLF换行，以保持检出后的发布文件与记录一致。若外部软件改写文件，需要重新生成这一验证集，不能忽略不匹配。数值源代码和配置另用LF归一化检查。

### 10.2 S17351规格书能补什么

用户提供的`K30-B60168_S17351 specification sheet.pdf`内部标题为S17351，日期2025-01-29，标为PRELIMINARY。已核对第1页表格及第2页脚注；结构化摘录和原PDF SHA-256在[`s17351_datasheet.json`](../configs/sii/s17351_datasheet.json)。原PDF保留为本地来源，运行程序不依赖重新分发该文档。

| 参数 | 规格书值与条件 | 程序处理 |
|---|---|---|
| 通道及微单元 | 8通道，每通道33792微单元；25 μm节距 | 总数270336与main一致 |
| PDE | 405 nm典型34%、450 nm典型37%，过压8.5 V | 用于核对，保留main完整曲线 |
| 暗计数 | 25℃每通道典型1.2 MHz、最大3 MHz | 独立计算1或8通道汇总情景；默认仍为0 |
| 瞬时串扰 | 25℃典型3%、最大5% | 保存数值，不自动叠加到实测电荷涨落上 |
| 后脉冲 | 典型5%，条件为−10℃ | 未给25℃值及延迟、电荷分布，保持关闭 |
| 雪崩增益 | 典型1.1×10⁶ | 不是mV电压增益，不替换实测SPE |
| 端电容 | 每通道750 pF，100 kHz测试 | 缺少淬灭电阻等信息，不能推出恢复时间或前端带宽 |
| 电压及温度系数 | 典型击穿52 V；工作电压为各器件击穿电压+8.5 V；54 mV/℃ | 记录条件，60.5 V不是所有器件统一适用的偏压 |
| 增益均匀性、稳定性 | 最大相对偏差典型1.7%、最大2.7%；24小时变化上限0.5% | 不是高斯RMS，不直接作为随机增益参数 |

PDF注明PDE不含串扰与后脉冲。现有电荷分布是否已包含这些事件，需要看原测量阈值、积分窗和事件筛选；不能因新增一个概率就再乘一次超额噪声。main在405 nm的曲线插值约35.56%，与表格典型34%不完全相同；这是不同输入的差异，不凭规格书典型点覆盖已有完整曲线。

```python
from sii_performance import datasheet_dark_scenario
# 假设25℃、过压8.5 V、8通道汇总，不表示已确认真实SII接线。
scenario = datasheet_dark_scenario(instrument, root, summed_channels=8, rating='typical')
assert scenario.dark_count_rate_hz == 9.6e6
# 背景改变后须按scenario重新标定，不沿用原缓存。
```

### 10.3 尚缺的电子学信息按六组收集

这些是六组数据，不是六个独立数字；部分需要曲线或波形。

| 组 | 需要的数据 | PDF是否补足 |
|---|---|---|
| ADC | 实际位数、输入满量程、偏置/削顶范围，必要时非线性标定 | 未提供；4 ns采样间隔已在main |
| 前端噪声与响应 | 遮光基线波形、噪声频谱与跨通道相关，确认SPE对应的前端设置 | 未提供；已有SPE可用，不重复加入同一脉冲响应 |
| 同步与时延 | 采样时钟抖动、通道间残余延迟和漂移 | 未提供 |
| 微单元恢复 | 双脉冲间隔—恢复电荷曲线、饱和行为 | 微单元数和端电容不足以确定恢复时间 |
| 相关雪崩 | 工作条件下串扰/后脉冲的概率、延迟与电荷分布，原电荷样本筛选规则 | 补了部分概率，未补时间结构及重复计数判定 |
| 运行条件与通道映射 | 实际温度、各器件偏压/过压、1或8通道如何组成SII通道，及该条件下暗计数/增益变化 | 只有厂家典型条件及上限，无LACT运行记录 |

优先提供实际温度/偏压及接线方式、遮光基线波形、ADC位数和满量程，即可先约束目前对灵敏度最直接的缺口。默认缺失项继续0/关闭；规格书情景与实测配置分开保存。
