"""SII的独立物理对照：不调用共享光电子对生成器来构造参考答案。"""
from pathlib import Path
import hashlib
import json

import numpy as np
from scipy.signal import fftconvolve
from scipy.stats import beta, chi2

from sii_unified import (
    Instrument, WaveformGLSCalibration, detected_star_rate_hz,
    load_measured_spe_template, load_empirical_charge_factors,
    load_optical_timing_mixture, waveform_instrument_signature,
    simulate_hbt_primary_pe, render_pe_waveform, waveform_cross_correlation,
)


def verify_main_parameters(root):
    """按已核对的main提交检查参数文件；后续实测数据更新时需要同步更新清单。"""
    root = Path(root)
    manifest = json.loads((root/'configs/sii/main_parameter_manifest.json').read_text(encoding='utf-8'))
    for entry in manifest['files']:
        content = (root/entry['path']).read_bytes().replace(b'\r\n', b'\n')
        if hashlib.sha256(content).hexdigest() != entry['sha256_lf']:
            raise ValueError(f"参数文件已改变，需要重新核对main来源：{entry['path']}")
    return manifest


def analytic_waveform_calibration(instrument: Instrument, magnitude=2.,
                                 block_duration_ns=20_000., max_lag_ns=200., fine_dt_ns=.05,
                                 residual_delay_ns=0.):
    """由连续SPE自相关、光学传递函数和Bartlett公式推导GLS模板与零信号协方差。

    这是线性探测、弱HBT、长于脉冲支持区间的块近似。电子学噪声可为白噪声，
    ADC剪裁、恢复、串扰和后脉冲不在此独立解析模型内，启用时必须另作标定。
    """
    if (instrument.adc_bits != 0 or instrument.microcell_recovery_time_ns != 0
            or instrument.sipm_crosstalk_probability != 0
            or instrument.sipm_afterpulse_probability != 0):
        raise ValueError('解析对照要求线性响应；请关闭非线性效应或使用波形标定')
    if (fine_dt_ns <= 0 or block_duration_ns <= 4*max_lag_ns
            or not np.isfinite(residual_delay_ns)):
        raise ValueError('细网格必须为正，标定块必须足够长')
    t, h = load_measured_spe_template(instrument.spe_template_path)
    fine_t = np.arange(t[0], t[-1]+fine_dt_ns/2, fine_dt_ns)
    pulse = np.interp(fine_t, t, h)
    autocorrelation = fftconvolve(pulse, pulse[::-1], mode='full')*fine_dt_ns
    correlation_time = np.arange(1-len(pulse), len(pulse))*fine_dt_ns
    pulse_square_area = float(autocorrelation[len(pulse)-1])
    charge = (load_empirical_charge_factors(instrument.charge_samples_path)
              if instrument.charge_samples_path else np.ones(1))
    second_moment = float(np.mean(charge**2))
    star = detected_star_rate_hz(magnitude, instrument)
    background = instrument.detected_nsb_rate_hz+instrument.dark_count_rate_hz
    total = star+background
    variance = total*1e-9*second_moment*pulse_square_area+instrument.electronic_noise_rms_mv**2
    # 频域乘|G|²等于把SPE自相关与两台望远镜的独立时间延迟差分布卷积。
    length = 1 << int(np.ceil(np.log2(4*len(pulse))))
    frequency = np.fft.rfftfreq(length, fine_dt_ns)
    transfer = np.ones(len(frequency), complex)
    if instrument.optical_timing_kernel_path:
        mixture = load_optical_timing_mixture(instrument.optical_timing_kernel_path)
        transfer = np.zeros(len(frequency), complex)
        for weight, mean, width in zip(mixture['weights'], mixture['mean_delay_ns'], mixture['std_delay_ns']):
            transfer += weight*np.exp(-2j*np.pi*frequency*mean-2*np.pi**2*frequency**2*width**2)
    transfer *= np.exp(-2*np.pi**2*frequency**2*instrument.intrinsic_time_jitter_ns**2)
    kernel = np.fft.fftshift(np.fft.irfft(
        abs(np.fft.rfft(pulse, length))**2*abs(transfer)**2, n=length))*fine_dt_ns
    kernel_time = (np.arange(length)-length//2)*fine_dt_ns
    dt = instrument.sample_width_ns
    lag_steps = np.arange(-int(max_lag_ns//dt), int(max_lag_ns//dt)+1)
    lags = lag_steps*dt
    span = int(np.ceil(np.ptp(t)/dt))+2
    steps = np.arange(-span, span+1)
    rho = total*1e-9*second_moment*np.interp(steps*dt, correlation_time, autocorrelation,
                                           left=0., right=0.)/variance
    rho[span] += instrument.electronic_noise_rms_mv**2/variance
    samples = block_duration_ns/dt
    # 块内减均值同时轻微降低信号和样本方差，保留其一阶有限块修正。
    normalization = 1-np.sum(np.maximum(samples-abs(steps), 0)*rho)/samples**2
    pair_rate = star**2*instrument.coherence_area_s
    # correlate(left,right)的峰位为-(右镜到达时间-左镜到达时间)。
    # 整数样本对齐后，把剩余分数时延留在模板中，避免对ADC做分数插值。
    sampled_signal = pair_rate*1e-9*np.interp(steps*dt+residual_delay_ns, kernel_time, kernel)/variance
    peak = (pair_rate*1e-9*np.interp(lags+residual_delay_ns, kernel_time, kernel)/variance
            -sampled_signal.sum()/samples)/normalization
    bartlett = fftconvolve(rho, rho[::-1], mode='full')/samples/normalization**2
    differences = lag_steps[:, None]-lag_steps[None, :]
    covariance = np.interp(differences, np.arange(1-len(rho), len(rho)), bartlett,
                           left=0., right=0.)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    eigenvalues = np.maximum(eigenvalues, eigenvalues.max()*1e-10)
    covariance = (eigenvectors*eigenvalues) @ eigenvectors.T
    calibration = WaveformGLSCalibration(
        lags, np.zeros(len(lags)), peak, covariance, block_duration_ns*1e-9,
        star, background, 1., 1., 0, 0, 0.,
        instrument_signature=waveform_instrument_signature(instrument), phase_model='linear_spe')
    diagnostics = dict(variance_mv2=variance, spe_area_mv_ns=float(np.trapezoid(pulse, fine_t)),
                       charge_second_moment=second_moment,
                       photon_degeneracy=star*instrument.coherence_area_s,
                       finite_block_normalization=normalization, fine_dt_ns=fine_dt_ns)
    return calibration, diagnostics


def thermal_mode_counts(rng, records, star_mean, background_mean, visibility2,
                        modes=8, polarizations=2):
    """独立复高斯热光模→光强→条件泊松计数；同时返回光强以验证Siegert关系。

    每条记录含modes×polarizations个独立模。此处是可精确验证的离散热光模型，
    不把纳秒采样窗口误称为单个光学相干模，也不把它用于整夜逐模仿真。
    """
    if (not isinstance(modes, int) or modes < 1 or polarizations not in (1, 2)
            or records < 1 or star_mean < 0 or background_mean < 0 or not 0 <= visibility2 <= 1):
        raise ValueError('热光模参数无效')
    count = modes*polarizations
    intensity = np.empty((records, 2))
    counts = np.empty((records, 2), dtype=np.int64)
    for first in range(0, records, 4096):
        last = min(first+4096, records)
        shape = (last-first, count)
        a = (rng.normal(size=shape)+1j*rng.normal(size=shape))/np.sqrt(2.)
        independent = (rng.normal(size=shape)+1j*rng.normal(size=shape))/np.sqrt(2.)
        b = np.sqrt(visibility2)*a+np.sqrt(1-visibility2)*independent
        intensity[first:last] = np.column_stack((np.mean(abs(a)**2, axis=1), np.mean(abs(b)**2, axis=1)))
        counts[first:last] = rng.poisson(star_mean*intensity[first:last]+background_mean)
    return counts, intensity


def proportion_interval(successes, total, confidence=.95):
    """二项覆盖率的Clopper–Pearson区间；报告区间，避免把有限重复比例当精确概率。"""
    tail = (1-confidence)/2
    return (0. if successes == 0 else float(beta.ppf(tail, successes, total-successes+1)),
            1. if successes == total else float(beta.ppf(1-tail, successes+1, total-successes)))


def profile_model_grid(uv, predictions, observations=None):
    """在预先计算的参数网格上剖面掉共享增益，返回每个模型的统计卡方加先验。

    predictions形状为(模型数,测量数)，可批量检验独立的观测实现。
    仅支持独立统计误差与一个共享高斯增益；不会把真值用作拟合初值。
    """
    prediction = np.asarray(predictions, float)
    measured = np.atleast_2d(uv.visibility_abs2 if observations is None else observations)
    if (uv.covariance is not None or prediction.ndim != 2
            or prediction.shape[1] != len(uv.sigma) or measured.shape[1] != len(uv.sigma)
            or np.any(~np.isfinite(prediction)) or np.any(~np.isfinite(measured))
            or np.any(~np.isfinite(uv.sigma)) or np.any(uv.sigma <= 0)
            or not np.isfinite(uv.calibration_relative_sigma) or uv.calibration_relative_sigma < 0):
        raise ValueError('参数网格需要有限模型、正的独立统计误差和有效共享标定先验')
    model = prediction/uv.sigma
    data = measured/uv.sigma
    model_square = np.sum(model**2, axis=1)
    cross = data @ model.T
    if uv.calibration_relative_sigma > 0:
        precision = 1/uv.calibration_relative_sigma**2
        gain = (cross+precision)/(model_square+precision)
    else:
        precision = 0.
        gain = np.ones_like(cross)
    statistic = (np.sum(data**2, axis=1)[:, None]-2*gain*cross
                 +gain**2*model_square+precision*(gain-1)**2)
    return statistic, gain


def profile_grid_interval(parameter, statistic, delta_chi2=3.841458820694124):
    """返回一维剖面区间，并报告截边或不连通；阈值是内部正则参数的渐近95%值。

    边界、低信噪比或模型选择问题不应直接套用此阈值，必须另做覆盖率校准。
    """
    parameter, statistic = np.asarray(parameter, float), np.asarray(statistic, float)
    if (parameter.ndim != 1 or statistic.shape != parameter.shape or len(parameter) < 3
            or np.any(np.diff(parameter) <= 0) or np.any(~np.isfinite(parameter+statistic))
            or not np.isfinite(delta_chi2) or delta_chi2 <= 0):
        raise ValueError('剖面区间需要递增有限参数网格和正阈值')
    difference = statistic-statistic.min()
    included = np.flatnonzero(difference <= delta_chi2)
    first, last = included[0], included[-1]
    low = (parameter[0] if first == 0 else np.interp(delta_chi2,
           difference[[first, first-1]], parameter[[first, first-1]]))
    high = (parameter[-1] if last == len(parameter)-1 else np.interp(delta_chi2,
            difference[[last, last+1]], parameter[[last, last+1]]))
    return dict(estimate=float(parameter[np.argmin(statistic)]), low=float(low), high=float(high),
                touches_boundary=bool(first == 0 or last == len(parameter)-1),
                disconnected=bool(np.any(np.diff(included) != 1)))


def standard_deviation_interval(values, confidence=.95):
    """正态零信号样本的标准差区间，仅适用于通过分布检查的块统计量。"""
    values = np.asarray(values)
    dof = len(values)-1
    tail = (1-confidence)/2
    return np.sqrt(dof*np.var(values, ddof=1)/chi2.ppf([1-tail, tail], dof))


def waveform_records(instrument, records, seed, visibility2=0., pair_scale=1.,
                     duration_ns=20_000., magnitude=2., max_lag_ns=200.):
    """产生专用于留出检验的波形相关记录；不拟合、选择或修正任何标定量。"""
    rng = np.random.default_rng(seed)
    template = load_measured_spe_template(instrument.spe_template_path)
    timing = (load_optical_timing_mixture(instrument.optical_timing_kernel_path)
              if instrument.optical_timing_kernel_path else None)
    padding = max(abs(template[0]))+5*(0. if timing is None else timing['rms_spread_ns'])
    padding += 5*instrument.intrinsic_time_jitter_ns
    correlations = []
    for _ in range(records):
        hits, _ = simulate_hbt_primary_pe(rng, duration_ns,
            detected_star_rate_hz(magnitude, instrument),
            instrument.detected_nsb_rate_hz+instrument.dark_count_rate_hz,
            visibility2, instrument, optical_timing_mixture=timing,
            padding_ns=padding, hbt_pair_rate_scale=pair_scale)
        signals = [render_pe_waveform(rng, hits.loc[hits.telescope_id == tel, 'time_ns'].to_numpy(),
                                      duration_ns, instrument, template=template)['adc_mv']
                   for tel in (0, 1)]
        lags, correlation = waveform_cross_correlation(*signals, instrument.sample_width_ns, max_lag_ns)
        correlations.append(correlation)
    return lags, np.asarray(correlations)
