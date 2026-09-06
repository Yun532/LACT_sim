"""多镜共享光子/波形、几何时延补偿与分块相关；短记录使用冻结几何。"""
from itertools import combinations

import numpy as np
from scipy.signal import correlate

from sii_unified import (
    C_M_S, source_direction_enu, load_measured_spe_template,
    load_optical_timing_mixture, sample_optical_delays_ns,
    render_pe_waveform, waveform_cross_correlation,
)


def _coherence_matrix(coherence):
    """复相干矩阵必须Hermitian、半正定且对角为1；任意镜对功率未必可联合实现。"""
    gamma = np.asarray(coherence, complex)
    if (gamma.ndim != 2 or gamma.shape[0] != gamma.shape[1] or len(gamma) < 2
            or not np.all(np.isfinite(gamma))
            or not np.allclose(gamma, gamma.conj().T, atol=1e-12, rtol=0)
            or not np.allclose(np.diag(gamma), 1., atol=1e-12, rtol=0)
            or np.linalg.eigvalsh(gamma).min() < -1e-12):
        raise ValueError('需要至少两镜的有效复相干矩阵：Hermitian、半正定、单位对角')
    return gamma


def _rates(value, size):
    rates = np.broadcast_to(np.asarray(value, float), (size,)).copy()
    if np.any(~np.isfinite(rates)) or np.any(rates < 0):
        raise ValueError('每镜光子率必须有限且非负')
    return rates


def geometric_arrival_delays_ns(positions_enu_m, hour_angle_rad, dec_rad, lat_rad,
                               reference=0):
    """源方向s指向天体，入射波传播方向为-s；到达时差为-(r_i-r_ref)·s/c。"""
    positions = np.asarray(positions_enu_m, float)
    if (positions.ndim != 2 or positions.shape[1] != 3 or len(positions) < 2
            or not isinstance(reference, (int, np.integer)) or not 0 <= reference < len(positions)
            or not np.all(np.isfinite(positions))
            or not np.all(np.isfinite([hour_angle_rad, dec_rad, lat_rad]))
            or abs(dec_rad) > np.pi/2 or abs(lat_rad) > np.pi/2):
        raise ValueError('需要有限ENU坐标、有效角度及参考镜索引')
    direction = source_direction_enu(hour_angle_rad, dec_rad, lat_rad)
    return -(positions-positions[reference]) @ direction / C_M_S*1e9


def simulate_array_photon_times(rng, duration_ns, star_rate_hz, background_rate_hz,
                               coherence, instrument, arrival_delays_ns=None,
                               pair_rate_scale=1., padding_ns=200.):
    """稀疏HBT对过程的多镜扩展：每镜只有一个事件流，由所有相关基线共享。

    镜对率为R_i R_j tau_eff |Gamma_ij|²。各镜独立星光率减去其参与的全部
    镜对率，保持边缘平均光子率。此近似不包含热光自聚束和三阶闭合相位项，
    放大注入仅用于处理链检验；不能把其高阶统计当作完整热光模型。
    """
    gamma = _coherence_matrix(coherence)
    size = len(gamma)
    star, background = _rates(star_rate_hz, size), _rates(background_rate_hz, size)
    delays = np.zeros(size) if arrival_delays_ns is None else np.asarray(arrival_delays_ns, float)
    if (delays.shape != (size,) or not np.all(np.isfinite(delays))
            or not np.all(np.isfinite([duration_ns, padding_ns, pair_rate_scale]))
            or duration_ns <= 0 or padding_ns < 0 or pair_rate_scale <= 0
            or not np.isfinite(instrument.coherence_area_s) or instrument.coherence_area_s < 0):
        raise ValueError('无效的时长、填充、到达时差或注入倍率')
    pair_rates = star[:, None]*star[None, :]*instrument.coherence_area_s*abs(gamma)**2
    np.fill_diagonal(pair_rates, 0.)
    injected = pair_rates*pair_rate_scale
    single_rates = star-injected.sum(axis=1)
    if np.any(single_rates < 0):
        raise ValueError('全部相关对率之和超过某镜星光率；必须降低注入倍率')
    # 扩大潜在事件区间后再延迟，避免记录首尾出现人工缺光。
    timing = (load_optical_timing_mixture(instrument.optical_timing_kernel_path)
              if instrument.optical_timing_kernel_path else None)
    timing_guard = (0. if timing is None else
                    float(np.max(abs(timing['mean_delay_ns'])+8*timing['std_delay_ns'])))
    guard = padding_ns+float(np.max(abs(delays)))+timing_guard
    start, stop = -guard, duration_ns+guard
    seconds = (stop-start)*1e-9
    streams = [[rng.uniform(start, stop, rng.poisson(rate*seconds))] for rate in single_rates]
    pair_counts = np.zeros((size, size), dtype=np.int64)
    for left, right in combinations(range(size), 2):
        count = rng.poisson(injected[left, right]*seconds)
        centers = rng.uniform(start, stop, count)
        streams[left].append(centers)
        streams[right].append(centers)
        pair_counts[left, right] = pair_counts[right, left] = count
    times = []
    for index in range(size):
        # 光学飞行时间独立作用于每个光子；同一事件在不同基线中不再重新抽样。
        stellar = np.concatenate(streams[index])+delays[index]
        if timing is not None:
            stellar += sample_optical_delays_ns(rng, len(stellar), timing)
        sky = rng.uniform(start, stop, rng.poisson(background[index]*seconds))
        times.append(np.sort(np.concatenate((stellar, sky))))
    return times, dict(duration_ns=float(duration_ns), guard_ns=guard,
                       arrival_delays_ns=delays, star_rate_hz=star,
                       background_rate_hz=background, physical_pair_rates_hz=pair_rates,
                       injected_pair_rates_hz=injected, single_star_rates_hz=single_rates,
                       pair_counts_with_padding=pair_counts, pair_rate_scale=float(pair_rate_scale),
                       model='weak_pair_clusters_no_thermal_higher_orders')


def simulate_array_waveforms(rng, duration_ns, star_rate_hz, background_rate_hz,
                             coherence, instrument, arrival_delays_ns=None,
                             pair_rate_scale=1.):
    """按main响应生成一次逐镜ADC；背景率应显式包含所需的NSB和暗计数。"""
    for name in ('telescope_gain_calibration_rms', 'per_night_gain_rms',
                 'baseline_zero_point_rms', 'residual_timing_rms_ns',
                 'transparency_fractional_rms', 'nsb_fractional_rms'):
        if getattr(instrument, name) != 0:
            raise NotImplementedError(f'联合短波形尚未实现时变状态：{name}')
    template = load_measured_spe_template(instrument.spe_template_path)
    padding = float(np.max(abs(template[0])))+8*instrument.intrinsic_time_jitter_ns
    times, metadata = simulate_array_photon_times(
        rng, duration_ns, star_rate_hz, background_rate_hz, coherence, instrument,
        arrival_delays_ns, pair_rate_scale, padding)
    waveforms = [render_pe_waveform(rng, events, duration_ns, instrument, template=template)
                 for events in times]
    return dict(adc_mv=np.stack([item['adc_mv'] for item in waveforms]),
                sample_time_ns=waveforms[0]['sample_time_ns'],
                sample_width_ns=instrument.sample_width_ns, metadata=metadata)


def align_waveforms(adc_mv, arrival_delays_ns, sample_width_ns, half_width=16):
    """取y_i(t)=x_i(t+d_i)，用有限Lanczos-sinc插值并裁到全镜共同有效区。

    不循环卷绕、不补零进入相关。整数延迟精确索引，分数延迟按带限近似插值。
    真实SPE未必严格带限，插值后的模板和噪声必须按同一处理重新标定。
    """
    values, delays = np.asarray(adc_mv, float), np.asarray(arrival_delays_ns, float)
    if (values.ndim != 2 or values.shape[0] < 2 or values.shape[1] < 2
            or delays.shape != (values.shape[0],) or not np.all(np.isfinite(values))
            or not np.all(np.isfinite(delays)) or not np.isfinite(sample_width_ns)
            or sample_width_ns <= 0 or not isinstance(half_width, int) or half_width < 2):
        raise ValueError('需要有限逐镜波形、每镜时延、正采样间隔及至少2点的插值半宽')
    offsets = delays/sample_width_ns
    kernels = []
    low, high = 0, values.shape[1]-1
    for offset in offsets:
        nearest = int(np.rint(offset))
        if abs(offset-nearest) < 1e-12:
            taps, weights = np.array([nearest]), np.ones(1)
        else:
            taps = int(np.floor(offset))+np.arange(1-half_width, half_width+1)
            distance = offset-taps
            weights = np.sinc(distance)*np.sinc(distance/half_width)
            weights /= weights.sum()  # 常量响应严格为1。
        kernels.append((taps, weights))
        low, high = max(low, -int(taps.min())), min(high, values.shape[1]-1-int(taps.max()))
    if high-low < 1:
        raise ValueError('时延与插值支持区间超过记录长度，没有共同有效数据')
    indices = np.arange(low, high+1)
    aligned = np.empty((len(values), len(indices)))
    for channel, (taps, weights) in enumerate(kernels):
        if len(taps) == 1:
            aligned[channel] = values[channel, indices+taps[0]]
        else:
            # 有限线性相关直接实现同一FIR求和，避免Python逐抽头循环；仍无周期卷绕。
            filtered = correlate(values[channel], weights, mode='valid', method='auto')
            aligned[channel] = filtered[indices+taps[0]]
    return dict(adc_mv=aligned, sample_time_ns=(indices+.5)*sample_width_ns,
                sample_width_ns=float(sample_width_ns), first_input_index=int(low),
                effective_duration_ns=float(len(indices)*sample_width_ns),
                discarded_samples=int(values.shape[1]-len(indices)), half_width=half_width)


def correlate_blocks(adc_mv, sample_width_ns, max_lag_ns=200., block_samples=None):
    """所有基线使用相同有效块；返回逐块相关及等曝光平均，尾部不足一块丢弃并计数。"""
    values = np.asarray(adc_mv, float)
    if (values.ndim != 2 or values.shape[0] < 2 or values.shape[1] < 2
            or not np.all(np.isfinite(values)) or not np.isfinite(sample_width_ns)
            or sample_width_ns <= 0 or not np.isfinite(max_lag_ns) or max_lag_ns < 0):
        raise ValueError('相关输入必须是有限逐镜波形和有效采样/时延范围')
    samples = values.shape[1] if block_samples is None else block_samples
    if (not isinstance(samples, (int, np.integer)) or samples < 2 or samples > values.shape[1]
            or max_lag_ns >= (samples-1)*sample_width_ns):
        raise ValueError('相关块过短或没有完整块')
    pairs = np.array(list(combinations(range(len(values)), 2)))
    count = values.shape[1]//samples
    correlations = []
    for block in range(count):
        record = values[:, block*samples:(block+1)*samples]
        if np.any(np.std(record, axis=1) == 0):
            raise ValueError('相关块含恒定通道，无法定义归一相关')
        curves = []
        for left, right in pairs:
            lags, curve = waveform_cross_correlation(record[left], record[right], sample_width_ns, max_lag_ns)
            curves.append(curve)
        correlations.append(curves)
    correlations = np.asarray(correlations)
    return dict(baselines=pairs, lags_ns=lags, block_correlations=correlations,
                mean_correlation=correlations.mean(axis=0), block_count=count,
                block_duration_s=samples*sample_width_ns*1e-9,
                effective_duration_s=count*samples*sample_width_ns*1e-9,
                discarded_tail_samples=values.shape[1]-count*samples)


def joint_thermal_mode_counts(rng, records, star_mean, background_mean, coherence,
                              modes=8, polarizations=2):
    """独立复高斯多镜场→光强→泊松计数，包含闭合相位；用于小规模高阶基准。"""
    gamma = _coherence_matrix(coherence)
    if (not isinstance(records, int) or records < 1 or not isinstance(modes, int)
            or modes < 1 or polarizations not in (1, 2)):
        raise ValueError('热光记录数和模数必须为正整数，偏振数为1或2')
    star, sky = _rates(star_mean, len(gamma)), _rates(background_mean, len(gamma))
    eigenvalues, vectors = np.linalg.eigh(gamma)
    root = vectors*np.sqrt(np.maximum(eigenvalues, 0.))
    intensities = np.empty((records, len(gamma)))
    counts = np.empty((records, len(gamma)), dtype=np.int64)
    for first in range(0, records, 4096):
        last = min(first+4096, records)
        shape = (last-first, modes*polarizations, len(gamma))
        independent = (rng.normal(size=shape)+1j*rng.normal(size=shape))/np.sqrt(2.)
        field = independent @ root.T
        intensities[first:last] = np.mean(abs(field)**2, axis=1)
        counts[first:last] = rng.poisson(star*intensities[first:last]+sky)
    return counts, intensities
