"""36镜角直径性能：采样相位模板、弱光误差界及均匀圆瞳对照。"""
from dataclasses import replace
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.linalg import cho_factor, cho_solve

import sii_unified as sii
from sii_validation import analytic_waveform_calibration


def datasheet_dark_scenario(instrument, root, summed_channels=8, rating='typical'):
    """仅启用25℃、过压8.5 V下的规格书暗计数假设；不推定真实接线或温度。"""
    data = json.loads((Path(root)/'configs/sii/s17351_datasheet.json').read_text(encoding='utf-8'))
    if (not isinstance(summed_channels, int) or not 1 <= summed_channels <= data['channels']
            or rating not in ('typical', 'max')):
        raise ValueError('需要1至8个汇总通道及typical/max额定类别')
    rate = summed_channels*data['dark_count_hz_per_channel'][rating]
    return replace(instrument, dark_count_rate_hz=float(rate))


def phase_template_bank(instrument, magnitude=2., block_duration_ns=24000., nodes=65,
                        fine_dt_ns=.05):
    """整数移位后，相关峰模板保留[-dt,dt]残差；所有相位共用零信号协方差。"""
    if not isinstance(nodes, int) or nodes < 5 or nodes % 2 != 1:
        raise ValueError('相位节点需要至少5个奇数节点')
    phase = np.linspace(-instrument.sample_width_ns, instrument.sample_width_ns, nodes)
    templates = []
    for delay in phase:
        calibration, _ = analytic_waveform_calibration(instrument, magnitude,
            block_duration_ns=block_duration_ns, residual_delay_ns=delay, fine_dt_ns=fine_dt_ns)
        templates.append(calibration.peak_per_visibility2)
    templates = np.asarray(templates)
    solved = cho_solve(cho_factor(calibration.covariance_per_block), templates.T).T
    information = np.einsum('ij,ij->i', templates, solved)
    return dict(phase_ns=phase, templates=templates, weights=solved/information[:, None],
                sigma_block=1/np.sqrt(information), block_duration_s=block_duration_ns*1e-9,
                lags_ns=calibration.lags_ns, covariance=calibration.covariance_per_block)


def interpolate_bank(bank, residual_ns, name):
    """逐列线性插值；禁止超出已计算的相位范围。"""
    residual = np.asarray(residual_ns, float)
    phase = bank['phase_ns']
    if np.any(~np.isfinite(residual)) or np.any(abs(residual) > phase[-1]+1e-10):
        raise ValueError('残余时延超出模板范围')
    values = bank[name]
    if values.ndim == 1:
        return np.interp(residual, phase, values)
    return np.stack([np.interp(residual, phase, column) for column in values.T], axis=-1)


def integer_align(adc_mv, delays_ns, dt_ns):
    """每镜仅平移整数样本，裁剪公共支持；不做循环移位或分数插值。"""
    values, delays = np.asarray(adc_mv, float), np.asarray(delays_ns, float)
    if (values.ndim != 2 or delays.shape != (len(values),) or not np.isfinite(dt_ns)
            or dt_ns <= 0 or not np.all(np.isfinite(values)) or not np.all(np.isfinite(delays))):
        raise ValueError('无效ADC、时差或采样间隔')
    shifts = np.rint(delays/dt_ns).astype(int)
    first, last = max(0, -shifts.min()), min(values.shape[1], values.shape[1]-shifts.max())
    if last-first < 2:
        raise ValueError('没有共同有效样本')
    aligned = np.stack([row[first+shift:last+shift] for row, shift in zip(values, shifts)])
    residual = delays-shifts*dt_ns
    pairs = np.array(list(combinations(range(len(values)), 2)))
    return aligned, residual[pairs[:, 1]]-residual[pairs[:, 0]], (last-first)*dt_ns*1e-9


def tracked_segment_precision(uvw, observation, instrument, bank, time_nodes=1200):
    """按段内实际几何相位积累逆方差；均匀P均值对应等时长块平均。

    使用平均块方差而非平均信息量，因此和既有等时间可见度模型严格对应。
    逐块模板跟随时角；1200个中点仅用于积分采样相位，不是ADC时间步长。
    """
    if not isinstance(time_nodes, int) or time_nodes < 2:
        raise ValueError('时延相位积分节点必须至少为2')
    result = np.empty(len(uvw))
    lat, dec = np.deg2rad([observation.site_lat_deg, observation.source_dec_deg])
    offsets = ((np.arange(time_nodes)+.5)/time_nodes-.5)*observation.segment_s
    for _, frame in uvw.groupby('segment', sort=False):
        h = frame.hour_angle_h.iloc[0]*np.pi/12+offsets*2*np.pi/sii.SIDEREAL_DAY_S
        direction = np.array([-np.cos(dec)*np.sin(h),
            np.sin(dec)*np.cos(lat)-np.cos(dec)*np.cos(h)*np.sin(lat),
            np.sin(dec)*np.sin(lat)+np.cos(dec)*np.cos(h)*np.cos(lat)])
        baselines = frame[['baseline_east_m','baseline_north_m','baseline_up_m']].to_numpy()
        delays = -baselines @ direction/sii.C_M_S*1e9
        # 相差一个样本仅平移无限滞后模板；200 ns窗口远大于剩余4 ns。
        residual = (delays+instrument.sample_width_ns/2)%instrument.sample_width_ns-instrument.sample_width_ns/2
        variance = interpolate_bank(bank, residual, 'sigma_block')**2
        result[frame.index] = np.sqrt(variance.mean(axis=1)*bank['block_duration_s']/
                                     (observation.segment_s*observation.nights))
    return result


def weak_light_covariance_budget(instrument, magnitude, telescopes=36, count_bin_ns=4.):
    """给出线性高斯Bartlett谱界，以及独立的未分辨多模热光计数精确四阶对照。

    前者含自聚束和共镜谱协方差：归一谱矩阵I+E满足||E||<=N*epsilon。
    后者是矩形时间箱模型，不冒充SPE滤波后四阶累积量的严格界。
    两者都按未分辨源|Gamma|=1计算；电子学非线性不在此近似内。
    """
    if telescopes < 2 or count_bin_ns <= 0:
        raise ValueError('无效镜数或计数窗口')
    star = sii.detected_star_rate_hz(magnitude, instrument)
    total = star+instrument.detected_nsb_rate_hz+instrument.dark_count_rate_hz
    tau = instrument.coherence_area_s
    epsilon = star**2*tau/(total*instrument.excess_noise_factor**2)
    spectral_bound = (1+telescopes*epsilon)**2-1
    duration = count_bin_ns*1e-9
    mu, mean = star*duration, total*duration
    modes = duration/tau
    a, b, c = mu**2/modes, mu**3/modes**2, mu**4/modes**3
    variance = mean+a
    diagonal = variance**2+a*a+a+4*b+6*c
    shared = variance*a+a*a+2*b+6*c
    disjoint = 2*a*a+6*c
    shared_count = 2*(telescopes-2)
    disjoint_count = (telescopes-2)*(telescopes-3)//2
    count_bound = (diagonal+shared_count*shared+disjoint_count*disjoint)/mean**2-1
    return dict(photon_degeneracy=star*tau, normalized_cross_spectrum_bound=epsilon,
        bartlett_covariance_operator_bound=spectral_bound,
        bartlett_sigma_inflation_bound=np.sqrt(1+spectral_bound)-1,
        thermal_count_covariance_operator_bound=count_bound, count_bin_ns=count_bin_ns,
        thermal_count_connected_diagonal_fraction=(a+4*b+6*c)/mean**2)


def pupil_difference_quadrature(radius_m=4., radial_nodes=12, angle_nodes=24):
    """两个独立均匀圆瞳入射点之差；圆瞳位于视线垂直平面。"""
    if radius_m < 0 or not np.isfinite(radius_m) or radial_nodes < 2 or angle_nodes < 4:
        raise ValueError('无效圆瞳积分参数')
    if radius_m == 0:
        return np.zeros((1, 2)), np.ones(1)
    x, weights = np.polynomial.legendre.leggauss(radial_nodes)
    separation = radius_m*(x+1)
    overlap = (2*radius_m**2*np.arccos(separation/(2*radius_m))
               -.5*separation*np.sqrt(4*radius_m**2-separation**2))
    radial_weights = weights*radius_m*2*separation*overlap/(np.pi*radius_m**4)
    if abs(radial_weights.sum()-1) > 2e-4:
        raise ValueError('圆瞳积分不收敛，请增加径向节点')
    angles = 2*np.pi*(np.arange(angle_nodes)+.5)/angle_nodes
    offsets = np.stack([separation[:,None]*np.cos(angles), separation[:,None]*np.sin(angles)], axis=-1)
    weight = np.repeat(radial_weights/angle_nodes, angle_nodes)
    return offsets.reshape(-1,2), weight/weight.sum()


def aperture_disk_power(uvw, observation, instrument, diameter_mas, radius_m=4.,
                        radial_nodes=12, angle_nodes=24):
    """均匀8 m瞳面情景的时间、波长、入瞳平均；不假装已知真实遮挡权重。"""
    u, v = sii.segment_uv_samples(uvw, observation, instrument)
    spectral = np.asarray(instrument.visibility_wavelength_nm or (instrument.wavelength_nm,))
    wavelengths = np.tile(spectral*1e-9, observation.visibility_subsamples_per_segment)
    points, weights = pupil_difference_quadrature(radius_m, radial_nodes, angle_nodes)
    result = np.zeros(len(uvw))
    sample_weights = sii.segment_sampling_weights(observation, instrument)
    for point, weight in zip(points, weights):
        radius = np.hypot(u+point[0]/wavelengths, v+point[1]/wavelengths)
        result += weight*(sii.uniform_disk_visibility(radius, diameter_mas)**2 @ sample_weights)
    return result


def disk_model_grid(uvw, observation, instrument, diameters):
    """拟合只用布局及候选直径，真值不进入网格选择。"""
    u, v = sii.segment_uv_samples(uvw, observation, instrument)
    radius = np.hypot(u, v)
    weights = sii.segment_sampling_weights(observation, instrument)
    return np.array([sii.uniform_disk_visibility(radius, d)**2 @ weights for d in diameters])


def compressed_profiles(models, sigma, truth, gain_sigma, trials, rng, anchors=81,
                        noise_relative_sigma=0., assumed_sigma_scale=1.):
    """用候选模型张成的子空间抽样高斯充分统计量，保留增益剖面。

    QR基仅由预定候选模型产生。逐模型重构误差随结果返回，超过阈值报错；
    正交补空间给所有候选相同的卡方常数，可从似然比中严格消去。
    """
    from scipy.linalg import qr
    if (gain_sigma < 0 or noise_relative_sigma < 0 or assumed_sigma_scale <= 0
            or not np.all(np.isfinite([gain_sigma,noise_relative_sigma,assumed_sigma_scale]))):
        raise ValueError('标定误差及拟合噪声倍率无效')
    fitted_sigma = np.asarray(sigma)*assumed_sigma_scale
    whitened = models/fitted_sigma
    q, _ = qr(whitened[np.linspace(0, len(models)-1, anchors).astype(int)].T, mode='economic')
    projected = whitened @ q
    remainder = whitened-projected @ q.T
    error = float(np.max(np.linalg.norm(remainder, axis=1)))
    if error > 1e-5:
        raise ValueError(f'候选模型子空间截断误差过大：{error}')
    target = np.asarray(truth)/fitted_sigma @ q
    gains = 1+rng.normal(0, gain_sigma, (trials,1))
    # 学习所得噪声尺度作为整次观测共享的不确定量；对数正态只描述标定误差，非天气。
    log_width = np.sqrt(np.log1p(noise_relative_sigma**2))
    noise_scale = np.exp(rng.normal(-.5*log_width**2,log_width,(trials,1)))
    observations = gains*target+noise_scale/assumed_sigma_scale*rng.normal(size=(trials, q.shape[1]))
    square = np.sum(projected**2, axis=1)
    cross = observations @ projected.T
    precision = 1/gain_sigma**2 if gain_sigma > 0 else 0.
    fit_gain = (cross+precision)/(square+precision) if gain_sigma > 0 else np.ones_like(cross)
    profile = -2*cross*fit_gain+square*fit_gain**2+precision*(fit_gain-1)**2
    fixed = -2*cross+square
    return profile, fixed, error


def refined_profile_intervals(grid, statistics, refinement=10):
    """细化平滑似然曲线，避免候选网格量化高信噪比直径；不依赖真值。"""
    from scipy.interpolate import CubicSpline
    from sii_validation import profile_grid_interval
    grid = np.asarray(grid)
    fine = np.linspace(grid[0],grid[-1],(len(grid)-1)*refinement+1)
    return [profile_grid_interval(fine,CubicSpline(grid,row-row.min())(fine)) for row in statistics]
