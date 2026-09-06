"""复现三镜共享波形、几何补偿、分块合并和独立热光高阶检查。"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'python'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import chi2, f

from sii_unified import (Instrument, detected_star_rate_hz, waveform_gls_weights,
                         waveform_instrument_signature, SIDEREAL_DAY_S)
from sii_validation import analytic_waveform_calibration, verify_main_parameters
from sii_layout import read_corsika_layout
from sii_observation import (geometric_arrival_delays_ns, simulate_array_waveforms,
                             align_waveforms, correlate_blocks, joint_thermal_mode_counts)


def sha256(path):
    return hashlib.sha256((ROOT/path).read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def mean_sem(values):
    return np.mean(values, axis=0), np.std(values, axis=0, ddof=1)/np.sqrt(len(values))


def common_interval_projections(raw, delays, dt, widths, weights):
    """所有插值宽度使用同一批样本；参考宽度是数值对照，不当作未知模拟真值。"""
    variants = [align_waveforms(raw, delays, dt, half_width=width) for width in widths]
    first = max(item['first_input_index'] for item in variants)
    last = min(item['first_input_index']+item['adc_mv'].shape[1] for item in variants)
    if last-first < 2:
        raise ValueError('核宽对照没有共同有效区间')
    projections = []
    for item in variants:
        offset = first-item['first_input_index']
        curves = correlate_blocks(item['adc_mv'][:, offset:offset+last-first], dt)
        projections.append(curves['mean_correlation'] @ weights)
    return np.asarray(projections), (last-first)*dt*1e-9


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--records', type=int, default=384)
    parser.add_argument('--calibration-records', type=int, default=256)
    parser.add_argument('--thermal-records', type=int, default=300_000)
    parser.add_argument('--seed', type=int, default=20260906)
    parser.add_argument('--output', type=Path, default=ROOT/'validation/sii_observation')
    args = parser.parse_args()
    if min(args.records, args.calibration_records) < 32 or args.thermal_records < 1000:
        parser.error('至少32条波形记录及1000条热光记录；正式结果建议采用默认值')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = verify_main_parameters(ROOT)
    instrument = Instrument.from_repository(ROOT)
    layout_path = 'configs/arrays/lact36_20260906.input'
    layout = read_corsika_layout(ROOT/layout_path).iloc[[0, 1, 5]]
    positions = layout[['east_m', 'north_m', 'up_m']].to_numpy()
    # 几何来自LACT布局；以下赤经对应的时角、赤纬及复相干矩阵是明确的检验场景。
    hour, dec, lat = .5, .3, np.deg2rad(29.3586)
    delays = geometric_arrival_delays_ns(positions, hour, dec, lat)
    duration = 24_000.
    gamma = np.array([[1., .7, .5], [.7, 1., .6], [.5, .6, 1.]])
    pairs = [(0, 1), (0, 2), (1, 2)]
    truth = np.array([abs(gamma[i, j])**2 for i, j in pairs])
    scale = 10_000.
    star = detected_star_rate_hz(2., instrument)
    background = instrument.detected_nsb_rate_hz+instrument.dark_count_rate_hz
    dt = instrument.sample_width_ns
    analytic, _ = analytic_waveform_calibration(instrument, block_duration_ns=duration)
    weights, _ = waveform_gls_weights(analytic, duration*1e-9)
    seed_sequences = np.random.SeedSequence(args.seed).spawn(6)
    series, curves = {}, {}
    widths = [16, 32, 64, 128, 256, 512]
    kernel_series = {}
    first_record = None
    for case_index, (name, count, coherence, injection) in enumerate([
            ('calibration', args.calibration_records, gamma, scale),
            ('injection_holdout', args.records, gamma, scale),
            ('physical_scale', args.records, gamma, 1.),
            ('null', args.records, np.eye(3), 1.)]):
        rng = np.random.default_rng(seed_sequences[case_index])
        projections, all_curves, width_records = [], [], []
        for record in range(count):
            raw = simulate_array_waveforms(rng, duration, star, background, coherence,
                                           instrument, delays, injection)
            aligned = align_waveforms(raw['adc_mv'], delays, dt)
            corrected = correlate_blocks(aligned['adc_mv'], dt)
            # 未补偿对照采用相同ADC索引区间，避免曝光长度差混入对照。
            first = aligned['first_input_index']
            length = aligned['adc_mv'].shape[1]
            uncorrected = correlate_blocks(raw['adc_mv'][:, first:first+length], dt)
            wide = align_waveforms(raw['adc_mv'], delays, dt, half_width=64)
            wide_curve = correlate_blocks(wide['adc_mv'], dt)
            blocks = correlate_blocks(aligned['adc_mv'], dt, block_samples=length//4)
            projections.append(np.stack([item['mean_correlation'] @ weights for item in
                                         (corrected, uncorrected, wide_curve, blocks)]))
            all_curves.append(np.stack((corrected['mean_correlation'], uncorrected['mean_correlation'])))
            if name in ('injection_holdout', 'null'):
                width_values, common_duration = common_interval_projections(
                    raw['adc_mv'], delays, dt, widths, weights)
                width_records.append(width_values)
            if first_record is None:
                first_record = dict(delays=delays, effective_duration_s=corrected['effective_duration_s'],
                                    wide_effective_duration_s=wide_curve['effective_duration_s'],
                                    block_effective_duration_s=blocks['effective_duration_s'],
                                    raw_duration_s=duration*1e-9,
                                    telescope_ids=layout.telescope_id.tolist(),
                                    injected_pair_fraction=(raw['metadata']['injected_pair_rates_hz'].sum(axis=1)/star))
            if (record+1) % 128 == 0:
                print(f'{name}: {record+1}/{count}', flush=True)
        series[name] = np.array(projections)
        curves[name] = np.array(all_curves)
        if width_records:
            kernel_series[name] = np.asarray(width_records)
    # 解析权重只定义固定投影；从独立训练样本标定处理后的响应，留出记录不参与训练。
    gain, gain_sem = mean_sem(series['calibration']/(scale*truth))
    rows, checks = [], []
    method_names = ['aligned_16', 'uncompensated', 'aligned_64', 'four_blocks']
    for name, values in series.items():
        for index, record in enumerate(values):
            for baseline, (left, right) in enumerate(pairs):
                row = dict(case=name, record=index, baseline=f'{left}-{right}')
                for method, label in enumerate(method_names):
                    row[label+'_reference_projection'] = record[method, baseline]
                rows.append(row)
    pd.DataFrame(rows).to_csv(output/'records.csv', index=False)
    # 未补偿曲线可能不包含峰，不用其接近零的训练响应去除噪声。
    calibrated = {name: value[:, [0, 2, 3], :]/gain[[0, 2, 3], :]
                  for name, value in series.items()}
    response_rows = []
    for method in [0, 2, 3]:
        for baseline in range(3):
            response_rows.append(dict(method=method_names[method], baseline=f'{pairs[baseline][0]}-{pairs[baseline][1]}',
                                      gain=gain[method, baseline], sem=gain_sem[method, baseline]))
    pd.DataFrame(response_rows).to_csv(output/'response.csv', index=False)
    holdout_mean, holdout_sem = mean_sem(calibrated['injection_holdout']/scale)
    total_sem = np.sqrt(holdout_sem**2+(truth*gain_sem[[0, 2, 3]]/gain[[0, 2, 3]])**2)
    response_z = (holdout_mean-truth)/total_sem
    checks.append(dict(name='independent_injection_response', max_abs_z=float(abs(response_z).max()),
                       passed=bool(np.all(abs(response_z) < 5))))
    null_mean, null_sem = mean_sem(calibrated['null'][:, 0])
    checks.append(dict(name='independent_null_mean', max_abs_z=float(abs(null_mean/null_sem).max()),
                       passed=bool(np.all(abs(null_mean/null_sem) < 5))))
    covariance_rows = []
    for name in ['null', 'physical_scale', 'injection_holdout']:
        values = calibrated[name][:, 0]
        covariance, correlation = np.cov(values.T), np.corrcoef(values.T)
        for i in range(3):
            for j in range(3):
                # Fisher区间是近似正态块估计量的有限重复诊断，不是全仪器系统误差区间。
                interval = ([1., 1.] if i == j else
                            np.tanh(np.arctanh(np.clip(correlation[i, j], -.999999, .999999))
                                    +np.array([-1., 1.])*1.96/np.sqrt(args.records-3)))
                covariance_rows.append(dict(case=name, row=i, column=j,
                                            covariance=covariance[i, j], correlation=correlation[i, j],
                                            correlation_low_95=interval[0], correlation_high_95=interval[1]))
    pd.DataFrame(covariance_rows).to_csv(output/'baseline_covariance.csv', index=False)
    # 同一条较长记录的整段估计与四块合并直接配对比较，保留二者相关性。
    difference = calibrated['injection_holdout'][:, 2]-calibrated['injection_holdout'][:, 0]
    diff_mean, diff_sem = mean_sem(difference/scale)
    train_difference = (series['calibration'][:, 3]/gain[3]-series['calibration'][:, 0]/gain[0])/scale
    _, train_diff_sem = mean_sem(train_difference)
    block_z = diff_mean/np.sqrt(diff_sem**2+train_diff_sem**2)
    checks.append(dict(name='whole_record_vs_four_blocks', max_abs_z=float(abs(block_z).max()),
                       passed=bool(np.all(abs(block_z) < 5))))
    width_difference = (series['calibration'][:, 2]-series['calibration'][:, 0])/(scale*truth)
    width_mean, width_sem = mean_sem(width_difference)
    null_std = np.std(calibrated['null'][:, 0], axis=0, ddof=1)
    dof = args.records-1
    null_std_interval = np.stack([null_std*np.sqrt(dof/chi2.ppf(.975, dof)),
                                 null_std*np.sqrt(dof/chi2.ppf(.025, dof))], axis=1)
    curve_mean, curve_sem = mean_sem(curves['injection_holdout'])
    pd.DataFrame([dict(baseline=f'{left}-{right}', lag_ns=lag,
                       aligned_mean=curve_mean[0, baseline, k], aligned_sem=curve_sem[0, baseline, k],
                       uncorrected_mean=curve_mean[1, baseline, k], uncorrected_sem=curve_sem[1, baseline, k])
                  for baseline, (left, right) in enumerate(pairs)
                  for k, lag in enumerate(analytic.lags_ns)]).to_csv(output/'correlation_curves.csv', index=False)

    # 复相位特意非零，以检验三镜闭合项；这里的模数是数学基准，不是LACT的纳秒模数。
    thermal_gamma = gamma.astype(complex)
    thermal_gamma[0, 1] *= np.exp(.4j)
    thermal_gamma[1, 0] = thermal_gamma[0, 1].conjugate()
    counts, intensity = joint_thermal_mode_counts(np.random.default_rng(seed_sequences[4]),
                                                 args.thermal_records, [4., 5., 6.], .3,
                                                 thermal_gamma, modes=2, polarizations=1)
    centered = intensity-1.
    thermal_rows = []
    for i, j in pairs:
        values = centered[:, i]*centered[:, j]
        mean, sem = mean_sem(values)
        thermal_rows.append(dict(moment=f'intensity_covariance_{i}{j}', expected=abs(thermal_gamma[i, j])**2/2,
                                 measured=mean, sem=sem))
    mean, sem = mean_sem(centered.prod(axis=1))
    thermal_rows.append(dict(moment='intensity_third_central',
                             expected=2*np.real(thermal_gamma[0, 1]*thermal_gamma[1, 2]*thermal_gamma[2, 0])/4,
                             measured=mean, sem=sem))
    for i, mu in enumerate([4., 5., 6.]):
        mean, sem = mean_sem((counts[:, i]-mu-.3)**2)
        thermal_rows.append(dict(moment=f'count_variance_{i}', expected=mu+.3+mu**2/2, measured=mean, sem=sem))
    thermal = pd.DataFrame(thermal_rows)
    thermal['z'] = (thermal.measured-thermal.expected)/thermal['sem']
    thermal.to_csv(output/'thermal_moments.csv', index=False)
    checks.append(dict(name='independent_thermal_moments', max_abs_z=float(abs(thermal.z).max()),
                       passed=bool(np.all(abs(thermal.z) < 5))))

    # 保留同记录配对差的不确定度；不同宽度分别定标会掩盖这里要测量的数值响应差。
    kernel_rows, kernel_records = [], []
    signal = kernel_series['injection_holdout']/(scale*truth)
    null = kernel_series['null']
    for index, width in enumerate(widths):
        difference = signal[:, index]-signal[:, -1]
        mean, sem = mean_sem(difference)
        response, response_sem = mean_sem(signal[:, index])
        noise_ratio = np.std(null[:, index], axis=0, ddof=1)/np.std(null[:, -1], axis=0, ddof=1)
        for baseline, (left, right) in enumerate(pairs):
            kernel_rows.append(dict(half_width=width, baseline=f'{left}-{right}',
                                    common_duration_s=common_duration,
                                    reference_half_width=widths[-1], response=response[baseline],
                                    response_sem=response_sem[baseline],
                                    response_minus_reference=mean[baseline], paired_sem=sem[baseline],
                                    noise_std_ratio_to_reference=noise_ratio[baseline]))
    for name, values in kernel_series.items():
        for record in range(args.records):
            for index, width in enumerate(widths):
                for baseline, (left, right) in enumerate(pairs):
                    kernel_records.append(dict(case=name, record=record, half_width=width,
                                               baseline=f'{left}-{right}', projection=values[record, index, baseline]))
    kernel_table = pd.DataFrame(kernel_rows)
    kernel_table.to_csv(output/'kernel_convergence.csv', index=False)
    pd.DataFrame(kernel_records).to_csv(output/'kernel_records.csv', index=False)

    # 用新随机流直接生成4倍长记录，不以短记录拼接或人为缩放噪声代替观测。
    rng = np.random.default_rng(seed_sequences[5])
    long_duration = 4*duration
    long_values = []
    for record in range(args.records):
        raw = simulate_array_waveforms(rng, long_duration, star, background, gamma, instrument, delays)
        aligned = align_waveforms(raw['adc_mv'], delays, dt)
        correlated = correlate_blocks(aligned['adc_mv'], dt)
        long_values.append(correlated['mean_correlation'] @ weights)
        if (record+1) % 128 == 0:
            print(f'long_physical_scale: {record+1}/{args.records}', flush=True)
    long_values = np.asarray(long_values)
    long_exposure = correlated['effective_duration_s']
    short_exposure = first_record['effective_duration_s']
    # 增益在比值中抵消；F区间要求块估计量近似正态，并不是极尾概率保证。
    ratio = (long_values.std(axis=0, ddof=1)/series['physical_scale'][:, 0].std(axis=0, ddof=1)
             *np.sqrt(long_exposure/short_exposure))
    ratio_low = ratio/np.sqrt(f.ppf(.975, args.records-1, args.records-1))
    ratio_high = ratio/np.sqrt(f.ppf(.025, args.records-1, args.records-1))
    scaling_rows = [dict(baseline=f'{left}-{right}', short_exposure_s=short_exposure,
                         long_exposure_s=long_exposure, scaled_sigma_ratio=ratio[index],
                         low_95=ratio_low[index], high_95=ratio_high[index])
                    for index, (left, right) in enumerate(pairs)]
    pd.DataFrame(scaling_rows).to_csv(output/'exposure_scaling.csv', index=False)
    pd.DataFrame([dict(record=record, baseline=f'{left}-{right}', projection=long_values[record, index])
                  for record in range(args.records) for index, (left, right) in enumerate(pairs)]).to_csv(
                      output/'long_records.csv', index=False)
    checks.append(dict(name='independent_long_record_noise_scaling',
                       min_scaled_sigma_ratio=float(ratio.min()), max_scaled_sigma_ratio=float(ratio.max()),
                       passed=bool(np.all((ratio/np.sqrt(f.ppf(.999, args.records-1, args.records-1)) < 1)
                                          & (ratio/np.sqrt(f.ppf(.001, args.records-1, args.records-1)) > 1)))))
    # 冻结几何的块内误差是实际转动速度下的上界检查，不把整夜当成固定时延。
    later_delays = geometric_arrival_delays_ns(positions, hour+2*np.pi*duration*1e-9/SIDEREAL_DAY_S, dec, lat)
    code_paths = ['python/sii_unified.py', 'python/sii_validation.py', 'python/sii_observation.py', 'python/sii_layout.py',
                  'tools/validate_sii_observation.py', 'tests/test_sii_observation.py']
    input_paths = [item['path'] for item in manifest['files']]+[layout_path, 'configs/arrays/lact36_20260906_coordinates.csv']
    input_paths += ['configs/optics/lact2_measured_single_pixel_400nm.csv',
                    'configs/optics/lact2_measured_single_pixel_400nm.provenance.json']
    summary = dict(seed=args.seed, records=args.records, calibration_records=args.calibration_records,
                   thermal_records=args.thermal_records, python=platform.python_version(),
                   numpy=np.__version__, scipy=scipy.__version__, main_commit=manifest['main_commit'],
                   instrument_signature=waveform_instrument_signature(instrument),
                   instrument_scenario=dict(magnitude_ab=2., star_rate_hz=star, background_rate_hz=background,
                                            sample_width_ns=dt, coherence_area_s=instrument.coherence_area_s,
                                            hour_angle_rad=hour, dec_rad=dec, lat_rad=lat,
                                            coherence=gamma, injection_scale=scale),
                   observation=first_record, frozen_delay_change_ns=later_delays-delays,
                   kernel_convergence=dict(widths=widths, common_duration_s=common_duration,
                                           reference_is_exact=False, paired_comparisons=kernel_rows),
                   exposure_scaling=scaling_rows, long_records=args.records,
                   holdout=dict(truth=truth, means=holdout_mean, total_sem=total_sem,
                                null_sigma=null_std, null_sigma_95_interval=null_std_interval,
                                block_minus_whole_mean=diff_mean, block_minus_whole_sem=diff_sem,
                                block_minus_whole_total_sem=np.sqrt(diff_sem**2+train_diff_sem**2),
                                wide_to_default_training_response=gain[2]/gain[0],
                                wide_minus_default_training_response=width_mean,
                                wide_minus_default_training_response_sem=width_sem),
                   checks=checks, code_sha256_lf={path: sha256(path) for path in code_paths},
                   input_sha256_lf={path: sha256(path) for path in sorted(set(input_paths))},
                   scope=dict(waveform_model='weak pair clusters, one stream per telescope',
                              thermal_model='independent discrete complex Gaussian modes with conditional Poisson counts',
                              physical_scale_covariance='finite Monte Carlo resolution; not a proof of independent baselines',
                              injection_covariance='artificial stress response; never extrapolate to physical exposure',
                              timing='constant within each 24 us record; geometry must be updated between observing blocks',
                              calibration='fixed analytic projection, independently trained processed response per baseline; gain SEM retained',
                              sigma_interval='conditional on fitted response; gain uncertainty reported separately in response.csv',
                              correlation_interval='Fisher approximate 95% interval for nearly normal block estimates',
                              interpolation='finite bandlimited approximation; no recovery of aliased analog frequencies',
                              long_exposure='direct 24/96 us comparison; hour-scale statistical branch still separate'))
    def serial(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(type(value).__name__)
    (output/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=serial,
                                                  allow_nan=False)+'\n', encoding='utf-8')
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for baseline, (left, right) in enumerate(pairs):
        axes[0, 0].plot(analytic.lags_ns, curve_mean[0, baseline], label=f'{left}-{right}, aligned')
    axes[0, 0].plot(analytic.lags_ns, curve_mean[1, 0], '--', color='0.5', label='0-1, uncorrected')
    axes[0, 0].set(xlabel='Lag (ns)', ylabel='Normalized correlation', title='Injected signal: shared raw waveforms')
    axes[0, 0].legend(fontsize=8)
    for index, label in enumerate(['16 taps/side', '64 taps/side', '4 blocks']):
        axes[0, 1].errorbar(np.arange(3)+(index-1)*.15, holdout_mean[index]/truth,
                           yerr=1.96*total_sem[index]/truth, fmt='o', capsize=3, label=label)
    axes[0, 1].axhline(1., color='0.5', ls='--')
    axes[0, 1].set(xticks=range(3), xticklabels=['0-1', '0-2', '1-2'], xlabel='Baseline',
                   ylabel='Recovered / injected power', title='Independent holdout; 95% normal error bars')
    axes[0, 1].legend(fontsize=8)
    for axis, name in zip(axes[1], ['physical_scale', 'injection_holdout']):
        matrix = np.corrcoef(calibrated[name][:, 0].T)
        plot = axis.imshow(matrix, vmin=-1, vmax=1, cmap='coolwarm')
        for i in range(3):
            for j in range(3):
                axis.text(j, i, f'{matrix[i, j]:.2f}', ha='center', va='center')
        axis.set(xticks=range(3), yticks=range(3), xticklabels=['0-1', '0-2', '1-2'],
                 yticklabels=['0-1', '0-2', '1-2'],
                 title=('Physical scale = 1' if name == 'physical_scale' else 'Artificial injection scale = 10000'))
    fig.colorbar(plot, ax=axes[1].tolist(), shrink=.8)
    fig.savefig(output/'observation_validation.png', dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for label, group in kernel_table.groupby('baseline', sort=False):
        axes[0].errorbar(group.half_width, group.response_minus_reference,
                         yerr=1.96*group.paired_sem, fmt='o-', capsize=3, label=label)
    axes[0].axhline(0., color='0.5', ls='--')
    axes[0].set(xscale='log', xlabel='Interpolation half width (samples)',
                ylabel='Response minus width 512', title='Same photons and common time interval')
    axes[0].legend(title='Baseline')
    axes[1].errorbar(np.arange(3), ratio, yerr=np.array([ratio-ratio_low, ratio_high-ratio]), fmt='o', capsize=4)
    axes[1].axhline(1., color='0.5', ls='--')
    axes[1].set(xticks=range(3), xticklabels=['0-1', '0-2', '1-2'], xlabel='Baseline',
                ylabel='SD(long) / SD(short) * sqrt(Tlong / Tshort)', title='Independent 24 / 96 us raw records; 95% CI')
    fig.savefig(output/'processing_convergence.png', dpi=160)
    plt.close(fig)
    print(json.dumps(dict(checks=checks, output=str(output)), ensure_ascii=False, indent=2))
    if not all(check['passed'] for check in checks):
        raise SystemExit('独立检查未通过；结果已保留，不能宣称验证通过')


if __name__ == '__main__':
    main()
