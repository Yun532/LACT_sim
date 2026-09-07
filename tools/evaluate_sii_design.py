"""复用连续波形GLS，评估仪器改动的条件收益；不改动实测配置。"""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import sys
import json
import hashlib
from dataclasses import replace
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'python'))
import sii_unified as sii
import sii_performance as perf
from sii_layout import read_corsika_layout
from sii_validation import analytic_waveform_calibration, verify_main_parameters


def diameter_information(uvw, observation, instrument, sigma, step=1e-5):
    """对直径与一个全阵列共同增益求局部Fisher矩阵，波长/曝光平均由主代码完成。"""
    minus, model, plus = perf.disk_model_grid(uvw, observation, instrument,
                                             [.16-step, .16, .16+step])
    jacobian = np.stack(((plus-minus)/(2*step), model))/sigma
    return jacobian @ jacobian.T


def diameter_sigma(information, gain_sigma):
    """共同增益先验只加入一次；零表示固定已知增益。返回mas。"""
    precision = information[0, 0]
    if gain_sigma > 0:
        precision -= information[0, 1]**2/(information[1, 1]+1/gain_sigma**2)
    return float(1/np.sqrt(precision))


def covariance_floor_probe(inst, bank, fine=.025):
    """独立重建未截断Bartlett矩阵，检查零电子噪声时数值正则化是否限制带宽。"""
    _, diag = analytic_waveform_calibration(inst, fine_dt_ns=fine, block_duration_ns=24000.)
    time, pulse = sii.load_measured_spe_template(inst.spe_template_path)
    h = np.interp(np.arange(time[0], time[-1]+fine/2, fine), time, pulse)
    auto = fftconvolve(h, h[::-1])*fine
    t_auto = np.arange(1-len(h), len(h))*fine
    dt = inst.sample_width_ns
    span = int(np.ceil(np.ptp(time)/dt))+2
    steps = np.arange(-span, span+1)
    total = sii.detected_star_rate_hz(2., inst)+inst.detected_nsb_rate_hz+inst.dark_count_rate_hz
    rho = total*1e-9*diag['charge_second_moment']*np.interp(steps*dt, t_auto, auto, left=0, right=0)/diag['variance_mv2']
    rho[span] += inst.electronic_noise_rms_mv**2/diag['variance_mv2']
    bartlett = fftconvolve(rho, rho[::-1])/(24000/dt)/diag['finite_block_normalization']**2
    lags = np.rint(bank['lags_ns']/dt).astype(int)
    cov = np.interp(lags[:, None]-lags[None, :], np.arange(1-len(rho),len(rho)), bartlett)
    values, vectors = np.linalg.eigh(cov)
    projection = bank['templates'] @ vectors
    results = []
    for floor in [1e-10, 1e-12, 1e-14]:
        sigma = 1/np.sqrt(np.sum(projection**2/np.maximum(values, values.max()*floor), axis=1))
        results.append((floor, int(np.sum(values < values.max()*floor)), dict(bank, sigma_block=sigma)))
    return results


def run():
    output = ROOT/'validation/sii_design'
    output.mkdir(parents=True, exist_ok=True)
    manifest = verify_main_parameters(ROOT)
    base = sii.Instrument.from_repository(ROOT)
    obs = sii.Observation()
    layout = read_corsika_layout(ROOT/'configs/arrays/lact36_20260906.input')
    uvw = sii.generate_uvw(layout, obs, base)
    gain_sigma = .0108040292195
    rows, banks, precisions, instruments = [], {}, {}, {}

    def evaluate(name, inst, magnitude=2., nodes=17, fine=.025):
        bank = perf.phase_template_bank(inst, magnitude, nodes=nodes, fine_dt_ns=fine)
        sigma = perf.tracked_segment_precision(uvw, obs, inst, bank, time_nodes=1200)
        information = diameter_information(uvw, obs, inst, sigma)
        cal, diagnostic = analytic_waveform_calibration(inst, magnitude, fine_dt_ns=fine)
        row = dict(case=name, magnitude=magnitude, sample_ns=inst.sample_width_ns,
            noise_mv=inst.electronic_noise_rms_mv, dark_hz=inst.dark_count_rate_hz,
            star_hz=sii.detected_star_rate_hz(magnitude, inst), nsb_hz=inst.detected_nsb_rate_hz,
            peak_C_per_P=float(cal.peak_per_visibility2.max()),
            variance_mv2=diagnostic['variance_mv2'], sigma_P_median=float(np.median(sigma)),
            diameter_sigma_mas=diameter_sigma(information, gain_sigma),
            diameter_sigma_fixed_gain_mas=diameter_sigma(information, 0))
        banks[name], precisions[name], instruments[name] = bank, sigma, inst
        print(name, row['sigma_P_median'], flush=True)
        return row, information

    for mag in [2., 6.]:
        for noise in [0., 1.]:
            group = f'm{int(mag)}_e{int(noise)}'
            cases = [('base', base)]
            for dt in [2., 1., .5]:
                cases.append((f'dt{dt:g}', replace(base, adc_sample_rate_hz=1e9/dt,
                                                   electronics_bandwidth_hz=.5e9/dt)))
            for scale in [1.2, 1.5, 2.]:
                cases.append((f'eff{scale:g}', replace(base, throughput=base.throughput*scale,
                    detected_nsb_rate_hz=base.detected_nsb_rate_hz*scale)))
            for dt in [4., 1.]:
                cases.append((f'iso_dt{dt:g}', replace(base, optical_timing_kernel_path=None,
                    adc_sample_rate_hz=1e9/dt, electronics_bandwidth_hz=.5e9/dt)))
            cases.append(('jitter1', replace(base, intrinsic_time_jitter_ns=1.)))
            # 理想硬件变更：压缩时间但保持脉冲面积；不是裁掉现有数据的尾部。
            t, h = sii.load_measured_spe_template(base.spe_template_path)
            for scale in [.5, .25]:
                path = output/f'spe_time_{scale:g}.csv'
                pd.DataFrame(dict(time_ns=t*scale, amplitude_mv=h/scale)).to_csv(path, index=False)
                cases.append((f'spe{scale:g}', replace(base, spe_template_path=str(path))))
            for suffix, inst in cases:
                # 同一模拟白噪声谱密度：采样带宽翻倍，样本RMS乘sqrt(2)。
                inst = replace(inst, electronic_noise_rms_mv=noise*np.sqrt(4/inst.sample_width_ns))
                row, _ = evaluate(f'{group}_{suffix}', inst, mag)
                row['reference'] = f'{group}_base'
                rows.append(row)

    frame = pd.DataFrame(rows).set_index('case', drop=False)
    frame['snr_gain'] = [frame.loc[r.reference, 'sigma_P_median']/r.sigma_P_median for r in frame.itertuples()]
    frame['diameter_precision_gain'] = [frame.loc[r.reference, 'diameter_sigma_mas']/r.diameter_sigma_mas for r in frame.itertuples()]
    frame.to_csv(output/'design.csv', index=False)

    # 固定399--401nm总通带分色；每通道独立探测器，积分节点不是物理通道。
    channels = []
    for count in [1, 2, 4, 8, 16]:
        for scenario, loss, noise, dark in [('ideal', 1., 0., 0.), ('loss_dark_noise', .8, 1., 9.6e6)]:
            total_info = np.zeros((2, 2))
            photon_rate = 0.
            for channel in range(count):
                low, high = 399+2*channel/count, 399+2*(channel+1)/count
                path = output/f'band_{count}_{channel}.csv'
                pd.DataFrame({'wavelength_nm': [low-1e-6, low, high, high+1e-6],
                              'transmission': [0., 1., 1., 0.]}).to_csv(path, index=False)
                inst = sii.Instrument.from_repository(ROOT, sii_bandpass_path=str(path))
                inst = replace(inst, throughput=inst.throughput*loss,
                    detected_nsb_rate_hz=inst.detected_nsb_rate_hz*loss,
                    electronic_noise_rms_mv=noise, dark_count_rate_hz=dark)
                row, information = evaluate(f'channel_{scenario}_{count}_{channel}', inst)
                total_info += information
                photon_rate += row['star_hz']
            channels.append(dict(channels=count, scenario=scenario, star_hz_total=photon_rate,
                diameter_sigma_mas=diameter_sigma(total_info, gain_sigma),
                diameter_sigma_fixed_gain_mas=diameter_sigma(total_info, 0)))
    channel_frame = pd.DataFrame(channels)
    channel_frame['gain_vs_ideal_single'] = frame.loc['m2_e0_base', 'diameter_sigma_mas']/channel_frame.diameter_sigma_mas
    channel_frame.to_csv(output/'channels.csv', index=False)

    # 恢复只给占据率和电荷矩，不能冒充非线性SII误差预测。
    rate = sii.detected_star_rate_hz(2., base)+base.detected_nsb_rate_hz
    recovery = []
    rng = np.random.default_rng(20260907)
    for cells in [base.microcells_per_pixel, 1000]:
        for tau in [10., 30., 100.]:
            for multiplier in [1., 100.]:
                x = rate*multiplier/cells*tau*1e-9
                f = -np.expm1(-rng.exponential(1/x, 500000))
                expected, second = 1/(1+x), 2/((1+x)*(2+x))
                fourth = 24/np.prod(np.arange(1., 5.)+x)
                assert abs(f.mean()-expected) < 6*np.sqrt((second-expected**2)/len(f))+1e-8
                assert abs(np.mean(f*f)-second) < 6*np.sqrt((fourth-second**2)/len(f))+1e-8
                recovery.append(dict(cells=cells, tau_ns=tau, rate_multiplier=multiplier, occupancy=x,
                    mean_charge=1/(1+x), charge_second_moment=2/((1+x)*(2+x)),
                    mc_mean=float(f.mean()), mc_second=float(np.mean(f*f))))
    pd.DataFrame(recovery).to_csv(output/'recovery.csv', index=False)

    checks = []
    for name in ['m2_e0_base', 'm2_e0_dt0.5', 'm2_e1_spe0.25']:
        inst = instruments[name]
        refined = perf.phase_template_bank(inst, nodes=33, fine_dt_ns=.0125)
        sigma = perf.tracked_segment_precision(uvw, obs, inst, refined, time_nodes=2400)
        relative = float(np.max(abs(sigma/precisions[name]-1)))
        checks.append(dict(case=name, refined_max_sigma_relative_change=relative))
        if relative > .02:
            raise ValueError(f'{name}: 相位/时间/细网格误差超过2%: {relative}')
    pd.DataFrame(checks).to_csv(output/'convergence.csv', index=False)
    nominal = diameter_information(uvw, obs, base, precisions['m2_e0_base'])
    refined = diameter_information(uvw, obs, base, precisions['m2_e0_base'], step=5e-6)
    np.testing.assert_allclose(nominal, refined, rtol=1e-6)
    assert abs(channel_frame.query("channels == 16 and scenario == 'ideal'").star_hz_total.iloc[0]/sii.detected_star_rate_hz(2., base)-1) < 1e-4
    for mag in [2, 6]:
        assert abs(frame.loc[f'm{mag}_e0_eff2', 'snr_gain']-2) < 1e-5

    floors = []
    for name in ['m2_e0_base', 'm2_e0_dt2', 'm2_e0_dt1', 'm2_e0_dt0.5', 'm2_e1_dt1']:
        for floor, modes, alternate in covariance_floor_probe(instruments[name], banks[name]):
            sigma = perf.tracked_segment_precision(uvw, obs, instruments[name], alternate)
            floors.append(dict(case=name, relative_eigenvalue_floor=floor, clipped_modes=modes,
                lag_modes=len(alternate['lags_ns']), sigma_P_median=float(np.median(sigma)),
                snr_vs_current_floor=float(np.median(precisions[name])/np.median(sigma))))
            if floor == 1e-10:
                np.testing.assert_allclose(sigma, precisions[name], rtol=2e-6)
    pd.DataFrame(floors).to_csv(output/'covariance_floor.csv', index=False)

    fig, panels = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
    t, h = sii.load_measured_spe_template(base.spe_template_path)
    for scale in [1., .5, .25]:
        panels[0].plot(t*scale, h/scale, label=f'time x{scale:g}, same area')
    panels[0].set(xlabel='Time (ns)', ylabel='Single PE voltage (mV)', xlim=(-10, 140))
    panels[0].legend(fontsize=8)
    ax = panels[1:]
    for noise in [0, 1]:
        for iso in [False, True]:
            suffixes = ['iso_dt4','iso_dt1'] if iso else ['base','dt2','dt1','dt0.5']
            part = frame.loc[[f'm2_e{noise}_{suffix}' for suffix in suffixes]]
            ax[0].plot(part.sample_ns, part.snr_gain, 'o--' if iso else 'o-',
                label=f'{noise} mV at 4 ns; '+('isochronous' if iso else 'current optics'))
    ax[0].set(xlabel='ADC interval (ns)', ylabel='SNR / same-noise 4 ns baseline', xscale='log')
    ax[0].invert_xaxis()
    ax[0].legend(fontsize=8)
    for scenario, part in channel_frame.groupby('scenario', sort=False):
        ax[1].plot(part.channels, part.gain_vs_ideal_single, 'o-', label=scenario)
    ax[1].plot([1,2,4,8,16], np.sqrt([1,2,4,8,16]), ':', label='sqrt(N), guide only')
    ax[1].set(xlabel='Independent spectral channels (fixed total 2 nm)', ylabel='Diameter precision / ideal single channel', xscale='log')
    ax[1].legend(fontsize=8)
    for a in panels:
        a.grid(alpha=.2)
    fig.savefig(output/'design.png', dpi=160)
    plt.close(fig)
    inputs = {entry['path'] for entry in manifest['files']}
    inputs.update(['tools/evaluate_sii_design.py', 'python/sii_unified.py', 'python/sii_performance.py',
        'python/sii_validation.py', 'python/sii_layout.py', 'configs/arrays/lact36_20260906.input'])
    digest = lambda p: hashlib.sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
    summary = dict(kind='conditional_linear_waveform_GLS_design_study', magnitude=2, diameter_mas=.16,
        telescopes=len(layout), rows=len(uvw), exposure_s=21600, gain_sigma=gain_sigma,
        phase_nodes=17, phase_time_nodes=1200, fine_dt_ns=.025,
        recovery_seed=20260907, recovery_trials=500000,
        input_sha256_lf={p: digest(ROOT/p) for p in sorted(inputs)},
        output_sha256_lf={p.name: digest(p) for p in sorted(output.glob('*.csv'))},
        limitations=['Fisher local precision, not image recovery or confidence coverage',
            'No new hardware calibration; isochrony holds collecting efficiency fixed',
            'Sampling scenarios assume current measured SPE interpolation and zero unknown jitter',
            'Spectral channels have independent noise; one shared gain prior',
            'Recovery charge moments assume uniform illumination and fixed trigger PDE'])
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    print('Design study complete', flush=True)


if __name__ == '__main__':
    run()
