"""复现36镜自转几何和七个时刻的三镜波形追踪；稀疏短记录不计作六小时曝光。"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
from dataclasses import asdict

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'python'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sii_layout import read_corsika_layout
from sii_unified import Instrument, BinarySource, detected_star_rate_hz, waveform_gls_weights
from sii_validation import analytic_waveform_calibration, verify_main_parameters
from sii_observation import (tracking_geometry, simulate_array_waveforms, align_waveforms,
                             correlate_blocks, source_coherence_spectrum)


def digest(path):
    return hashlib.sha256((ROOT/path).read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--records', type=int, default=128)
    parser.add_argument('--calibration-records', type=int, default=96)
    parser.add_argument('--seed', type=int, default=20260907)
    parser.add_argument('--source-case', choices=['fixed','single_disk'], default='fixed')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if min(args.records, args.calibration_records) < 32:
        parser.error('训练和留出至少各32条记录')
    output = args.output or ROOT/'validation'/('sii_tracking' if args.source_case=='fixed' else 'sii_source_tracking')
    output.mkdir(parents=True, exist_ok=True)
    manifest = verify_main_parameters(ROOT)
    instrument = Instrument.from_repository(ROOT)
    input_path = 'configs/arrays/lact36_20260906.input'
    layout = read_corsika_layout(ROOT/input_path)
    positions = layout[['east_m', 'north_m', 'up_m']].to_numpy()
    array_positions = positions.copy()
    dec, lat = .3, np.deg2rad(29.3586)
    epochs = np.linspace(-10800., 10800., 7)
    initial = tracking_geometry(positions, 0., dec, lat, epochs[0])
    # 一分钟几何表包含全部镜；波形只在七个时刻对Tel.1/2/6独立采样。
    geometry_rows = []
    for elapsed in np.linspace(-10800., 10800., 361):
        state = tracking_geometry(positions, 0., dec, lat, elapsed)
        for channel in range(36):
            geometry_rows.append(dict(elapsed_s=elapsed, telescope_id=channel+1,
                delay_ns=state['arrival_delays_ns'][channel],
                delay_rate_ns_per_s=state['arrival_delay_rates_ns_per_s'][channel],
                initial_delay_residual_ns=state['arrival_delays_ns'][channel]-initial['arrival_delays_ns'][channel],
                curvature_bound_ns_per_s2=state['curvature_bound_ns_per_s2'][channel]))
    geometry = pd.DataFrame(geometry_rows)
    geometry.to_csv(output/'geometry.csv', index=False)
    selected = [0, 1, 5]
    positions = positions[selected]
    initial_delays = initial['arrival_delays_ns'][selected]
    duration, scale = 24000., 10000.
    dt = instrument.sample_width_ns
    star = detected_star_rate_hz(2., instrument)
    background = instrument.detected_nsb_rate_hz+instrument.dark_count_rate_hz
    gamma = np.array([[1., .7, .5], [.7, 1., .6], [.5, .6, 1.]])
    truth = np.array([.49, .25, .36])
    source = BinarySource(primary_diameter_mas=.16)
    calibration_truth = truth if args.source_case=='fixed' else np.ones(3)
    spectral_weights = None
    analytic, _ = analytic_waveform_calibration(instrument, block_duration_ns=duration)
    weights, _ = waveform_gls_weights(analytic, duration*1e-9)
    rows, result_rows, checks, epoch_rows, source_rows = [], [], [], [], []
    seeds = np.random.SeedSequence(args.seed).spawn(len(epochs)*3)
    for epoch_index, elapsed in enumerate(epochs):
        state = tracking_geometry(positions, 0., dec, lat, elapsed)
        delays, rates = state['arrival_delays_ns'], state['arrival_delay_rates_ns_per_s']
        source_change = 0.
        if args.source_case != 'fixed':
            spectrum = source_coherence_spectrum(array_positions,state['hour_angle_rad'],dec,lat,
                                                 source,instrument,args.source_case)
            gamma = spectrum['coherence'][:,selected][:,:,selected]
            spectral_weights = spectrum['spectral_weights']
            truth = spectrum['pair_visibility2'][[0,0,1],[1,5,5]]
            later = tracking_geometry(array_positions,0.,dec,lat,elapsed+duration*1e-9)
            later_source = source_coherence_spectrum(array_positions,later['hour_angle_rad'],dec,lat,
                                                     source,instrument,args.source_case)
            source_change = float(abs(later_source['pair_visibility2']-spectrum['pair_visibility2']).max())
            for left in range(36):
                for right in range(left+1,36):
                    source_rows.append(dict(epoch=epoch_index,elapsed_s=elapsed,telescope_i=left+1,
                        telescope_j=right+1,visibility2=spectrum['pair_visibility2'][left,right]))
        calibrator = gamma if args.source_case=='fixed' else np.ones_like(gamma)
        null_coherence = np.eye(3) if args.source_case=='fixed' else np.broadcast_to(np.eye(3),gamma.shape)
        series = {}
        for case_index, (case, count, coherence, injection) in enumerate([
                ('calibration', args.calibration_records, calibrator, scale),
                ('holdout', args.records, gamma, scale),
                ('null', args.records, null_coherence, 1.)]):
            rng = np.random.default_rng(seeds[epoch_index*3+case_index])
            projections = []
            for record in range(count):
                raw = simulate_array_waveforms(rng, duration, star, background, coherence,
                    instrument, delays, injection, arrival_delay_rates_ns_per_s=rates,
                    spectral_weights=spectral_weights)
                variants = [align_waveforms(raw['adc_mv'], delays, dt,
                                           arrival_delay_rates_ns_per_s=rates),
                            align_waveforms(raw['adc_mv'], initial_delays, dt)]
                first = max(item['first_input_index'] for item in variants)
                last = min(item['first_input_index']+item['adc_mv'].shape[1] for item in variants)
                projected = []
                for item in variants:
                    start = first-item['first_input_index']
                    curves = correlate_blocks(item['adc_mv'][:, start:start+last-first], dt)
                    projected.append(curves['mean_correlation'] @ weights)
                projections.append(projected)
                for baseline in range(3):
                    rows.append(dict(epoch=epoch_index, elapsed_s=elapsed, case=case, record=record,
                        baseline=['0-1','0-2','1-2'][baseline], effective_duration_s=(last-first)*dt*1e-9,
                        tracked_projection=projected[0][baseline], initial_delay_projection=projected[1][baseline]))
            series[case] = np.asarray(projections)
        train = series['calibration'][:, 0]/(scale*calibration_truth)
        gain, gain_sem = train.mean(axis=0), train.std(axis=0, ddof=1)/np.sqrt(len(train))
        if np.any(gain <= 0):
            raise RuntimeError('追踪响应训练得到非正增益，不能继续归一化')
        holdout = series['holdout']/(scale*gain)
        mean = holdout.mean(axis=0)
        sem = holdout.std(axis=0, ddof=1)/np.sqrt(args.records)
        total_sem = np.sqrt(sem[0]**2+(truth*gain_sem/gain)**2)
        control_sem = np.sqrt(sem[1]**2+(mean[1]*gain_sem/gain)**2)
        null = series['null'][:, 0]/gain
        null_sem = null.std(axis=0, ddof=1)/np.sqrt(args.records)
        response_z, null_z = (mean[0]-truth)/total_sem, null.mean(axis=0)/null_sem
        checks.append(dict(epoch=epoch_index, max_response_abs_z=float(abs(response_z).max()),
                           max_null_abs_z=float(abs(null_z).max()),
                           passed=bool(np.all(abs(response_z)<5) and np.all(abs(null_z)<5))))
        for baseline in range(3):
            result_rows.append(dict(epoch=epoch_index, elapsed_s=elapsed,
                baseline=['0-1','0-2','1-2'][baseline], truth=truth[baseline],
                calibration_truth=calibration_truth[baseline],
                gain=gain[baseline], gain_sem=gain_sem[baseline],
                tracked_power=mean[0,baseline], tracked_total_sem=total_sem[baseline],
                initial_delay_power=mean[1,baseline], initial_delay_holdout_sem=sem[1,baseline],
                initial_delay_total_sem=control_sem[baseline],
                null_mean=null[:,baseline].mean(), null_sem=null_sem[baseline]))
        epoch_rows.append(dict(elapsed_s=elapsed, effective_record_duration_s=(last-first)*dt*1e-9,
            delay_rate_max_ns_per_s=float(abs(rates).max()),
            source_power_change_over_record=source_change,
            local_linear_error_bound_ns=float(.5*state['curvature_bound_ns_per_s2'].max()*(duration*1e-9)**2)))
        print(f'追踪时刻 {epoch_index+1}/{len(epochs)}：独立响应/零信号检查 {checks[-1]["passed"]}', flush=True)
    records, results = pd.DataFrame(rows), pd.DataFrame(result_rows)
    records.to_csv(output/'records.csv', index=False)
    results.to_csv(output/'response.csv', index=False)
    pd.DataFrame(epoch_rows).to_csv(output/'epochs.csv', index=False)
    if source_rows:
        pd.DataFrame(source_rows).to_csv(output/'source_predictions.csv', index=False)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for telescope, group in geometry.groupby('telescope_id'):
        axes[0].plot(group.elapsed_s/3600, group.initial_delay_residual_ns, lw=.7, alpha=.65)
    axes[0].set(xlabel='Time from transit [h]', ylabel='Delay minus initial value [ns]',
                title='36 telescope delays relative to Tel.1')
    for baseline, group in results.groupby('baseline'):
        divisor = group.truth if args.source_case=='fixed' else 1.
        line = axes[1].errorbar(group.elapsed_s/3600, group.tracked_power/divisor,
            yerr=1.96*group.tracked_total_sem/divisor, marker='o', ms=3, label=baseline+' tracked')
        if args.source_case=='fixed':
            axes[1].errorbar(group.elapsed_s/3600, group.initial_delay_power/divisor,
                yerr=1.96*group.initial_delay_total_sem/divisor, fmt='--',
                color=line[0].get_color(), label=baseline+' initial delay')
        else:
            axes[1].plot(group.elapsed_s/3600,group.truth,':',color=line[0].get_color(),label=baseline+' model')
    if args.source_case=='fixed':
        axes[1].axhline(1., color='k', lw=.7)
    axes[1].set(xlabel='Time from transit [h]', ylabel='Recovered / injected power' if args.source_case=='fixed' else 'Squared visibility',
                title='Independent short-record holdout')
    axes[1].legend(fontsize=7, ncol=2)
    for axis in axes:
        axis.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(output/'tracking_validation.png', dpi=180)
    plt.close(fig)
    paths = ['python/sii_unified.py','python/sii_observation.py','python/sii_layout.py',
             'python/sii_validation.py','tools/validate_sii_tracking.py','tests/test_sii_observation.py',
             'tests/test_sii_source_coherence.py']
    inputs = [item['path'] for item in manifest['files']]+[input_path,
        'configs/optics/lact2_measured_single_pixel_400nm.csv',
        'configs/optics/lact2_measured_single_pixel_400nm.provenance.json']
    summary = dict(seed=args.seed, records=args.records, calibration_records=args.calibration_records,
        source_case=args.source_case, source_parameters=asdict(source) if args.source_case!='fixed' else None,
        telescopes=36, raw_telescope_ids=[1,2,6], epochs=epoch_rows, checks=checks,
        code_sha256_lf={path:digest(path) for path in paths},
        input_sha256_lf={path:digest(path) for path in inputs},
        max_array_delay_rate_ns_per_s=float(abs(geometry.delay_rate_ns_per_s).max()),
        max_array_pair_delay_rate_ns_per_s=float(geometry.groupby('elapsed_s').delay_rate_ns_per_s.agg(lambda x:x.max()-x.min()).max()),
        max_array_24us_linear_error_bound_ns=float(.5*geometry.curvature_bound_ns_per_s2.max()*(duration*1e-9)**2),
        max_array_stale_pair_residual_ns=float(geometry.groupby('elapsed_s').initial_delay_residual_ns.agg(lambda x:x.max()-x.min()).max()),
        actual_raw_records=len(epochs)*(args.calibration_records+2*args.records),
        actual_raw_duration_s=len(epochs)*(args.calibration_records+2*args.records)*duration*1e-9,
        scope='Seven independent short-record epochs; not continuous six-hour ADC or full-array covariance. '+
              ('Fixed test coherence.' if args.source_case=='fixed' else 'Static uniform disk; rotating baselines and HBT spectral power average. Independent unresolved calibrator.'),
        timing_model='SI sidereal rotation at each epoch, local affine arrival-time map before optical response; receiver stellar rate includes Jacobian.',
        geometry_dec_rad=dec, geometry_lat_rad=lat, duration_ns=duration, pair_rate_scale=scale,
        star_rate_hz=star, background_rate_hz=background,
        fixed_coherence=gamma.tolist() if args.source_case=='fixed' else None,
        spectral_wavelength_nm=list(instrument.visibility_wavelength_nm),
        spectral_weights=list(instrument.visibility_spectral_weights),
        calibration='Fixed analytic projection; independent per-epoch training response. Null intervals conditional on training gain.',
        initial_delay_control='Same raw record and common reference-time interval, divided by tracked gain; control is not calibrated as a separate detector.')
    (output/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    if not all(item['passed'] for item in checks):
        raise RuntimeError('追踪验证未通过；保留全部结果供检查')


if __name__ == '__main__':
    main()
