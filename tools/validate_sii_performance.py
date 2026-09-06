"""固定36镜SII角直径验证集；所有随机种子、场景和失败结果均保留。"""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
import hashlib
import json
import sys
from dataclasses import replace

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'python'))
import sii_unified as sii
import sii_performance as perf
from sii_layout import read_corsika_layout
from sii_observation import (tracking_geometry, source_coherence_spectrum,
                             simulate_array_waveforms, correlate_blocks)
from sii_validation import profile_grid_interval, proportion_interval, verify_main_parameters


def digest(path):
    return hashlib.sha256((ROOT/path).read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def raw_array_validation(instrument, layout, output, seed, records):
    """全部36镜共用ADC，点源训练、独立圆盘留出及零信号噪声校验。"""
    positions = layout[['east_m','north_m','up_m']].to_numpy()
    star = sii.detected_star_rate_hz(2., instrument)
    background = instrument.detected_nsb_rate_hz
    duration, scale, dt = 24000., 300., instrument.sample_width_ns
    source = sii.BinarySource(ab_magnitude=2., primary_diameter_mas=.16)
    arrays, rows, phase_rows, epoch_rows = {}, [], [], []
    seeds = np.random.SeedSequence(seed).spawn(9)
    for epoch, hour in enumerate([-.5, 0., .5]):
        state = tracking_geometry(positions, hour, np.deg2rad(22.), np.deg2rad(29.36))
        delay = state['arrival_delays_ns']
        _, residual, exposure = perf.integer_align(np.zeros((36,6000)), delay, dt)
        bank = perf.phase_template_bank(instrument, block_duration_ns=exposure*1e9)
        weights = perf.interpolate_bank(bank, residual, 'weights')
        sigma = perf.interpolate_bank(bank, residual, 'sigma_block')
        spectrum = source_coherence_spectrum(positions, hour, np.deg2rad(22.), np.deg2rad(29.36),
                                             source, instrument, 'single_disk')
        pairs = np.array(list(__import__('itertools').combinations(range(36),2)))
        truth = spectrum['pair_visibility2'][pairs[:,0],pairs[:,1]]
        data = {}
        for case_index, case in enumerate(['train','null','holdout']):
            rng = np.random.default_rng(seeds[epoch*3+case_index])
            count = records if case != 'null' else 2*records
            gamma = (np.ones((36,36)) if case=='train' else np.eye(36) if case=='null' else spectrum['coherence'])
            spectral = spectrum['spectral_weights'] if case=='holdout' else None
            data[case] = []
            for record in range(count):
                waveform = simulate_array_waveforms(rng, duration, star, background, gamma, instrument,
                    arrival_delays_ns=delay, pair_rate_scale=1. if case=='null' else scale,
                    arrival_delay_rates_ns_per_s=state['arrival_delay_rates_ns_per_s'], spectral_weights=spectral)
                aligned, _, actual = perf.integer_align(waveform['adc_mv'], delay, dt)
                curves = correlate_blocks(aligned, dt)['mean_correlation']
                data[case].append(np.einsum('ij,ij->i', curves, weights))
            data[case] = np.asarray(data[case])
            arrays[f'epoch_{epoch}_{case}'] = data[case]
            print(f'36镜共享ADC {epoch+1}/3 {case}: {count}条完成', flush=True)
        # 每条记录内630条基线共享波形；用“逐记录阵列均值”计算标定SEM，不能当独立630倍。
        train = (data['train']/scale).mean(axis=1)
        gain, gain_sem = train.mean(), train.std(ddof=1)/np.sqrt(records)
        null_energy = np.mean((data['null']/sigma)**2, axis=1)
        noise_scale = np.sqrt(null_energy.mean())
        noise_sem = null_energy.std(ddof=1)/np.sqrt(len(null_energy))/(2*noise_scale)
        arrays[f'epoch_{epoch}_truth'] = truth
        arrays[f'epoch_{epoch}_sigma'] = sigma
        arrays[f'epoch_{epoch}_residual_ns'] = residual
        epoch_rows.append(dict(epoch=epoch, hour_angle_rad=hour, gain=gain, gain_sem=gain_sem,
            null_sigma_scale=noise_scale, null_sigma_scale_sem=noise_sem,
            effective_record_s=exposure, raw_records=4*records,
            max_record_delay_drift_ns=float(np.ptp(state['arrival_delay_rates_ns_per_s'])*duration*1e-9)))
        for baseline, (i,j) in enumerate(pairs):
            hold = data['holdout'][:,baseline]/scale/gain
            null = data['null'][:,baseline]/sigma[baseline]
            sem = np.sqrt(hold.var(ddof=1)/records+(truth[baseline]*gain_sem/gain)**2)
            rows.append(dict(epoch=epoch, telescope_i=int(i+1), telescope_j=int(j+1),
                residual_ns=residual[baseline], truth=truth[baseline], mean=hold.mean(), sem=sem,
                null_mean_sigma=null.mean(), null_std_ratio=null.std(ddof=1)))
        # 相位分组留出只作验证，不用目标真值修正训练增益。
        bins = np.floor((residual+dt)/(2*dt)*8).astype(int).clip(0,7)
        for group in np.unique(bins):
            chosen = bins==group
            expected = truth[chosen].mean()
            measured = (data['holdout'][:,chosen]/scale/gain).mean(axis=1)
            sem = np.sqrt(measured.var(ddof=1)/records+(expected*gain_sem/gain)**2)
            null_group=(data['null'][:,chosen]/sigma[chosen]).mean(axis=1)
            null_z=null_group.mean()/(null_group.std(ddof=1)/np.sqrt(len(null_group)))
            phase_rows.append(dict(epoch=epoch, phase_bin=int(group), baselines=int(chosen.sum()),
                truth=expected, mean=measured.mean(), sem=sem, z=(measured.mean()-expected)/sem,
                null_z=null_z))
    np.savez_compressed(output/'raw_projections.npz', **arrays)
    pd.DataFrame(rows).to_csv(output/'raw_baselines.csv',index=False)
    pd.DataFrame(phase_rows).to_csv(output/'raw_phase_holdout.csv',index=False)
    epochs = pd.DataFrame(epoch_rows)
    epochs.to_csv(output/'raw_epochs.csv',index=False)
    # 一个统一响应由独立训练决定；跨时刻离散作为额外标定不确定度保守保留。
    inverse = 1/epochs.gain_sem.to_numpy()**2
    gain = float(np.average(epochs.gain, weights=inverse))
    gain_sem = float(np.sqrt(1/inverse.sum()+np.var(epochs.gain,ddof=1)))
    noise_scale = float(np.sqrt(np.mean(epochs.null_sigma_scale**2)))
    noise_sem = float(np.sqrt(np.mean(epochs.null_sigma_scale_sem**2)+np.var(epochs.null_sigma_scale,ddof=1)))
    return dict(gain=gain, gain_relative_sigma=gain_sem/gain, noise_scale=noise_scale,
        noise_relative_sigma=noise_sem/noise_scale, raw_records=12*records,
        raw_duration_s=12*records*duration*1e-9,
        effective_duration_s=float(np.sum(epochs.effective_record_s*epochs.raw_records)),
        phase_holdout_max_abs_z=float(np.max(np.abs(pd.DataFrame(phase_rows).z))),
        phase_null_max_abs_z=float(np.max(np.abs(pd.DataFrame(phase_rows).null_z))),
        scope='36 shared ADC channels, 3 short-record epochs; weak pair process. Full night uses sufficient statistics.')


def run(args):
    output = ROOT/args.output
    output.mkdir(parents=True,exist_ok=True)
    instrument = sii.Instrument.from_repository(ROOT)
    layout = read_corsika_layout(ROOT/'configs/arrays/lact36_20260906.input')
    manifest = verify_main_parameters(ROOT)
    raw = raw_array_validation(instrument, layout, output, args.seed, args.records)
    print('36镜训练结果', raw, flush=True)
    bank = perf.phase_template_bank(instrument)
    pd.DataFrame({'residual_ns':bank['phase_ns'], 'sigma_24us':bank['sigma_block']}).to_csv(output/'phase_templates.csv',index=False)
    # 0至0.64 mas预先固定，三个真值都是网格节点；不围着真值窄化搜索区间。
    grid = np.linspace(0., .64, 1281)
    magnitudes, diameters, hours = [2.,4.,6.], [.08,.16,.32], [1.,3.,6.]
    result_rows, budget_rows, approximation_rows = [], [], []
    rng = np.random.default_rng(args.seed+100)
    star_reference = sii.detected_star_rate_hz(2.,instrument)
    total_reference = star_reference+instrument.detected_nsb_rate_hz
    raw_scale = raw['noise_scale']/raw['gain']
    for magnitude in magnitudes:
        budget_rows.append(dict(magnitude=magnitude, **perf.weak_light_covariance_budget(instrument,magnitude)))
    for hour in hours:
        observation = sii.Observation(hours_per_night=hour, segment_s=1200.)
        uvw = sii.generate_uvw(layout, observation, instrument)
        sigma_base = perf.tracked_segment_precision(uvw, observation, instrument, bank)*raw_scale
        sigma_fine = perf.tracked_segment_precision(uvw, observation, instrument, bank, time_nodes=2400)*raw_scale
        approximation_rows.append(dict(hours=hour, effect='phase_time_1200_to_2400',
            max_abs_power_delta=0., max_relative_sigma_delta=float(np.max(abs(sigma_base/sigma_fine-1))),
            diameter_mas=0., fitted_bias_mas=0.))
        print(f'计算{hour:g}小时的候选直径模型',flush=True)
        models = perf.disk_model_grid(uvw, observation, instrument, grid)
        # 预定33个非网格中点的真实模型对照；检查似然细化所依赖的平滑插值。
        from scipy.interpolate import CubicSpline
        midpoint=grid[np.linspace(0,len(grid)-2,33).astype(int)]+.00025
        midpoint_models=perf.disk_model_grid(uvw,observation,instrument,midpoint)
        interpolated=CubicSpline(grid,models,axis=0)(midpoint)
        interpolation_error=float(np.max(np.linalg.norm((midpoint_models-interpolated)/sigma_base,axis=1)))
        approximation_rows.append(dict(hours=hour,effect='diameter_model_interpolation_sigma',diameter_mas=0.,
            max_abs_power_delta=float(np.max(abs(midpoint_models-interpolated))),
            max_relative_sigma_delta=interpolation_error,fitted_bias_mas=0.))
        for magnitude in magnitudes:
            star = sii.detected_star_rate_hz(magnitude,instrument)
            rate_scale = (star+instrument.detected_nsb_rate_hz)/star**2/(total_reference/star_reference**2)
            bound = perf.weak_light_covariance_budget(instrument,magnitude)['bartlett_covariance_operator_bound']
            # 把共镜Bartlett谱界作为保守方差放大，覆盖界内任何协方差方向。
            sigma = sigma_base*rate_scale
            assumed_scale = np.sqrt(1+bound)*(1+1.96*raw['noise_relative_sigma'])
            for diameter in diameters:
                target = models[int(round(diameter/.0005))]
                profile, fixed, error = perf.compressed_profiles(models,sigma,target,
                    raw['gain_relative_sigma'],args.trials,rng,
                    noise_relative_sigma=raw['noise_relative_sigma'],assumed_sigma_scale=assumed_scale)
                for name, statistics in [('profile_gain',profile),('fixed_gain',fixed)]:
                    intervals = perf.refined_profile_intervals(grid,statistics)
                    estimates = np.array([item['estimate'] for item in intervals])
                    covered = sum(item['low']<=diameter<=item['high'] for item in intervals)
                    lower,upper = proportion_interval(covered,len(intervals))
                    result_rows.append(dict(hours=hour,magnitude=magnitude,diameter_mas=diameter,method=name,
                        trials=args.trials,bias_mas=float(estimates.mean()-diameter),
                        rmse_mas=float(np.sqrt(np.mean((estimates-diameter)**2))),
                        std_mas=float(estimates.std(ddof=1)), median_width_mas=float(np.median([i['high']-i['low'] for i in intervals])),
                        coverage=covered/args.trials,covered=covered,coverage_ci_low=lower,coverage_ci_high=upper,
                        boundary_intervals=sum(i['touches_boundary'] for i in intervals),
                        disconnected_intervals=sum(i['disconnected'] for i in intervals),
                        compression_error_sigma=error,sigma_visibility2_median=float(np.median(sigma*assumed_scale)),
                        assumed_noise_inflation=assumed_scale))
                print(f'性能 {hour:g}h m={magnitude:g} d={diameter:g} mas 完成',flush=True)
        if hour==6.:
            for diameter in diameters:
                target = models[int(round(diameter/.0005))]
                aperture = perf.aperture_disk_power(uvw,observation,instrument,diameter)
                fine_aperture = perf.aperture_disk_power(uvw,observation,instrument,diameter,radial_nodes=18,angle_nodes=36)
                # 系统偏差用真值附近的局部线性响应，不让0.0005 mas网格掩盖小偏差。
                step = 1e-5
                derivative_models = perf.disk_model_grid(uvw,observation,instrument,[diameter+step,diameter-step])
                derivative = (derivative_models[0]-derivative_models[1])/(2*step)
                design = np.column_stack([derivative,target])/sigma_base[:,None]
                normal = design.T@design
                normal[1,1] += 1/raw['gain_relative_sigma']**2
                bias = np.linalg.solve(normal, design.T@((aperture-target)/sigma_base))[0]
                approximation_rows.append(dict(hours=hour,effect='uniform_pupil_radius_4m',diameter_mas=diameter,
                    max_abs_power_delta=float(np.max(abs(aperture-target))),max_relative_sigma_delta=0.,fitted_bias_mas=float(bias)))
                approximation_rows.append(dict(hours=hour,effect='pupil_quadrature_12x24_to_18x36',diameter_mas=diameter,
                    max_abs_power_delta=float(np.max(abs(aperture-fine_aperture))),max_relative_sigma_delta=0.,fitted_bias_mas=0.))
            # 光学核作为已计算响应；去掉它是有方向的对照，不是其真实误差分布。
            no_optical = replace(instrument,optical_timing_kernel_path=None)
            no_optical_bank = perf.phase_template_bank(no_optical)
            sigma_no_optical = perf.tracked_segment_precision(uvw,observation,no_optical,no_optical_bank)
            approximation_rows.append(dict(hours=hour,effect='omit_optical_timing_kernel',diameter_mas=0.,
                max_abs_power_delta=0.,max_relative_sigma_delta=float(np.max(abs(sigma_no_optical/(sigma_base/raw_scale)-1))),fitted_bias_mas=0.))
    # 规格书暗计数只作为单独25℃情景。PDE/SPE保持main；串扰不重复计入电荷分布。
    detector_rows=[]
    for magnitude in magnitudes:
        star=sii.detected_star_rate_hz(magnitude,instrument)
        for channels in [1,8]:
            for rating in ['typical','max']:
                scenario=perf.datasheet_dark_scenario(instrument,ROOT,channels,rating)
                detector_rows.append(dict(magnitude=magnitude,summed_channels=channels,rating=rating,
                    dark_rate_mhz=scenario.dark_count_rate_hz/1e6,
                    sigma_ratio=(star+instrument.detected_nsb_rate_hz+scenario.dark_count_rate_hz)/(star+instrument.detected_nsb_rate_hz),
                    temperature_c=25,overvoltage_v=8.5))
    performance=pd.DataFrame(result_rows)
    performance.to_csv(output/'performance.csv',index=False)
    pd.DataFrame(budget_rows).to_csv(output/'covariance_budget.csv',index=False)
    pd.DataFrame(approximation_rows).to_csv(output/'approximation_budget.csv',index=False)
    detector=pd.DataFrame(detector_rows)
    detector.to_csv(output/'datasheet_scenarios.csv',index=False)
    fig,axes=plt.subplots(1,3,figsize=(13,3.8))
    for diameter in diameters:
        group=performance[(performance.method=='profile_gain')&(performance.hours==6)&(performance.diameter_mas==diameter)]
        axes[0].plot(group.magnitude,group.rmse_mas/diameter,'o-',label=f'{diameter:g} mas')
        axes[1].errorbar(group.magnitude,group.coverage,
            yerr=np.array([group.coverage-group.coverage_ci_low,group.coverage_ci_high-group.coverage]),fmt='o-',label=f'{diameter:g} mas')
    axes[0].set(xlabel='AB magnitude',ylabel='Relative diameter RMSE',yscale='log',title='36 telescopes, 6 h')
    axes[1].axhline(.95,color='k',ls='--');axes[1].set(xlabel='AB magnitude',ylabel='95% interval coverage',ylim=(0,1.02))
    for rating,group in detector[detector.summed_channels==8].groupby('rating'):
        axes[2].plot(group.magnitude,group.sigma_ratio,'o-',label=rating)
    axes[2].set(xlabel='AB magnitude',ylabel='Statistical sigma / dark-free sigma',title='25 C, 8 channels summed')
    for axis in axes:
        axis.grid(alpha=.2);axis.legend(fontsize=8)
    fig.tight_layout();fig.savefig(output/'performance.png',dpi=180);plt.close(fig)
    paths=['python/sii_performance.py','python/sii_validation.py','python/sii_unified.py','python/sii_observation.py',
           'python/sii_layout.py','tools/validate_sii_performance.py','tests/test_sii_performance.py']
    inputs=[entry['path'] for entry in manifest['files']]+['configs/sii/s17351_datasheet.json',
        'configs/arrays/lact36_20260906.input','configs/optics/lact2_measured_single_pixel_400nm.csv',
        'configs/optics/lact2_measured_single_pixel_400nm.provenance.json']
    summary=dict(seed=args.seed,trials=args.trials,raw_records_per_case=args.records,raw=raw,
        code_sha256_lf={p:digest(p) for p in paths},input_sha256_lf={p:digest(p) for p in inputs},
        science_scope='36-telescope sensitivity and uniform-disk diameter, conditional on main response and disabled unknown electronics',
        scenarios=27,methods=['profile_gain','fixed_gain'],magnitudes=magnitudes,diameters_mas=diameters,hours=hours,
        diameter_grid_mas=[0.,.64,.0005],likelihood_refinement_step_mas=.00005,
        pupil='Uniform circular radius 4m comparison; actual pupil weights unavailable',
        covariance='Bartlett weak-light spectral operator bound propagated as variance inflation; separate exact thermal count fourth-moment check',
        default_electronics_changed=False,
        instrument_signature=sii.waveform_instrument_signature(instrument),
        artifact_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir()
                         if p.suffix in ('.csv','.npz','.png')},
        checks=dict(raw_phase_holdout=raw['phase_holdout_max_abs_z']<5,
            raw_phase_null=raw['phase_null_max_abs_z']<5,
            phase_quadrature=max(r['max_relative_sigma_delta'] for r in approximation_rows if r['effect']=='phase_time_1200_to_2400')<.002,
            pupil_quadrature=max(r['max_abs_power_delta'] for r in approximation_rows if r['effect']=='pupil_quadrature_12x24_to_18x36')<1e-6,
            model_compression=float(performance.compression_error_sigma.max())<1e-5,
            diameter_interpolation=max(r['max_relative_sigma_delta'] for r in approximation_rows if r['effect']=='diameter_model_interpolation_sigma')<.001,
            finite_results=bool(np.all(np.isfinite(performance.select_dtypes('number')))) ))
    (output/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    if not all(summary['checks'].values()):
        raise RuntimeError('固定验证集有未通过项，结果已经保存，不能宣称全部完成')
    print('固定36镜性能验证完成',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='validation/sii_performance')
    parser.add_argument('--seed',type=int,default=20260909)
    parser.add_argument('--records',type=int,default=96)
    parser.add_argument('--trials',type=int,default=500)
    args=parser.parse_args()
    if args.records<8 or args.trials<20:
        parser.error('至少8条波形记录及20次参数重复')
    run(args)
