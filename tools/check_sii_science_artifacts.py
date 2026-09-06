"""核对两份中文说明、执行Notebook及结果记录的可追溯性；不替代科学误差检验。"""
from pathlib import Path
import hashlib
import json
import math
import re
import sys

import nbformat

ROOT=Path(__file__).resolve().parents[1]


def main():
    for name in ['SII_PHYSICS_ZH.md','SII_IMPLEMENTATION_ZH.md']:
        path=ROOT/'docs'/name
        content=path.read_text(encoding='utf-8')
        assert content.count('$$') % 2 == 0, name
        assert content.replace('$$','').count('$') % 2 == 0, name
        table_columns=0
        for line in content.splitlines():
            if line.startswith('|'):
                columns=line.count('|')
                assert not table_columns or columns==table_columns, (name,line)
                table_columns=columns
            else:
                table_columns=0
        for target in re.findall(r'\]\(([^)]+)\)',content):
            if target.startswith(('http:','https:','#')):
                continue
            assert (path.parent/target.split('#')[0]).exists(), (name,target)
    notebook=nbformat.read(ROOT/'notebooks/sii_complete_waveform_report.ipynb',as_version=4)
    nbformat.validate(notebook)
    code=[cell for cell in notebook.cells if cell.cell_type=='code']
    assert all(cell.execution_count is not None for cell in code)
    assert [cell.execution_count for cell in code] == list(range(1,len(code)+1))
    assert not any(out.output_type=='error' for cell in code for out in cell.outputs)
    summary=json.loads((ROOT/'validation/sii_science/summary.json').read_text(encoding='utf-8'))
    for section in ['code_sha256_lf','derived_input_sha256_lf']:
        for path,expected in summary[section].items():
            content=(ROOT/path).read_bytes().replace(b'\r\n',b'\n')
            assert hashlib.sha256(content).hexdigest()==expected, path
    assert summary['scope']['notebook_execution']=='all cells executed in this run'
    assert all(value==0 for value in summary['missing_effects_disabled'].values())
    print(f'说明链接、公式分隔符和输入哈希通过；Notebook共{len(code)}个代码单元已顺序执行且无错误。')
    observation=ROOT/'validation/sii_observation'
    joint=json.loads((observation/'summary.json').read_text(encoding='utf-8'))
    for section in ['code_sha256_lf','input_sha256_lf']:
        for path,expected in joint[section].items():
            content=(ROOT/path).read_bytes().replace(b'\r\n',b'\n')
            assert hashlib.sha256(content).hexdigest()==expected, path
    assert len(joint['checks'])==5 and all(check['passed'] for check in joint['checks'])
    assert summary['joint_observation']['execution']=='fresh waveform simulation executed by this notebook'
    assert hashlib.sha256((observation/'summary.json').read_bytes().replace(b'\r\n',b'\n')).hexdigest()==summary['joint_observation']['summary_sha256_lf']
    sys.path.insert(0,str(ROOT/'python'))
    from sii_unified import Instrument, waveform_instrument_signature, detected_star_rate_hz
    instrument=Instrument.from_repository(ROOT)
    # 教学主线必须是同一基准源、同一条物理倍率记录，而不只是摘要数字相等。
    import numpy as np
    import sii_unified as sii
    from sii_observation import source_coherence_spectrum, tracking_geometry, simulate_array_photon_times
    from sii_performance import integer_align
    from sii_validation import analytic_waveform_calibration
    from sii_layout import read_corsika_layout
    import pandas as pd
    layout=read_corsika_layout(ROOT/'configs/arrays/lact36_20260906.input')
    pd.testing.assert_frame_equal(layout,pd.read_csv(ROOT/'configs/arrays/lact36_20260906_coordinates.csv'),
                                  check_exact=False,atol=1e-10,rtol=0)
    assert len(layout)==36
    lesson_dir=ROOT/'validation/sii_science'
    lesson=json.loads((lesson_dir/'walkthrough_parameters.json').read_text(encoding='utf-8'))
    pair=json.loads((lesson_dir/'walkthrough_pair_result.json').read_text(encoding='utf-8'))
    assert lesson==summary['walkthrough'] and pair==summary['walkthrough_pair']
    assert lesson['source_case']=='single_disk' and lesson['raw_pair_scale']==1.
    assert lesson['diameter_mas']==.16 and lesson['magnitude']==2. and lesson['pair']==[1,2]
    assert lesson['hours']==summary['observation']['hours']==6.
    assert lesson['segment_s']==summary['observation']['integration_s_per_measurement']==1200.
    assert lesson['wavelength_nm']==instrument.wavelength_nm and lesson['width_nm']==instrument.optical_width_nm
    observation=sii.Observation()
    source=sii.BinarySource(ab_magnitude=lesson['magnitude'],primary_diameter_mas=lesson['diameter_mas'])
    positions=layout[['east_m','north_m','up_m']].to_numpy()[:2]
    dec,lat=np.deg2rad([lesson['source_dec_deg'],lesson['site_lat_deg']])
    geometry=tracking_geometry(positions,lesson['hour_angle_rad'],dec,lat)
    spectrum=source_coherence_spectrum(positions,lesson['hour_angle_rad'],dec,lat,source,instrument,'single_disk')
    assert math.isclose(lesson['pair_power'],spectrum['pair_visibility2'][0,1],rel_tol=1e-12)
    assert math.isclose(lesson['star_rate_hz'],detected_star_rate_hz(2.,instrument),rel_tol=1e-12)
    with np.load(lesson_dir/'walkthrough_record.npz',allow_pickle=False) as record:
        np.testing.assert_allclose(record['arrival_delays_ns'],geometry['arrival_delays_ns'],atol=1e-12)
        np.testing.assert_allclose(record['time_ns'],(np.arange(6000)+.5)*instrument.sample_width_ns)
        assert record['adc_mv'].shape==(2,6000)
        # 独立重放公开种子，核对未放大对率的事件与实测SPE所产生的ADC。
        rng=np.random.default_rng(int(record['seed']))
        template=sii.load_measured_spe_template(instrument.spe_template_path)
        events,_=simulate_array_photon_times(rng,lesson['raw_record_ns'],lesson['star_rate_hz'],
            lesson['nsb_rate_hz'],spectrum['coherence'],instrument,geometry['arrival_delays_ns'],
            pair_rate_scale=1.,padding_ns=float(np.max(abs(template[0]))),
            spectral_weights=spectrum['spectral_weights'])
        for i,event in enumerate(events):
            np.testing.assert_allclose(event,record[f'telescope_{i+1}_events_ns'],atol=1e-12)
            rendered=sii.render_pe_waveform(rng,event,lesson['raw_record_ns'],instrument,template=template)
            np.testing.assert_allclose(rendered['adc_mv'],record['adc_mv'][i],rtol=1e-12,atol=1e-12)
        aligned,residual,exposure=integer_align(record['adc_mv'],record['arrival_delays_ns'],instrument.sample_width_ns)
    calibration,_=analytic_waveform_calibration(instrument,block_duration_ns=exposure*1e9,
                                               residual_delay_ns=float(residual[0]))
    lag,correlation=sii.waveform_cross_correlation(*aligned,instrument.sample_width_ns,200.)
    saved=pd.read_csv(lesson_dir/'walkthrough_pair_correlation.csv')
    np.testing.assert_allclose(saved.lag_ns,lag)
    np.testing.assert_allclose(saved.measured_C,correlation,rtol=1e-12,atol=1e-15)
    np.testing.assert_allclose(saved.expected_C,lesson['pair_power']*calibration.peak_per_visibility2,rtol=1e-12,atol=1e-15)
    estimate,sigma=sii.estimate_visibility2_gls(correlation[None,:],calibration)
    # 固定ADC的线性代数在不同BLAS线程下可有末位差异；容差远低于统计误差。
    assert math.isclose(pair['estimated_power'],estimate[0],rel_tol=1e-9)
    assert math.isclose(pair['sigma_short'],sigma,rel_tol=1e-10)
    assert math.isclose(pair['effective_record_s'],exposure,rel_tol=1e-12)
    assert math.isclose(pair['sigma_1200s'],sii.waveform_gls_weights(calibration,1200.)[1],rel_tol=1e-10)
    disk=pd.read_csv(lesson_dir/'walkthrough_disk_measurements.csv')
    expected=sii.segment_averaged_visibility2(disk,source,observation,instrument,'single_disk')
    np.testing.assert_allclose(disk.visibility2_true,expected,rtol=1e-11,atol=1e-13)
    assert len(disk)==11340 and (disk.visibility2_measured<0).any()
    ideal=pd.read_csv(lesson_dir/'walkthrough_ideal_hbt.csv')
    assert math.isclose(ideal.loc[ideal.lag_ps.abs().idxmin(),'g2_zero'],1.5,rel_tol=1e-12)
    print('教学主线：物理倍率事件和ADC可重放；同一记录的相关/GLS、基准圆盘的全部UV期望一致。')
    observation=ROOT/'validation/sii_observation'
    assert summary['observation']['telescopes']==len(layout)
    assert summary['observation']['uv_measurements']==len(layout)*(len(layout)-1)//2*18
    baselines=pd.read_csv(ROOT/'validation/sii_science/array_baselines.csv')
    assert len(baselines)==len(layout)*(len(layout)-1)//2
    assert not baselines.duplicated(['telescope_i','telescope_j']).any()
    assert summary['array_geometry']['baselines']==len(baselines)
    for bound in ['min','max']:
        assert math.isclose(summary['array_geometry'][f'{bound}_baseline_m'],
                            getattr(baselines.baseline_m,bound)(),rel_tol=1e-12)
    print('36台原始input与ENU坐标一致；630条基线和11340个UV测量的结果计数通过。')
    assert waveform_instrument_signature(instrument)==joint['instrument_signature']
    scenario=joint['instrument_scenario']
    assert math.isclose(detected_star_rate_hz(scenario['magnitude_ab'],instrument),
                        scenario['star_rate_hz'],rel_tol=1e-12)
    assert math.isclose(instrument.detected_nsb_rate_hz+instrument.dark_count_rate_hz,
                        scenario['background_rate_hz'],rel_tol=1e-12)
    records=pd.read_csv(observation/'records.csv')
    assert len(records)==3*(joint['calibration_records']+3*joint['records'])
    assert len(pd.read_csv(observation/'baseline_covariance.csv'))==27
    assert len(pd.read_csv(observation/'kernel_convergence.csv'))==18
    assert len(pd.read_csv(observation/'kernel_records.csv'))==2*joint['records']*6*3
    assert len(pd.read_csv(observation/'long_records.csv'))==3*joint['long_records']
    assert (observation/'observation_validation.png').is_file()
    print('三镜观测验证的输入/代码哈希、独立检查状态及输出条数通过。')
    for directory,summary_key in [('sii_tracking','tracking_observation'),
                                  ('sii_source_tracking','source_observation')]:
        tracking=ROOT/'validation'/directory
        tracked=json.loads((tracking/'summary.json').read_text(encoding='utf-8'))
        for section in ['code_sha256_lf','input_sha256_lf']:
            for path,expected in tracked[section].items():
                assert hashlib.sha256((ROOT/path).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==expected, path
        assert summary[summary_key]['execution']=='fresh sparse epoch waveforms executed by this notebook'
        assert hashlib.sha256((tracking/'summary.json').read_bytes().replace(b'\r\n',b'\n')).hexdigest()==summary[summary_key]['summary_sha256_lf']
        assert len(tracked['checks'])==7 and all(item['passed'] for item in tracked['checks'])
        assert len(pd.read_csv(tracking/'geometry.csv'))==361*36
        responses=pd.read_csv(tracking/'response.csv')
        samples=pd.read_csv(tracking/'records.csv')
        assert len(responses)==7*3 and not responses.duplicated(['epoch','baseline']).any()
        assert len(samples)==3*tracked['actual_raw_records']
        assert not samples.duplicated(['epoch','case','record','baseline']).any()
        assert samples.effective_duration_s.between(0,tracked['duration_ns']*1e-9,inclusive='right').all()
        # 从逐记录投影重算训练和留出均值，避免只核对摘要、漏掉表之间的不一致。
        for row in responses.itertuples():
            group=samples[(samples.epoch==row.epoch)&(samples.baseline==row.baseline)]
            training=group[group.case=='calibration'].tracked_projection
            holdout=group[group.case=='holdout']
            assert len(training)==tracked['calibration_records'] and len(holdout)==tracked['records']
            assert math.isclose(training.mean()/(tracked['pair_rate_scale']*row.calibration_truth),row.gain,rel_tol=1e-12)
            assert math.isclose(holdout.tracked_projection.mean()/(tracked['pair_rate_scale']*row.gain),
                                row.tracked_power,rel_tol=1e-12)
            assert math.isclose(holdout.initial_delay_projection.mean()/(tracked['pair_rate_scale']*row.gain),
                                row.initial_delay_power,rel_tol=1e-12,abs_tol=1e-14)
        assert tracked['actual_raw_records']==7*(tracked['calibration_records']+2*tracked['records'])
        assert math.isclose(tracked['actual_raw_duration_s'],tracked['actual_raw_records']*tracked['duration_ns']*1e-9)
        if directory=='sii_source_tracking':
            assert tracked['source_case']=='single_disk'
            assert (responses.calibration_truth==1.).all()
            predictions=pd.read_csv(tracking/'source_predictions.csv')
            assert len(predictions)==7*630
            assert not predictions.duplicated(['epoch','telescope_i','telescope_j']).any()
            for row in responses.itertuples():
                left,right=[tracked['raw_telescope_ids'][int(value)] for value in row.baseline.split('-')]
                prediction=predictions[(predictions.epoch==row.epoch)&(predictions.telescope_i==left)
                                       &(predictions.telescope_j==right)].visibility2.item()
                assert math.isclose(prediction,row.truth,rel_tol=1e-12)
        print(f'{directory}：七时刻的输入/代码/结果哈希、留出检查及实际曝光计数通过。')

    import numpy as np
    performance_dir=ROOT/'validation/sii_performance'
    performance=json.loads((performance_dir/'summary.json').read_text(encoding='utf-8'))
    assert all(performance['checks'].values())
    assert performance['instrument_signature']==waveform_instrument_signature(instrument)
    for section in ['code_sha256_lf','input_sha256_lf']:
        for path,expected in performance[section].items():
            assert hashlib.sha256((ROOT/path).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==expected,path
    for path,expected in performance['artifact_sha256'].items():
        assert hashlib.sha256((performance_dir/path).read_bytes()).hexdigest()==expected,path
    assert hashlib.sha256((performance_dir/'summary.json').read_bytes().replace(b'\r\n',b'\n')).hexdigest()==summary['performance']['summary_sha256_lf']
    table=pd.read_csv(performance_dir/'performance.csv')
    assert len(table)==54 and not table.duplicated(['hours','magnitude','diameter_mas','method']).any()
    assert (table.trials==500).all() and performance['scenarios']==27
    np.testing.assert_allclose(table.coverage,table.covered/table.trials,rtol=1e-12)
    assert ((table.coverage_ci_low<=table.coverage)&(table.coverage<=table.coverage_ci_high)).all()
    epochs=pd.read_csv(performance_dir/'raw_epochs.csv')
    baselines=pd.read_csv(performance_dir/'raw_baselines.csv')
    assert len(epochs)==3 and len(baselines)==3*630
    assert not baselines.duplicated(['epoch','telescope_i','telescope_j']).any()
    records=performance['raw_records_per_case']
    with np.load(performance_dir/'raw_projections.npz') as raw:
        for epoch in epochs.itertuples():
            train=raw[f'epoch_{epoch.epoch}_train']
            null=raw[f'epoch_{epoch.epoch}_null']
            holdout=raw[f'epoch_{epoch.epoch}_holdout']
            assert train.shape==holdout.shape==(records,630) and null.shape==(2*records,630)
            per_record=(train/300.).mean(axis=1)
            assert math.isclose(per_record.mean(),epoch.gain,rel_tol=1e-12)
            assert math.isclose(per_record.std(ddof=1)/math.sqrt(records),epoch.gain_sem,rel_tol=1e-12)
            selected=baselines[baselines.epoch==epoch.epoch]
            np.testing.assert_allclose(holdout.mean(axis=0)/300./epoch.gain,selected['mean'],atol=1e-13)
            np.testing.assert_allclose(raw[f'epoch_{epoch.epoch}_truth'],selected.truth,atol=1e-13)
            sigma=raw[f'epoch_{epoch.epoch}_sigma']
            assert math.isclose(float(np.sqrt(np.mean((null/sigma)**2))),epoch.null_sigma_scale,rel_tol=1e-12)
    assert performance['raw']['raw_records']==12*records
    assert math.isclose(performance['raw']['effective_duration_s'],float(np.sum(epochs.raw_records*epochs.effective_record_s)),rel_tol=1e-12)
    datasheet=json.loads((ROOT/'configs/sii/s17351_datasheet.json').read_text(encoding='utf-8'))
    assert datasheet['channels']*datasheet['microcells_per_channel']==instrument.microcells_per_pixel
    pdf=ROOT/datasheet['source']
    if pdf.exists():
        assert hashlib.sha256(pdf.read_bytes()).hexdigest()==datasheet['source_sha256']
    assert instrument.dark_count_rate_hz==0 and not datasheet['apply_to_default_instrument']
    print('36镜性能：54行结果、1152条共享ADC记录、逐记录标定、PDF来源及全部哈希通过。')


if __name__=='__main__':
    main()
