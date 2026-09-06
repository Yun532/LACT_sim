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
    from sii_layout import read_corsika_layout
    import pandas as pd
    layout=read_corsika_layout(ROOT/'configs/arrays/lact36_20260906.input')
    pd.testing.assert_frame_equal(layout,pd.read_csv(ROOT/'configs/arrays/lact36_20260906_coordinates.csv'),
                                  check_exact=False,atol=1e-10,rtol=0)
    assert len(layout)==36
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
    tracking=ROOT/'validation/sii_tracking'
    tracked=json.loads((tracking/'summary.json').read_text(encoding='utf-8'))
    for section in ['code_sha256_lf','input_sha256_lf']:
        for path,expected in tracked[section].items():
            assert hashlib.sha256((ROOT/path).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==expected, path
    assert summary['tracking_observation']['execution']=='fresh sparse epoch waveforms executed by this notebook'
    assert hashlib.sha256((tracking/'summary.json').read_bytes().replace(b'\r\n',b'\n')).hexdigest()==summary['tracking_observation']['summary_sha256_lf']
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
        assert math.isclose(training.mean()/(tracked['pair_rate_scale']*row.truth),row.gain,rel_tol=1e-12)
        assert math.isclose(holdout.tracked_projection.mean()/(tracked['pair_rate_scale']*row.gain),
                            row.tracked_power,rel_tol=1e-12)
        assert math.isclose(holdout.initial_delay_projection.mean()/(tracked['pair_rate_scale']*row.gain),
                            row.initial_delay_power,rel_tol=1e-12,abs_tol=1e-14)
    assert tracked['actual_raw_records']==7*(tracked['calibration_records']+2*tracked['records'])
    assert math.isclose(tracked['actual_raw_duration_s'],tracked['actual_raw_records']*tracked['duration_ns']*1e-9)
    print('七时刻追踪验证的输入/代码/结果哈希、留出检查及实际曝光计数通过。')


if __name__=='__main__':
    main()
