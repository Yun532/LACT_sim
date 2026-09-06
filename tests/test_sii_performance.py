"""相位响应、热光四阶误差及圆瞳平均的独立极限检验。"""
from pathlib import Path
import sys
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'python'))
import sii_unified as sii
import sii_performance as perf
from sii_validation import analytic_waveform_calibration


def test_phase_shift_sign_and_integer_alignment():
    x=np.zeros((2,100));x[0,30]=1;x[1,37]=1
    aligned,residual,exposure=perf.integer_align(x,[0,29.],4.)
    assert np.argmax(aligned[0])==np.argmax(aligned[1])
    assert residual[0]==1 and exposure==pytest.approx(93*4e-9)
    instrument=sii.Instrument.from_repository(ROOT)
    minus,_=analytic_waveform_calibration(instrument,residual_delay_ns=-1.)
    plus,_=analytic_waveform_calibration(instrument,residual_delay_ns=1.)
    np.testing.assert_allclose(minus.peak_per_visibility2,plus.peak_per_visibility2[::-1],rtol=1e-7,atol=1e-14)
    assert plus.peak_per_visibility2[49]>plus.peak_per_visibility2[51]


def test_uniform_pupil_against_independent_disk_pairs():
    points,weights=perf.pupil_difference_quadrature(4.,18,36)
    assert weights.sum()==pytest.approx(1.)
    # 两个均匀圆盘差向量的E|d|²=R²；用正亮度圆盘可见度与直接随机点对比较。
    assert np.sum(weights*np.sum(points**2,axis=1))==pytest.approx(16.,rel=1e-5)
    rng=np.random.default_rng(810)
    radius=4*np.sqrt(rng.random((2,200000)))
    angle=rng.uniform(0,2*np.pi,(2,200000))
    delta=np.array([np.diff(radius*np.cos(angle),axis=0)[0],np.diff(radius*np.sin(angle),axis=0)[0]]).T
    wave=400e-9
    direct=sii.uniform_disk_visibility(np.hypot(100+delta[:,0],delta[:,1])/wave,.32)**2
    integral=weights@sii.uniform_disk_visibility(np.hypot(100+points[:,0],points[:,1])/wave,.32)**2
    assert abs(direct.mean()-integral)<5*direct.std(ddof=1)/np.sqrt(len(direct))
    zero,w=perf.pupil_difference_quadrature(0.)
    np.testing.assert_array_equal(zero,[[0.,0.]])


def test_thermal_covariance_count_formula_and_spectral_bound():
    # 未分辨单模场的正光强Gamma分布可给出全部阶矩；直接数值积分条件Poisson的中心矩。
    from scipy.special import roots_genlaguerre, gamma
    instrument=sii.Instrument.from_repository(ROOT)
    budget=perf.weak_light_covariance_budget(instrument,2.)
    assert 0<budget['bartlett_covariance_operator_bound']<.01
    # 选择人为强占据度来检验四阶公式，避免真实微弱效应被MC噪声淹没。
    mu,background,modes=.7,.2,3.
    nodes,weights=roots_genlaguerre(24,modes-1)
    intensity=mu*nodes/modes
    shift=intensity-mu
    conditional_variance=intensity+background
    v=mu+background+mu**2/modes
    a,b,c=mu**2/modes,mu**3/modes**2,mu**4/modes**3
    direct=weights@((conditional_variance+shift**2)**2)/gamma(modes)-a*a
    assert direct==pytest.approx(v*v+a*a+a+4*b+6*c)
    shared=weights@((conditional_variance+shift**2)*shift**2)/gamma(modes)-a*a
    assert shared==pytest.approx(v*a+a*a+2*b+6*c)


def test_datasheet_conditions_and_default_isolation():
    instrument=sii.Instrument.from_repository(ROOT)
    assert instrument.dark_count_rate_hz==0
    scenario=perf.datasheet_dark_scenario(instrument,ROOT,8)
    assert scenario.dark_count_rate_hz==9.6e6
    assert scenario.sipm_crosstalk_probability==0
    assert scenario.sipm_afterpulse_probability==0
    assert scenario.spe_template_path==instrument.spe_template_path
    with pytest.raises(ValueError):
        perf.datasheet_dark_scenario(instrument,ROOT,9)


def test_compressed_likelihood_retains_model_geometry():
    # 线性两参数模型的QR子空间应精确保留所有模型间欧氏距离。
    x=np.linspace(0,1,100)
    grid=np.linspace(0,1,25)
    models=np.array([np.sin(x)+v*np.cos(x) for v in grid])
    profile,fixed,error=perf.compressed_profiles(models,np.ones(100),models[12],.02,100,np.random.default_rng(5),anchors=5)
    assert error<1e-12
    assert profile.shape==fixed.shape==(100,25)


def test_likelihood_refinement_recovers_subgrid_minimum():
    grid=np.linspace(0,1,21)
    row=((grid-.513)/.07)**2
    result=perf.refined_profile_intervals(grid,[row])[0]
    assert abs(result['estimate']-.513)<.0026
    assert result['low']==pytest.approx(.513-1.95996398454*.07,abs=1e-4)
