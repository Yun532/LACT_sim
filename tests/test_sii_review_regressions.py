"""跨入口、非默认仪器与共享误差检查，防止只验证默认局部公式。"""
from dataclasses import replace
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'python'))
import sii_unified as sii
import sii_reconstruction as reco
from sii_validation import analytic_waveform_calibration
from sii_performance import phase_template_bank, interpolate_bank


def small_observation(instrument, hours=1., segment=1200.):
    layout = pd.DataFrame({'east_m':[0.,100.,30.], 'north_m':[0.,0.,120.], 'up_m':[0.,0.,0.]})
    observation = sii.Observation(hours_per_night=hours, segment_s=segment)
    return sii.generate_uvw(layout, observation, instrument), observation


def test_phase_response_and_main_long_exposure():
    instrument = sii.Instrument.from_repository(ROOT)
    cal, _ = analytic_waveform_calibration(instrument)
    bank = phase_template_bank(instrument, block_duration_ns=20000., calibration=cal)
    # 独立连续积分的2 ns信号不能用零相位模板；匹配模板应恢复单位可见度。
    shifted, _ = analytic_waveform_calibration(instrument, residual_delay_ns=2.)
    zero_weights, _ = sii.waveform_gls_weights(cal,1200.)
    assert zero_weights@shifted.peak_per_visibility2 < .8
    assert np.isclose(interpolate_bank(bank,2.,'weights')@shifted.peak_per_visibility2,1,rtol=1e-8)
    uvw, observation = small_observation(instrument)
    frame, _ = sii.simulate_uv_observation(uvw, sii.BinarySource(), observation, instrument,
        estimator='waveform_gls', waveform_calibration=cal)
    zero_sigma = sii.waveform_gls_weights(cal,1200)[1]
    assert np.all(frame.sigma_visibility2_stat > zero_sigma*1.02)
    assert np.ptp(frame.sigma_visibility2_stat) > 1e-7


def test_full_uv_roundtrip_keeps_integral_and_shared_prior(tmp_path):
    instrument = sii.Instrument.from_repository(ROOT)
    uvw, observation = small_observation(instrument,6.,21600.)
    frame, _ = sii.simulate_uv_observation(uvw, sii.BinarySource(), observation, instrument)
    frame['calibration_relative_sigma'] = .1
    frame['visibility2_measured'] = [-.2,.3,.8]
    uv = sii.prepare_reconstruction_uv(frame, cell_mlambda=None)
    path = tmp_path/'observation.npz'
    reco.write_uv_data(path, uv)
    loaded = reco.read_uv_measurements(path, 'unused')
    for key in vars(uv):
        if key == 'sampling':
            for first, second in zip(uv.sampling, loaded.sampling):
                np.testing.assert_array_equal(first,second)
        elif getattr(uv,key) is not None:
            np.testing.assert_array_equal(getattr(uv,key),getattr(loaded,key))
    image = np.zeros((12,12)); image[5:7,5:7] = .25
    kernel = reco.power_sampling_kernel(uv,12,.7)
    np.testing.assert_array_equal(kernel,reco.power_sampling_kernel(loaded,12,.7))
    point_kernel = reco.power_sampling_kernel(replace(uv,sampling=None),12,.7)
    assert np.max(abs(kernel-point_kernel)) > .1
    prediction = reco.power_from_image(kernel,image)
    assert not np.allclose(prediction,reco.power_from_image(point_kernel,image),rtol=1e-8,atol=1e-8)
    assert reco.statistical_loss(uv,prediction-uv.visibility_abs2)[0] == reco.statistical_loss(
        loaded,prediction-loaded.visibility_abs2)[0]
    frame.to_csv(tmp_path/'unsafe.csv',index=False)
    with pytest.raises(ValueError,match='NPZ'):
        reco.read_uv_measurements(tmp_path/'unsafe.csv','visibility2_measured')


def test_baseline_zero_covariance_survives_grouping_and_disk(tmp_path):
    n = 12
    frame = pd.DataFrame(dict(u_lambda=np.repeat([1e6,2e6],n//2),v_lambda=0.,
        visibility2_measured=.3,sigma_visibility2=np.hypot(.01,.05),
        sigma_visibility2_stat=.01,baseline_zero_point_sigma=.05,
        telescope_i='A',telescope_j='B'),index=np.arange(n))
    uv = sii.prepare_reconstruction_uv(frame,cell_mlambda=.1)
    expected = np.diag(np.full(2,.01**2/(n//2)))+np.full((2,2),.05**2)
    np.testing.assert_allclose(uv.covariance,expected)
    rng = np.random.default_rng(120)
    trials = .3+rng.normal(0,.05,(20000,1))+rng.normal(0,.01,(20000,n))
    merged = trials.reshape(20000,2,n//2).mean(axis=2)
    np.testing.assert_allclose(np.cov(merged.T),expected,rtol=.025)
    reco.write_uv_data(tmp_path/'shared.npz',uv)
    np.testing.assert_array_equal(reco.read_uv_data(tmp_path/'shared.npz').covariance,uv.covariance)


def test_matched_noise_is_only_counted_once():
    base = sii.Instrument.from_repository(ROOT,electronic_noise_rms_mv=5.)
    instrument = sii.with_matched_effective_bandwidth(base)
    uvw, observation = small_observation(instrument)
    frame, _ = sii.simulate_uv_observation(uvw,sii.BinarySource(),observation,instrument)
    star = sii.detected_star_rate_hz(2.,instrument)
    total = star+instrument.detected_nsb_rate_hz+instrument.dark_count_rate_hz
    expected = total*instrument.excess_noise_factor**2/(star**2*instrument.coherence_area_s
                *np.sqrt(2*instrument.electronics_bandwidth_hz*observation.segment_s))
    np.testing.assert_allclose(frame.sigma_visibility2_stat,expected,rtol=2e-5)


def test_matched_bandwidth_responds_to_intrinsic_jitter():
    instrument = sii.Instrument.from_repository(ROOT)
    jitter = replace(instrument,intrinsic_time_jitter_ns=1.)
    assert sii.matched_effective_bandwidth_hz(jitter) < .85*sii.matched_effective_bandwidth_hz(instrument)
    first,_ = analytic_waveform_calibration(instrument)
    second,_ = analytic_waveform_calibration(jitter)
    assert sii.waveform_gls_weights(second,1200)[1] > sii.waveform_gls_weights(first,1200)[1]*1.1


def test_narrow_off_node_bandpass_uses_same_positive_spectrum(tmp_path):
    path = tmp_path/'narrow.csv'
    path.write_text('wavelength_nm,transmission\n399.19,0\n399.2,1\n399.21,0\n')
    instrument = sii.Instrument.from_repository(ROOT,sii_bandpass_path=str(path))
    nodes = np.array(instrument.visibility_wavelength_nm)
    weights = np.array(instrument.visibility_spectral_weights)
    assert np.all(np.isfinite(weights)) and np.isclose(weights.sum(),1)
    assert np.all((nodes > 399.19)&(nodes < 399.21))
    assert abs(weights@nodes-399.2) < 1e-5
    assert sii.detected_star_rate_hz(2.,instrument) > 0
    # 三角透射的平方积分/一阶积分平方给出4/(3*全宽)。
    expected = instrument.polarization_factor*(399.2e-9)**2/sii.C_M_S*4/(3*.02e-9)
    assert np.isclose(instrument.coherence_area_s,expected,rtol=1e-5)
