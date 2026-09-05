"""独立物理恒等式、共同标定误差和参数来源的回归测试。"""
from dataclasses import replace
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'python'))
import sii_unified as sii
from sii_reconstruction import UvData, statistical_loss, profile_calibration_gain
from sii_validation import (analytic_waveform_calibration, thermal_mode_counts,
                            verify_main_parameters, profile_model_grid, profile_grid_interval)


def test_main_parameters_and_missing_effects_zero():
    assert len(verify_main_parameters(ROOT)['files']) >= 10
    instrument = sii.Instrument.from_repository(ROOT)
    for name in ('microcell_recovery_time_ns', 'intrinsic_time_jitter_ns',
                 'electronic_noise_rms_mv', 'sipm_crosstalk_probability',
                 'sipm_afterpulse_probability', 'dark_count_rate_hz', 'adc_bits'):
        assert getattr(instrument, name) == 0, name
    assert np.array_equal(sii.apply_exponential_microcell_recovery([0., 0., 1.], [1, 1, 1], 0.), np.ones(3))
    assert sii.mean_recovery_fraction(1e10, 1, 0.) == 1.


def test_hour_angle_points_west_after_transit():
    hour, dec, lat = .4, .2, .5
    assert sii.source_direction_enu(hour, dec, lat)[0] < 0
    epsilon = 1e-6
    derivative = (sii.source_direction_enu(hour+epsilon, dec, lat)
                  -sii.source_direction_enu(hour-epsilon, dec, lat))/(2*epsilon)
    u, _, _ = sii.celestial_tangent_axes_enu(hour, dec, lat)
    assert np.allclose(-derivative/np.cos(dec), u, atol=1e-9)
    layout=pd.DataFrame({'east_m':[0.,100.],'north_m':[0.,0.],'up_m':[0.,0.]})
    uv=sii.generate_uvw(layout,sii.Observation(hours_per_night=2.,segment_s=3600.))
    assert np.isclose(np.diff(uv.hour_angle_h)[0],86400./sii.SIDEREAL_DAY_S)


def test_invalid_exposure_and_invisible_source_fail():
    layout = pd.DataFrame({'east_m':[0.,100.], 'north_m':[0.,0.], 'up_m':[0.,0.]})
    with pytest.raises(ValueError, match='integer number'):
        sii.generate_uvw(layout, sii.Observation(hours_per_night=1, segment_s=1000))
    with pytest.raises(ValueError, match='horizon'):
        sii.generate_uvw(layout, sii.Observation(source_dec_deg=-85.))


def test_profiled_gain_matches_explicit_fit_and_gradient():
    observed = np.array([.8, .4, .2, .5])
    sigma = np.full(4, .03)
    uv = UvData(np.zeros(4), np.zeros(4), observed, sigma, 1/sigma**2,
                np.ones(4), 4, 4, 0, calibration_relative_sigma=.1)
    prediction = np.array([.72, .36, .17, .48])
    gain, uncertainty = profile_calibration_gain(uv, prediction)
    explicit = minimize_scalar(lambda g: .5*np.sum(((g*prediction-observed)/sigma)**2)
                                +.5*((g-1)/.1)**2)
    value, gradient = statistical_loss(uv, prediction-observed)
    assert np.isclose(gain, explicit.x)
    assert np.isclose(value, explicit.fun)
    assert 0 < uncertainty < .1
    direction = np.array([.1, -.3, .2, .5])
    epsilon = 1e-6
    finite_difference = (statistical_loss(uv, prediction+epsilon*direction-observed)[0]
        -statistical_loss(uv, prediction-epsilon*direction-observed)[0])/(2*epsilon)
    assert np.isclose(gradient @ direction, finite_difference, rtol=1e-7)


def test_uv_grouping_preserves_shared_calibration_prior():
    data = pd.DataFrame({'u_lambda':[1e8,1e8], 'v_lambda':[0.,0.],
        'visibility2_measured':[.5,.6], 'sigma_visibility2':[.1,.1],
        'calibration_relative_sigma':[.07,.07], 'calibration_id':['same','same']})
    uv = sii.prepare_reconstruction_uv(data)
    assert np.isclose(uv.sigma[0], .1/np.sqrt(2))
    assert uv.calibration_relative_sigma == .07
    data['visibility2_true'] = [-999., np.nan]
    assert sii.prepare_reconstruction_uv(data).calibration_relative_sigma == .07
    data.loc[1, 'calibration_id'] = 'different'
    with pytest.raises(ValueError, match='independent calibrations'):
        sii.prepare_reconstruction_uv(data)


def test_thermal_fields_siegert_and_photon_variance():
    counts, intensity = thermal_mode_counts(np.random.default_rng(44), 150_000,
                                            star_mean=6., background_mean=2., visibility2=.6)
    assert np.allclose(intensity.mean(axis=0), 1., atol=.006)
    assert np.isclose(np.cov(intensity.T)[0, 1], .6/16, rtol=.03)
    assert np.allclose(counts.mean(axis=0), 8., atol=.04)
    assert np.allclose(counts.var(axis=0), 8.+36/16, rtol=.025)
    assert np.isclose(np.cov(counts.T)[0, 1], 36*.6/16, rtol=.06)


def test_analytic_waveform_fine_grid_convergence():
    instrument = sii.Instrument.from_repository(ROOT)
    first, diagnostics = analytic_waveform_calibration(instrument, fine_dt_ns=.05)
    second, _ = analytic_waveform_calibration(instrument, fine_dt_ns=.025)
    assert diagnostics['photon_degeneracy'] < 1e-4
    assert np.allclose(first.peak_per_visibility2, second.peak_per_visibility2, atol=1e-9)
    assert np.isclose(sii.waveform_gls_weights(first, 1200)[1],
                      sii.waveform_gls_weights(second, 1200)[1], rtol=.01)


def test_common_gain_generates_cross_measurement_covariance():
    instrument = sii.Instrument.from_repository(ROOT)
    cal, _ = analytic_waveform_calibration(instrument)
    cal = replace(cal, response_relative_uncertainty=.1)
    rng = np.random.default_rng(55)
    truth = np.array([.7,.4])
    observations = np.array([sii.sample_waveform_gls_visibility2(truth, cal, 1e8, rng)[0]
                             for _ in range(2000)])
    assert np.isclose(np.cov(observations.T)[0, 1], .1**2*.7*.4, rtol=.1)


def test_profile_grid_matches_explicit_likelihood_and_interval():
    sigma=np.array([.1,.2,.1])
    uv=UvData(np.zeros(3),np.zeros(3),np.array([.8,.4,.2]),sigma,1/sigma**2,
              np.ones(3),3,3,0,calibration_relative_sigma=.05)
    predictions=np.array([[.8,.4,.2],[.7,.4,.15]])
    statistic,gain=profile_model_grid(uv,predictions)
    for i,model in enumerate(predictions):
        assert np.isclose(statistic[0,i],2*statistical_loss(uv,model-uv.visibility_abs2)[0],atol=1e-10)
        assert np.isclose(gain[0,i],profile_calibration_gain(uv,model)[0])
    parameter=np.linspace(-5,5,1001)
    interval=profile_grid_interval(parameter,parameter**2)
    assert np.isclose(interval['low'],-1.9599639845,atol=1e-5)
    assert np.isclose(interval['high'],1.9599639845,atol=1e-5)
    assert not interval['touches_boundary'] and not interval['disconnected']
    assert profile_grid_interval(parameter,parameter)['touches_boundary']


def test_single_peak_diagnostic_does_not_invent_second_source():
    from sii_reconstruction import _peak_diagnostic
    axis=np.linspace(-1.,1.,32)
    xx,yy=np.meshgrid(axis,axis)
    metrics=_peak_diagnostic(np.exp(-(xx*xx+yy*yy)/.1),axis,.2)
    assert len(metrics['peaks']) == 1
    assert 'two_peak_separation_mas' not in metrics
