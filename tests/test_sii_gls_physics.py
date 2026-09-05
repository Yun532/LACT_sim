from dataclasses import replace
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'python'))
import sii_unified as sii


def test_shared_passband_rates_and_coherence(tmp_path):
    narrow = sii.Instrument.from_repository(ROOT, optical_width_nm=1.)
    wide = sii.Instrument.from_repository(ROOT, optical_width_nm=2.)
    assert np.isclose(wide.detected_nsb_rate_hz/narrow.detected_nsb_rate_hz, 2, rtol=.03)
    assert np.isclose(sii.detected_star_rate_hz(2, wide)/sii.detected_star_rate_hz(2, narrow), 2, rtol=.03)
    assert np.isclose(wide.coherence_area_s/narrow.coherence_area_s, .5, rtol=.01)
    triangle = tmp_path/'triangle.csv'
    triangle.write_text('wavelength_nm,transmission\n399,0\n400,1\n401,0\n')
    shaped = sii.Instrument.from_repository(ROOT, sii_bandpass_path=str(triangle))
    # 三角形透过率使两种光子率减半，相干面积增大为原来的4/3。
    assert np.isclose(shaped.detected_nsb_rate_hz/wide.detected_nsb_rate_hz, .5, rtol=.03)
    assert np.isclose(sii.detected_star_rate_hz(2, shaped)/sii.detected_star_rate_hz(2, wide), .5, rtol=.03)
    assert np.isclose(shaped.coherence_area_s/wide.coherence_area_s, 4/3, rtol=.01)


@pytest.mark.parametrize('duration', [32., 33., 34., 35.])
def test_continuous_spe_sampling_matches_direct_sum(duration):
    instrument = sii.Instrument(adc_sample_rate_hz=250e6, electronics_bandwidth_hz=125e6)
    template = sii.make_fast_spe_template(rise_ns=.4, fall_ns=2., support_ns=10., step_ns=.02)
    times = np.array([-5., .1, 8.1, 11.9, 25.2, 33.1])
    record = sii.render_pe_waveform(np.random.default_rng(31), times, duration,
                                    instrument, template=template)
    expected = sii.convolve_pe_times(times, record['recovery_fraction']*record['charge_factor'],
                                     record['sample_time_ns'], *template)
    assert np.allclose(record['analog_mv'], expected, atol=1e-13)
    assert len(record['adc_mv']) == len(record['sample_time_ns'])
    assert np.all(record['sample_time_ns'] < duration)
    outputs = [sii.render_pe_waveform(np.random.default_rng(1), [t], duration,
                                     instrument, template=template)['analog_mv'] for t in (8.1, 11.9)]
    assert not np.allclose(*outputs)


@pytest.mark.parametrize('field,value', [('sipm_crosstalk_probability', .1),
    ('sipm_afterpulse_probability', .1), ('electronics_bandwidth_hz', 1e6)])
def test_unsupported_waveform_settings_fail(field, value):
    instrument = replace(sii.Instrument(), **{field: value})
    with pytest.raises((ValueError, NotImplementedError)):
        sii.render_pe_waveform(np.random.default_rng(1), [1.], 32., instrument,
                               template=sii.make_fast_spe_template())


def test_calibration_signature_tracks_contents_and_actual_template(tmp_path):
    path = tmp_path/'spe.csv'
    path.write_text('time_ns,amplitude_mv\n0,0\n1,1\n2,0\n')
    instrument = sii.Instrument(spe_template_path=str(path))
    first = sii.waveform_instrument_signature(instrument)
    copy = tmp_path/'copy.csv'
    copy.write_bytes(path.read_bytes())
    assert sii.waveform_instrument_signature(replace(instrument, spe_template_path=str(copy))) == first
    path.write_text('time_ns,amplitude_mv\n0,0\n1,2\n2,0\n')
    assert sii.waveform_instrument_signature(instrument) != first
    assert sii.waveform_instrument_signature(instrument, sii.make_fast_spe_template()) != first


def test_calibration_splits_and_rejects_different_instrument():
    instrument = sii.Instrument.from_repository(ROOT)
    cal, diagnostics = sii.simulate_waveform_gls_calibration(instrument,
        block_duration_ns=2000, max_lag_ns=8, null_records=32, signal_records=32, seed=81)
    assert diagnostics['null_scale_start'] > diagnostics['null_training_records']
    assert diagnostics['signal_test_start'] > diagnostics['signal_training_records']
    assert len(diagnostics['shrinkage_selection']) == 4
    assert cal.response_relative_uncertainty > 0
    assert cal.sigma_relative_uncertainty > 0
    observation = sii.Observation(hours_per_night=1/3, segment_s=1200)
    layout = pd.DataFrame({'east_m':[0.,100.], 'north_m':[0.,0.], 'up_m':[0.,0.]})
    uvw = sii.generate_uvw(layout, observation, instrument)
    frame, _ = sii.simulate_uv_observation(uvw, sii.BinarySource(), observation, instrument,
        estimator='waveform_gls', waveform_calibration=cal)
    assert frame.sigma_visibility2_calibration.notna().all()
    for changed in (replace(instrument, adc_bits=2, adc_full_scale_mv=.001),
                    replace(instrument, electronic_noise_rms_mv=100),
                    replace(instrument, optical_timing_kernel_path=None)):
        with pytest.raises(ValueError, match='signature mismatch'):
            sii.simulate_uv_observation(uvw, sii.BinarySource(), observation, changed,
                estimator='waveform_gls', waveform_calibration=cal)


def test_calibration_nonintegral_duration_fails_before_simulation():
    instrument = sii.Instrument.from_repository(ROOT)
    with pytest.raises(ValueError, match='integer number'):
        sii.simulate_waveform_gls_calibration(instrument, block_duration_ns=33.)
