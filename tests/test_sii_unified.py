from pathlib import Path
import json
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

import sii_unified as sii
from config_io import expand_component_config


def test_uvw_preserves_baseline_norm():
    baseline = np.array([421.0, -317.0, 12.0])
    uvw = sii.uvw_from_enu(baseline, 0.37, 0.22, 0.51)
    assert np.isclose(uvw @ uvw, baseline @ baseline, rtol=1e-12)
    assert np.allclose(
        sii.uvw_from_enu(-baseline, 0.37, 0.22, 0.51), -uvw)


def test_celestial_tangent_axes_are_orthonormal_and_right_handed():
    u_axis, v_axis, w_axis = sii.celestial_tangent_axes_enu(0.37, 0.22, 0.51)
    assert np.allclose([u_axis @ u_axis, v_axis @ v_axis, w_axis @ w_axis], 1.0)
    assert np.allclose([u_axis @ v_axis, u_axis @ w_axis, v_axis @ w_axis], 0.0)
    assert np.allclose(np.cross(u_axis, v_axis), w_axis)


def test_zero_adc_interface_preserves_analog_waveform():
    waveform = np.array([-1.2, 0.0, 3.4])
    assert np.array_equal(sii.digitize_adc(waveform, 0, 0.0), waveform)


def test_binary_visibility_is_normalized_and_power_is_even():
    u = np.array([0.0, 1.2e9, -2.4e9])
    v = np.array([0.0, -0.7e9, 0.4e9])
    visibility = sii.binary_visibility(u, v)
    reverse = sii.binary_visibility(-u, -v)
    assert np.isclose(abs(visibility[0]) ** 2, 1.0)
    assert np.allclose(abs(visibility) ** 2, abs(reverse) ** 2)


def test_exponential_recovery_at_one_time_constant():
    fractions = sii.apply_exponential_microcell_recovery(
        [0.0, 0.0, 10.0], [7, 7, 7], 10.0)
    assert np.allclose(fractions, [1.0, 0.0, 1.0 - np.exp(-1.0)])


def test_many_microcells_make_recovery_loss_tiny_at_reference_rate():
    fraction = sii.mean_recovery_fraction(200e6, 270_336, 10.0)
    assert fraction > 0.99999


def test_adc_quantization_clips_and_uses_finite_levels():
    digitized = sii.digitize_adc([-1000, -0.2, 0.2, 1000], 8, 200.0)
    assert digitized.min() >= -100.0
    assert digitized.max() < 100.0
    assert np.all(np.isfinite(digitized))


def test_empirical_charge_loader_is_positive_and_mean_one():
    factors = sii.load_empirical_charge_factors(
        ROOT / "configs" / "electronics" / "parameters"
        / "spe_charge_samples_measured.csv")
    assert len(factors) == 537
    assert np.all(factors > 0.0)
    assert np.isclose(np.mean(factors), 1.0, atol=1e-14)


def test_optical_timing_kernel_is_normalized_and_broadens_arrivals():
    mixture = sii.load_optical_timing_mixture(
        ROOT / "configs" / "optics" / "lact_1229_onaxis_timing_kernel.csv")
    assert np.isclose(np.sum(mixture["weights"]), 1.0)
    assert 0.5 < mixture["rms_spread_ns"] < 0.7
    assert abs(np.sum(mixture["weights"] * mixture["mean_delay_ns"])) < 1e-12


def test_optical_timing_transfer_and_residual_delay_are_physical():
    mixture = sii.load_optical_timing_mixture(
        ROOT / "configs" / "optics" / "lact_1229_onaxis_timing_kernel.csv")
    efficiency = sii.optical_timing_transfer_efficiency(mixture, 200e6)
    response = sii.residual_delay_response(mixture, [0.0, 0.2, 1.0], 200e6)
    assert 0.0 < efficiency < 1.0
    assert np.isclose(response[0], 1.0, atol=1e-12)
    assert 1.0 >= response[1] >= response[2]


def test_short_waveform_accepts_optical_timing_kernel():
    mixture = sii.load_optical_timing_mixture(
        ROOT / "configs" / "optics" / "lact_1229_onaxis_timing_kernel.csv")
    template_t = np.array([-2.0, 0.0, 2.0, 4.0])
    template_v = np.array([0.0, 1.0, 0.2, 0.0])
    record = sii.simulate_short_pair_waveforms(
        32.0, 20e6, 5e6, 0.5, sii.Instrument(), template_t, template_v,
        optical_timing_mixture=mixture, seed=7)
    assert np.isclose(record["optical_timing_rms_ns"],
                      mixture["rms_spread_ns"])


def test_two_microsecond_waveform_has_padding_for_full_spe_tail():
    template_t = np.array([-40.0, 0.0, 170.0])
    template_v = np.array([0.0, 1.0, 0.0])
    record = sii.simulate_short_pair_waveforms(
        2000.0, 100e6, 20e6, 0.5, sii.Instrument(),
        template_t, template_v, seed=17)
    assert len(record["sample_time_ns"]) == 1250
    assert record["simulated_padding_each_side_ns"] >= 170.0
    assert np.min(record["pe_times_a_ns"]) < 0.0
    assert np.max(record["pe_times_a_ns"]) > 2000.0


def test_repository_instrument_follows_main_configs():
    instrument = sii.Instrument.from_repository(ROOT)
    assert np.isclose(instrument.effective_area_m2, 24.576860)
    assert np.isclose(
        instrument.electronics_bandwidth_hz,
        instrument.adc_sample_rate_hz/2.0)
    assert np.isclose(instrument.detected_nsb_rate_hz, 0.573214e6, rtol=2e-4)
    assert instrument.microcells_per_pixel == 270_336
    assert np.isclose(instrument.adc_sample_rate_hz, 250e6)
    assert instrument.adc_bits == 0
    assert np.isclose(instrument.adc_full_scale_mv, 0.0)
    assert np.isclose(instrument.electronic_noise_rms_mv, 0.0)
    assert np.isclose(instrument.residual_timing_rms_ns, 0.0)
    assert Path(instrument.spe_template_path).is_file()
    assert Path(instrument.microcell_device_path).is_file()
    full = ROOT / "configs" / "optics" / "lact2_measured_single_pixel_400nm.csv"
    fallback = ROOT / "configs" / "optics" / "lact_1229_onaxis_timing_kernel.csv"
    if full.exists():
        provenance = json.loads(
            full.with_suffix(".provenance.json").read_text(encoding="utf-8"))
        expected_throughput = (
            provenance["central_pixel_effective_detection_area_m2"]
            / instrument.effective_area_m2 * 0.7836336)
        assert np.isclose(instrument.throughput, expected_throughput, rtol=.02)
        assert provenance["input_photons"] == 1_000_000
        assert provenance["illuminated_pixel_count"] == 1
        assert np.isclose(sum(provenance["pixel_signal_fractions"].values()), 1.0)
        assert np.isclose(pd.read_csv(full).weight.sum(), 1.0)
    else:
        assert np.isclose(instrument.throughput, 0.20, rtol=2e-4)
    assert Path(instrument.optical_timing_kernel_path) == (
        full if full.exists() else fallback).resolve()


def test_repository_response_sets_matched_effective_bandwidth():
    instrument = sii.Instrument.from_repository(ROOT)
    bandwidth_hz = sii.matched_effective_bandwidth_hz(instrument)
    assert np.isclose(bandwidth_hz, 110.910077e6, rtol=2e-6)
    calibrated = sii.with_matched_effective_bandwidth(instrument)
    assert calibrated.optical_timing_in_effective_bandwidth
    assert np.isclose(calibrated.electronics_bandwidth_hz, bandwidth_hz)


def test_matched_bandwidth_does_not_double_count_optical_timing():
    instrument = sii.with_matched_effective_bandwidth(
        sii.Instrument.from_repository(ROOT))
    layout = pd.DataFrame({
        "name": ["A", "B"],
        "east_m": [0.0, 100.0],
        "north_m": [0.0, 0.0],
        "up_m": [0.0, 0.0],
    })
    observation = sii.Observation(hours_per_night=1.0/3.0, segment_s=1200.0)
    uvw = sii.generate_uvw(layout, observation, instrument)
    _, metadata = sii.simulate_uv_observation(
        uvw, sii.BinarySource(), observation, instrument, seed=17)
    assert metadata["optical_timing_efficiency"] == 1.0


def test_uv_measurement_uses_segment_averaged_visibility():
    layout = pd.DataFrame({
        "name": ["A", "B"],
        "east_m": [0.0, 1000.0],
        "north_m": [0.0, 0.0],
        "up_m": [0.0, 0.0],
    })
    observation = sii.Observation(
        hours_per_night=1.0/3.0, segment_s=1200.0,
        visibility_subsamples_per_segment=9)
    dense_observation = sii.Observation(
        hours_per_night=1.0/3.0, segment_s=1200.0,
        visibility_subsamples_per_segment=1001)
    source = sii.BinarySource(separation_mas=1.0)
    instrument = sii.Instrument(detected_nsb_rate_hz=0.0)
    uvw = sii.generate_uvw(layout, observation, instrument)
    averaged = sii.segment_averaged_visibility2(
        uvw, source, observation, instrument)
    dense = sii.segment_averaged_visibility2(
        uvw, source, dense_observation, instrument)
    center = np.abs(sii.binary_visibility(
        uvw.u_lambda, uvw.v_lambda, source))**2
    assert np.allclose(averaged, dense, atol=5e-4)
    assert np.max(np.abs(averaged-center)) > 1e-3

    measurements, metadata = sii.simulate_uv_observation(
        uvw, source, observation, instrument, seed=17,
        electronics_case="ideal")
    assert np.allclose(measurements.visibility2_true, averaged)
    assert np.allclose(
        measurements.segment_time_smearing_delta, averaged-center)
    assert measurements.baseline_integration_s.eq(1200.0).all()
    assert metadata["visibility_subsamples_per_segment"] == 9


def test_single_pixel_electronics_entry_uses_recovery_config():
    config = ROOT / "configs" / "sii" / "single_pixel_electronics.cfg"
    values, component_paths = expand_component_config(config)
    assert component_paths["electronics"].is_file()
    assert values["electronics.n_pixels"] == "1"
    assert values["electronics.microcell.recovery_enabled"] == "true"
    assert values["electronics.microcell.recovery_time_ns"] == "10.0"


def test_dark_counts_are_sampled_per_segment_and_telescope(monkeypatch):
    inner = np.random.default_rng(17)
    poisson_sizes = []

    class RecordingGenerator:
        def poisson(self, lam, size=None):
            poisson_sizes.append(size)
            return inner.poisson(lam, size)

        def __getattr__(self, name):
            return getattr(inner, name)

    recorder = RecordingGenerator()
    monkeypatch.setattr(sii.np.random, "default_rng", lambda seed: recorder)
    layout = pd.DataFrame({
        "telescope_id": [1, 2], "name": ["A", "B"],
        "east_m": [0.0, 100.0], "north_m": [0.0, 0.0],
        "up_m": [0.0, 0.0],
    })
    observation = sii.Observation(
        hours_per_night=2.0/3600.0, segment_s=1.0)
    instrument = sii.Instrument(
        detected_nsb_rate_hz=0.0, dark_count_rate_hz=1.0e3)
    uvw = sii.generate_uvw(layout, observation, instrument)
    sii.simulate_uv_observation(
        uvw, sii.BinarySource(), observation, instrument, seed=17)
    assert (2, 2) in poisson_sizes


def test_short_pair_rate_and_long_exposure_snr_use_same_coherence_area():
    instrument = sii.Instrument(
        optical_width_nm=2.0, electronics_bandwidth_hz=110.0e6,
        detected_nsb_rate_hz=70.0e6, polarization_factor=0.5,
        spectral_shape_factor=0.842)
    integration_s = 1200.0
    star = sii.detected_star_rate_hz(2.0, instrument)
    total = star + instrument.detected_nsb_rate_hz
    pair_rate_at_unit_visibility = sii.hbt_correlated_pair_rate_hz(
        star, star, instrument.coherence_area_s, 1.0)
    snr_implied_by_short_pairs = (
        pair_rate_at_unit_visibility/total
        * np.sqrt(2.0*instrument.electronics_bandwidth_hz*integration_s)
        / instrument.excess_noise_factor**2)
    assert np.isclose(
        sii.unit_visibility_snr(2.0, integration_s, instrument),
        snr_implied_by_short_pairs)


def test_one_sided_nyquist_bandwidth_counts_all_real_samples():
    """B=fs/2 时，2BT 必须等于时域独立样点数 fs*T。"""
    sample_rate_hz = 250.0e6
    integration_s = 10.0
    bandwidth_hz = sample_rate_hz/2.0
    assert np.isclose(
        2.0*bandwidth_hz*integration_s,
        sample_rate_hz*integration_s)


def test_full_optical_response_uses_brightest_pixel_and_weight(tmp_path):
    hits = pd.DataFrame({
        "mirror_id": [1, 2, 3], "time_ns": [80.0, 82.0, 90.0],
        "u_m": [0.0, 0.001, 0.02], "v_m": [0.0, 0.0, 0.0],
        "hit_camera": [1, 1, 1], "accepted": [1, 1, 1],
        "pixel_id": [7, 7, 8], "weight": [1.0, 1.0, 1.0],
        "relative_efficiency": [0.5, 0.4, 0.1],
    })
    input_csv = tmp_path / "hits.csv"
    output_csv = tmp_path / "kernel.csv"
    provenance = tmp_path / "kernel.provenance.json"
    hits.to_csv(input_csv, index=False)
    subprocess.run([
        sys.executable, str(ROOT/"tools"/"derive_full_optical_response.py"),
        str(input_csv), str(output_csv), "--input-photons", "10",
        "--sampling-radius-m", "1", "--provenance-json", str(provenance),
    ], check=True)
    kernel = pd.read_csv(output_csv)
    info = json.loads(provenance.read_text(encoding="utf-8"))
    assert set(kernel.mirror_id) == {1, 2}
    assert np.isclose(kernel.weight.sum(), 1.0)
    assert info["central_pixel_id"] == 7
    assert np.isclose(info["central_pixel_detection_probability"], 0.09)
    assert np.isclose(info["central_pixel_effective_detection_area_m2"],
                      np.pi*0.09)


def test_hbt_pair_rate_and_main_compatible_primary_stream(tmp_path):
    instrument = sii.Instrument(optical_width_nm=1e-5)
    expected = sii.hbt_correlated_pair_rate_hz(
        1e6, 1e6, instrument.coherence_area_s, 0.5)
    hits, metadata = sii.simulate_hbt_primary_pe(
        np.random.default_rng(12), 1e6, 1e6, 0.0, 0.5,
        instrument, padding_ns=0.0)
    assert np.isclose(metadata["hbt_pair_rate_hz"], expected)
    paired = hits[hits.hbt_pair_id >= 0]
    assert paired.groupby("hbt_pair_id").telescope_id.nunique().eq(2).all()
    output = sii.write_main_primary_pe_csv(hits, tmp_path / "primary.csv")
    header = output.read_text(encoding="utf-8").splitlines()[0]
    assert header == ("event_id,telescope_id,pixel_id,time_ns,sensor_x_m,"
                      "sensor_y_m,primary_pe,wavelength_nm,origin")


def test_fast_waveform_renderer_and_fft_cross_correlation():
    instrument = sii.Instrument(
        charge_samples_path=None, electronic_noise_rms_mv=0.0,
        adc_sample_rate_hz=4e9, electronics_bandwidth_hz=2e9)
    template = sii.make_fast_spe_template()
    times = np.array([10.0, 20.0, 51.0, 77.0])
    left = sii.render_pe_waveform(
        np.random.default_rng(1), times, 100.0, instrument,
        template=template, electronic_noise_rms_mv=0.0)
    right = sii.render_pe_waveform(
        np.random.default_rng(2), times, 100.0, instrument,
        template=template, electronic_noise_rms_mv=0.0)
    lags, correlation = sii.waveform_cross_correlation(
        left["analog_mv"], right["analog_mv"], 0.25, 5.0)
    assert np.isclose(lags[np.argmax(correlation)], 0.0)
    assert np.isclose(np.max(correlation), 1.0)
    assert len(left["sample_time_ns"]) == 400


def test_waveform_gls_recovers_injected_visibility_and_time_scaling():
    rng = np.random.default_rng(14)
    lags = np.arange(-2.0, 3.0)
    peak = np.exp(-0.5*(lags/0.8)**2)*2.0e-4
    covariance = 2.0e-4*np.fromfunction(
        lambda i, j: 0.35**np.abs(i-j), (len(lags), len(lags)))
    null = rng.multivariate_normal(np.zeros(len(lags)), covariance, 800)
    signal = rng.multivariate_normal(0.7*peak, covariance, 800)
    calibration = sii.calibrate_waveform_gls(
        lags, null, signal, 1.0, 10e6, 2e6,
        calibration_visibility2=0.7, covariance_shrinkage=0.05)

    noiseless = calibration.null_mean+0.4*calibration.peak_per_visibility2
    estimate, sigma_1s = sii.estimate_visibility2_gls(noiseless, calibration)
    _, sigma_4s = sii.waveform_gls_weights(calibration, 4.0)
    assert np.isclose(estimate, 0.4)
    assert np.isclose(sigma_4s, sigma_1s/2.0)


def test_ellipse_and_transit_visibility_are_normalized_and_hermitian():
    u = np.array([0.0, 1.1e8, -0.7e8])
    v = np.array([0.0, -0.4e8, 1.3e8])
    for case, source in (
            ("ellipse", sii.EllipseSource()),
            ("transit", sii.TransitSource())):
        visibility = sii.source_visibility(u, v, source, case)
        mirrored = sii.source_visibility(-u, -v, source, case)
        assert np.isclose(visibility[0], 1.0)
        assert np.allclose(mirrored, np.conj(visibility))


def test_transit_requires_planet_inside_stellar_disk():
    source = sii.TransitSource(
        stellar_diameter_mas=0.2, planet_diameter_mas=0.1,
        planet_east_offset_mas=0.08, planet_north_offset_mas=0.0)
    with pytest.raises(ValueError, match="completely inside"):
        sii.transit_visibility(np.array([0.0]), np.array([0.0]), source)


def test_waveform_gls_long_uv_path_uses_calibrated_peak_covariance():
    instrument = sii.Instrument(detected_nsb_rate_hz=3.0e6)
    star_rate = sii.detected_star_rate_hz(2.0, instrument)
    calibration = sii.WaveformGLSCalibration(
        lags_ns=np.array([-1.0, 0.0, 1.0]),
        null_mean=np.array([0.01, -0.02, 0.01]),
        peak_per_visibility2=np.array([0.2, 1.0, 0.2]),
        covariance_per_block=np.diag([0.5, 1.0, 0.5]),
        block_duration_s=1.0, star_rate_hz=star_rate,
        background_rate_hz=instrument.detected_nsb_rate_hz,
        calibration_visibility2=1.0, hbt_pair_rate_scale=10.0,
        null_records=20, signal_records=20, covariance_shrinkage=0.05,
        instrument_signature=sii.waveform_instrument_signature(instrument))
    layout = pd.DataFrame({
        "name": ["A", "B"], "east_m": [0.0, 100.0],
        "north_m": [0.0, 0.0], "up_m": [0.0, 0.0]})
    observation = sii.Observation(
        hours_per_night=1.0/3.0, segment_s=1200.0)
    uvw = sii.generate_uvw(layout, observation, instrument)
    measured, metadata = sii.simulate_uv_observation(
        uvw, sii.BinarySource(), observation, instrument, seed=7,
        estimator="waveform_gls", waveform_calibration=calibration)
    _, expected_sigma = sii.waveform_gls_weights(calibration, 1200.0)
    assert metadata["estimator"] == "waveform_gls"
    assert np.allclose(measured.sigma_visibility2, expected_sigma)
    assert measured.visibility2_measured.notna().all()


def test_hbt_calibration_scale_does_not_change_reported_physical_rate():
    instrument = sii.Instrument(optical_width_nm=0.01)
    rng = np.random.default_rng(7)
    _, reference = sii.simulate_hbt_primary_pe(
        rng, 10.0, 1e6, 0.0, 0.5, instrument,
        padding_ns=0.0, hbt_pair_rate_scale=1.0)
    _, enhanced = sii.simulate_hbt_primary_pe(
        rng, 10.0, 1e6, 0.0, 0.5, instrument,
        padding_ns=0.0, hbt_pair_rate_scale=10.0)
    assert enhanced["hbt_pair_rate_hz"] == reference["hbt_pair_rate_hz"]
    assert enhanced["injected_hbt_pair_rate_hz"] == 10*reference["hbt_pair_rate_hz"]


def test_complete_pipeline_without_reconstruction():
    layout = np.array([
        (1, "A", 0.0, 0.0, 0.0),
        (2, "B", 100.0, 0.0, 0.0),
        (3, "C", 0.0, 100.0, 0.0),
    ], dtype=[("telescope_id", "i4"), ("name", "U2"),
              ("east_m", "f8"), ("north_m", "f8"), ("up_m", "f8")])
    import pandas as pd
    observation = sii.Observation(hours_per_night=1.0, segment_s=1200.0)
    result = sii.run_sii_pipeline(
        pd.DataFrame(layout), observation=observation,
        instrument=sii.Instrument(), do_reconstruction=False, seed=3,
        electronics_case="ideal")
    assert len(result.uvw) == 3 * 3  # 3 baselines × 3 time segments
    assert len(result.measurements) == len(result.uvw)
    assert result.measurements.visibility2_measured.notna().all()
    assert result.reconstruction is None


def test_reconstruction_import_preserves_loaded_pyplot_backend():
    import importlib
    import matplotlib
    import matplotlib.pyplot  # noqa: F401 - establishes the active backend.

    before = matplotlib.get_backend()
    sys.modules.pop("sii_reconstruction", None)
    importlib.import_module("sii_reconstruction")
    assert matplotlib.get_backend() == before


def test_exact_averaged_power_and_gradient():
    from sii_reconstruction import (UvData, power_sampling_kernel, power_from_image,
                                    _power_gradient, MAS_TO_RAD)
    rng = np.random.default_rng(45)
    u, v = rng.normal(0, 1e9, (2, 12))
    group = np.repeat(np.arange(4), 3)
    fraction = np.tile([0.2, 0.3, 0.5], 4)
    uv = UvData(u[:4], v[:4], np.zeros(4), np.ones(4), np.ones(4),
                np.ones(4), 4, 4, 0, sampling=(u, v, group, fraction))
    image = rng.uniform(size=(12, 12))
    image /= image.sum()
    theta = np.linspace(-.35, .35, 12)*MAS_TO_RAD
    yy, xx = np.meshgrid(theta, theta, indexing="ij")
    field = np.exp(-2j*np.pi*(u[:, None]*xx.ravel()+v[:, None]*yy.ravel())) @ image.ravel()
    direct = np.bincount(group, weights=fraction*abs(field)**2)
    kernel = power_sampling_kernel(uv, 12, .7)
    assert np.allclose(power_from_image(kernel, image), direct, atol=1e-13)
    influence = rng.normal(size=4)
    gradient = _power_gradient(kernel, influence, image)
    direction = rng.normal(size=image.shape)
    epsilon = 1e-7
    finite_difference = (power_from_image(kernel, image+epsilon*direction)
                         - power_from_image(kernel, image-epsilon*direction)) @ influence/(2*epsilon)
    assert np.isclose(finite_difference, np.sum(gradient*direction), rtol=1e-7)


def test_absolute_likelihood_sigma_scaling_and_covariance():
    from dataclasses import replace
    from sii_reconstruction import UvData, statistical_loss
    uv = UvData(np.arange(3.), np.zeros(3), np.zeros(3), np.array([.1, .2, .3]),
                np.ones(3), np.ones(3), 3, 3, 0)
    residual = np.array([.02, -.3, .1])
    loss, gradient = statistical_loss(uv, residual)
    scaled_loss, scaled_gradient = statistical_loss(replace(uv, sigma=10*uv.sigma), residual)
    assert np.isclose(scaled_loss, loss/100)
    assert np.allclose(scaled_gradient, gradient/100)
    covariance = np.diag(uv.sigma**2)
    covariance[0, 1] = covariance[1, 0] = .005
    full_loss, full_gradient = statistical_loss(replace(uv, covariance=covariance), residual)
    assert np.allclose(full_gradient, np.linalg.solve(covariance, residual))
    assert np.isclose(full_loss, residual @ full_gradient/2)
    with pytest.raises(ValueError, match="diagonal"):
        statistical_loss(replace(uv, sigma=uv.sigma*1e-8,
                                 covariance=2*covariance*1e-16), residual)
    with pytest.raises(ValueError, match="positive sigma"):
        statistical_loss(replace(uv, sigma=np.zeros(3)), residual)


def test_uv_grouping_keeps_forward_weights_without_artificial_dc():
    data = pd.DataFrame({"u_lambda": [1e8, 1.1e8], "v_lambda": [0., 0.],
                         "visibility2_measured": [.4, .8], "sigma_visibility2": [.1, .2],
                         "uv_samples_u": [(8e7, 1.2e8), (1e8, 1.2e8)],
                         "uv_samples_v": [(0., 0.), (0., 0.)]})
    uv = sii.prepare_reconstruction_uv(data, cell_mlambda=500)
    assert len(uv.u_lambda) == 1
    assert np.isclose(uv.visibility_abs2[0], .48)
    assert np.isclose(uv.sigma[0], np.sqrt(1/125))
    assert np.allclose(uv.sampling[3], [.4, .4, .1, .1])
    assert len(sii.prepare_reconstruction_uv(data, cell_mlambda=None).u_lambda) == 2
    # 真值列即使被污染，也不能改变重建输入；只允许观测和采样几何进入。
    data["visibility2_true"] = [np.nan, -999.0]
    poisoned = sii.prepare_reconstruction_uv(data, cell_mlambda=500)
    assert np.array_equal(poisoned.visibility_abs2, uv.visibility_abs2)
    assert all(np.array_equal(a, b) for a, b in zip(poisoned.sampling, uv.sampling))


def test_statistical_reconstruction_cv_and_returned_prediction():
    from sii_reconstruction import (UvData, power_sampling_kernel, power_from_image,
                                    reconstruct_uv_data)
    rng = np.random.default_rng(91)
    u, v = rng.normal(0, 3e8, (2, 24))
    uv = UvData(u, v, np.zeros(24), np.full(24, .03), np.full(24, 1/.03**2),
                np.ones(24), 24, 24, 0)
    image = np.zeros((12, 12))
    image[5:7, 5:7] = .25
    kernel = power_sampling_kernel(uv, 12, .7)
    uv.visibility_abs2 = power_from_image(kernel, image)+rng.normal(0, .03, 24)
    result = reconstruct_uv_data(uv, grid_size=12, fov_mas=.7, starts=1,
        max_iter=100, smoothness="cv", smoothness_candidates=(0., .01), seed=5)
    assert len(result.metrics["smoothness_selection"]) == 2
    assert np.allclose(result.predicted_visibility_abs2, power_from_image(kernel, result.image))
    residual = (result.predicted_visibility_abs2-uv.visibility_abs2)/uv.sigma
    assert np.isclose(result.metrics["chi2"], np.sum(residual**2))
    assert np.isclose(result.image.sum(), 1)


def test_flux_parameterization_gradient_and_zero_vector_guard(monkeypatch):
    import sii_reconstruction as reconstruction
    rng = np.random.default_rng(81)
    uv = reconstruction.UvData(rng.normal(0, 1e9, 20), rng.normal(0, 1e9, 20),
        np.full(20, .2), np.full(20, .1), np.full(20, 100.), np.ones(20), 20, 20, 0)
    original = reconstruction.minimize

    def checked(fun, initial, **kwargs):
        direction = initial*rng.normal(size=len(initial))
        epsilon = 1e-4
        value, gradient = fun(initial)
        finite_difference = (fun(initial+epsilon*direction)[0]
                             - fun(initial-epsilon*direction)[0])/(2*epsilon)
        assert np.isclose(finite_difference, gradient @ direction, rtol=1e-5, atol=1e-7)
        invalid, escape = fun(np.zeros_like(initial))
        assert np.isinf(invalid) and np.all(np.isfinite(escape))
        return original(fun, initial, **kwargs)

    monkeypatch.setattr(reconstruction, 'minimize', checked)
    fit = reconstruction.reconstruct_uv_data(uv, grid_size=12, fov_mas=.7,
        starts=1, max_iter=200, smoothness=.001, parameterization='flux')
    assert np.all(fit.image >= 0) and np.isclose(fit.image.sum(), 1)
    assert fit.metrics['simplex_stationarity_gap'] >= -1e-8
    assert fit.metrics['parameterization'] == 'flux'
