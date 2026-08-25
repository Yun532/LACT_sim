from pathlib import Path
import json
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

import sii_unified as sii


def test_uvw_preserves_baseline_norm():
    baseline = np.array([421.0, -317.0, 12.0])
    uvw = sii.uvw_from_enu(baseline, 0.37, 0.22, 0.51)
    assert np.isclose(uvw @ uvw, baseline @ baseline, rtol=1e-12)


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
    assert np.isclose(instrument.detected_nsb_rate_hz, 70.527e6, rtol=2e-4)
    assert instrument.microcells_per_pixel == 270_336
    assert np.isclose(instrument.adc_sample_rate_hz, 250e6)
    assert Path(instrument.spe_template_path).is_file()
    assert Path(instrument.microcell_device_path).is_file()
    full = ROOT / "configs" / "optics" / "lact2_measured_full_response_400nm.csv"
    fallback = ROOT / "configs" / "optics" / "lact_1229_onaxis_timing_kernel.csv"
    if full.exists():
        provenance = json.loads(
            full.with_suffix(".provenance.json").read_text(encoding="utf-8"))
        expected_throughput = (
            provenance["central_pixel_effective_detection_area_m2"]
            / instrument.effective_area_m2 * 0.7836336)
        assert np.isclose(instrument.throughput, expected_throughput)
        assert provenance["input_photons"] == 1_000_000
        assert provenance["illuminated_pixel_count"] == 4
        assert np.isclose(sum(provenance["pixel_signal_fractions"].values()), 1.0)
        assert np.isclose(pd.read_csv(full).weight.sum(), 1.0)
    else:
        assert np.isclose(instrument.throughput, 0.20, rtol=2e-4)
    assert Path(instrument.optical_timing_kernel_path) == (
        full if full.exists() else fallback).resolve()


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
        adc_sample_rate_hz=4e9)
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
