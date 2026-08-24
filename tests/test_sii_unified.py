from pathlib import Path
import sys

import numpy as np

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
