#!/usr/bin/env python3
"""Unified stellar-intensity-interferometry simulation primitives.

The module deliberately separates three scales:

1. optical coherence and source visibility;
2. short digitized electronics waveforms used to validate delays, SPE shaping,
   ADC effects, noise, and SiPM recovery;
3. long-exposure |V|^2 estimators used for an observing-night simulation.

It is not scientifically meaningful to allocate a full night at 625 MS/s in
memory.  The short waveform calibrates response effects, while the sufficient
long-exposure statistic is simulated with its explicit uncertainty.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math

import numpy as np

try:
    from scipy.special import j1
except ImportError:  # pragma: no cover
    j1 = None


C_M_S = 299_792_458.0
H_J_S = 6.626_070_15e-34
JY_W_M2_HZ = 1.0e-26
MAS_TO_RAD = math.pi / (180.0 * 3600.0 * 1000.0)


@dataclass(frozen=True)
class Instrument:
    wavelength_nm: float = 400.0
    optical_width_nm: float = 2.0
    effective_area_m2: float = 24.576860
    throughput: float = 0.20
    electronics_bandwidth_hz: float = 200.0e6
    adc_sample_rate_hz: float = 625.0e6
    adc_bits: int = 8
    adc_full_scale_mv: float = 200.0
    detected_nsb_rate_hz: float = 70.527e6
    electronic_noise_rms_mv: float = 0.35
    excess_noise_factor: float = 1.016142
    polarization_factor: float = 0.5
    spectral_shape_factor: float = 0.842
    microcells_per_pixel: int = 270_336
    microcell_recovery_time_ns: float = 10.0

    @property
    def wavelength_m(self) -> float:
        return self.wavelength_nm * 1.0e-9

    @property
    def optical_width_m(self) -> float:
        return self.optical_width_nm * 1.0e-9

    @property
    def optical_bandwidth_hz(self) -> float:
        return C_M_S * self.optical_width_m / self.wavelength_m**2

    @property
    def sample_width_ns(self) -> float:
        return 1.0e9 / self.adc_sample_rate_hz

    @property
    def coherence_area_s(self) -> float:
        return (self.polarization_factor * self.spectral_shape_factor
                * self.wavelength_m**2
                / (C_M_S * self.optical_width_m))


@dataclass(frozen=True)
class BinarySource:
    separation_mas: float = 0.20
    position_angle_deg: float = 35.0
    flux_ratio_secondary_to_primary: float = 0.55
    primary_diameter_mas: float = 0.060
    secondary_diameter_mas: float = 0.040


def source_direction_enu(hour_angle_rad: float, dec_rad: float,
                         lat_rad: float) -> np.ndarray:
    """Topocentric source unit vector in East, North, Up coordinates."""
    vector = np.array([
        math.cos(dec_rad) * math.sin(hour_angle_rad),
        math.sin(dec_rad) * math.cos(lat_rad)
        - math.cos(dec_rad) * math.cos(hour_angle_rad) * math.sin(lat_rad),
        math.sin(dec_rad) * math.sin(lat_rad)
        + math.cos(dec_rad) * math.cos(hour_angle_rad) * math.cos(lat_rad),
    ])
    return vector / np.linalg.norm(vector)


def uvw_from_enu(baseline_enu_m, hour_angle_rad: float, dec_rad: float,
                 lat_rad: float) -> np.ndarray:
    """Project a local ENU baseline onto east/north/source sky axes."""
    source = source_direction_enu(hour_angle_rad, dec_rad, lat_rad)
    u_axis = np.cross(np.array([0.0, 0.0, 1.0]), source)
    if np.linalg.norm(u_axis) < 1.0e-14:
        u_axis = np.array([1.0, 0.0, 0.0])
    else:
        u_axis /= np.linalg.norm(u_axis)
    v_axis = np.cross(source, u_axis)
    baseline = np.asarray(baseline_enu_m, dtype=float)
    return np.array([baseline @ u_axis, baseline @ v_axis,
                     baseline @ source])


def uniform_disk_visibility(q_lambda, diameter_mas: float) -> np.ndarray:
    """Real visibility of a circular uniform disk."""
    if j1 is None:
        raise RuntimeError("SciPy is required for uniform-disk visibility")
    x = np.pi * diameter_mas * MAS_TO_RAD * np.asarray(q_lambda, float)
    output = np.ones_like(x)
    nonzero = np.abs(x) > 1.0e-12
    output[nonzero] = 2.0 * j1(x[nonzero]) / x[nonzero]
    return output


def binary_visibility(u_lambda, v_lambda,
                      source: BinarySource = BinarySource()) -> np.ndarray:
    """Normalized complex visibility of two uniform-disk components."""
    u = np.asarray(u_lambda, dtype=float)
    v = np.asarray(v_lambda, dtype=float)
    q = np.hypot(u, v)
    pa = math.radians(source.position_angle_deg)
    dx = 0.5 * source.separation_mas * math.sin(pa) * MAS_TO_RAD
    dy = 0.5 * source.separation_mas * math.cos(pa) * MAS_TO_RAD
    first = (uniform_disk_visibility(q, source.primary_diameter_mas)
             * np.exp(-2j * np.pi * (-u * dx - v * dy)))
    second = (uniform_disk_visibility(q, source.secondary_diameter_mas)
              * np.exp(-2j * np.pi * (u * dx + v * dy)))
    ratio = source.flux_ratio_secondary_to_primary
    return (first + ratio * second) / (1.0 + ratio)


def ab_photon_spectral_density(magnitude: float, wavelength_m: float) -> float:
    """AB-source photon spectral density in photons s^-1 m^-2 Hz^-1."""
    f_nu = 3631.0 * JY_W_M2_HZ * 10.0 ** (-0.4 * magnitude)
    photon_energy = H_J_S * C_M_S / wavelength_m
    return f_nu / photon_energy


def detected_star_rate_hz(magnitude: float,
                          instrument: Instrument = Instrument()) -> float:
    return (instrument.effective_area_m2 * instrument.throughput
            * ab_photon_spectral_density(magnitude, instrument.wavelength_m)
            * instrument.optical_bandwidth_hz)


def unit_visibility_snr(magnitude: float, integration_s: float,
                        instrument: Instrument = Instrument(),
                        spectral_channels: int = 1,
                        nsb_rate_hz: float | None = None) -> float:
    """Shot-noise approximation for an unresolved equal-telescope pair."""
    nsb = (instrument.detected_nsb_rate_hz if nsb_rate_hz is None
           else nsb_rate_hz)
    star = detected_star_rate_hz(magnitude, instrument)
    total = star + nsb
    one_channel = (
        star**2 / (total * instrument.optical_bandwidth_hz)
        * math.sqrt(instrument.electronics_bandwidth_hz * integration_s / 2.0)
        / instrument.excess_noise_factor
    )
    return one_channel * math.sqrt(spectral_channels)


def mean_recovery_fraction(rate_hz: float, microcells: int,
                           recovery_time_ns: float) -> float:
    """Mean charge fraction for uniform random illumination of many cells.

    Each cell sees a Poisson rate ``rate_hz/microcells``.  Averaging
    ``1-exp(-dt/tau)`` over the exponential inter-arrival distribution gives
    exactly ``1/(1 + rate_per_cell*tau)``.
    """
    if rate_hz < 0 or microcells <= 0 or recovery_time_ns <= 0:
        raise ValueError("rate, microcell count, and recovery time are invalid")
    occupancy = rate_hz * recovery_time_ns * 1.0e-9 / microcells
    return 1.0 / (1.0 + occupancy)


def apply_exponential_microcell_recovery(times_ns, cell_ids,
                                         recovery_time_ns: float) -> np.ndarray:
    """Return charge fraction for a time-ordered sequence of cell hits."""
    times = np.asarray(times_ns, dtype=float)
    cells = np.asarray(cell_ids, dtype=np.int64)
    if times.shape != cells.shape or recovery_time_ns <= 0:
        raise ValueError("times/cells or recovery time are invalid")
    order = np.argsort(times, kind="stable")
    fractions = np.empty_like(times)
    last: dict[int, float] = {}
    for index in order:
        cell = int(cells[index])
        previous = last.get(cell)
        if previous is None:
            fraction = 1.0
        else:
            dt = max(0.0, float(times[index] - previous))
            fraction = -math.expm1(-dt / recovery_time_ns)
        fractions[index] = fraction
        if fraction > 0.0:
            last[cell] = float(times[index])
    return fractions


def load_measured_spe_template(path) -> tuple[np.ndarray, np.ndarray]:
    """Load the repository two-column time/amplitude measured SPE CSV."""
    data = np.genfromtxt(Path(path), delimiter=",", names=True)
    names = list(data.dtype.names or ())
    if len(names) < 2:
        raise ValueError(f"{path} does not contain two numeric columns")
    time_ns = np.asarray(data[names[0]], dtype=float)
    amplitude_mv = np.asarray(data[names[1]], dtype=float)
    finite = np.isfinite(time_ns) & np.isfinite(amplitude_mv)
    return time_ns[finite], amplitude_mv[finite]


def convolve_pe_times(times_ns, amplitudes, sample_times_ns,
                      template_time_ns, template_amplitude_mv) -> np.ndarray:
    """Directly sum shifted calibrated SPE templates on an ADC time grid."""
    output = np.zeros_like(np.asarray(sample_times_ns, dtype=float))
    for time_ns, amplitude in zip(times_ns, amplitudes):
        output += float(amplitude) * np.interp(
            output * 0.0 + sample_times_ns - time_ns,
            template_time_ns, template_amplitude_mv, left=0.0, right=0.0)
    return output


def digitize_adc(waveform_mv, bits: int, full_scale_mv: float) -> np.ndarray:
    """Symmetric signed ADC quantization in millivolts."""
    if bits < 2 or full_scale_mv <= 0:
        raise ValueError("invalid ADC specification")
    half = full_scale_mv / 2.0
    levels = 2**bits
    step = full_scale_mv / levels
    clipped = np.clip(np.asarray(waveform_mv, float), -half, half - step)
    codes = np.rint((clipped + half) / step)
    return codes * step - half


def simulate_short_pair_waveforms(
        duration_ns: float,
        star_rate_hz: float,
        nsb_rate_hz: float,
        visibility2: float,
        instrument: Instrument,
        template_time_ns,
        template_amplitude_mv,
        delay_ns: float = 0.0,
        recovery_time_ns: float | None = None,
        seed: int = 20260824) -> dict[str, np.ndarray | float]:
    """Simulate a representative two-telescope ADC record.

    Correlated star counts use a shared-Poisson construction whose covariance
    per ADC bin equals ``r_star^2 * coherence_area * dt * |V|^2``.  NSB and
    remaining star counts are independent.  This reproduces the required
    second moment without pretending to resolve the femtosecond optical field.
    """
    rng = np.random.default_rng(seed)
    dt_ns = instrument.sample_width_ns
    edges = np.arange(0.0, duration_ns + dt_ns, dt_ns)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bins = len(centers)
    dt_s = dt_ns * 1.0e-9
    shared_mean = (star_rate_hz**2 * instrument.coherence_area_s * dt_s
                   * max(0.0, float(visibility2)))
    star_mean = star_rate_hz * dt_s
    nsb_mean = nsb_rate_hz * dt_s
    if shared_mean > star_mean:
        raise ValueError("shared thermal component exceeds the star count")
    common = rng.poisson(shared_mean, bins)
    counts_a = common + rng.poisson(star_mean - shared_mean, bins)
    counts_b = common + rng.poisson(star_mean - shared_mean, bins)
    counts_a += rng.poisson(nsb_mean, bins)
    counts_b += rng.poisson(nsb_mean, bins)

    def expand(counts, time_shift_ns):
        indices = np.repeat(np.arange(bins), counts)
        times = edges[indices] + rng.random(len(indices)) * dt_ns + time_shift_ns
        cells = rng.integers(0, instrument.microcells_per_pixel, len(times))
        tau = (instrument.microcell_recovery_time_ns if recovery_time_ns is None
               else recovery_time_ns)
        amplitudes = apply_exponential_microcell_recovery(times, cells, tau)
        return times, amplitudes

    times_a, amplitudes_a = expand(counts_a, 0.0)
    times_b, amplitudes_b = expand(counts_b, delay_ns)
    waveform_a = convolve_pe_times(
        times_a, amplitudes_a, centers, template_time_ns,
        template_amplitude_mv)
    waveform_b = convolve_pe_times(
        times_b, amplitudes_b, centers, template_time_ns,
        template_amplitude_mv)
    waveform_a += rng.normal(0.0, instrument.electronic_noise_rms_mv, bins)
    waveform_b += rng.normal(0.0, instrument.electronic_noise_rms_mv, bins)
    adc_a = digitize_adc(waveform_a, instrument.adc_bits,
                         instrument.adc_full_scale_mv)
    adc_b = digitize_adc(waveform_b, instrument.adc_bits,
                         instrument.adc_full_scale_mv)
    return {
        "sample_time_ns": centers,
        "adc_a_mv": adc_a,
        "adc_b_mv": adc_b,
        "analog_a_mv": waveform_a,
        "analog_b_mv": waveform_b,
        "pe_times_a_ns": times_a,
        "pe_times_b_ns": times_b,
        "recovery_a": amplitudes_a,
        "recovery_b": amplitudes_b,
        "shared_count_mean_per_bin": shared_mean,
        "expected_sample_contrast": (
            instrument.coherence_area_s / dt_s * visibility2),
    }


def normalized_cross_correlation(left, right) -> tuple[np.ndarray, np.ndarray]:
    left = np.asarray(left, float) - np.mean(left)
    right = np.asarray(right, float) - np.mean(right)
    scale = np.std(left) * np.std(right) * len(left)
    correlation = np.correlate(left, right, mode="full") / scale
    lags = np.arange(-len(left) + 1, len(left))
    return lags, correlation
