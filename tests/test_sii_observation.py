"""用独立解析信号及热光矩关系检查共享数据和到达时间约定。"""
from dataclasses import replace
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'python'))
from sii_unified import Instrument, C_M_S
from sii_observation import (geometric_arrival_delays_ns, align_waveforms,
                             correlate_blocks, simulate_array_photon_times,
                             simulate_array_waveforms, joint_thermal_mode_counts, tracking_geometry)


def test_tracking_derivative_and_curvature_against_rotating_direction():
    from sii_layout import read_corsika_layout
    positions = read_corsika_layout(ROOT/'configs/arrays/lact36_20260906.input')[
        ['east_m', 'north_m', 'up_m']].to_numpy()
    for elapsed in [-10800., 0., 10800.]:
        state = tracking_geometry(positions, .1, .3, .5, elapsed)
        before = tracking_geometry(positions, .1, .3, .5, elapsed-.1)['arrival_delays_ns']
        after = tracking_geometry(positions, .1, .3, .5, elapsed+.1)['arrival_delays_ns']
        np.testing.assert_allclose((after-before)/.2, state['arrival_delay_rates_ns_per_s'], atol=2e-11)
        exact = tracking_geometry(positions, .1, .3, .5, elapsed+10.)['arrival_delays_ns']
        linear = state['arrival_delays_ns']+10*state['arrival_delay_rates_ns_per_s']
        assert np.all(abs(exact-linear) <= .5*state['curvature_bound_ns_per_s2']*100+1e-10)


def test_varying_delay_tracks_independent_analytic_signal():
    dt = 4.
    time = (np.arange(12000)+.5)*dt
    delays, rates = np.array([0., 13.3, -19.1]), np.array([0., 2e6, -3e6])
    # 人工加大时间伸缩，只验证符号、半样本坐标及变化读取；不是地球自转场景。
    def signal(t):
        return 2.+np.cos(.021*t)+.4*np.sin(.113*t)
    adc = np.array([signal((time-delay)/(1+rate*1e-9)) for delay, rate in zip(delays, rates)])
    result = align_waveforms(adc, delays, dt, half_width=64, arrival_delay_rates_ns_per_s=rates)
    expected = np.tile(signal(result['sample_time_ns']), (3, 1))
    np.testing.assert_allclose(result['adc_mv'], expected, atol=3e-5, rtol=0)
    frozen = align_waveforms(adc, delays, dt, half_width=64)
    assert np.sqrt(np.mean((frozen['adc_mv'][1]-signal(frozen['sample_time_ns']))**2)) > .1
    np.testing.assert_array_equal(align_waveforms(adc, delays, dt)['adc_mv'],
                                  align_waveforms(adc, delays, dt, arrival_delay_rates_ns_per_s=[0.,0.,0.])['adc_mv'])


def test_time_mapped_photons_keep_shared_wavefronts_and_receiver_rate():
    instrument = replace(Instrument.from_repository(ROOT), optical_timing_kernel_path=None)
    # 两镜全相关且独立余项为零；接收时刻逆映射必须还原同一组真实共享事件。
    rate = 1e9
    scale = 1/(rate*instrument.coherence_area_s)
    factors = np.array([1., 1.05])
    delays = np.array([0., 17.])
    times, meta = simulate_array_photon_times(np.random.default_rng(25), 1e5, rate, 0.,
        np.ones((2,2)), instrument, delays, scale, padding_ns=0.,
        arrival_delay_rates_ns_per_s=(factors-1)*1e9)
    np.testing.assert_allclose(times[0], (times[1]-delays[1])/factors[1], atol=1e-10, rtol=0)
    for index in range(2):
        observed = np.count_nonzero((times[index]>=0)&(times[index]<1e5))
        expected = rate/factors[index]*1e-4
        assert abs(observed-expected) < 6*np.sqrt(expected)
    np.testing.assert_allclose(meta['received_star_rate_hz'], rate/factors)


@pytest.mark.parametrize('rates', [[0., np.nan], [0., -1e9], [0.]])
def test_invalid_time_mapping_is_rejected(rates):
    with pytest.raises(ValueError, match='时延率'):
        align_waveforms(np.ones((2, 200)), [0., 0.], 4., arrival_delay_rates_ns_per_s=rates)


def test_wavefront_arrival_sign_and_reference():
    positions = [[0, 0, 0], [0, 0, 10], [10, 0, 0]]
    # 赤道天顶：高10米的望远镜先接收波面，提前10/c；东向没有时差。
    delays = geometric_arrival_delays_ns(positions, 0., 0., 0.)
    np.testing.assert_allclose(delays, [0., -10/C_M_S*1e9, 0.], atol=1e-12)
    shifted = geometric_arrival_delays_ns(np.asarray(positions)+[7, 8, 9], 0., 0., 0., 1)
    np.testing.assert_allclose(shifted, delays-delays[1], atol=1e-12)


@pytest.mark.parametrize('delays', [[0., 12., -20.], [0., 13.3, -19.1]])
def test_delayed_bandlimited_signal_is_aligned_without_wrap(delays):
    dt = 4.
    time = (np.arange(2048)+.5)*dt
    def signal(t):
        return 2.+np.cos(.021*t)+.4*np.sin(.113*t)
    adc = np.array([signal(time-delay) for delay in delays])
    result = align_waveforms(adc, delays, dt)
    np.testing.assert_allclose(result['adc_mv'], np.tile(signal(result['sample_time_ns']), (3, 1)),
                               atol=4e-4, rtol=0)
    fine = align_waveforms(adc, delays, dt, half_width=64)
    np.testing.assert_allclose(fine['adc_mv'], np.tile(signal(fine['sample_time_ns']), (3, 1)),
                               atol=3e-5, rtol=0)
    assert result['effective_duration_ns'] == result['adc_mv'].shape[1]*dt
    assert result['discarded_samples'] > 0
    # 边缘脉冲只允许裁掉，不能从另一端重新进入。
    impulse = np.zeros((2, 64))
    impulse[1, 0] = 1.
    assert not np.any(align_waveforms(impulse, [0., 4.], 4.)['adc_mv'])


def test_block_exposure_and_reused_channel():
    rng = np.random.default_rng(33)
    a, b = rng.normal(size=(2, 1030))
    result = correlate_blocks(np.array([a, b, b]), 4., 40., block_samples=256)
    # 镜1和镜2完全相同，所以基线01与02必须逐值相同。
    np.testing.assert_array_equal(result['block_correlations'][:, 0], result['block_correlations'][:, 1])
    np.testing.assert_allclose(result['block_correlations'][:, 2, 10], 1., atol=1e-14)
    assert result['block_count'] == 4 and result['discarded_tail_samples'] == 6
    assert result['effective_duration_s'] == 1024*4e-9
    with pytest.raises(ValueError, match='恒定'):
        correlate_blocks(np.zeros((3, 100)), 4., 40.)


def test_joint_pair_stream_marginals_and_shared_events():
    instrument = replace(Instrument.from_repository(ROOT), optical_timing_kernel_path=None)
    rates = np.array([1e9, .9e9, .8e9])
    gamma = np.full((3, 3), .6)+np.eye(3)*.4
    times, meta = simulate_array_photon_times(np.random.default_rng(5), 1e5, rates, 0.,
                                              gamma, instrument, pair_rate_scale=1000., padding_ns=0.)
    for index in range(3):
        assert abs(len(times[index])-rates[index]*1e-4) < 6*np.sqrt(rates[index]*1e-4)
    for left, right in ((0, 1), (0, 2), (1, 2)):
        assert len(np.intersect1d(times[left], times[right])) == meta['pair_counts_with_padding'][left, right]
    np.testing.assert_allclose(meta['single_star_rates_hz']+meta['injected_pair_rates_hz'].sum(axis=1), rates)
    with pytest.raises(ValueError, match='全部相关对率'):
        simulate_array_photon_times(np.random.default_rng(5), 100., rates, 0., np.ones((3, 3)),
                                    instrument, pair_rate_scale=10000.)


def test_thermal_joint_second_and_third_moments():
    gamma = np.full((3, 3), .7, complex)+np.eye(3)*.3
    gamma[0, 1] *= np.exp(.4j)
    gamma[1, 0] = gamma[0, 1].conjugate()
    means = np.array([4., 5., 6.])
    counts, intensities = joint_thermal_mode_counts(np.random.default_rng(15), 240_000,
                                                    means, .3, gamma, modes=2, polarizations=1)
    centered = intensities-1.
    for left, right in ((0, 1), (0, 2), (1, 2)):
        product = centered[:, left]*centered[:, right]
        expected = abs(gamma[left, right])**2/2
        assert abs(product.mean()-expected) < 6*product.std(ddof=1)/np.sqrt(len(product))
    triple = centered.prod(axis=1)
    expected = 2*np.real(gamma[0, 1]*gamma[1, 2]*gamma[2, 0])/2**2
    assert abs(triple.mean()-expected) < 6*triple.std(ddof=1)/np.sqrt(len(triple))
    np.testing.assert_allclose(counts.mean(axis=0), means+.3, rtol=.01)
    np.testing.assert_allclose(counts.var(axis=0), means+.3+means**2/2, rtol=.03)


@pytest.mark.parametrize('gamma', [np.array([[1, .9, .9], [.9, 1, -.9], [.9, -.9, 1]]),
                                   np.array([[1, .5j], [.5j, 1]]), np.ones((2, 3))])
def test_invalid_joint_coherence_is_rejected(gamma):
    with pytest.raises(ValueError, match='复相干矩阵'):
        joint_thermal_mode_counts(np.random.default_rng(0), 10, 1., 0., gamma)


def test_measured_array_waveform_smoke_and_empty_overlap():
    instrument = Instrument.from_repository(ROOT)
    record = simulate_array_waveforms(np.random.default_rng(9), 4000., 2e8, 5e5,
                                      np.eye(3), instrument, [0., 23.5, -46.2])
    aligned = align_waveforms(record['adc_mv'], [0., 23.5, -46.2], instrument.sample_width_ns)
    curves = correlate_blocks(aligned['adc_mv'], instrument.sample_width_ns, 40.)
    assert curves['mean_correlation'].shape == (3, 21)
    with pytest.raises(ValueError, match='共同有效'):
        align_waveforms(record['adc_mv'], [0., 1e5, 0.], instrument.sample_width_ns)


def test_unimplemented_instrument_state_is_not_silently_ignored():
    instrument = replace(Instrument.from_repository(ROOT), transparency_fractional_rms=.1)
    with pytest.raises(NotImplementedError, match='transparency'):
        simulate_array_waveforms(np.random.default_rng(0), 4000., 2e8, 5e5, np.eye(3), instrument)
