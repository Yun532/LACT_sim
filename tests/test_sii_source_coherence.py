"""以点源Gram矩阵、非相干谱叠加和现有时间/光谱平均检查天体到联合波形接口。"""
from dataclasses import replace
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'python'))
import sii_unified as sii
from sii_layout import read_corsika_layout
from sii_observation import (source_coherence_spectrum, band_averaged_coherence_power,
                             simulate_array_photon_times)


def test_binary_point_source_matches_independent_field_gram():
    positions = np.array([[0.,0.,0.],[50.,0.,0.],[0.,80.,0.],[20.,40.,10.]])
    source = sii.BinarySource(primary_diameter_mas=0., secondary_diameter_mas=0.,
                              separation_mas=.13, position_angle_deg=90., flux_ratio_secondary_to_primary=.4)
    instrument = replace(sii.Instrument(), visibility_wavelength_nm=(390.,410.),
                         visibility_spectral_weights=(.3,.7))
    result = source_coherence_spectrum(positions,0.,0.,0.,source,instrument)
    # H=dec=lat=0时天球东西轴即ENU东；两点在东西轴±separation/2。
    for index, wave in enumerate(result['wavelength_nm']):
        phase = 2*np.pi*positions[:,0]/(wave*1e-9)*source.separation_mas/2*sii.MAS_TO_RAD
        fields = np.column_stack((np.exp(-1j*phase),np.sqrt(.4)*np.exp(1j*phase)))
        expected = fields @ fields.conj().T/1.4
        np.testing.assert_allclose(result['coherence'][index], expected, atol=1e-14)


def test_spectral_power_is_not_coherently_averaged_amplitude():
    gamma = np.array([[[1.,1.],[1.,1.]], [[1.,-1.],[-1.,1.]]])
    power = band_averaged_coherence_power(gamma,[.5,.5])
    assert power[0,1] == 1.
    assert abs(gamma.mean(axis=0)[0,1])**2 == 0.
    instrument = replace(sii.Instrument.from_repository(ROOT),optical_timing_kernel_path=None)
    _, meta = simulate_array_photon_times(np.random.default_rng(2),10000.,2e8,0.,gamma,
        instrument,spectral_weights=[.5,.5])
    assert meta['physical_pair_rates_hz'][0,1] == (2e8)**2*instrument.coherence_area_s


def test_uniform_disk_against_direct_positive_brightness_integral():
    positions=np.array([[0.,0.,0.],[100.,0.,0.],[50.,60.,0.]])
    source=sii.BinarySource(primary_diameter_mas=.16)
    result=source_coherence_spectrum(positions,0.,0.,0.,source,sii.Instrument(),'single_disk')
    axis=np.linspace(-.08,.08,401)
    x,y=np.meshgrid(axis,axis)
    inside=x*x+y*y <= .08**2
    phase=2*np.pi*(positions[:,0,None]*x[inside]+positions[:,1,None]*y[inside])*sii.MAS_TO_RAD/400e-9
    fields=np.exp(1j*phase)
    expected=fields @ fields.conj().T/inside.sum()
    np.testing.assert_allclose(result['coherence'][0],expected,atol=2e-4,rtol=0)


@pytest.mark.parametrize('case,source', [('binary',sii.BinarySource()),
    ('single_disk',sii.BinarySource(primary_diameter_mas=.16)),
    ('ellipse',sii.EllipseSource()),('transit',sii.TransitSource())])
def test_full_array_source_power_matches_existing_uv_spectral_average(case,source):
    layout=read_corsika_layout(ROOT/'configs/arrays/lact36_20260906.input')
    positions=layout[['east_m','north_m','up_m']].to_numpy()
    instrument=sii.Instrument.from_repository(ROOT)
    observation=sii.Observation(hours_per_night=1.,segment_s=1200,visibility_subsamples_per_segment=1)
    uv=sii.generate_uvw(layout,observation,instrument)
    uv=uv[uv.segment==1]
    state=source_coherence_spectrum(positions,0.,np.deg2rad(observation.source_dec_deg),
        np.deg2rad(observation.site_lat_deg),source,instrument,case)
    expected=sii.segment_averaged_visibility2(uv,source,observation,instrument,case)
    actual=state['pair_visibility2'][uv.telescope_i_index.to_numpy(),uv.telescope_j_index.to_numpy()]
    np.testing.assert_allclose(actual,expected,atol=2e-14)
    assert np.linalg.eigvalsh(state['coherence']).min() >= -1e-12


@pytest.mark.parametrize('weights', [[1.], [np.nan,0.], [-.1,1.1], [.4,.4]])
def test_invalid_spectral_weights_are_rejected(weights):
    with pytest.raises(ValueError,match='权重'):
        band_averaged_coherence_power(np.stack([np.eye(2)]*2),weights)


def test_each_spectral_matrix_must_be_physical():
    bad=np.array([[1.,.9,.9],[.9,1.,-.9],[.9,-.9,1.]])
    with pytest.raises(ValueError,match='复相干矩阵'):
        band_averaged_coherence_power(np.stack([np.eye(3),bad]),[.99,.01])
