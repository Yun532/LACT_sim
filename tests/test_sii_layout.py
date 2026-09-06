"""用手算坐标、距离不变性和实际UV条数检验input接口。"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'python'))
from sii_layout import read_corsika_layout
from sii_unified import generate_uvw, Observation, Instrument


def test_new_layout_conversion_and_baselines():
    layout = read_corsika_layout(ROOT/'configs/arrays/lact36_20260906.input')
    assert layout.name.tolist() == [f'Tel.{i}' for i in range(1, 37)]
    np.testing.assert_allclose(layout.loc[[0, 35], ['east_m', 'north_m', 'up_m']].to_numpy(float),
                               [[-434.92063, 417.46032, 0.], [266.66667, 136.50794, 0.]], atol=1e-10)
    np.testing.assert_array_equal(layout.radius_m, np.full(36, 4.))
    csv = pd.read_csv(ROOT/'configs/arrays/lact36_20260906_coordinates.csv')
    pd.testing.assert_frame_equal(layout, csv, check_exact=False, atol=1e-10, rtol=0)
    positions = layout[['east_m', 'north_m', 'up_m']].to_numpy()
    original = layout[['corsika_north_cm', 'corsika_west_cm', 'corsika_up_cm']].to_numpy()
    np.testing.assert_allclose(np.linalg.norm(positions[:, None]-positions, axis=2),
                               np.linalg.norm(original[:, None]-original, axis=2)*.01, atol=1e-10)
    uv = generate_uvw(layout, Observation(), Instrument.from_repository(ROOT))
    assert len(uv[['telescope_i', 'telescope_j']].drop_duplicates()) == 630
    assert len(uv) == 630*18
    pd.testing.assert_frame_equal(uv, generate_uvw(ROOT/'configs/arrays/lact36_20260906.input',
                                                  Observation(), Instrument.from_repository(ROOT)))


def test_input_height_order_and_other_cards(tmp_path):
    path = tmp_path/'array.input'
    path.write_text('RUNNR 1\nTELESCOPE 100 200 300 400 Z\nTELESCOPE -2D2 100 -50 400 A # note\n')
    layout = read_corsika_layout(path)
    assert layout.name.tolist() == ['Z', 'A']
    assert layout.source_input_line.tolist() == [2, 3]
    np.testing.assert_allclose(layout[['east_m', 'north_m', 'up_m']], [[-2, 1, 3], [-1, -2, -.5]])


@pytest.mark.parametrize('bad', ['TELESCOPE nan 0 0 400 B', 'TELESCOPE 0 0 0 -1 B',
                                  'TELESCOPE 0 0 0 400 A', 'TELESCOPE 0 0'])
def test_invalid_telescope_input_is_rejected(tmp_path, bad):
    path = tmp_path/'bad.input'
    path.write_text('TELESCOPE 0 0 0 400 A\n'+bad+'\n')
    with pytest.raises(ValueError):
        read_corsika_layout(path)


def test_nonfinite_dataframe_layout_is_rejected():
    with pytest.raises(ValueError, match='有限'):
        generate_uvw(pd.DataFrame({'east_m':[0., np.nan], 'north_m':[0., 1.], 'up_m':[0., 0.]}),
                     Observation())
