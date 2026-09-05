"""执行独立波形标定检查和小规模重建检查。"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'python'))
import sii_unified as sii


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--records', type=int, default=1024)
    parser.add_argument('--seed', type=int, default=20260824)
    args = parser.parse_args()
    instrument = sii.Instrument.from_repository(ROOT)
    options = dict(null_records=args.records, signal_records=max(16, args.records//4))
    calibration, diagnostics = sii.simulate_waveform_gls_calibration(
        instrument, seed=args.seed, **options)
    print('Reference calibration complete', flush=True)
    _, independent = sii.simulate_waveform_gls_calibration(
        instrument, seed=args.seed+1, **options)
    null, sigma = sii.estimate_visibility2_gls(independent['null_correlations'], calibration)
    signal, _ = sii.estimate_visibility2_gls(independent['signal_correlations'], calibration)
    signal = signal/(calibration.calibration_visibility2*calibration.hbt_pair_rate_scale)
    ratio = float(null.std(ddof=1)/sigma)
    response = float(signal.mean())
    response_sem = float(signal.std(ddof=1)/np.sqrt(len(signal)))
    ratio_uncertainty = np.hypot(calibration.sigma_relative_uncertainty,
                                1/np.sqrt(2*(len(null)-1)))
    assert abs(ratio-1) < 5*ratio_uncertainty, (ratio, ratio_uncertainty)
    assert abs(response-1) < 5*np.hypot(response_sem, calibration.response_relative_uncertainty)
    assert abs(null.mean()/sigma) < 5/np.sqrt(len(null))

    layout = pd.DataFrame({'east_m':[0.,100.,0.,100.], 'north_m':[0.,0.,100.,100.],
                           'up_m':np.zeros(4)})
    pipeline = sii.run_sii_pipeline(layout, instrument=instrument,
        observation=sii.Observation(hours_per_night=1, segment_s=1200),
        estimator='waveform_gls', waveform_calibration=calibration,
        reconstruction_kwargs=dict(grid_size=12, starts=1, max_iter=200, smoothness=.001))
    assert np.isclose(pipeline.reconstruction.image.sum(), 1.)
    assert np.all(np.isfinite(pipeline.reconstruction.predicted_visibility_abs2))

    results = {
        'seed':args.seed, 'null_records':args.records,
        'star_rate_hz':calibration.star_rate_hz,
        'nsb_rate_hz':instrument.detected_nsb_rate_hz,
        'coherence_area_s':instrument.coherence_area_s,
        'sigma_visibility2_1200_s':sii.waveform_gls_weights(calibration, 1200)[1],
        'analytic_matched_sigma_1200_s':1/sii.unit_visibility_snr(
            2, 1200, sii.with_matched_effective_bandwidth(instrument)),
        'selected_shrinkage':calibration.covariance_shrinkage,
        'shrinkage_selection':diagnostics['shrinkage_selection'],
        'independent_null_std_over_sigma':ratio,
        'independent_null_mean_over_sigma':float(null.mean()/sigma),
        'independent_response':response, 'independent_response_sem':response_sem,
        'response_relative_calibration_uncertainty':calibration.response_relative_uncertainty,
        'sigma_relative_calibration_uncertainty':calibration.sigma_relative_uncertainty,
        'reconstruction_measurements':len(pipeline.measurements),
        'reconstruction_image_flux':float(pipeline.reconstruction.image.sum()),
        'reconstruction_optimizer_success':pipeline.reconstruction.metrics['optimizer_success'],
        'uncertainty_model':pipeline.metadata['uncertainty_model'],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
