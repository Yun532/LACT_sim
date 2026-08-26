#!/usr/bin/env python3
"""Plot the two NSB2 interface modes with pyLAST and measure the pedestal."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from pylast.visualize import plot_event_quicklook, plot_lact_waveforms


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expectation-root", type=Path, required=True)
    parser.add_argument("--waveform-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--baseline-samples", type=int, default=16)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    expectation = plot_event_quicklook(
        args.expectation_root,
        output_dir=args.output_dir,
        image_level="dl0",
        plot_gathered=False,
        plot_telescopes=False,
        plot_event=True,
        include_non_triggered=True,
        show=False,
    )
    event_a = expectation["event"]
    tel_a = expectation["source"].get_readout_tels(event_a)[0]
    image_a = np.asarray(event_a.dl0.tels[tel_a].image, dtype=float)

    raw = plot_lact_waveforms(
        args.waveform_root,
        level="raw",
        max_pixels=4,
        output_path=args.output_dir / "mode_b_raw_mv.png",
        show=False,
    )
    pixels = raw["pixel_indices"]
    r1 = plot_lact_waveforms(
        args.waveform_root,
        level="r1",
        pixel_indices=pixels,
        output_path=args.output_dir / "mode_b_r1_no_baseline.png",
        show=False,
    )
    corrected = plot_lact_waveforms(
        args.waveform_root,
        level="r1",
        baseline_samples=args.baseline_samples,
        pixel_indices=pixels,
        output_path=args.output_dir / "mode_b_r1_baseline_subtracted.png",
        show=False,
    )

    count = min(args.baseline_samples, raw["waveform"].shape[1])
    raw_baselines = raw["waveform"][:, :count].mean(axis=1)
    r1_baselines = r1["waveform"][:, :count].mean(axis=1)
    corrected_baselines = corrected["waveform"][:, :count].mean(axis=1)
    summary = {
        "mode_a": {
            "exposure_ns": 1.0e9,
            "camera_total_pe": float(image_a.sum()),
            "camera_mean_pe_per_pixel": float(image_a.mean()),
            "camera_nonzero_pixels": int(np.count_nonzero(image_a)),
        },
        "mode_b": {
            "integration_ns": 128.0,
            "baseline_samples": count,
            "raw_baseline_mean_mv": float(raw_baselines.mean()),
            "raw_baseline_median_mv": float(np.median(raw_baselines)),
            "raw_baseline_pixel_rms_mv": float(raw_baselines.std()),
            "r1_baseline_mean_pe_per_bin_before": float(r1_baselines.mean()),
            "r1_baseline_mean_pe_per_bin_after": float(corrected_baselines.mean()),
            "r1_baseline_abs_max_pe_per_bin_after": float(
                np.abs(corrected_baselines).max()
            ),
            "plotted_pixel_indices": pixels,
        },
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

