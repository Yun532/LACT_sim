#!/usr/bin/env python3
"""Plot the two NSB2 interface modes with pyLAST and measure the pedestal."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np

from pylast.helper import Calibrator, LactEventSource
from pylast.visualize import (
    EventVisualizer,
    plot_event_cameras,
    plot_event_quicklook,
    plot_lact_waveforms,
)


def calibrate_full_waveform(root_file: Path, baseline_samples: int):
    """Return one pyLAST-calibrated event and its integrated camera image."""

    source = LactEventSource(
        str(root_file), baseline_samples=int(baseline_samples)
    )
    event = source[0]
    Calibrator(
        source.subarray,
        config_str=json.dumps({"image_extractor_type": "FullWaveFormExtractor"}),
    )(event)
    tel_id = source.get_readout_tels(event.event_id)[0]
    image = np.asarray(event.dl0.tels[tel_id].image, dtype=float)
    return source, event, tel_id, image


def plot_signed_camera(source, tel_id: int, image: np.ndarray, output_path: Path):
    """Draw a baseline-subtracted camera while preserving signed residuals."""

    visualizer = EventVisualizer(source, enable_secondary_axes=True)
    scale = max(float(np.max(np.abs(image))), 1.0e-12)
    norm = TwoSlopeNorm(vmin=-scale, vcenter=0.0, vmax=scale)
    figure, axis = plt.subplots(figsize=(7.2, 6.2))
    visualizer._draw_camera_image(
        axis,
        visualizer.tel_geoms[tel_id],
        image,
        norm,
        plt.get_cmap("RdBu_r"),
    )
    axis.set_title(
        "Mode B: pyLAST FullWaveFormExtractor\n"
        "per-pixel baseline subtracted; signed residual"
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(figure)


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

    source_b, event_b, tel_b, image_b = calibrate_full_waveform(
        args.waveform_root, baseline_samples=0
    )
    plot_event_cameras(
        event_b,
        source=source_b,
        image_level="dl0",
        output_path=args.output_dir / "mode_b_camera_integrated_pe.png",
        include_non_triggered=True,
        show=False,
    )
    source_b_corrected, _, tel_b_corrected, image_b_corrected = (
        calibrate_full_waveform(
            args.waveform_root, baseline_samples=args.baseline_samples
        )
    )
    plot_signed_camera(
        source_b_corrected,
        tel_b_corrected,
        image_b_corrected,
        args.output_dir / "mode_b_camera_baseline_subtracted_residual_pe.png",
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
            "camera_total_integrated_pe_before": float(image_b.sum()),
            "camera_mean_integrated_pe_per_pixel_before": float(image_b.mean()),
            "camera_total_residual_pe_after": float(image_b_corrected.sum()),
            "camera_mean_residual_pe_per_pixel_after": float(
                image_b_corrected.mean()
            ),
            "camera_residual_rms_pe_after": float(image_b_corrected.std()),
            "camera_residual_min_pe_after": float(image_b_corrected.min()),
            "camera_residual_max_pe_after": float(image_b_corrected.max()),
            "plotted_pixel_indices": pixels,
        },
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
