#!/usr/bin/env python3
"""把 main 的逐光线输出压缩成可复用、可审计的单光子响应。"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("hits_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--provenance-json", type=Path, required=True)
    parser.add_argument("--input-photons", type=int, required=True)
    parser.add_argument("--sampling-radius-m", type=float, default=4.0)
    parser.add_argument("--source-config", type=Path)
    args = parser.parse_args()

    columns = ["mirror_id", "time_ns", "u_m", "v_m", "hit_camera",
               "accepted", "pixel_id", "weight", "relative_efficiency"]
    frame = pd.read_csv(args.hits_csv, usecols=columns)
    signal_weight = (frame.weight * frame.relative_efficiency).clip(lower=0.0)
    detected = frame[(frame.hit_camera != 0) & (signal_weight > 0.0)].copy()
    detected["signal_weight"] = signal_weight.loc[detected.index]
    if detected.empty:
        raise ValueError("main output contains no detected camera response")

    pixel_signal = detected.groupby("pixel_id").signal_weight.sum().sort_values(
        ascending=False)
    central_pixel = int(pixel_signal.index[0])
    selected = detected.loc[detected.pixel_id == central_pixel].copy()
    grouped = selected.groupby("mirror_id", sort=True)
    rows = []
    for mirror_id, group in grouped:
        weights = group.signal_weight.to_numpy(float)
        times = group.time_ns.to_numpy(float)
        mean = np.average(times, weights=weights)
        variance = np.average((times-mean)**2, weights=weights)
        rows.append({
            "mirror_id": int(mirror_id), "count": len(group),
            "signal_weight": weights.sum(), "mean_time_ns": mean,
            "std_time_ns": np.sqrt(max(variance, 0.0)),
            "min_time_ns": times.min(), "max_time_ns": times.max(),
        })
    output = pd.DataFrame(rows)
    output["weight"] = output.signal_weight / output.signal_weight.sum()
    output = output[["mirror_id", "count", "signal_weight", "weight",
                     "mean_time_ns", "std_time_ns", "min_time_ns",
                     "max_time_ns"]]
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output_csv, index=False, float_format="%.12g")

    weights = selected.signal_weight.to_numpy(float)
    all_detected_weight = float(detected.signal_weight.sum())
    central_weight = float(weights.sum())
    mean_time = np.average(selected.time_ns, weights=weights)
    sampling_area = math.pi*args.sampling_radius_m**2
    provenance = {
        "description": "main full optical response: measured LACT2 mirror, errors, obstruction, camera, collector and spectral efficiencies",
        "source_hits_csv": str(args.hits_csv.resolve()),
        "source_hits_sha256": sha256(args.hits_csv),
        "input_photons": args.input_photons,
        "surface_rows": int(len(frame)),
        "camera_rows": int((frame.hit_camera != 0).sum()),
        "detected_signal_weight_all_pixels": all_detected_weight,
        "total_detection_probability_all_pixels": float(
            all_detected_weight/args.input_photons),
        "central_pixel_id": central_pixel,
        "illuminated_pixel_count": int(len(pixel_signal)),
        "pixel_signal_fractions": {
            str(int(pixel)): float(value/pixel_signal.sum())
            for pixel, value in pixel_signal.items()
        },
        "central_pixel_detection_probability": float(central_weight/args.input_photons),
        "source_sampling_radius_m": args.sampling_radius_m,
        "source_sampling_area_m2": sampling_area,
        "effective_detection_area_all_pixels_m2": float(
            sampling_area*all_detected_weight/args.input_photons),
        "central_pixel_effective_detection_area_m2": float(
            sampling_area*central_weight/args.input_photons),
        "mean_arrival_time_ns": float(mean_time),
        "arrival_time_rms_ns": float(np.sqrt(np.average(
            (selected.time_ns-mean_time)**2, weights=weights))),
        "psf_rms_m": float(np.sqrt(np.average(
            detected.u_m**2+detected.v_m**2,
            weights=detected.signal_weight))),
        "components": int(len(output)),
    }
    if args.source_config:
        provenance["source_config"] = str(args.source_config.resolve())
        provenance["source_config_sha256"] = sha256(args.source_config)
    args.provenance_json.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_json.write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2)+"\n",
        encoding="utf-8")


if __name__ == "__main__":
    main()
