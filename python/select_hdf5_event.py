#!/usr/bin/env python3
"""Select a bright common event from one or more LACT_sim HDF5 camera files."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Select a common event_id with enough camera signal from HDF5 files."
    )
    parser.add_argument("h5", nargs="+", help="HDF5 files to intersect by event_id")
    parser.add_argument("--quantity", default="pe", choices=("pe", "signal", "photon_count"))
    parser.add_argument("--min-images", type=int, default=1)
    parser.add_argument("--output-env", default=None, help="write shell variables for sourcing")
    parser.add_argument("--output-summary", default=None, help="write a human-readable summary")
    return parser.parse_args()


def event_scores(path: str, quantity: str) -> dict[int, dict[str, float]]:
    with h5py.File(path, "r") as h5:
        if "images/index" not in h5:
            raise SystemExit(f"{path} has no images/index dataset")
        rows = h5["images/index"][:]
        scores: dict[int, dict[str, float]] = {}
        for row in rows:
            event_id = int(row["event_id"])
            entry = scores.setdefault(
                event_id,
                {"total": 0.0, "images": 0, "telescopes": set()},
            )
            if quantity == "pe":
                value = float(row["total_pe"])
            elif quantity == "signal":
                value = float(row["total_signal"])
            else:
                value = float(row["count"])
            entry["total"] += value
            entry["images"] += 1
            entry["telescopes"].add(int(row["telescope_id"]))
        for entry in scores.values():
            entry["n_telescopes"] = len(entry["telescopes"])
            del entry["telescopes"]
        return scores


def event_metadata(path: str, event_id: int) -> dict[str, float | int] | None:
    with h5py.File(path, "r") as h5:
        if "events/corsika" not in h5:
            return None
        rows = h5["events/corsika"][:]
        matches = rows[rows["event_id"] == event_id]
        if len(matches) == 0:
            return None
        row = matches[0]
        result = {}
        for name in row.dtype.names or ():
            value = row[name]
            if isinstance(value, np.generic):
                value = value.item()
            result[name] = value
        return result


def main() -> None:
    args = parse_args()
    all_scores = [event_scores(path, args.quantity) for path in args.h5]
    common = set(all_scores[0])
    for scores in all_scores[1:]:
        common &= set(scores)
    if not common:
        raise SystemExit("No common event_id appears in all requested HDF5 files.")

    candidates = []
    for event_id in sorted(common):
        min_images = min(scores[event_id]["images"] for scores in all_scores)
        if min_images < args.min_images:
            continue
        combined_total = sum(scores[event_id]["total"] for scores in all_scores)
        baseline_total = all_scores[0][event_id]["total"]
        min_telescopes = min(scores[event_id]["n_telescopes"] for scores in all_scores)
        candidates.append((combined_total, baseline_total, min_telescopes, event_id))
    if not candidates:
        raise SystemExit(
            f"No common event_id has at least {args.min_images} image rows in every HDF5 file."
        )
    candidates.sort(reverse=True)
    combined_total, baseline_total, min_telescopes, event_id = candidates[0]
    meta = event_metadata(args.h5[0], event_id) or {}

    summary_lines = [
        "Selected official CORSIKA plot event",
        "====================================",
        f"event_id: {event_id}",
        f"quantity: {args.quantity}",
        f"combined_total: {combined_total:.9g}",
        f"baseline_total: {baseline_total:.9g}",
        f"min_telescopes_with_images: {min_telescopes}",
        f"n_common_events: {len(common)}",
    ]
    for key in (
        "shower_event_id",
        "array_id",
        "core_x_north_m",
        "core_y_west_m",
        "azimuth_north_to_east_deg",
        "energy_gev",
    ):
        if key in meta:
            summary_lines.append(f"{key}: {meta[key]}")

    text = "\n".join(summary_lines) + "\n"
    print(text, end="")
    if args.output_summary:
        out = Path(args.output_summary)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)

    if args.output_env:
        out = Path(args.output_env)
        out.parent.mkdir(parents=True, exist_ok=True)
        array_id = int(meta.get("array_id", event_id % 100))
        shower_event_id = int(meta.get("shower_event_id", event_id // 100))
        env = [
            f"LACT_SELECTED_EVENT_ID={event_id}",
            f"LACT_SELECTED_ARRAY_ID={array_id}",
            f"LACT_SELECTED_SHOWER_EVENT_ID={shower_event_id}",
        ]
        out.write_text("\n".join(env) + "\n")


if __name__ == "__main__":
    main()
