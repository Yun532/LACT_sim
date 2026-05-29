#!/usr/bin/env python3
"""Plot input-plane x/y distributions for CORSIKA whiteboard hits."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_xy(path: Path, event_id: int | None, telescope_id: int | None) -> tuple[np.ndarray, np.ndarray]:
    xs: list[float] = []
    ys: list[float] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if "input_x_m" not in (reader.fieldnames or ()) or "input_y_m" not in (reader.fieldnames or ()):
            raise SystemExit(
                f"{path} does not contain input_x_m/input_y_m. "
                "Re-run run_corsika_trace with the updated whiteboard CSV writer."
            )
        for row in reader:
            if event_id is not None and int(row["event_id"]) != event_id:
                continue
            if telescope_id is not None and int(row["telescope_id"]) != telescope_id:
                continue
            xs.append(float(row["input_x_m"]))
            ys.append(float(row["input_y_m"]))
    return np.asarray(xs), np.asarray(ys)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare input x/y distributions of photons that hit the whiteboard."
    )
    parser.add_argument("csv", nargs="+", help="whiteboard hit CSV files")
    parser.add_argument("--label", action="append", help="legend label; repeat once per CSV")
    parser.add_argument("--event-id", type=int)
    parser.add_argument("--telescope-id", type=int, help="0-based telescope id")
    parser.add_argument("--output", required=True)
    parser.add_argument("--bins", type=int, default=180)
    parser.add_argument("--max-points", type=int, default=250000)
    args = parser.parse_args()

    paths = [Path(p) for p in args.csv]
    labels = args.label or [p.stem for p in paths]
    if len(labels) != len(paths):
        raise SystemExit("--label must be repeated once per CSV")

    datasets: list[tuple[str, np.ndarray, np.ndarray]] = []
    for label, path in zip(labels, paths):
        x, y = read_xy(path, args.event_id, args.telescope_id)
        if len(x) == 0:
            raise SystemExit(f"No rows selected from {path}")
        datasets.append((label, x, y))

    all_x = np.concatenate([x for _, x, _ in datasets])
    all_y = np.concatenate([y for _, _, y in datasets])
    xpad = max(0.5, 0.04 * (float(all_x.max()) - float(all_x.min()) + 1e-9))
    ypad = max(0.5, 0.04 * (float(all_y.max()) - float(all_y.min()) + 1e-9))
    extent = [
        float(all_x.min()) - xpad,
        float(all_x.max()) + xpad,
        float(all_y.min()) - ypad,
        float(all_y.max()) + ypad,
    ]

    fig, axes = plt.subplots(1, len(datasets), figsize=(6.2 * len(datasets), 5.6), squeeze=False)
    for ax, (label, x, y) in zip(axes[0], datasets):
        hist, xedges, yedges = np.histogram2d(x, y, bins=args.bins, range=[[extent[0], extent[1]], [extent[2], extent[3]]])
        cmap = plt.get_cmap("viridis").copy()
        cmap.set_under((1, 1, 1, 0))
        im = ax.imshow(
            hist.T,
            origin="lower",
            extent=extent,
            cmap=cmap,
            vmin=0.5,
            aspect="equal",
            interpolation="nearest",
        )
        if args.max_points > 0 and len(x) > 0:
            step = max(1, len(x) // args.max_points)
            ax.scatter(x[::step], y[::step], s=0.15, c="black", alpha=0.10, linewidths=0)
        ax.set_title(f"{label}\nN={len(x):,}")
        ax.set_xlabel("input x (m)")
        ax.set_ylabel("input y (m)")
        ax.grid(alpha=0.18, linewidth=0.5)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.035)
        cbar.set_label("count")

    subtitle = []
    if args.event_id is not None:
        subtitle.append(f"event_id={args.event_id}")
    if args.telescope_id is not None:
        subtitle.append(f"telescope_id={args.telescope_id}")
    if subtitle:
        fig.suptitle("Whiteboard-hit photons: input-plane x/y distribution (" + ", ".join(subtitle) + ")")
    else:
        fig.suptitle("Whiteboard-hit photons: input-plane x/y distribution")
    fig.tight_layout()
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
