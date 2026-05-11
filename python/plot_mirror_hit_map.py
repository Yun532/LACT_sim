#!/usr/bin/env python3
"""Plot mirror hit positions, optionally overlaid with projected facet outlines."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm

from config_io import (
    apply_telescope_frame_to_facets,
    expand_component_config,
    load_facets_from_config,
    telescope_frame_from_config,
)
from plot_optical_layout_3d import aperture_polygon


def projected_facet_outlines(config_path: Path) -> list[np.ndarray]:
    cfg, _ = expand_component_config(config_path)
    frame = telescope_frame_from_config(cfg)
    facets = apply_telescope_frame_to_facets(load_facets_from_config(config_path, cfg), frame)
    outlines = []
    for facet in facets:
        polygon = aperture_polygon(facet)[:, :2] * 1000.0
        outlines.append(np.vstack([polygon, polygon[0]]))
    return outlines


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hits_csv", help="whiteboard/camera hit CSV from run_optical_sim")
    parser.add_argument("--config", help="optical config used to run the simulation")
    parser.add_argument("--output", required=True, help="output PNG/PDF path")
    parser.add_argument("--require-surface", action="store_true", help="keep only hit_surface=1 rows")
    parser.add_argument("--overlay-facets", action="store_true", help="draw mirror facet x-y projected outlines")
    parser.add_argument("--bins", type=int, default=380, help="2D histogram bins")
    parser.add_argument("--dpi", type=int, default=350, help="output DPI")
    parser.add_argument("--title", default="Mirror hit points", help="plot title")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    hits_path = Path(args.hits_csv)
    df = pd.read_csv(hits_path)
    required = {"hit_mirror", "mirror_x", "mirror_y"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"{hits_path} is missing required columns: {', '.join(sorted(missing))}")

    df = df[df["hit_mirror"] == 1].copy()
    if args.require_surface:
        if "hit_surface" not in df.columns:
            raise SystemExit("--require-surface needs a hit_surface column")
        df = df[df["hit_surface"] == 1].copy()
    if df.empty:
        raise SystemExit("no mirror hits selected")

    x = df["mirror_x"].to_numpy(dtype=float) * 1000.0
    y = df["mirror_y"].to_numpy(dtype=float) * 1000.0
    hist, xedges, yedges = np.histogram2d(x, y, bins=max(1, args.bins))
    masked = np.ma.masked_where(hist.T <= 0, hist.T)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad((1.0, 1.0, 1.0, 0.0))
    positive = hist[hist > 0]
    norm = LogNorm(vmin=max(1, positive.min()), vmax=positive.max()) if positive.size else None

    fig, ax = plt.subplots(figsize=(7.4, 6.5), dpi=args.dpi)
    image = ax.imshow(
        masked,
        origin="lower",
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
    )

    if args.overlay_facets:
        if not args.config:
            raise SystemExit("--overlay-facets requires --config")
        outlines = projected_facet_outlines(Path(args.config).resolve())
        ax.add_collection(LineCollection(outlines, colors="black", linewidths=0.75, alpha=0.9))

    mean_x = float(np.mean(x))
    mean_y = float(np.mean(y))
    ax.scatter(mean_x, mean_y, marker="x", s=85, linewidths=2.6, color="yellow", label="Mean", zorder=5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Mirror hit x [mm]")
    ax.set_ylabel("Mirror hit y [mm]")
    ax.set_title(args.title)
    ax.grid(alpha=0.16, linewidth=0.6)
    ax.legend(loc="upper right", frameon=True, framealpha=0.95, edgecolor="black")
    cbar = fig.colorbar(image, ax=ax, pad=0.02)
    cbar.set_label("Count / bin")

    notes = [f"N = {len(df)}", f"mean = ({mean_x:.1f}, {mean_y:.1f}) mm"]
    if args.overlay_facets:
        notes.insert(1, "facet outlines: x-y projection")
    ax.text(
        0.02,
        0.02,
        "\n".join(notes),
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        fontsize=9,
        bbox=dict(facecolor="white", edgecolor="0.75", alpha=0.9, boxstyle="round,pad=0.28"),
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    print(f"Saved mirror hit map = {output}")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
