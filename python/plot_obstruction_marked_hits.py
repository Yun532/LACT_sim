#!/usr/bin/env python3
"""Plot optical hits with obstruction-marked photons highlighted."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm

from plot_mirror_hit_map import mirror_display_basis, project_points_3d, projected_facet_outlines


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hits_csv", help="CSV written by run_optical_sim")
    parser.add_argument("--space", choices=["spot", "mirror"], required=True)
    parser.add_argument("--config", help="optical cfg; needed for --space mirror --overlay-facets")
    parser.add_argument("--overlay-facets", action="store_true")
    parser.add_argument("--output", required=True)
    parser.add_argument("--bins", type=int, default=420)
    parser.add_argument("--dpi", type=int, default=350)
    parser.add_argument("--title", default=None)
    parser.add_argument(
        "--sky-up",
        action="store_true",
        help="mirror-space diagnostic only: project global +z/up onto the mirror plane",
    )
    return parser.parse_args()


def require_columns(df: pd.DataFrame, columns: set[str], path: Path) -> None:
    missing = columns - set(df.columns)
    if missing:
        raise SystemExit(f"{path} is missing required columns: {', '.join(sorted(missing))}")


def plot_hist_with_marked(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    blocked: np.ndarray,
    xlabel: str,
    ylabel: str,
    bins: int,
) -> None:
    unblocked = ~blocked
    if not np.any(unblocked):
        raise SystemExit("no unblocked hits selected")

    hist, xedges, yedges = np.histogram2d(x[unblocked], y[unblocked], bins=max(1, bins))
    masked = np.ma.masked_where(hist.T <= 0, hist.T)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad((1.0, 1.0, 1.0, 0.0))
    positive = hist[hist > 0]
    norm = LogNorm(vmin=max(1.0, float(positive.min())), vmax=float(positive.max())) if positive.size else None
    image = ax.imshow(
        masked,
        origin="lower",
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        aspect="auto",
    )

    if np.any(blocked):
        ax.scatter(
            x[blocked],
            y[blocked],
            s=1.1,
            c="#d7301f",
            marker=".",
            alpha=0.48,
            linewidths=0,
            rasterized=True,
            label="Marked obstructed",
            zorder=4,
        )

    ax.scatter(
        float(np.mean(x[unblocked])),
        float(np.mean(y[unblocked])),
        marker="x",
        s=85,
        linewidths=2.5,
        color="yellow",
        label="Unblocked centroid",
        zorder=5,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.16, linewidth=0.6)
    cbar = ax.figure.colorbar(image, ax=ax, pad=0.02)
    cbar.set_label("Unblocked count / bin")


def main() -> None:
    args = parse_args()
    hits_path = Path(args.hits_csv)
    df = pd.read_csv(hits_path)
    require_columns(
        df,
        {"hit_mirror", "hit_surface", "obstruction_blocked"},
        hits_path,
    )
    df = df[(df["hit_mirror"] == 1) & (df["hit_surface"] == 1)].copy()
    if df.empty:
        raise SystemExit("no rows with hit_mirror=1 and hit_surface=1")

    blocked = df["obstruction_blocked"].to_numpy(dtype=int) != 0
    mirror_basis = None
    oriented = False
    if args.sky_up:
        if not args.config:
            raise SystemExit("--sky-up requires --config")
        if args.space == "spot":
            raise SystemExit("whiteboard plots always use stored physical (u, v); --sky-up is mirror-only")
    if args.space == "spot":
        require_columns(df, {"u_m", "v_m"}, hits_path)
        x = df["u_m"].to_numpy(dtype=float) * 1000.0
        y = df["v_m"].to_numpy(dtype=float) * 1000.0
        xlabel = "LACT focal-plane u [mm]"
        ylabel = "LACT focal-plane v [mm]"
        default_title = "Outer-ring parallel spot with marked obstruction"
    else:
        required = {"mirror_x", "mirror_y"}
        if args.sky_up:
            required.add("mirror_z")
        require_columns(df, required, hits_path)
        if args.sky_up:
            mirror_basis = mirror_display_basis(Path(args.config).resolve())
            projected = project_points_3d(
                df[["mirror_x", "mirror_y", "mirror_z"]].to_numpy(dtype=float),
                mirror_basis,
            )
            x = projected[:, 0] * 1000.0
            y = projected[:, 1] * 1000.0
            oriented = True
        else:
            x = df["mirror_x"].to_numpy(dtype=float) * 1000.0
            y = df["mirror_y"].to_numpy(dtype=float) * 1000.0
        xlabel = "Mirror display x [mm]" if oriented else "Mirror hit x [mm]"
        ylabel = "Mirror display y [mm] (global +z/up)" if oriented else "Mirror hit y [mm]"
        default_title = "Outer-ring mirror hit distribution with marked obstruction"

    fig, ax = plt.subplots(figsize=(7.4, 6.5), dpi=args.dpi)
    plot_hist_with_marked(ax, x, y, blocked, xlabel, ylabel, args.bins)
    ax.set_title(args.title or default_title)
    ax.set_aspect("equal", adjustable="box")

    if args.space == "mirror" and args.overlay_facets:
        if not args.config:
            raise SystemExit("--overlay-facets requires --config")
        outlines = projected_facet_outlines(Path(args.config).resolve(), display_basis=mirror_basis)
        ax.add_collection(LineCollection(outlines, colors="black", linewidths=0.75, alpha=0.9))

    n_total = len(df)
    n_blocked = int(np.count_nonzero(blocked))
    notes = [
        f"saved surface hits = {n_total}",
        f"marked obstructed = {n_blocked} ({n_blocked / n_total:.3%})",
        f"unblocked = {n_total - n_blocked}",
    ]
    if oriented:
        notes.append("display +y = global +z/up")
    note_x, note_y = (0.34, 0.86) if args.space == "mirror" else (0.02, 0.02)
    note_va = "top" if args.space == "mirror" else "bottom"
    ax.text(
        note_x,
        note_y,
        "\n".join(notes),
        transform=ax.transAxes,
        va=note_va,
        ha="left",
        fontsize=9,
        bbox=dict(facecolor="white", edgecolor="0.75", alpha=0.92, boxstyle="round,pad=0.28"),
    )
    ax.legend(loc="upper right", frameon=True, framealpha=0.95, edgecolor="black")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    print(f"Saved {args.space} obstruction-marked plot = {output}")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
