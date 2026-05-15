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


def _normalize(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm <= 0:
        raise ValueError("cannot normalize a zero-length vector")
    return vector / norm


def mirror_display_basis(config_path: Path) -> dict[str, np.ndarray | bool]:
    """Return a 3D mirror-plane display basis with global +z shown upward.

    The mirror/facet coordinates are 3D after telescope pointing is applied.
    For a tilted telescope, simply plotting global x-y discards the component
    that actually corresponds to sky/up. This basis projects global +z onto the
    mirror transverse plane and uses that as display +y.
    """

    cfg, _ = expand_component_config(config_path)
    frame = telescope_frame_from_config(cfg)
    normal = _normalize(np.asarray(frame["z_axis"], dtype=float))
    global_up = np.array([0.0, 0.0, 1.0])
    display_y = global_up - float(np.dot(global_up, normal)) * normal
    sky_up_projected = True
    if np.linalg.norm(display_y) < 1.0e-10:
        # When the telescope points at zenith, global up is parallel to the
        # optical axis, so its projection onto the mirror plane is undefined.
        # Fall back to the native mirror local +y direction.
        display_y = np.asarray(frame["y_axis"], dtype=float)
        sky_up_projected = False
    display_y = _normalize(display_y)
    display_x = _normalize(np.cross(display_y, normal))
    return {
        "origin": np.asarray(frame["origin"], dtype=float),
        "x_axis": display_x,
        "y_axis": display_y,
        "sky_up_projected": sky_up_projected,
    }


def project_points_3d(points_m: np.ndarray, display_basis: dict[str, np.ndarray | bool]) -> np.ndarray:
    origin = np.asarray(display_basis["origin"], dtype=float)
    x_axis = np.asarray(display_basis["x_axis"], dtype=float)
    y_axis = np.asarray(display_basis["y_axis"], dtype=float)
    rel = points_m - origin
    return np.column_stack((rel @ x_axis, rel @ y_axis))


def projected_facet_outlines(config_path: Path, display_basis=None) -> list[np.ndarray]:
    cfg, _ = expand_component_config(config_path)
    frame = telescope_frame_from_config(cfg)
    facets = apply_telescope_frame_to_facets(load_facets_from_config(config_path, cfg), frame)
    outlines = []
    for facet in facets:
        polygon_3d = aperture_polygon(facet)
        if display_basis is not None:
            polygon = project_points_3d(polygon_3d, display_basis) * 1000.0
        else:
            polygon = polygon_3d[:, :2] * 1000.0
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
    parser.add_argument(
        "--sky-up",
        action="store_true",
        help="with --config, rotate the plot so display +y is global +z/up projected onto the mirror plane",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    hits_path = Path(args.hits_csv)
    df = pd.read_csv(hits_path)
    required = {"hit_mirror", "mirror_x", "mirror_y"}
    if args.sky_up:
        required.add("mirror_z")
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

    display_basis = None
    sky_up_projected = False
    if args.sky_up:
        if not args.config:
            raise SystemExit("--sky-up requires --config")
        display_basis = mirror_display_basis(Path(args.config).resolve())
        projected = project_points_3d(df[["mirror_x", "mirror_y", "mirror_z"]].to_numpy(dtype=float), display_basis)
        x = projected[:, 0] * 1000.0
        y = projected[:, 1] * 1000.0
        sky_up_projected = bool(display_basis["sky_up_projected"])
    else:
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
        outlines = projected_facet_outlines(Path(args.config).resolve(), display_basis=display_basis)
        ax.add_collection(LineCollection(outlines, colors="black", linewidths=0.75, alpha=0.9))

    mean_x = float(np.mean(x))
    mean_y = float(np.mean(y))
    ax.scatter(mean_x, mean_y, marker="x", s=85, linewidths=2.6, color="yellow", label="Mean", zorder=5)
    ax.set_aspect("equal", adjustable="box")
    if args.sky_up:
        ax.set_xlabel("Mirror display x [mm]")
        if sky_up_projected:
            ax.set_ylabel("Mirror display y [mm] (global +z/up)")
        else:
            ax.set_ylabel("Mirror display y [mm] (native mirror +y)")
    else:
        ax.set_xlabel("Mirror hit x [mm]")
        ax.set_ylabel("Mirror hit y [mm]")
    ax.set_title(args.title)
    ax.grid(alpha=0.16, linewidth=0.6)
    ax.legend(loc="upper right", frameon=True, framealpha=0.95, edgecolor="black")
    cbar = fig.colorbar(image, ax=ax, pad=0.02)
    cbar.set_label("Count / bin")

    notes = [f"N = {len(df)}", f"mean = ({mean_x:.1f}, {mean_y:.1f}) mm"]
    if args.overlay_facets:
        notes.insert(1, "facet outlines: mirror-plane projection" if args.sky_up else "facet outlines: x-y projection")
    if args.sky_up:
        if sky_up_projected:
            notes.append("display +y = global +z/up projected on mirror plane")
        else:
            notes.append("display +y = native mirror +y; global +z projection is degenerate")
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
