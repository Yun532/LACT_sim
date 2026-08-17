#!/usr/bin/env python3
import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

from config_io import (
    apply_telescope_frame_to_facets,
    expand_component_config,
    load_facets_from_config,
    point_to_global,
    telescope_frame_from_config,
)


def parse_vec3(text, fallback):
    if text is None:
        return np.array(fallback, dtype=float)
    return np.array([float(x.strip()) for x in text.split(",")], dtype=float)


def local_frame(normal):
    n = normal / np.linalg.norm(normal)
    ref = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(ref, n)
    u = u / np.linalg.norm(u)
    v = np.cross(n, u)
    v = v / np.linalg.norm(v)
    return u, v


def facet_polygon_xy(row):
    center = np.array(row["center"], dtype=float)
    normal = np.array(row["normal"], dtype=float)
    normal = normal / np.linalg.norm(normal)
    if "u_axis" in row and "v_axis" in row:
        u = np.array(row["u_axis"], dtype=float)
        v = np.array(row["v_axis"], dtype=float)
        u = u / np.linalg.norm(u)
        v = v / np.linalg.norm(v)
    else:
        u, v = local_frame(normal)

    size = float(row["size1"])
    shape = row["shape"].strip().lower()
    rotation = float(row.get("rotation", 0.0) or 0.0)

    if shape == "hexagon":
        a = 0.5 * size
        local = np.array([
            [a, a / math.sqrt(3.0)],
            [0.0, 2.0 * a / math.sqrt(3.0)],
            [-a, a / math.sqrt(3.0)],
            [-a, -a / math.sqrt(3.0)],
            [0.0, -2.0 * a / math.sqrt(3.0)],
            [a, -a / math.sqrt(3.0)],
        ])
    elif shape == "square":
        h = 0.5 * size
        local = np.array([[h, h], [-h, h], [-h, -h], [h, -h]])
    else:
        angles = np.linspace(0, 2.0 * math.pi, 72, endpoint=False)
        local = np.column_stack([size * np.cos(angles), size * np.sin(angles)])

    if rotation != 0.0:
        c = math.cos(rotation)
        s = math.sin(rotation)
        R = np.array([[c, -s], [s, c]])
        local = local @ R.T

    points = np.array([center + x * u + y * v for x, y in local])
    return points[:, :2], center[:2]


def main():
    parser = argparse.ArgumentParser(description="Publication-style top view of mirror facets.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--elevation-deg", type=float, default=None, help="override telescope/mirror elevation state")
    parser.add_argument(
        "--frame",
        choices=("local", "global"),
        default="local",
        help="plot in telescope-local coordinates for mirror checks, or global array coordinates",
    )
    parser.add_argument("--output", default="mirror_top_view.png")
    parser.add_argument("--dpi", type=int, default=350)
    args = parser.parse_args()

    cfg_path = Path(args.config).resolve()
    cfg, _ = expand_component_config(cfg_path)
    if args.elevation_deg is not None:
        cfg["telescope.pointing_el_deg"] = str(args.elevation_deg)
        cfg["mirror.series_elevation_deg"] = str(args.elevation_deg)
    frame = telescope_frame_from_config(cfg)
    plane_point = parse_vec3(cfg.get("output.plane_point"), [0, 0, 0])
    facets = load_facets_from_config(cfg_path, cfg)
    if args.frame == "global":
        plane_point = point_to_global(plane_point, frame)
        facets = apply_telescope_frame_to_facets(facets, frame)

    patches = []
    centers = []
    for row in facets:
        poly_xy, center_xy = facet_polygon_xy(row)
        patches.append(Polygon(poly_xy, closed=True))
        centers.append(center_xy)

    centers = np.array(centers)
    radial = np.linalg.norm(centers, axis=1)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "figure.dpi": args.dpi,
        "savefig.dpi": args.dpi,
    })

    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    collection = PatchCollection(
        patches,
        cmap="viridis",
        edgecolor=(0.12, 0.12, 0.12, 0.75),
        linewidth=0.55,
        alpha=0.94,
    )
    collection.set_array(radial)
    ax.add_collection(collection)

    ax.scatter(centers[:, 0], centers[:, 1], s=5, c="black", alpha=0.45, zorder=3)
    ax.scatter([plane_point[0]], [plane_point[1]], s=32, c="#B2182B", zorder=4, label="Output plane center")

    all_xy = np.concatenate([p.get_xy() for p in patches])
    min_xy = all_xy.min(axis=0)
    max_xy = all_xy.max(axis=0)
    span = max(max_xy - min_xy)
    center = 0.5 * (min_xy + max_xy)
    half = 0.56 * span
    ax.set_xlim(center[0] - half, center[0] + half)
    ax.set_ylim(center[1] - half, center[1] + half)
    ax.set_aspect("equal", adjustable="box")

    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(f"Mirror facet layout, {args.frame} top view")
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.legend(loc="upper right", frameon=True, framealpha=0.92)

    cbar = fig.colorbar(collection, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("facet radial position [m]")

    ax.text(
        0.02,
        0.02,
        (
            f"N facets = {len(patches)}\n"
            f"hex flat-to-flat = 0.80 m\n"
            f"elevation = {float(cfg.get('mirror.series_elevation_deg', cfg.get('telescope.pointing_el_deg', 90.0))):.1f} deg\n"
            f"frame = {args.frame}"
        ),
        transform=ax.transAxes,
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.82", alpha=0.9),
    )

    fig.tight_layout()
    fig.savefig(args.output, bbox_inches="tight")
    print(f"Saved top view = {args.output}")

    if "agg" not in plt.get_backend().lower():
        plt.show()


if __name__ == "__main__":
    main()
