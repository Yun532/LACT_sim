#!/usr/bin/env python3
"""Plot only a LACT obstruction-primitives CSV as a standalone 3D layout."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

from config_io import expand_component_config, telescope_frame_from_config
from plot_optical_layout_3d import (
    box_corners,
    box_edges,
    box_faces,
    load_obstruction_primitives,
    regular_prism_edges,
    regular_prism_faces,
    set_equal_3d,
    transform_points,
)


def collect_primitive_points(primitives, frame):
    points = []
    for primitive in primitives:
        points.extend(transform_points([primitive["p0"], primitive["p1"]], frame))
        ptype = primitive["type"]
        if ptype in {"box", "aabb"}:
            points.extend(transform_points(box_corners(primitive["center"], primitive["half"]), frame))
            if primitive["hole_radius"] > 0 and primitive["hole_sides"] >= 3:
                for edge in regular_prism_edges(
                    primitive["center"],
                    primitive["hole_radius"],
                    primitive["half"][2] * 2.0,
                    primitive["hole_rotation"],
                    primitive["hole_sides"],
                ):
                    points.extend(transform_points(edge, frame))
        elif ptype == "polygon_prism":
            for face in regular_prism_faces(
                primitive["center"],
                primitive["radius"],
                primitive["height"],
                primitive["rotation"],
                primitive["sides"],
            ):
                points.extend(transform_points(face, frame))
    if not points:
        raise SystemExit("no drawable primitive points found")
    return np.asarray(points, dtype=float)


def draw_primitives(ax, primitives, frame, label_names=False):
    cylinder_segments = []
    cylinder_widths = []
    solid_faces = []
    wire_segments = []
    wire_colors = []
    wire_widths = []
    label_points = []

    for primitive in primitives:
        ptype = primitive["type"]
        name = primitive["name"]
        role = primitive["role"]
        if ptype == "cylinder":
            segment = transform_points([primitive["p0"], primitive["p1"]], frame)
            cylinder_segments.append(segment)
            cylinder_widths.append(max(0.8, primitive["radius"] * 24.0))
            label_points.append((name, np.mean(segment, axis=0)))
        elif ptype in {"box", "aabb"}:
            faces = [transform_points(face, frame) for face in box_faces(primitive["center"], primitive["half"])]
            solid_faces.extend(faces)
            for edge in box_edges(primitive["center"], primitive["half"]):
                wire_segments.append(transform_points(edge, frame))
                wire_colors.append((0.10, 0.44, 0.22, 0.95))
                wire_widths.append(1.2)
            if primitive["hole_radius"] > 0 and primitive["hole_sides"] >= 3:
                for edge in regular_prism_edges(
                    primitive["center"],
                    primitive["hole_radius"],
                    primitive["half"][2] * 2.0,
                    primitive["hole_rotation"],
                    primitive["hole_sides"],
                ):
                    wire_segments.append(transform_points(edge, frame))
                    wire_colors.append((0.78, 0.05, 0.14, 1.0))
                    wire_widths.append(1.5)
            label_points.append((name, transform_points([primitive["center"]], frame)[0]))
        elif ptype == "polygon_prism":
            faces = [
                transform_points(face, frame)
                for face in regular_prism_faces(
                    primitive["center"],
                    primitive["radius"],
                    primitive["height"],
                    primitive["rotation"],
                    primitive["sides"],
                )
            ]
            solid_faces.extend(faces)
            for edge in regular_prism_edges(
                primitive["center"],
                primitive["radius"],
                primitive["height"],
                primitive["rotation"],
                primitive["sides"],
            ):
                wire_segments.append(transform_points(edge, frame))
                wire_colors.append((0.0, 0.58, 0.68, 1.0))
                wire_widths.append(1.3)
            label_points.append((name, transform_points([primitive["center"]], frame)[0]))
        else:
            print(f"warning: unsupported primitive type skipped: {ptype} ({name})", file=sys.stderr)
            continue

        if role:
            ax.scatter([], [], [], label=role)

    if cylinder_segments:
        ax.add_collection3d(
            Line3DCollection(
                cylinder_segments,
                colors=(0.03, 0.03, 0.03, 0.76),
                linewidths=cylinder_widths,
                label="support cylinders",
            )
        )
    if solid_faces:
        ax.add_collection3d(
            Poly3DCollection(
                solid_faces,
                facecolors=(0.08, 0.08, 0.08, 0.14),
                edgecolors=(0.03, 0.03, 0.03, 0.32),
                linewidths=0.35,
            )
        )
    if wire_segments:
        ax.add_collection3d(
            Line3DCollection(
                wire_segments,
                colors=wire_colors,
                linewidths=wire_widths,
                label="box/prism edges",
            )
        )
    if label_names:
        for name, point in label_points:
            ax.text(point[0], point[1], point[2], name, fontsize=6.5, color="0.12")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("primitives_csv", help="obstruction primitives CSV")
    parser.add_argument("--config", help="optional cfg; if provided, plot after telescope pointing transform")
    parser.add_argument("--output", required=True, help="output PNG/PDF path")
    parser.add_argument("--dpi", type=int, default=320)
    parser.add_argument("--view", default="30,-58", help="matplotlib 3D view as elev,azim")
    parser.add_argument("--label-names", action="store_true", help="label each primitive name")
    parser.add_argument("--title", default="Obstruction primitives")
    return parser.parse_args()


def main():
    args = parse_args()
    frame = None
    if args.config:
        cfg, _ = expand_component_config(Path(args.config).resolve())
        frame = telescope_frame_from_config(cfg)

    primitives = load_obstruction_primitives(Path(args.primitives_csv).resolve())
    if not primitives:
        raise SystemExit("no obstruction primitives loaded")

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 12,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    })
    fig = plt.figure(figsize=(7.2, 6.4), dpi=args.dpi)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("white")

    draw_primitives(ax, primitives, frame, label_names=args.label_names)
    set_equal_3d(ax, collect_primitive_points(primitives, frame), pad=0.10)
    elev, azim = [float(x) for x in args.view.split(",", 1)]
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    ax.set_title(args.title)
    ax.grid(True, alpha=0.18)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
        axis.pane.set_edgecolor((0.82, 0.82, 0.82, 1.0))

    role_counts = {}
    for primitive in primitives:
        role = primitive["role"] or primitive["type"]
        role_counts[role] = role_counts.get(role, 0) + 1
    ax.text2D(
        0.02,
        0.02,
        "\n".join([f"N primitives = {len(primitives)}"] + [f"{k}: {v}" for k, v in sorted(role_counts.items())]),
        transform=ax.transAxes,
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.82", alpha=0.88),
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    print(f"Saved obstruction 3D layout = {output}")


if __name__ == "__main__":
    main()
