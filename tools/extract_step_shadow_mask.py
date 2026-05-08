#!/usr/bin/env python3
"""Extract a lightweight 2D obstruction mask from a STEP support model.

This is intentionally a projection tool, not a full CAD kernel. It samples
`CARTESIAN_POINT` entities from the STEP file, maps the support-model axes to
the LACT local aperture frame, rasterizes the projected points, and dilates the
result by a configurable number of cells.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np


POINT_RE = re.compile(r"CARTESIAN_POINT\('',\(([^)]*)\)\)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", type=Path, help="input STEP/STP file")
    parser.add_argument("--output", type=Path, required=True, help="output mask CSV")
    parser.add_argument("--preview", type=Path, help="optional preview PNG")
    parser.add_argument("--extent-m", type=float, default=4.5, help="half-width of mask grid")
    parser.add_argument("--cell-size-m", type=float, default=0.02, help="mask cell size")
    parser.add_argument("--dilate-cells", type=int, default=2, help="binary dilation cells")
    parser.add_argument("--plane-z-m", type=float, default=-16.0, help="LACT local z plane")
    return parser.parse_args()


def parse_step_point(line: str) -> tuple[float, float, float] | None:
    match = POINT_RE.search(line)
    if not match:
        return None
    vals = [float(v.replace("D", "E")) for v in match.group(1).split(",")]
    if len(vals) != 3:
        return None
    return vals[0], vals[1], vals[2]


def dilate_mask(mask: np.ndarray, iterations: int) -> np.ndarray:
    out = mask.copy()
    for _ in range(max(0, iterations)):
        padded = np.pad(out, 1, mode="constant", constant_values=False)
        expanded = np.zeros_like(out)
        for dy in range(3):
            for dx in range(3):
                expanded |= padded[dy : dy + out.shape[0], dx : dx + out.shape[1]]
        out = expanded
    return out


def main() -> None:
    args = parse_args()
    half = float(args.extent_m)
    cell = float(args.cell_size_m)
    if half <= 0 or cell <= 0:
        raise SystemExit("--extent-m and --cell-size-m must be positive")

    x_min = -half
    y_min = -half
    nx = int(round(2 * half / cell))
    ny = nx
    mask = np.zeros((ny, nx), dtype=bool)
    n_points = 0
    n_used = 0

    with args.step.open(errors="ignore") as handle:
        for line in handle:
            point = parse_step_point(line)
            if point is None:
                continue
            n_points += 1
            step_x_mm, _step_y_mm, step_z_mm = point
            # STEP model axes for this support:
            #   STEP X -> LACT local x
            #   STEP Z -> LACT local y
            #   STEP Y -> LACT local z, projected onto plane_z_m
            x = step_x_mm * 0.001
            y = step_z_mm * 0.001
            if not (x_min <= x < -x_min and y_min <= y < -y_min):
                continue
            ix = int((x - x_min) / cell)
            iy = int((y - y_min) / cell)
            mask[iy, ix] = True
            n_used += 1

    if args.dilate_cells > 0:
        mask = dilate_mask(mask, args.dilate_cells)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        handle.write("# format=LACT_obstruction_mask_v1\n")
        handle.write("# source_step=%s\n" % args.step)
        handle.write("# axis_mapping=STEP_X_to_local_x,STEP_Z_to_local_y,STEP_Y_to_local_z\n")
        handle.write("# x_min_m=%.10g\n" % x_min)
        handle.write("# y_min_m=%.10g\n" % y_min)
        handle.write("# cell_size_m=%.10g\n" % cell)
        handle.write("# nx=%d\n" % nx)
        handle.write("# ny=%d\n" % ny)
        handle.write("# plane_z_m=%.10g\n" % args.plane_z_m)
        handle.write("# input_cartesian_points=%d\n" % n_points)
        handle.write("# projected_points_in_grid=%d\n" % n_used)
        handle.write("# dilate_cells=%d\n" % args.dilate_cells)
        writer = csv.writer(handle)
        writer.writerow(["ix", "iy"])
        ys, xs = np.nonzero(mask)
        for ix, iy in zip(xs, ys):
            writer.writerow([int(ix), int(iy)])

    if args.preview:
        import matplotlib.pyplot as plt

        args.preview.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(6.0, 6.0))
        ax.imshow(
            mask,
            origin="lower",
            extent=[x_min, -x_min, y_min, -y_min],
            cmap="gray_r",
            interpolation="nearest",
        )
        circle = plt.Circle((0, 0), 4.0, fill=False, color="tab:red", lw=1.2)
        ax.add_patch(circle)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Local x (m)")
        ax.set_ylabel("Local y (m)")
        ax.set_title("Projected Support Obstruction Mask")
        ax.text(
            0.02,
            0.02,
            f"blocked cells={int(mask.sum())}",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=9,
        )
        fig.tight_layout()
        fig.savefig(args.preview, dpi=250)
        plt.close(fig)

    aperture = np.fromfunction(
        lambda j, i: (x_min + (i + 0.5) * cell) ** 2
        + (y_min + (j + 0.5) * cell) ** 2
        <= 4.0**2,
        (ny, nx),
    )
    frac = float(np.logical_and(mask, aperture).sum()) / float(aperture.sum())
    print(f"wrote {args.output}")
    print(f"input_cartesian_points={n_points}")
    print(f"projected_points_in_grid={n_used}")
    print(f"blocked_cells={int(mask.sum())}")
    print(f"aperture_blocked_fraction_estimate={frac:.6f}")


if __name__ == "__main__":
    main()
