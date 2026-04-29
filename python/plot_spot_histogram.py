#!/usr/bin/env python3
import argparse
import math

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def containment_radii(center, points, fractions=(0.68, 0.80)):
    r = np.linalg.norm(points - center, axis=1)
    return [float(np.quantile(r, f)) for f in fractions]


def weighted_centroid(x, y, w):
    sw = np.sum(w)
    if sw <= 0.0:
        return float(np.mean(x)), float(np.mean(y))
    return float(np.sum(w * x) / sw), float(np.sum(w * y) / sw)


def main():
    parser = argparse.ArgumentParser(description="Plot optical output-plane spot as a 2D histogram.")
    parser.add_argument("csv", help="surface hit CSV from run_optical_sim")
    parser.add_argument("--output", default="spot_histogram.png", help="output PNG/PDF path")
    parser.add_argument(
        "--binsize-mm",
        type=float,
        default=None,
        help="histogram bin size in mm; default chooses an adaptive value from data range and hit count",
    )
    parser.add_argument(
        "--max-bins",
        type=int,
        default=320,
        help="upper cap for automatically chosen bin count along the longest axis",
    )
    parser.add_argument("--focal-length-m", type=float, default=8.0, help="plate-scale focal length")
    parser.add_argument("--dpi", type=int, default=350)
    parser.add_argument("--title", default="Optical output-plane spot")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    df = df[df["hit_surface"] == 1].copy()
    if len(df) == 0:
        raise SystemExit("No output-plane hits in file.")

    x_mm = df["u_m"].to_numpy(float) * 1000.0
    y_mm = df["v_m"].to_numpy(float) * 1000.0
    w = df["weight"].to_numpy(float) * df["relative_efficiency"].to_numpy(float)

    cx_mm, cy_mm = weighted_centroid(x_mm, y_mm, w)
    points_mm = np.column_stack([x_mm, y_mm])
    r68_mm, r80_mm = containment_radii(np.array([cx_mm, cy_mm]), points_mm)
    rms_mm = math.sqrt(np.sum(w * ((x_mm - cx_mm) ** 2 + (y_mm - cy_mm) ** 2)) / np.sum(w))

    xmin, xmax = float(np.min(x_mm)), float(np.max(x_mm))
    ymin, ymax = float(np.min(y_mm)), float(np.max(y_mm))
    span = max(xmax - xmin, ymax - ymin, 1e-9)
    if args.binsize_mm is None:
        target_bins = int(np.clip(np.sqrt(len(df)) / 2.0, 140, args.max_bins))
        binsize_mm = span / target_bins
        bins = target_bins
    else:
        binsize_mm = args.binsize_mm
        bins = max(20, int(math.ceil(span / binsize_mm)))

    pad = 0.08 * span
    x_range = (0.5 * (xmin + xmax) - 0.5 * span - pad,
               0.5 * (xmin + xmax) + 0.5 * span + pad)
    y_range = (0.5 * (ymin + ymax) - 0.5 * span - pad,
               0.5 * (ymin + ymax) + 0.5 * span + pad)

    hist, xedges, yedges = np.histogram2d(x_mm, y_mm, bins=bins, range=[x_range, y_range])

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_under("white")

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

    fig, ax = plt.subplots(figsize=(6.5, 5.8))
    image = ax.imshow(
        hist.T,
        interpolation="nearest",
        cmap=cmap,
        origin="lower",
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        vmin=0.01,
        aspect="equal",
    )

    cbar = fig.colorbar(image, ax=ax, pad=0.03)
    cbar.set_label("Count / bin")

    ax.scatter(cx_mm, cy_mm, color="yellow", marker="x", s=64, linewidths=2.0, label="Centroid")
    circle68 = patches.Circle(
        (cx_mm, cy_mm),
        r68_mm,
        fill=False,
        edgecolor="orange",
        linestyle="-",
        linewidth=2.2,
        label=f"R68 = {r68_mm:.2f} mm",
    )
    circle80 = patches.Circle(
        (cx_mm, cy_mm),
        r80_mm,
        fill=False,
        edgecolor="white",
        linestyle="--",
        linewidth=1.6,
        label=f"R80 = {r80_mm:.2f} mm",
    )
    ax.add_patch(circle68)
    ax.add_patch(circle80)

    def mm_to_deg(mm):
        return np.degrees(np.asarray(mm, dtype=float) / (args.focal_length_m * 1000.0))

    def deg_to_mm(deg):
        return np.radians(np.asarray(deg, dtype=float)) * args.focal_length_m * 1000.0

    ax.set_xlabel("u [mm]")
    ax.set_ylabel("v [mm]")
    ax.set_title(args.title)
    secax_x = ax.secondary_xaxis("top", functions=(mm_to_deg, deg_to_mm))
    secax_x.set_xlabel("u [deg]")
    secax_y = ax.secondary_yaxis("right", functions=(mm_to_deg, deg_to_mm))
    secax_y.set_ylabel("v [deg]")

    ax.grid(True, color="white", alpha=0.22, linewidth=0.5)
    ax.legend(frameon=True, framealpha=1.0, edgecolor="black", loc="best")

    text = (
        f"N = {len(df)}\n"
        f"Centroid = ({cx_mm:.2f}, {cy_mm:.2f}) mm\n"
        f"RMS = {rms_mm:.2f} mm"
    )
    ax.text(
        0.02,
        0.02,
        text,
        transform=ax.transAxes,
        fontsize=8,
        va="bottom",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.82", alpha=0.92),
    )

    fig.tight_layout()
    fig.savefig(args.output, bbox_inches="tight")

    print(f"N hits = {len(df)}")
    print(f"Histogram span [mm] = {span:.6f}")
    print(f"Histogram bins = {bins}")
    print(f"Histogram bin size [mm] = {binsize_mm:.6f}")
    print(f"Centroid [mm] = ({cx_mm:.6f}, {cy_mm:.6f})")
    print(f"R68 [mm] = {r68_mm:.6f}")
    print(f"R80 [mm] = {r80_mm:.6f}")
    print(f"Weighted RMS about centroid [mm] = {rms_mm:.6f}")
    print(f"Saved plot = {args.output}")

    if "agg" not in plt.get_backend().lower():
        plt.show()


if __name__ == "__main__":
    main()
