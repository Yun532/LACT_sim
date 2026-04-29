#!/usr/bin/env python3
import argparse
import csv
import math
import subprocess
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config_io import expand_component_config


def parse_float_list(text):
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def weighted_centroid(x, y, w):
    sw = np.sum(w)
    if sw <= 0.0:
        return float(np.mean(x)), float(np.mean(y))
    return float(np.sum(w * x) / sw), float(np.sum(w * y) / sw)


def containment_radii(center, points, fractions=(0.68, 0.80)):
    r = np.linalg.norm(points - center, axis=1)
    return [float(np.quantile(r, f)) for f in fractions]


def compute_metrics(csv_path):
    df = pd.read_csv(csv_path)
    df = df[df["hit_surface"] == 1].copy()
    if len(df) == 0:
        raise RuntimeError(f"No output-plane hits in {csv_path}")

    x_mm = df["u_m"].to_numpy(float) * 1000.0
    y_mm = df["v_m"].to_numpy(float) * 1000.0
    w = df["weight"].to_numpy(float) * df["relative_efficiency"].to_numpy(float)

    cx_mm, cy_mm = weighted_centroid(x_mm, y_mm, w)
    pts = np.column_stack([x_mm, y_mm])
    centered = pts - np.array([cx_mm, cy_mm])
    r68_mm, r80_mm = containment_radii(np.array([cx_mm, cy_mm]), pts)
    rms_mm = math.sqrt(np.sum(w * (centered[:, 0] ** 2 + centered[:, 1] ** 2)) / np.sum(w))

    return {
        "n_hits": int(len(df)),
        "centroid_x_mm": cx_mm,
        "centroid_y_mm": cy_mm,
        "centroid_r_mm": math.hypot(cx_mm, cy_mm),
        "r68_mm": r68_mm,
        "r80_mm": r80_mm,
        "rms_mm": rms_mm,
        "x_mm": x_mm,
        "y_mm": y_mm,
        "w": w,
        "centered": centered,
    }


def make_spot_grid(results, output_path, focal_length_m, dpi):
    all_centered = np.concatenate([r["centered"] for r in results.values()], axis=0)
    max_abs = max(np.max(np.abs(all_centered[:, 0])), np.max(np.abs(all_centered[:, 1])), 1.0)
    span = 2.2 * max_abs
    bins = 220
    rng = [[-span / 2.0, span / 2.0], [-span / 2.0, span / 2.0]]

    fig, axes = plt.subplots(2, 5, figsize=(15.2, 6.4), constrained_layout=True)
    axes = axes.ravel()
    vmax = 0.0
    hists = []
    for elevation in sorted(results):
        centered = results[elevation]["centered"]
        hist, xedges, yedges = np.histogram2d(centered[:, 0], centered[:, 1], bins=bins, range=rng)
        hists.append((elevation, hist, xedges, yedges))
        vmax = max(vmax, float(np.max(hist)))

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_under("white")

    for ax, (elevation, hist, xedges, yedges) in zip(axes, hists):
        image = ax.imshow(
            hist.T,
            interpolation="nearest",
            cmap=cmap,
            origin="lower",
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
            vmin=0.01,
            vmax=vmax,
            aspect="equal",
        )
        r = results[elevation]
        circle68 = patches.Circle((0.0, 0.0), r["r68_mm"], fill=False, edgecolor="orange", linewidth=1.4)
        ax.add_patch(circle68)
        ax.scatter([0.0], [0.0], color="yellow", marker="x", s=28, linewidths=1.2)
        ax.set_title(f"El = {int(elevation)} deg", fontsize=10)
        ax.grid(True, color="white", alpha=0.15, linewidth=0.45)
        ax.text(
            0.03,
            0.03,
            (
                f"N={r['n_hits']}\n"
                f"dC=({r['centroid_x_mm']:.1f},{r['centroid_y_mm']:.1f}) mm\n"
                f"RMS={r['rms_mm']:.2f} mm"
            ),
            transform=ax.transAxes,
            fontsize=7,
            va="bottom",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="0.82", alpha=0.92),
        )
        ax.set_xlabel("u - centroid [mm]")
        ax.set_ylabel("v - centroid [mm]")

    cbar = fig.colorbar(image, ax=axes.tolist(), shrink=0.94, pad=0.02)
    cbar.set_label("Count / bin")
    fig.suptitle(
        f"1229 elevation deformation: parallel-light spots (f = {focal_length_m:.2f} m)",
        fontsize=12,
        y=1.02,
    )
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def make_metrics_plot(results, output_path, dpi):
    elevations = np.array(sorted(results), dtype=float)
    centroid_r = np.array([results[e]["centroid_r_mm"] for e in elevations])
    rms = np.array([results[e]["rms_mm"] for e in elevations])
    r68 = np.array([results[e]["r68_mm"] for e in elevations])

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.1), constrained_layout=True)
    axes[0].plot(elevations, centroid_r, marker="o", lw=1.8)
    axes[0].set_xlabel("Telescope elevation [deg]")
    axes[0].set_ylabel("Centroid offset [mm]")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(elevations, rms, marker="o", lw=1.8, label="RMS")
    axes[1].plot(elevations, r68, marker="s", lw=1.4, label="R68")
    axes[1].set_xlabel("Telescope elevation [deg]")
    axes[1].set_ylabel("Spot size [mm]")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(frameon=True)

    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Run a parallel-light elevation scan and plot spot summaries.")
    parser.add_argument("--config", default="configs/official_tests/deformation_parallel_whiteboard.cfg")
    parser.add_argument("--run-binary", default="build_check/run_optical_sim")
    parser.add_argument("--elevations", default="0,10,20,30,40,50,60,70,80,90")
    parser.add_argument("--n-bunches", type=int, default=100000)
    parser.add_argument("--beam-theta-deg", type=float, default=0.0)
    parser.add_argument("--beam-phi-deg", type=float, default=0.0)
    parser.add_argument("--output-dir", default="run_logs/elevation_parallel_scan")
    parser.add_argument("--dpi", type=int, default=250)
    args = parser.parse_args()

    base_cfg_path = Path(args.config).resolve()
    expanded_cfg, component_paths = expand_component_config(base_cfg_path)
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    telescope_cfg = component_paths.get("telescope")
    mirror_cfg = component_paths.get("mirror")
    source_cfg = component_paths.get("source")
    output_cfg = component_paths.get("output")
    efficiency_cfg = component_paths.get("efficiency")
    error_cfg = component_paths.get("error")

    if telescope_cfg is None or mirror_cfg is None or source_cfg is None:
        raise SystemExit("The scan config must resolve telescope.config, mirror.config, and source.config")

    focal_length_m = float(expanded_cfg.get("telescope.focal_length_m", "8.0"))
    repo_root = Path(__file__).resolve().parents[1]
    run_binary = Path(args.run_binary).resolve()
    elevations = parse_float_list(args.elevations)

    metrics_rows = []
    results = {}

    for elevation in elevations:
        tag = format(int(elevation)) if float(elevation).is_integer() else str(elevation).replace(".", "p")
        temp_cfg = output_dir / f"run_el_{tag}.cfg"
        hit_csv = output_dir / f"hits_el_{tag}.csv"
        log_path = output_dir / f"run_el_{tag}.log"

        lines = [
            f"telescope.config={telescope_cfg}",
            f"mirror.config={mirror_cfg}",
            f"source.config={source_cfg}",
            f"telescope.pointing_el_deg={elevation}",
            f"source.n_bunches={args.n_bunches}",
            f"output.csv={hit_csv}",
        ]
        if output_cfg is not None:
            lines.append(f"output.config={output_cfg}")
        if efficiency_cfg is not None:
            lines.append(f"efficiency.config={efficiency_cfg}")
        if error_cfg is not None:
            lines.append(f"error.config={error_cfg}")
        if "propagation.speed_of_light_m_per_ns" in expanded_cfg:
            lines.append(
                f"propagation.speed_of_light_m_per_ns={expanded_cfg['propagation.speed_of_light_m_per_ns']}"
            )
        if abs(args.beam_theta_deg) > 0.0 or abs(args.beam_phi_deg) > 0.0:
            lines.append(f"source.beam_theta_deg={args.beam_theta_deg}")
            lines.append(f"source.beam_phi_deg={args.beam_phi_deg}")

        temp_cfg.write_text("\n".join(lines) + "\n")

        completed = subprocess.run(
            [str(run_binary), str(temp_cfg)],
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=True,
        )
        log_path.write_text(completed.stdout + completed.stderr)

        metrics = compute_metrics(hit_csv)
        results[elevation] = metrics
        metrics_rows.append(
            {
                "elevation_deg": elevation,
                "n_hits": metrics["n_hits"],
                "centroid_x_mm": metrics["centroid_x_mm"],
                "centroid_y_mm": metrics["centroid_y_mm"],
                "centroid_r_mm": metrics["centroid_r_mm"],
                "r68_mm": metrics["r68_mm"],
                "r80_mm": metrics["r80_mm"],
                "rms_mm": metrics["rms_mm"],
                "hit_csv": str(hit_csv),
                "log": str(log_path),
            }
        )
        print(
            f"elevation={elevation:5.1f} deg  hits={metrics['n_hits']:7d}  "
            f"centroid_r={metrics['centroid_r_mm']:8.2f} mm  rms={metrics['rms_mm']:7.2f} mm"
        )

    summary_csv = output_dir / "elevation_scan_metrics.csv"
    with summary_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(metrics_rows[0].keys()))
        writer.writeheader()
        writer.writerows(metrics_rows)

    grid_png = output_dir / "elevation_spot_grid.png"
    metrics_png = output_dir / "elevation_metrics.png"
    make_spot_grid(results, grid_png, focal_length_m, args.dpi)
    make_metrics_plot(results, metrics_png, args.dpi)

    print(f"Saved summary CSV = {summary_csv}")
    print(f"Saved spot grid = {grid_png}")
    print(f"Saved metrics plot = {metrics_png}")


if __name__ == "__main__":
    main()
