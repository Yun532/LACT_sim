#!/usr/bin/env python3
"""Plot array-run summaries and per-telescope output-plane spots."""

import argparse
import math
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def as_number(series, default=np.nan):
    return pd.to_numeric(series, errors="coerce").fillna(default)


def weighted_centroid(x, y, w):
    sw = np.sum(w)
    if sw <= 0.0:
        return float(np.mean(x)), float(np.mean(y))
    return float(np.sum(w * x) / sw), float(np.sum(w * y) / sw)


def spot_metrics(hit_csv):
    df = pd.read_csv(hit_csv)
    if "hit_surface" in df.columns:
        df = df[df["hit_surface"] == 1].copy()
    if len(df) == 0:
        return None

    x_mm = df["u_m"].to_numpy(float) * 1000.0
    y_mm = df["v_m"].to_numpy(float) * 1000.0
    if "relative_efficiency" in df.columns:
        rel = df["relative_efficiency"].to_numpy(float)
    else:
        rel = np.ones(len(df), dtype=float)
    if "weight" in df.columns:
        w = df["weight"].to_numpy(float) * rel
    else:
        w = rel

    cx, cy = weighted_centroid(x_mm, y_mm, w)
    r = np.sqrt((x_mm - cx) ** 2 + (y_mm - cy) ** 2)
    rms = math.sqrt(np.sum(w * r * r) / np.sum(w)) if np.sum(w) > 0.0 else float("nan")
    return {
        "n": len(df),
        "x_mm": x_mm,
        "y_mm": y_mm,
        "w": w,
        "cx_mm": cx,
        "cy_mm": cy,
        "r68_mm": float(np.quantile(r, 0.68)),
        "r80_mm": float(np.quantile(r, 0.80)),
        "rms_mm": rms,
    }


def plot_metrics(summary, output, dpi):
    has_event = "event_id" in summary.columns and summary["event_id"].notna().any()
    labels = []
    for row in summary.itertuples(index=False):
        event_text = ""
        if has_event:
            event = getattr(row, "event_id", "")
            if str(event) and str(event).lower() != "nan":
                event_text = f"event {int(float(event))}\n"
        labels.append(f"{event_text}tel {int(row.telescope_id)}\n{row.name}")
    total = as_number(summary["total_photons"])
    hit = as_number(summary["hit_output_plane"])
    rms_mm = as_number(summary["weighted_spot_rms_m"]) * 1000.0

    frac = np.divide(hit, total, out=np.zeros_like(hit, dtype=float), where=total > 0)
    x = np.arange(len(summary))

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "figure.dpi": dpi,
        "savefig.dpi": dpi,
    })

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2), constrained_layout=True)
    axes[0].bar(x, frac, color="#4c78a8", edgecolor="black", linewidth=0.5)
    axes[0].set_ylabel("output-plane hit fraction")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels)
    axes[0].set_ylim(0.0, max(1.0, float(np.nanmax(frac)) * 1.12))
    axes[0].grid(axis="y", alpha=0.25)

    axes[1].bar(x, rms_mm, color="#f58518", edgecolor="black", linewidth=0.5)
    axes[1].set_ylabel("weighted RMS [mm]")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels)
    if np.isfinite(rms_mm).any():
        axes[1].set_ylim(0.0, float(np.nanmax(rms_mm)) * 1.18)
    axes[1].grid(axis="y", alpha=0.25)

    fig.suptitle("Array optical run summary")
    fig.savefig(output, bbox_inches="tight")
    print(f"Saved metrics plot = {output}")


def plot_spot_grid(summary, output, dpi, max_bins, focal_length_m):
    metrics = []
    for row in summary.itertuples(index=False):
        path = Path(str(row.output_csv))
        if path.exists():
            item = spot_metrics(path)
        else:
            item = None
        metrics.append(item)

    n_tel = len(summary)
    ncols = min(3, n_tel)
    nrows = int(math.ceil(n_tel / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(4.2 * ncols, 3.9 * nrows), squeeze=False,
        constrained_layout=True,
    )

    finite_points = [m for m in metrics if m is not None]
    if finite_points:
        max_abs = max(
            float(np.max(np.abs(np.concatenate([m["x_mm"], m["y_mm"]]))))
            for m in finite_points
        )
        lim = max(max_abs * 1.08, 1.0)
    else:
        lim = 1.0

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_under((1, 1, 1, 0))

    for idx, row in enumerate(summary.itertuples(index=False)):
        ax = axes[idx // ncols][idx % ncols]
        m = metrics[idx]
        title = f"tel {int(row.telescope_id)}: {row.name}"
        if "event_id" in summary.columns:
            event = getattr(row, "event_id", "")
            if str(event) and str(event).lower() != "nan":
                title = f"event {int(float(event))}, {title}"
        if m is None:
            ax.text(0.5, 0.5, "no hits", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title)
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_aspect("equal")
            continue

        span = max(
            float(np.max(m["x_mm"]) - np.min(m["x_mm"])),
            float(np.max(m["y_mm"]) - np.min(m["y_mm"])),
            1e-9,
        )
        bins = int(np.clip(np.sqrt(m["n"]) / 2.0, 80, max_bins))
        hist, xedges, yedges = np.histogram2d(
            m["x_mm"], m["y_mm"], bins=bins, range=[[-lim, lim], [-lim, lim]]
        )
        image = ax.imshow(
            hist.T,
            origin="lower",
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
            interpolation="nearest",
            cmap=cmap,
            vmin=0.01,
            aspect="equal",
        )
        ax.scatter(m["cx_mm"], m["cy_mm"], color="yellow", marker="x", s=42, linewidths=1.7)
        ax.add_patch(
            patches.Circle(
                (m["cx_mm"], m["cy_mm"]),
                m["r68_mm"],
                fill=False,
                edgecolor="orange",
                linewidth=1.7,
            )
        )
        ax.set_title(f"{title}\nRMS={m['rms_mm']:.2f} mm, R68={m['r68_mm']:.2f} mm")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("u [mm]")
        ax.set_ylabel("v [mm]")
        ax.grid(color="white", alpha=0.22, linewidth=0.5)
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03, label="Count / bin")

        def mm_to_deg(mm):
            return np.degrees(np.asarray(mm, dtype=float) / (focal_length_m * 1000.0))

        def deg_to_mm(deg):
            return np.radians(np.asarray(deg, dtype=float)) * focal_length_m * 1000.0

        secx = ax.secondary_xaxis("top", functions=(mm_to_deg, deg_to_mm))
        secy = ax.secondary_yaxis("right", functions=(mm_to_deg, deg_to_mm))
        secx.set_xlabel("u [deg]")
        secy.set_ylabel("v [deg]")

    for idx in range(n_tel, nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.suptitle("Array output-plane spots")
    fig.savefig(output, bbox_inches="tight")
    print(f"Saved spot grid = {output}")


def main():
    parser = argparse.ArgumentParser(description="Plot outputs from python/run_array_sim.py.")
    parser.add_argument("--summary-csv", required=True)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--dpi", type=int, default=350)
    parser.add_argument("--max-bins", type=int, default=260)
    parser.add_argument("--focal-length-m", type=float, default=8.0)
    args = parser.parse_args()

    summary_path = Path(args.summary_csv).resolve()
    output_dir = Path(args.output_dir).resolve() if args.output_dir else summary_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    summary = pd.read_csv(summary_path)
    if len(summary) == 0:
        raise SystemExit("array summary has no rows")
    if "output_csv" not in summary.columns:
        raise SystemExit("array summary is missing output_csv")

    plot_metrics(summary, output_dir / "array_metrics.png", args.dpi)
    plot_spot_grid(
        summary,
        output_dir / "array_spot_grid.png",
        args.dpi,
        args.max_bins,
        args.focal_length_m,
    )


if __name__ == "__main__":
    main()
