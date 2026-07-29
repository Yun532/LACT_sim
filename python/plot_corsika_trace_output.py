#!/usr/bin/env python3
"""Plot one event from run_corsika_trace output.

The script reads only the C++ trace output file. It supports both:
  - pixel aggregate CSV from camera mode
  - whiteboard photon-hit CSV from whiteboard mode
"""

import argparse
import csv
import math
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


def read_camera_pixels(path):
    pixels = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            pixels.append({
                "id": int(row["id"]),
                "x_m": float(row["x_m"]),
                "y_m": float(row["y_m"]),
                "size_m": float(row["size_m"]),
            })
    return pixels


def setup_style(dpi):
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


def plot_pixel_image(df, pixels, event_label, telescope_id, output):
    df = df.groupby("pixel_id", as_index=False).agg(
        photon_count=("photon_count", "sum"),
        pe=("pe", "sum"),
        signal=("signal", "sum"),
    )
    values_by_id = {
        int(row.pixel_id): float(row.signal)
        for row in df.itertuples(index=False)
    }
    patch_list = []
    values = []
    for pixel in pixels:
        size = pixel["size_m"]
        x = pixel["x_m"]
        y = pixel["y_m"]
        patch_list.append(Rectangle((x - 0.5 * size, y - 0.5 * size), size, size))
        values.append(values_by_id.get(pixel["id"], 0.0))

    fig, ax = plt.subplots(figsize=(6.8, 6.2))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_under((1.0, 1.0, 1.0, 0.0))
    positive = [v for v in values if v > 0.0]

    collection = PatchCollection(
        patch_list,
        cmap=cmap,
        edgecolor=(0.15, 0.15, 0.15, 0.25),
        linewidth=0.25,
    )
    collection.set_array(np.asarray(values, dtype=float))
    collection.set_clim(vmin=0.5, vmax=max(positive) if positive else 1.0)
    ax.add_collection(collection)
    ax.set_aspect("equal", adjustable="box")
    ax.autoscale_view()
    ax.set_xlabel("LACT focal-plane u [m]")
    ax.set_ylabel("LACT focal-plane v [m]")
    ax.set_title(f"CORSIKA {event_label} telescope {telescope_id} camera")
    ax.grid(True, alpha=0.18, linewidth=0.5)
    cbar = fig.colorbar(collection, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("weighted photons")
    ax.text(
        0.02,
        0.02,
        f"filled pixels = {len(df)}\nsignal = {df['signal'].sum():.1f}\n"
        "stored coordinates = physical (u, v)",
        transform=ax.transAxes,
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.82", alpha=0.9),
    )
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def weighted_centroid(x, y, w):
    sw = np.sum(w)
    if sw <= 0.0:
        return float(np.mean(x)), float(np.mean(y))
    return float(np.sum(w * x) / sw), float(np.sum(w * y) / sw)


def plot_whiteboard(df, event_label, telescope_id, output, max_bins):
    x_mm = df["u_m"].to_numpy(float) * 1000.0
    y_mm = df["v_m"].to_numpy(float) * 1000.0
    w = df["signal_weight"].to_numpy(float)
    cx, cy = weighted_centroid(x_mm, y_mm, w)
    r = np.sqrt((x_mm - cx) ** 2 + (y_mm - cy) ** 2)
    r68 = float(np.quantile(r, 0.68))

    span = max(float(np.max(x_mm) - np.min(x_mm)),
               float(np.max(y_mm) - np.min(y_mm)),
               1e-9)
    bins = int(np.clip(np.sqrt(len(df)) / 2.0, 140, max_bins))
    pad = 0.08 * span
    xmid = 0.5 * (float(np.max(x_mm)) + float(np.min(x_mm)))
    ymid = 0.5 * (float(np.max(y_mm)) + float(np.min(y_mm)))
    hist, xedges, yedges = np.histogram2d(
        x_mm,
        y_mm,
        bins=bins,
        range=[
            (xmid - 0.5 * span - pad, xmid + 0.5 * span + pad),
            (ymid - 0.5 * span - pad, ymid + 0.5 * span + pad),
        ],
    )

    fig, ax = plt.subplots(figsize=(6.5, 5.8))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_under((1.0, 1.0, 1.0, 0.0))
    image = ax.imshow(
        hist.T,
        interpolation="nearest",
        cmap=cmap,
        origin="lower",
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        vmin=0.01,
        aspect="equal",
    )
    fig.colorbar(image, ax=ax, pad=0.03).set_label("Count / bin")
    ax.scatter(cx, cy, color="yellow", marker="x", s=64, linewidths=2.0, label="Centroid")
    ax.add_patch(patches.Circle(
        (cx, cy),
        r68,
        fill=False,
        edgecolor="orange",
        linewidth=2.2,
        label=f"R68 = {r68:.2f} mm",
    ))
    ax.set_xlabel("LACT focal-plane u [mm]")
    ax.set_ylabel("LACT focal-plane v [mm]")
    ax.set_title(f"CORSIKA {event_label} telescope {telescope_id} whiteboard")
    ax.grid(True, color="white", alpha=0.22, linewidth=0.5)
    ax.legend(frameon=True, framealpha=1.0, edgecolor="black", loc="best")
    ax.text(
        0.02,
        0.02,
        f"N = {len(df)}\nRMS = {math.sqrt(np.average(r * r, weights=w)):.2f} mm\n"
        "stored coordinates = physical (u, v)",
        transform=ax.transAxes,
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.82", alpha=0.92),
    )
    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def default_summary_path(trace_path):
    candidates = [
        trace_path.with_name("corsika_trace_summary.csv"),
        trace_path.with_name("camera_summary.csv"),
        trace_path.with_name("whiteboard_summary.csv"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def shower_event_from_event_id(event_id):
    return int(event_id) // 100


def select_shower_event_from_summary(summary_csv, shower_event_number):
    summary_path = Path(summary_csv)
    if not summary_path.exists():
        raise SystemExit(
            "--shower-event-number requires a summary CSV with the full event list. "
            f"Expected {summary_path}; pass --summary-csv explicitly if it is elsewhere."
        )

    summary = pd.read_csv(summary_path, usecols=["event_id"])
    shower_ids = []
    seen = set()
    for event_id in summary["event_id"].astype(np.int64).to_numpy():
        shower_id = int(event_id // 100)
        if shower_id not in seen:
            seen.add(shower_id)
            shower_ids.append(shower_id)
    if shower_event_number < 1 or shower_event_number > len(shower_ids):
        raise SystemExit(
            f"--shower-event-number must be in 1..{len(shower_ids)} for {summary_path}."
        )
    return shower_ids[shower_event_number - 1], len(shower_ids)


def main():
    parser = argparse.ArgumentParser(description="Plot run_corsika_trace output.")
    parser.add_argument("trace_csv", help="pixel aggregate CSV or whiteboard hit CSV")
    parser.add_argument("--event-id", type=int, default=None)
    parser.add_argument(
        "--event-number",
        type=int,
        default=None,
        help="1-based index in the sorted event_id list contained in the trace CSV",
    )
    parser.add_argument(
        "--shower-event-id",
        type=int,
        default=None,
        help="Original CORSIKA shower event id for event_array100 outputs",
    )
    parser.add_argument(
        "--shower-event-number",
        type=int,
        default=None,
        help="1-based index in the sorted original shower-event list for event_array100 outputs",
    )
    parser.add_argument(
        "--array-id",
        type=int,
        default=None,
        help="Optional array id when selecting by --shower-event-id/--shower-event-number",
    )
    parser.add_argument("--telescope-id", type=int, default=None)
    parser.add_argument(
        "--summary-csv",
        default=None,
        help="run_corsika_trace summary CSV; used to resolve --shower-event-number",
    )
    parser.add_argument("--camera-csv", default="configs/cameras/new_camera_pixels.csv")
    parser.add_argument(
        "--config",
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--dpi", type=int, default=350)
    parser.add_argument("--max-bins", type=int, default=520)
    args = parser.parse_args()

    setup_style(args.dpi)
    trace_path = Path(args.trace_csv)
    df = pd.read_csv(trace_path)
    is_pixel = {"pixel_id", "signal", "photon_count"}.issubset(df.columns)
    is_whiteboard = {"u_m", "v_m", "signal_weight"}.issubset(df.columns)
    if not is_pixel and not is_whiteboard:
        raise SystemExit("Input is neither run_corsika_trace pixel CSV nor whiteboard CSV.")

    selection_count = sum(
        x is not None
        for x in (
            args.event_id,
            args.event_number,
            args.shower_event_id,
            args.shower_event_number,
        )
    )
    if selection_count == 0:
        raise SystemExit(
            "Please set one of --event-id, --event-number, "
            "--shower-event-id, or --shower-event-number."
        )
    if selection_count > 1:
        raise SystemExit(
            "Use only one of --event-id, --event-number, "
            "--shower-event-id, or --shower-event-number."
        )
    if args.array_id is not None and args.shower_event_id is None and args.shower_event_number is None:
        raise SystemExit("--array-id is only valid with a shower-event selection.")

    event_label = None
    if args.event_number is not None:
        event_ids = sorted(int(x) for x in df["event_id"].unique())
        if args.event_number < 1 or args.event_number > len(event_ids):
            raise SystemExit(
                f"--event-number must be in 1..{len(event_ids)} for this file."
            )
        args.event_id = event_ids[args.event_number - 1]
        print(f"Selected event_id={args.event_id} from event_number={args.event_number}")
    if args.event_id is not None:
        df = df[df["event_id"] == args.event_id].copy()
        event_label = f"event_id {args.event_id}"
    else:
        if args.shower_event_number is not None:
            summary_csv = args.summary_csv or default_summary_path(trace_path)
            args.shower_event_id, n_shower = select_shower_event_from_summary(
                summary_csv,
                args.shower_event_number,
            )
            print(
                "Selected shower_event_id="
                f"{args.shower_event_id} from shower_event_number={args.shower_event_number} "
                f"using {summary_csv} ({n_shower} shower events)"
            )
        if args.array_id is None:
            df = df[(df["event_id"].astype(np.int64) // 100) == args.shower_event_id].copy()
            event_label = f"shower_event {args.shower_event_id} all arrays"
        else:
            selected_event_id = args.shower_event_id * 100 + args.array_id
            df = df[df["event_id"] == selected_event_id].copy()
            event_label = (
                f"shower_event {args.shower_event_id} array {args.array_id} "
                f"(event_id {selected_event_id})"
            )

    if args.telescope_id is not None:
        df = df[df["telescope_id"] == args.telescope_id].copy()
    if df.empty:
        raise SystemExit("No rows match the requested event/telescope selection.")

    default_tag = event_label.replace(" ", "_").replace("(", "").replace(")", "")
    outdir = Path(args.output_dir) if args.output_dir else trace_path.parent / "plots" / default_tag
    outdir.mkdir(parents=True, exist_ok=True)
    telescope_ids = sorted(int(x) for x in df["telescope_id"].unique())
    pixels = read_camera_pixels(args.camera_csv) if is_pixel else None
    for tel in telescope_ids:
        one = df[df["telescope_id"] == tel].copy()
        output = outdir / f"tel{tel:03d}_{'camera' if is_pixel else 'whiteboard'}.png"
        if is_pixel:
            plot_pixel_image(one, pixels, event_label, tel, output)
        else:
            plot_whiteboard(one, event_label, tel, output, args.max_bins)
        print(f"Saved {output}")


if __name__ == "__main__":
    main()
