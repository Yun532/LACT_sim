#!/usr/bin/env python3
"""Plot a camera image directly from a LACT_sim HDF5 file."""

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


def event_id_from_shower_number(h5, shower_event_number, array_id):
    if "events/table" not in h5:
        raise SystemExit("This HDF5 file has no events/table; use --event-id or --image-index.")
    events = h5["events/table"][:]
    if "shower_event_id" not in events.dtype.names:
        if shower_event_number < 1 or shower_event_number > len(events):
            raise SystemExit(
                f"--shower-event-number must be 1..{len(events)} for this file."
            )
        return int(events[shower_event_number - 1]["event_id"])
    shower_ids = []
    seen = set()
    for row in events:
        shower = int(row["shower_event_id"])
        if shower not in seen:
            seen.add(shower)
            shower_ids.append(shower)
    if shower_event_number < 1 or shower_event_number > len(shower_ids):
        raise SystemExit(
            f"--shower-event-number must be 1..{len(shower_ids)} for this file."
        )
    shower = shower_ids[shower_event_number - 1]
    matches = events[
        (events["shower_event_id"] == shower) & (events["array_id"] == array_id)
    ]
    if len(matches) == 0:
        raise SystemExit(
            f"No event found for shower_event_number={shower_event_number} "
            f"(shower_event_id={shower}) and array_id={array_id}."
        )
    if len(matches) > 1:
        raise SystemExit("Internal error: multiple event rows match shower/array selection.")
    return int(matches[0]["event_id"])


def find_image(index, event_id=None, telescope_id=None, image_index=None):
    if image_index is not None:
        rows = index[index["image_index"] == image_index]
    else:
        rows = index[index["event_id"] == event_id]
        if telescope_id is not None:
            rows = rows[rows["telescope_id"] == telescope_id]
    if len(rows) == 0:
        raise SystemExit("No image matches the requested selection.")
    if len(rows) > 1:
        choices = ", ".join(
            f"image_index={int(r['image_index'])},event_id={int(r['event_id'])},tel={int(r['telescope_id'])}"
            for r in rows[:10]
        )
        raise SystemExit(f"Selection matches multiple images. Add --telescope-id or --image-index. Choices: {choices}")
    return rows[0]


def main():
    parser = argparse.ArgumentParser(description="Plot one camera image from LACT_sim HDF5.")
    parser.add_argument("h5", help="HDF5 file from run_corsika_trace or export_trace_hdf5.py")
    parser.add_argument("--event-id", type=int, default=None)
    parser.add_argument(
        "--shower-event-number",
        type=int,
        default=None,
        help="1-based event order inside events/table; legacy files use shower-event order",
    )
    parser.add_argument("--array-id", type=int, default=0)
    parser.add_argument("--telescope-id", type=int, default=None)
    parser.add_argument("--image-index", type=int, default=None)
    parser.add_argument("--quantity", choices=("signal", "pe", "photon_count"), default="signal")
    parser.add_argument("--output", default="hdf5_camera.png")
    parser.add_argument("--dpi", type=int, default=350)
    args = parser.parse_args()

    with h5py.File(args.h5, "r") as h5:
        if args.event_id is not None and args.shower_event_number is not None:
            raise SystemExit("Use either --event-id or --shower-event-number, not both.")
        event_id = args.event_id
        if args.shower_event_number is not None:
            event_id = event_id_from_shower_number(h5, args.shower_event_number, args.array_id)
        if args.image_index is None and event_id is None:
            raise SystemExit("Set --image-index, --event-id, or --shower-event-number.")

        camera = h5["camera/pixels"][:]
        index = h5["images/index"][:]
        image = find_image(index, event_id, args.telescope_id, args.image_index)
        values_by_pixel = {}

        if "dense" in h5["images"] and args.quantity in h5["images/dense"]:
            pixel_axis = h5["images/dense/pixel_id_axis"][:]
            values = h5[f"images/dense/{args.quantity}"][int(image["image_index"]), :]
            values_by_pixel = {int(pid): float(v) for pid, v in zip(pixel_axis, values)}
        else:
            start = int(image["start"])
            count = int(image["count"])
            rows = h5["images/sparse/pixels"][start:start + count]
            values_by_pixel = {int(r["pixel_id"]): float(r[args.quantity]) for r in rows}

    patches = []
    values = []
    for p in camera:
        size = float(p["size_m"])
        patches.append(Rectangle(
            (float(p["x_m"]) - 0.5 * size, float(p["y_m"]) - 0.5 * size),
            size,
            size,
        ))
        values.append(values_by_pixel.get(int(p["pixel_id"]), 0.0))

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
    fig, ax = plt.subplots(figsize=(6.8, 6.2))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_under((1, 1, 1, 0))
    positive = [v for v in values if v > 0]
    collection = PatchCollection(
        patches,
        cmap=cmap,
        edgecolor=(0.15, 0.15, 0.15, 0.25),
        linewidth=0.25,
    )
    collection.set_array(np.asarray(values, dtype=float))
    collection.set_clim(vmin=0.5, vmax=max(positive) if positive else 1.0)
    ax.add_collection(collection)
    ax.set_aspect("equal", adjustable="box")
    ax.autoscale_view()
    ax.set_xlabel("camera x [m]")
    ax.set_ylabel("camera y [m]")
    ax.set_title(
        f"HDF5 camera event {int(image['event_id'])}, telescope {int(image['telescope_id'])}"
    )
    ax.grid(True, alpha=0.18, linewidth=0.5)
    cbar = fig.colorbar(collection, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(args.quantity)
    ax.text(
        0.02,
        0.02,
        f"filled = {int(image['count'])}\n"
        f"signal = {float(image['total_signal']):.1f}\n"
        f"pe = {float(image['total_pe']):.1f}",
        transform=ax.transAxes,
        fontsize=8,
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.82", alpha=0.9),
    )
    fig.tight_layout()
    out = Path(args.output)
    if out.parent:
        out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
