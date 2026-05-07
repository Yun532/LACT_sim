#!/usr/bin/env python3
"""Plot a camera image directly from a LACT_sim HDF5 file."""

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


def event_table(h5):
    if "events/corsika" in h5:
        return h5["events/corsika"][:]
    if "events/table" in h5:
        return h5["events/table"][:]
    raise SystemExit("This HDF5 file has no events table; use --image-index.")


def event_id_from_shower_id(h5, shower_event_id, array_id):
    events = event_table(h5)
    names = events.dtype.names or ()
    if "shower_event_id" not in names or "array_id" not in names:
        raise SystemExit(
            "--shower-event-id requires events/corsika metadata. "
            "Use --event-id for this file."
        )
    matches = events[
        (events["shower_event_id"] == shower_event_id) & (events["array_id"] == array_id)
    ]
    if len(matches) == 0:
        raise SystemExit(
            f"No event found for shower_event_id={shower_event_id} and array_id={array_id}."
        )
    if len(matches) > 1:
        raise SystemExit("Internal error: multiple event rows match shower/array selection.")
    return int(matches[0]["event_id"])


def event_id_from_shower_number(h5, shower_event_number, array_id):
    events = event_table(h5)
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


def resolve_event_id(h5, event_id=None, shower_event_id=None, shower_event_number=None, array_id=0):
    choices = [
        event_id is not None,
        shower_event_id is not None,
        shower_event_number is not None,
    ]
    if sum(choices) > 1:
        raise SystemExit(
            "Use only one event selector: --event-id, --shower-event-id, "
            "or --shower-event-number."
        )
    if event_id is not None:
        return int(event_id)
    if shower_event_id is not None:
        return event_id_from_shower_id(h5, shower_event_id, array_id)
    if shower_event_number is not None:
        return event_id_from_shower_number(h5, shower_event_number, array_id)
    return None


def find_images(index, event_id=None, telescope_id=None, image_index=None):
    if image_index is not None:
        rows = index[index["image_index"] == image_index]
    else:
        rows = index[index["event_id"] == event_id]
        if telescope_id is not None:
            rows = rows[rows["telescope_id"] == telescope_id]
    if len(rows) == 0:
        raise SystemExit("No image matches the requested selection.")
    if image_index is not None and len(rows) > 1:
        choices = ", ".join(
            f"image_index={int(r['image_index'])},event_id={int(r['event_id'])},tel={int(r['telescope_id'])}"
            for r in rows[:10]
        )
        raise SystemExit(f"Selection matches multiple images. Choices: {choices}")
    return rows


def values_for_image(h5, image, quantity):
    if "dense" in h5["images"] and quantity in h5["images/dense"]:
        pixel_axis = h5["images/dense/pixel_id_axis"][:]
        values = h5[f"images/dense/{quantity}"][int(image["image_index"]), :]
        return {int(pid): float(v) for pid, v in zip(pixel_axis, values)}

    if quantity in ("cherenkov_pe", "nsb_pe"):
        raise SystemExit(
            f"{quantity} is only available for dense HDF5 files written with "
            "output.hdf5_write_components=true."
        )
    start = int(image["start"])
    count = int(image["count"])
    rows = h5["images/sparse/pixels"][start:start + count]
    return {int(r["pixel_id"]): float(r[quantity]) for r in rows}


def draw_camera(camera, image, values_by_pixel, quantity, dpi):
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
        "figure.dpi": dpi,
        "savefig.dpi": dpi,
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
        f"HDF5 camera event {int(image['event_id'])}, telescope {int(image['telescope_id']) + 1}"
    )
    ax.grid(True, alpha=0.18, linewidth=0.5)
    cbar = fig.colorbar(collection, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(quantity)
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
    return fig


def output_path(base_output, image, quantity, multiple):
    out = Path(base_output)
    if not multiple:
        if out.parent:
            out.parent.mkdir(parents=True, exist_ok=True)
        return out
    out_dir = out if out.suffix == "" else out.with_suffix("")
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / (
        f"event_{int(image['event_id'])}_tel_{int(image['telescope_id']) + 1:02d}_{quantity}.png"
    )


def main():
    parser = argparse.ArgumentParser(description="Plot one camera image from LACT_sim HDF5.")
    parser.add_argument("h5", help="HDF5 file from run_corsika_trace or export_trace_hdf5.py")
    parser.add_argument("--event-id", type=int, default=None)
    parser.add_argument(
        "--shower-event-number",
        "--event-number",
        dest="shower_event_number",
        type=int,
        default=None,
        help="1-based original CORSIKA shower-event order in the HDF5 file.",
    )
    parser.add_argument(
        "--shower-event-id",
        type=int,
        default=None,
        help="Original CORSIKA shower event id. Use with --array-id to select one output event.",
    )
    parser.add_argument(
        "--array-id",
        type=int,
        default=0,
        help=(
            "CORSIKA CSCAT/MC_TELOFF array-use index for the selected shower event; "
            "not the telescope ID."
        ),
    )
    parser.add_argument("--telescope-id", type=int, default=None)
    parser.add_argument("--image-index", type=int, default=None)
    parser.add_argument(
        "--quantity",
        choices=("signal", "pe", "photon_count", "cherenkov_pe", "nsb_pe"),
        default="signal",
        help="Dense image quantity to plot. Component quantities require output.hdf5_write_components=true.",
    )
    parser.add_argument("--output", default="hdf5_camera.png")
    parser.add_argument("--dpi", type=int, default=350)
    args = parser.parse_args()

    with h5py.File(args.h5, "r") as h5:
        event_id = resolve_event_id(
            h5,
            event_id=args.event_id,
            shower_event_id=args.shower_event_id,
            shower_event_number=args.shower_event_number,
            array_id=args.array_id,
        )
        if args.image_index is None and event_id is None:
            raise SystemExit(
                "Set --image-index, --event-id, --shower-event-id, or --shower-event-number."
            )

        camera = h5["camera/pixels"][:]
        index = h5["images/index"][:]
        images = find_images(index, event_id, args.telescope_id, args.image_index)
        multiple = len(images) > 1

        saved = []
        for image in images:
            values_by_pixel = values_for_image(h5, image, args.quantity)
            fig = draw_camera(camera, image, values_by_pixel, args.quantity, args.dpi)
            out = output_path(args.output, image, args.quantity, multiple)
            fig.savefig(out, bbox_inches="tight")
            plt.close(fig)
            saved.append(out)

    if len(saved) == 1:
        print(f"Saved {saved[0]}")
    else:
        print(f"Saved {len(saved)} camera images to {saved[0].parent}")


if __name__ == "__main__":
    main()
