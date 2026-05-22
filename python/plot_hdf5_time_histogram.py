#!/usr/bin/env python3
"""Plot arrival-time histograms from LACT_sim HDF5 waveform output."""

import argparse
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

from plot_hdf5_camera import find_images, resolve_event_id


def decode_attr(value, default=""):
    if value is None:
        return default
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def output_path(base_output, image, multiple):
    out = Path(base_output)
    if not multiple:
        if out.parent:
            out.parent.mkdir(parents=True, exist_ok=True)
        return out
    out_dir = out if out.suffix == "" else out.with_suffix("")
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / (
        f"event_{int(image['event_id'])}_tel_{int(image['telescope_id']) + 1:02d}_time_hist.png"
    )


def quantity_label(quantity):
    labels = {
        "photon_count": "Photon count",
        "cherenkov_pe": "Cherenkov p.e.",
        "nsb_pe": "NSB p.e.",
        "pe": "p.e.",
    }
    return labels.get(quantity, quantity)


def main():
    parser = argparse.ArgumentParser(
        description="Plot camera arrival-time histograms from HDF5 /waveforms datasets."
    )
    parser.add_argument("h5", help="HDF5 file written by run_corsika_trace")
    parser.add_argument("--event-id", type=int, default=None)
    parser.add_argument("--shower-event-number", "--event-number", type=int, default=None)
    parser.add_argument("--shower-event-id", type=int, default=None)
    parser.add_argument("--array-id", type=int, default=0)
    parser.add_argument("--telescope-id", type=int, default=None)
    parser.add_argument("--image-index", type=int, default=None)
    parser.add_argument(
        "--quantity",
        choices=("photon_count", "cherenkov_pe", "nsb_pe", "pe"),
        default="pe",
    )
    parser.add_argument(
        "--mode",
        choices=("combined", "per-telescope"),
        default="combined",
        help=(
            "combined sums all selected telescopes into one histogram; "
            "per-telescope writes one histogram per selected telescope."
        ),
    )
    parser.add_argument("--output", default="time_hist.png")
    parser.add_argument("--dpi", type=int, default=350)
    args = parser.parse_args()

    with h5py.File(args.h5, "r") as h5:
        if "waveforms" not in h5:
            raise SystemExit("This HDF5 file has no /waveforms group.")
        dataset_name = f"waveforms/{args.quantity}"
        if dataset_name not in h5:
            available = ", ".join(h5["waveforms"].keys())
            raise SystemExit(f"{dataset_name} is missing. Available waveform datasets: {available}")

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

        index = h5["images/index"][:]
        images = find_images(index, event_id, args.telescope_id, args.image_index)
        time_edges = h5["waveforms/time_edges_ns"][:]
        dataset = h5[dataset_name]
        time_reference = decode_attr(h5["waveforms"].attrs.get("time_reference"), "absolute")

        histograms = []
        for image in images:
            image_index = int(image["image_index"])
            if image_index < 0 or image_index >= dataset.shape[0]:
                raise SystemExit(f"image_index={image_index} is outside waveform dataset shape.")
            hist = np.asarray(dataset[image_index, :, :], dtype=float).sum(axis=1)
            histograms.append((image, hist))

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

    def draw_one(output, rows, title_suffix):
        fig, ax = plt.subplots(figsize=(7.6, 4.4))
        total = np.zeros_like(rows[0][1], dtype=float)
        for _, hist in rows:
            total += hist
        ax.stairs(total, time_edges, color="#1f77b4", linewidth=1.8, fill=True, alpha=0.28)
        ax.stairs(total, time_edges, color="#1f77b4", linewidth=1.8)
        if time_reference == "image_first":
            xlabel = "time - T0 [ns]"
        elif time_reference == "image_mean":
            xlabel = "time - image mean [ns]"
        else:
            xlabel = "time [ns]"
        ax.set_xlabel(xlabel)
        ax.set_ylabel(quantity_label(args.quantity))
        event = int(rows[0][0]["event_id"])
        ax.set_title(f"event {event} arrival-time histogram{title_suffix}")
        ax.grid(True, alpha=0.22, linewidth=0.5)
        ax.text(
            0.98,
            0.92,
            f"sum = {float(total.sum()):.3g}",
            ha="right",
            va="top",
            transform=ax.transAxes,
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="0.82", alpha=0.92),
        )
        fig.tight_layout()
        fig.savefig(output, bbox_inches="tight")
        plt.close(fig)
        return output

    saved = []
    if args.mode == "combined":
        out = Path(args.output)
        if out.parent:
            out.parent.mkdir(parents=True, exist_ok=True)
        saved.append(draw_one(out, histograms, ""))
    else:
        multiple = len(histograms) > 1
        for image, hist in histograms:
            out = output_path(args.output, image, multiple)
            saved.append(draw_one(out, [(image, hist)], f", telescope {int(image['telescope_id']) + 1}"))

    if len(saved) == 1:
        print(f"Saved {saved[0]}")
    else:
        print(f"Saved {len(saved)} time histograms to {saved[0].parent}")


if __name__ == "__main__":
    main()
