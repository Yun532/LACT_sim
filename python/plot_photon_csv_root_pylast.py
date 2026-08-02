#!/usr/bin/env python3
"""Plot one LACT ROOT event with pyLAST's native EventVisualizer."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLBACKEND", "Agg")


def select_event(source, event_id: int | None, event_index: int):
    if event_id is None:
        return source[event_index]
    for index in range(len(source)):
        event = source[index]
        if int(event.event_id) == event_id:
            return event
    raise RuntimeError(f"event {event_id} is not present in the ROOT file")


def plot_root_dl0_with_pylast(root_file: Path, event_id: int, output: Path) -> None:
    """Render observations.image_pe with pyLAST's native camera primitive."""
    import matplotlib.pyplot as plt
    import uproot
    from pylast.visualize.visualize import plot_camera_image

    with uproot.open(root_file) as root:
        camera = root["camera_pixels"].arrays(library="np")
        observations = root["observations"].arrays(library="ak")
        matches = np.nonzero(
            np.asarray(observations.event_id, dtype=int) == int(event_id)
        )[0]
        if len(matches) == 0:
            raise RuntimeError(f"event {event_id} has no ROOT observation")
        row = int(matches[0])
        stored_by_id = {
            int(pixel_id): float(value)
            for pixel_id, value in zip(
                observations.pixel_id[row], observations.image_pe[row]
            )
        }
        pixel_ids = np.asarray(camera["pixel_id"], dtype=int)
        image = np.asarray(
            [stored_by_id.get(int(pixel_id), 0.0) for pixel_id in pixel_ids],
            dtype=float,
        )
        focal_length_m = float(
            root["optics"]["effective_focal_length_m"].array(library="np")[0]
        )
        # LactEventSource's canonical boundary mapping is
        # pix_x=-LACT v, pix_y=-LACT u. plot_camera_image then applies the
        # established pyLAST display-axis convention.
        pix_x_deg = np.degrees(
            np.arctan2(-np.asarray(camera["y_m"], dtype=float), focal_length_m)
        )
        pix_y_deg = np.degrees(
            np.arctan2(-np.asarray(camera["x_m"], dtype=float), focal_length_m)
        )
        pixel_size_deg = float(
            np.median(
                np.degrees(
                    np.arctan2(
                        np.asarray(camera["size_m"], dtype=float), focal_length_m
                    )
                )
            )
        )

    axis = plot_camera_image(
        pix_x_deg,
        pix_y_deg,
        pixel_size_deg,
        image,
        mask=image != 0.0,
        vmin=0.0,
        vmax=float(np.max(image)) if image.size else 1.0,
        title=f"LACT event {event_id} DL0 image",
        pixel_shape="square",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    axis.figure.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(axis.figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_file", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--event-id", type=int)
    selection.add_argument("--event-index", type=int, default=0)
    args = parser.parse_args()

    from pylast.io import LactEventSource
    from pylast.visualize import EventVisualizer

    source = LactEventSource(str(args.root_file), max_events=-1)
    event = select_event(source, args.event_id, args.event_index)
    visualizer = EventVisualizer(source)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure, _axes = visualizer.plot_event(
        event,
        output_path=str(args.output),
        # LACT_sim's observations.image_pe is the final detector image and is
        # exposed by pyLAST at DL0.  Simulation truth is intentionally kept
        # separate from this readout-level README example.
        image_level="dl0",
        show_hillas=False,
        # This README example deliberately disables the trigger, so the
        # camera image is valid even though its trigger flag is false.
        include_non_triggered=True,
        show=False,
    )
    if figure is None:
        # Older installed pyLAST native readers may not yet populate DL0 from
        # the current ROOT schema. Keep the same pyLAST rendering style while
        # reading the canonical observations.image_pe branch directly.
        plot_root_dl0_with_pylast(
            args.root_file, int(event.event_id), args.output
        )

    print(f"event_id={int(event.event_id)}; plotted with native pyLAST camera style")
    print(args.output)


if __name__ == "__main__":
    main()
