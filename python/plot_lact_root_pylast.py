#!/usr/bin/env python3
"""Read an image_pe LACT ROOT file and plot it with native pyLAST."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

from pylast.io import LactEventSource
from pylast.visualize import EventVisualizer


def select_event(source, event_id: int | None, event_index: int):
    if event_id is None:
        return source[event_index]
    for index in range(len(source)):
        event = source[index]
        if int(event.event_id) == event_id:
            return event
    raise RuntimeError(f"event {event_id} is not present in the ROOT file")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_file", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--event-id", type=int)
    selection.add_argument("--event-index", type=int, default=0)
    args = parser.parse_args()

    source = LactEventSource(str(args.root_file), max_events=-1)
    event = select_event(source, args.event_id, args.event_index)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    visualizer = EventVisualizer(source)
    figure, _axes = visualizer.plot_event(
        event,
        output_path=str(args.output),
        image_level="dl0",
        show_hillas=False,
        include_non_triggered=False,
        show=False,
    )
    if figure is None:
        raise RuntimeError(f"event {int(event.event_id)} has no camera images to plot")
    print(f"event_id={int(event.event_id)}; plotted with native pyLAST EventVisualizer")
    print(args.output)


if __name__ == "__main__":
    main()
