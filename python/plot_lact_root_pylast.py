#!/usr/bin/env python3
"""Read a timeseries_pe LACT ROOT file through pyLAST and plot DL0 cameras."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np

from pylast.calib import Calibrator
from pylast.io import LactEventSource
from pylast.visualize import EventVisualizer, plot_raw_images


CALIBRATOR_CONFIG = """{
  "Calibrator": {
    "image_extractor_type": "LocalPeakExtractor",
    "LocalPeakExtractor": {
      "window_width": 7,
      "window_shift": 3,
      "apply_correction": false
    }
  }
}
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root_file", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--event-index", type=int, default=0)
    args = parser.parse_args()

    source = LactEventSource(str(args.root_file), max_events=-1)
    event = source[args.event_index]
    calibrator = Calibrator(source.subarray, config_str=CALIBRATOR_CONFIG)
    calibrator(event)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    visualizer = EventVisualizer(source)
    plot_raw_images(
        event,
        visualizer=visualizer,
        output_path=str(args.output),
        include_non_triggered=False,
        show=False,
    )

    sums = ", ".join(
        f"tel {tel_id}: {float(np.sum(tel.image)):.0f} p.e."
        for tel_id, tel in sorted(event.dl0.tels.items())
    )
    print(f"event_id={int(event.event_id)}; {sums}")
    print(args.output)


if __name__ == "__main__":
    main()
