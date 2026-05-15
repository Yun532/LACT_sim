#!/usr/bin/env python3
"""Plot a static 3D coordinate-system check for a LACT optical cfg.

The plot includes mirror facets, the output/camera plane, configured camera
+u/+v image axes, camera normal, global +z/up, and obstruction primitives when
the cfg enables them.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from config_io import expand_component_config
from plot_optical_layout_3d import main as plot_layout_main


def _truthy(text: str | None) -> bool:
    return str(text or "").strip().lower() in {"1", "true", "yes", "on"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="LACT run cfg to inspect")
    parser.add_argument(
        "--output",
        default="coordinate_system_3d.png",
        help="output PNG/PDF path",
    )
    parser.add_argument("--dpi", type=int, default=350, help="output DPI")
    parser.add_argument("--view", default="32,-58", help="matplotlib 3D view: elev,azim")
    parser.add_argument(
        "--elevation-deg",
        type=float,
        default=None,
        help="optional telescope/mirror elevation override",
    )
    parser.add_argument(
        "--ray-stride",
        type=int,
        default=12,
        help="draw every Nth center ray",
    )
    parser.add_argument(
        "--normal-scale",
        type=float,
        default=0.35,
        help="mirror-normal arrow length in m",
    )
    parser.add_argument(
        "--no-obstruction",
        action="store_true",
        help="do not draw obstruction primitives even if the cfg enables them",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_path = Path(args.config).resolve()
    cfg, _ = expand_component_config(config_path)

    plot_args = [
        "plot_optical_layout_3d.py",
        "--config",
        str(config_path),
        "--output",
        args.output,
        "--dpi",
        str(args.dpi),
        "--view",
        args.view,
        "--ray-stride",
        str(args.ray_stride),
        "--normal-scale",
        str(args.normal_scale),
        "--show-camera-axes",
    ]
    if args.elevation_deg is not None:
        plot_args.extend(["--elevation-deg", str(args.elevation_deg)])
    if not args.no_obstruction and _truthy(cfg.get("obstruction.enabled")):
        plot_args.append("--show-obstruction")

    old_argv = sys.argv
    try:
        sys.argv = plot_args
        plot_layout_main()
    finally:
        sys.argv = old_argv


if __name__ == "__main__":
    main()
