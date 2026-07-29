#!/usr/bin/env python3
"""Run the baseline and 0..90 degree deformation cases with run_optical_sim."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def replace_value(text: str, key: str, value: str) -> str:
    prefix = f"{key}="
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if line.strip().startswith(prefix):
            lines[index] = prefix + value
            break
    else:
        lines.append(prefix + value)
    return "\n".join(lines) + "\n"


def run_case(binary: Path, config: Path, log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [str(binary), str(config)],
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="build/run_optical_sim")
    parser.add_argument("--output-root", default="run_logs/coordinate_refresh_main")
    args = parser.parse_args()

    repo = Path.cwd()
    binary = (repo / args.binary).resolve()
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = repo / output_root
    output_root.mkdir(parents=True, exist_ok=True)

    config_dir = repo / "configs" / "examples"
    baseline_source = config_dir / "coordinate_parallel_camera.cfg"
    deformation_source = config_dir / "coordinate_parallel_deformation_camera.cfg"

    baseline_text = baseline_source.read_text(encoding="utf-8")
    baseline_text = replace_value(
        baseline_text, "output.csv", (output_root / "parallel_baseline_hits.csv").as_posix()
    )
    baseline_text = replace_value(
        baseline_text,
        "output.pixel_csv",
        (output_root / "parallel_baseline_camera.csv").as_posix(),
    )
    generated_baseline = config_dir / "coordinate_generated_baseline.cfg"
    generated_baseline.write_text(baseline_text, encoding="utf-8")
    run_case(binary, generated_baseline, output_root / "parallel_baseline.log")

    deformation_template = deformation_source.read_text(encoding="utf-8")
    elevation_root = output_root / "elevation"
    elevation_root.mkdir(parents=True, exist_ok=True)
    for elevation in range(0, 91, 10):
        tag = str(elevation)
        text = replace_value(
            deformation_template, "telescope.pointing_el_deg", tag
        )
        text = replace_value(
            text, "output.csv", (elevation_root / f"hits_el_{tag}.csv").as_posix()
        )
        text = replace_value(
            text,
            "output.pixel_csv",
            (elevation_root / f"camera_el_{tag}.csv").as_posix(),
        )
        generated = config_dir / f"coordinate_generated_el_{tag}.cfg"
        generated.write_text(text, encoding="utf-8")
        run_case(binary, generated, elevation_root / f"run_el_{tag}.log")


if __name__ == "__main__":
    main()
