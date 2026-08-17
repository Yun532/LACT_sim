"""Reconcile direct-PDE baseline and pitch-saturation ROOT outputs."""

from pathlib import Path
import sys

import numpy as np
import uproot


def load_decisions(run_dir: Path) -> dict[str, np.ndarray]:
    with uproot.open(run_dir / "lact_events.root") as root_file:
        return root_file["microcell_decisions"].arrays(library="np")


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(
            "usage: summarize_s17351_validation.py BASELINE_DIR SATURATION_DIR"
        )
    baseline = load_decisions(Path(sys.argv[1]))
    saturation = load_decisions(Path(sys.argv[2]))
    identity_fields = [
        "event_id",
        "telescope_id",
        "pixel_id",
        "time_ns",
        "sensor_x_m",
        "sensor_y_m",
        "grid_column",
        "grid_row",
        "channel_id",
        "microcell_id",
        "origin",
    ]
    primary_identity = all(
        np.array_equal(baseline[field], saturation[field])
        for field in identity_fields
    )
    primary_count = len(baseline["event_id"])
    saturation_primary_count = len(saturation["event_id"])
    print(f"baseline_primary={primary_count}")
    print(f"baseline_fired={int(baseline['fired'].sum())}")
    print(f"saturation_primary={saturation_primary_count}")
    print(f"saturation_fired={int(saturation['fired'].sum())}")
    print(
        "saturation_rejected="
        f"{int(saturation['saturation_rejected'].sum())}"
    )
    print(f"primary_identity={str(primary_identity).lower()}")
    if not primary_identity:
        raise SystemExit("primary p.e. sequences differ")


if __name__ == "__main__":
    main()
