#!/usr/bin/env python3
"""Convert hessioxxx read_iact text output into LACT_sim PhotonCsv.

This is a diagnostic bridge, not the final high-performance EventIO reader.
It parses lines printed by `read_iact -n N file.corsika.gz`.
"""

import argparse
import csv
import math
import re
from pathlib import Path


OUTPUT_COLUMNS = [
    "x_m",
    "y_m",
    "z_m",
    "dir_x",
    "dir_y",
    "dir_z",
    "time_ns",
    "wavelength_nm",
    "weight",
    "multiplicity",
    "event_id",
    "telescope_id",
]

TEL_RE = re.compile(
    r"Telescope no\.\s+(?P<tel>\d+)\s+in array\s+(?P<array>\d+)\s+gets\s+"
    r"(?P<photons>[-+0-9.eE]+)\s+photons in\s+(?P<nbunches>\d+)\s+bunches"
)

EVENT_RE = re.compile(r"Event(?:\s+no\.?| number)?\s+(?P<event>\d+)", re.IGNORECASE)

NORMAL_RE = re.compile(
    r"Bunch at\s+(?P<x>[-+0-9.eE]+)\s+m,\s+(?P<y>[-+0-9.eE]+)\s+m\s+"
    r"in direction\s+(?P<cx>[-+0-9.eE]+)\s*,\s*(?P<cy>[-+0-9.eE]+),\s+"
    r"arrival time\s+(?P<t>[-+0-9.eE]+)\s+ns,\s+"
    r"emission at\s+(?P<zem>[-+0-9.eE]+)\s+m height,\s+with\s+"
    r"(?P<photons>[-+0-9.eE]+)\s+photons of wavelength\s+"
    r"(?P<wlen>[-+0-9.eE]+)\s+nm\."
)

THREED_RE = re.compile(
    r"Bunch at\s+(?P<x>[-+0-9.eE]+)\s+m,\s+(?P<y>[-+0-9.eE]+)\s+m,\s+"
    r"(?P<z>[-+0-9.eE]+)\s+m\s+in direction\s+"
    r"(?P<cx>[-+0-9.eE]+)\s*,\s*(?P<cy>[-+0-9.eE]+)\s*,\s*(?P<cz>[-+0-9.eE]+),\s+"
    r"arrival time\s+(?P<t>[-+0-9.eE]+)\s+ns,\s+"
    r"emission at\s+(?P<dist>[-+0-9.eE]+)\s+m distance,\s+with\s+"
    r"(?P<photons>[-+0-9.eE]+)\s+photons of wavelength\s+"
    r"(?P<wlen>[-+0-9.eE]+)\s+nm\."
)

COMPACT_RE = re.compile(
    r"Bunch \(compact\) at\s+(?P<x_cm>[-+0-9.eE]+)\s*,\s*(?P<y_cm>[-+0-9.eE]+)\s+"
    r"direction\s+(?P<cx>[-+0-9.eE]+)\s*,\s*(?P<cy>[-+0-9.eE]+),\s+"
    r"arrival time\s+(?P<t>[-+0-9.eE]+),\s+"
    r"emission at\s+(?P<zem_cm>[-+0-9.eE]+),\s+with\s+"
    r"(?P<photons>[-+0-9.eE]+)\s+photons of wavelength\s+"
    r"(?P<wlen>[-+0-9.eE]+)\s+nm\."
)


def direction_from_2d(cx, cy):
    cz2 = max(0.0, 1.0 - cx * cx - cy * cy)
    return cx, cy, -math.sqrt(cz2)


def output_event_id(event_id, array_id, event_id_mode):
    if event_id_mode == "event_array100":
        return event_id * 100 + array_id
    return event_id


def parse_text(path, default_event_id, event_id_mode):
    current_tel = None
    current_array = 0
    current_event = default_event_id
    rows = []
    with open(path) as f:
        for line in f:
            event_match = EVENT_RE.search(line)
            if event_match:
                current_event = int(event_match.group("event"))

            tel_match = TEL_RE.search(line)
            if tel_match:
                current_tel = int(tel_match.group("tel"))
                current_array = int(tel_match.group("array"))
                continue

            match = THREED_RE.search(line)
            if match and current_tel is not None:
                g = {k: float(v) for k, v in match.groupdict().items()}
                rows.append({
                    "x_m": g["x"],
                    "y_m": g["y"],
                    "z_m": g["z"],
                    "dir_x": g["cx"],
                    "dir_y": g["cy"],
                    "dir_z": g["cz"],
                    "time_ns": g["t"],
                    "wavelength_nm": g["wlen"],
                    "weight": 1.0,
                    "multiplicity": g["photons"],
                    "event_id": output_event_id(current_event, current_array, event_id_mode),
                    "telescope_id": current_tel,
                })
                continue

            match = NORMAL_RE.search(line)
            if match and current_tel is not None:
                g = {k: float(v) for k, v in match.groupdict().items()}
                dx, dy, dz = direction_from_2d(g["cx"], g["cy"])
                rows.append({
                    "x_m": g["x"],
                    "y_m": g["y"],
                    "z_m": 0.0,
                    "dir_x": dx,
                    "dir_y": dy,
                    "dir_z": dz,
                    "time_ns": g["t"],
                    "wavelength_nm": g["wlen"],
                    "weight": 1.0,
                    "multiplicity": g["photons"],
                    "event_id": output_event_id(current_event, current_array, event_id_mode),
                    "telescope_id": current_tel,
                })
                continue

            match = COMPACT_RE.search(line)
            if match and current_tel is not None:
                g = {k: float(v) for k, v in match.groupdict().items()}
                dx, dy, dz = direction_from_2d(g["cx"], g["cy"])
                # compact read_iact output prints core coordinates in cm.
                rows.append({
                    "x_m": g["x_cm"] / 100.0,
                    "y_m": g["y_cm"] / 100.0,
                    "z_m": 0.0,
                    "dir_x": dx,
                    "dir_y": dy,
                    "dir_z": dz,
                    "time_ns": g["t"],
                    "wavelength_nm": g["wlen"],
                    "weight": 1.0,
                    "multiplicity": g["photons"],
                    "event_id": output_event_id(current_event, current_array, event_id_mode),
                    "telescope_id": current_tel,
                })
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Parse hessioxxx read_iact text output into LACT_sim PhotonCsv."
    )
    parser.add_argument("input_txt")
    parser.add_argument("output_csv")
    parser.add_argument("--default-event-id", type=int, default=0)
    parser.add_argument(
        "--event-id-mode",
        choices=["event", "event_array100"],
        default="event",
        help=(
            "event: keep the CORSIKA shower event number; "
            "event_array100: encode read_iact array blocks as event*100+array, "
            "matching the provided ROOT runid convention."
        ),
    )
    args = parser.parse_args()

    rows = parse_text(args.input_txt, args.default_event_id, args.event_id_mode)
    if not rows:
        raise SystemExit("no photon bunch rows parsed from read_iact text")

    out_path = Path(args.output_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {out_path} with {len(rows)} photon bunches")


if __name__ == "__main__":
    main()
