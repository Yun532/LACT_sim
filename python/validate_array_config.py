#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path

from config_io import read_key_value_config, resolve_workspace_path


REQUIRED_COLUMNS = {
    "telescope_id",
    "name",
    "config_path",
    "position_x_m",
    "position_y_m",
    "position_z_m",
    "pointing_az_deg",
    "pointing_el_deg",
}


def main():
    parser = argparse.ArgumentParser(description="Validate a LACT_sim array manifest.")
    parser.add_argument("--config", required=True, help="array config with mode=csv and csv_path=...")
    parser.add_argument("--summary-csv", default=None)
    args = parser.parse_args()

    cfg_path = Path(args.config).resolve()
    cfg = read_key_value_config(cfg_path)
    mode = cfg.get("mode", "csv").strip().lower()
    if mode != "csv":
        raise SystemExit(f"unsupported array mode: {mode}")
    csv_value = cfg.get("csv_path")
    if not csv_value:
        raise SystemExit("array config must define csv_path")
    csv_path = resolve_workspace_path(cfg_path, csv_value)

    rows = []
    ids = set()
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        missing = REQUIRED_COLUMNS.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"array CSV missing columns: {sorted(missing)}")
        for line_no, row in enumerate(reader, start=2):
            tid = int(row["telescope_id"])
            if tid in ids:
                raise SystemExit(f"duplicate telescope_id={tid} at line {line_no}")
            ids.add(tid)

            tel_cfg = resolve_workspace_path(cfg_path, row["config_path"])
            if not tel_cfg.exists():
                raise SystemExit(f"missing telescope config at line {line_no}: {tel_cfg}")

            position = [
                float(row["position_x_m"]),
                float(row["position_y_m"]),
                float(row["position_z_m"]),
            ]
            az = float(row["pointing_az_deg"])
            el = float(row["pointing_el_deg"])
            if not all(math.isfinite(v) for v in position + [az, el]):
                raise SystemExit(f"non-finite coordinate/pointing at line {line_no}")
            if el < -90.0 or el > 90.0:
                raise SystemExit(f"pointing_el_deg out of range at line {line_no}: {el}")

            rows.append({
                "telescope_id": tid,
                "name": row["name"],
                "config_path": str(tel_cfg),
                "position_x_m": position[0],
                "position_y_m": position[1],
                "position_z_m": position[2],
                "pointing_az_deg": az,
                "pointing_el_deg": el,
            })

    if not rows:
        raise SystemExit("array CSV has no telescope rows")

    if args.summary_csv:
        out_path = Path(args.summary_csv)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    print(f"Array telescopes = {len(rows)}")
    for row in rows:
        print(
            f"  id={row['telescope_id']} name={row['name']} "
            f"pos=({row['position_x_m']:.3f},{row['position_y_m']:.3f},{row['position_z_m']:.3f}) m "
            f"az/el=({row['pointing_az_deg']:.3f},{row['pointing_el_deg']:.3f}) deg"
        )
    if args.summary_csv:
        print(f"Saved summary = {args.summary_csv}")


if __name__ == "__main__":
    main()
