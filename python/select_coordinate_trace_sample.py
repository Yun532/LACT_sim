#!/usr/bin/env python3
"""Select a deterministic one-telescope sample from a real trace CSV."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def record_path(path: Path) -> str:
    try:
        return path.relative_to(Path.cwd().resolve()).as_posix()
    except ValueError:
        return path.name


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv")
    parser.add_argument("output_csv")
    parser.add_argument("--provenance", required=True)
    parser.add_argument("--eventio", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--cpp-source", required=True)
    parser.add_argument("--summary-csv", default="")
    parser.add_argument("--summary-output", default="")
    parser.add_argument("--rows", type=int, default=64)
    parser.add_argument("--command", default="")
    parser.add_argument("--shower-event-id", type=int)
    parser.add_argument("--event-energy-gev", type=float)
    parser.add_argument("--event-arrival-az-deg", type=float)
    parser.add_argument("--event-altitude-deg", type=float)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_csv = Path(args.input_csv).resolve()
    output_csv = Path(args.output_csv).resolve()
    provenance_path = Path(args.provenance).resolve()
    eventio_path = Path(args.eventio).resolve()
    config_path = Path(args.config).resolve()
    cpp_path = Path(args.cpp_source).resolve()
    summary_path = Path(args.summary_csv).resolve() if args.summary_csv else None
    summary_output = Path(args.summary_output).resolve() if args.summary_output else None
    if bool(summary_path) != bool(summary_output):
        raise SystemExit("--summary-csv and --summary-output must be used together")
    if args.rows <= 0:
        raise SystemExit("--rows must be positive")

    group_counts: Counter[tuple[str, str]] = Counter()
    total_rows = 0
    with input_csv.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        required = {
            "event_id", "telescope_id", "photon_index",
            "input_x_m", "input_y_m", "input_z_m",
            "input_dir_x", "input_dir_y", "input_dir_z",
            "mirror_x_m", "mirror_y_m", "mirror_z_m",
            "surface_x_m", "surface_y_m", "surface_z_m", "u_m", "v_m",
        }
        missing = sorted(required - set(fieldnames))
        if missing:
            raise SystemExit("input CSV is missing: " + ", ".join(missing))
        for row in reader:
            group_counts[(row["event_id"], row["telescope_id"])] += 1
            total_rows += 1

    if not group_counts:
        raise SystemExit("input trace CSV has no rows")
    selected_group, selected_count = max(
        group_counts.items(), key=lambda item: (item[1], item[0])
    )
    sample_count = min(args.rows, selected_count)
    if sample_count == 1:
        selected_indices = {0}
    else:
        selected_indices = {
            round(index * (selected_count - 1) / (sample_count - 1))
            for index in range(sample_count)
        }

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    selected_index = 0
    written = 0
    with (
        input_csv.open(newline="") as source,
        output_csv.open("w", newline="", encoding="utf-8") as target,
    ):
        reader = csv.DictReader(source)
        writer = csv.DictWriter(target, fieldnames=reader.fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in reader:
            if (row["event_id"], row["telescope_id"]) != selected_group:
                continue
            if selected_index in selected_indices:
                writer.writerow(row)
                written += 1
            selected_index += 1

    summary_record = None
    if summary_path and summary_output:
        with summary_path.open(newline="") as source:
            reader = csv.DictReader(source)
            matching = [
                row for row in reader
                if (row.get("event_id"), row.get("telescope_id")) == selected_group
            ]
            if len(matching) != 1:
                raise SystemExit("expected exactly one matching summary row")
            summary_output.parent.mkdir(parents=True, exist_ok=True)
            with summary_output.open("w", newline="", encoding="utf-8") as target:
                writer = csv.DictWriter(
                    target, fieldnames=reader.fieldnames, lineterminator="\n"
                )
                writer.writeheader()
                writer.writerow(matching[0])
        summary_record = {
            "file": summary_output.name,
            "source_sha256": sha256(summary_path),
            "sample_sha256": sha256(summary_output),
        }

    provenance = {
        "description": "Deterministic subset of actual run_corsika_trace whiteboard output",
        "run_command": args.command,
        "source_eventio": {
            "file": eventio_path.name,
            "size_bytes": eventio_path.stat().st_size,
            "sha256": sha256(eventio_path),
        },
        "program_inputs": {
            "config": record_path(config_path),
            "config_sha256": sha256(config_path),
            "cpp_source": record_path(cpp_path),
            "cpp_source_sha256": sha256(cpp_path),
        },
        "full_trace": {
            "file": input_csv.name,
            "rows": total_rows,
            "sha256": sha256(input_csv),
        },
        "selection": {
            "event_id": selected_group[0],
            "telescope_id": selected_group[1],
            "available_rows": selected_count,
            "sample_rows": written,
            "method": "uniform indices across the event/telescope group with most saved output hits",
        },
        "sample": {
            "file": output_csv.name,
            "sha256": sha256(output_csv),
        },
    }
    if args.shower_event_id is not None:
        provenance["corsika_event"] = {
            "shower_event_id": args.shower_event_id,
            "energy_gev": args.event_energy_gev,
            "arrival_azimuth_north_to_east_deg": args.event_arrival_az_deg,
            "altitude_deg": args.event_altitude_deg,
        }
    if summary_record:
        provenance["summary"] = summary_record
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(
        f"selected event={selected_group[0]} telescope={selected_group[1]} "
        f"rows={written}/{selected_count} from full_rows={total_rows}"
    )


if __name__ == "__main__":
    main()
