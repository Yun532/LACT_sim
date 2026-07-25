#!/usr/bin/env python3
"""Summarize and sanity-check a LACT_sim trust-validation HDF5 result."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Any

import h5py
import numpy as np


def scalar(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        return value.item()
    return value


def stats(values: np.ndarray) -> dict[str, Any]:
    data = np.asarray(values, dtype=np.float64)
    finite = np.isfinite(data)
    valid = data[finite]
    result: dict[str, Any] = {
        "count": int(data.size),
        "finite": int(np.count_nonzero(finite)),
        "nonfinite": int(data.size - np.count_nonzero(finite)),
        "negative": int(np.count_nonzero(valid < 0.0)),
    }
    if valid.size:
        result.update(
            {
                "sum": float(np.sum(valid, dtype=np.float64)),
                "min": float(np.min(valid)),
                "max": float(np.max(valid)),
                "mean": float(np.mean(valid)),
                "median": float(np.median(valid)),
            }
        )
    return result


def dataset_shapes(handle: h5py.File) -> dict[str, list[int]]:
    shapes: dict[str, list[int]] = {}

    def visit(name: str, obj: Any) -> None:
        if isinstance(obj, h5py.Dataset):
            shapes[f"/{name}"] = [int(v) for v in obj.shape]

    handle.visititems(visit)
    return shapes


def structured_field(handle: h5py.File, path: str, name: str) -> np.ndarray:
    data = handle[path][()]
    if data.dtype.names is None or name not in data.dtype.names:
        raise KeyError(f"{path} has no field {name}")
    return np.asarray(data[name])


def string_attribute(group: h5py.Group, name: str) -> str:
    return str(scalar(group.attrs[name]))


def verify_hdf5_trigger_from_waveform(handle: h5py.File) -> dict[str, Any]:
    required = (
        "/waveforms/samples",
        "/waveforms/pixel_id_axis",
        "/waveforms/reference_time_ns",
        "/images/index",
        "/trigger/telescope",
    )
    if any(path not in handle for path in required):
        return {"available": False}

    samples = handle["/waveforms/samples"][()]
    if samples.dtype.names is None or "pe" not in samples.dtype.names:
        return {"available": False}
    images = handle["/images/index"][()]
    trigger_rows = handle["/trigger/telescope"][()]
    pixel_ids = np.asarray(
        handle["/waveforms/pixel_id_axis"][()], dtype=np.int64
    )
    reference_times = np.asarray(
        handle["/waveforms/reference_time_ns"][()], dtype=np.float64
    )
    n_pixels = int(pixel_ids.size)
    n_bins = int(handle["/waveforms/time_centers_ns"].shape[0])
    trigger_meta = handle["/metadata/trigger"]
    waveform_meta = handle["/metadata/waveform"]
    threshold = float(string_attribute(trigger_meta, "pixel_threshold_pe"))
    multiplicity = int(string_attribute(trigger_meta, "camera_multiplicity"))
    coincidence_ns = float(
        string_attribute(trigger_meta, "camera_coincidence_window_ns")
    )
    bin_width_ns = float(string_attribute(waveform_meta, "time_bin_width_ns"))
    window_start_ns = float(
        string_attribute(waveform_meta, "time_window_start_ns")
    )
    window_bins = min(
        n_bins, max(1, int(math.ceil(coincidence_ns / bin_width_ns)))
    )
    pixel_to_col = {int(pixel_id): col for col, pixel_id in enumerate(pixel_ids)}
    samples_by_image: dict[int, list[Any]] = {}
    for row in samples:
        samples_by_image.setdefault(int(row["image_index"]), []).append(row)
    trigger_by_key = {
        (int(row["event_id"]), int(row["telescope_id"])): row
        for row in trigger_rows
    }

    failures: list[str] = []
    checked = 0
    for row_index, image in enumerate(images):
        image_index = int(image["image_index"])
        grid = np.zeros((n_pixels, n_bins), dtype=np.float64)
        for sample in samples_by_image.get(image_index, []):
            col = pixel_to_col.get(int(sample["pixel_id"]))
            bin_index = int(sample["time_bin"])
            if col is not None and 0 <= bin_index < n_bins:
                grid[col, bin_index] += float(sample["pe"])

        multiplicities = np.zeros(n_bins, dtype=np.int64)
        for col in range(n_pixels):
            window_pe = float(np.sum(grid[col, :window_bins]))
            if window_pe >= threshold:
                multiplicities[0] += 1
            for start in range(1, n_bins):
                window_pe -= grid[col, start - 1]
                add_bin = start + window_bins - 1
                if add_bin < n_bins:
                    window_pe += grid[col, add_bin]
                if window_pe >= threshold:
                    multiplicities[start] += 1

        triggered_starts = np.flatnonzero(multiplicities >= multiplicity)
        triggered = triggered_starts.size > 0
        first_start = int(triggered_starts[0]) if triggered else 0
        max_start = int(np.argmax(multiplicities))

        def window_center(start: int) -> float:
            center_offset = min(n_bins - 1 - start, window_bins // 2)
            first_bin_center = (
                reference_times[row_index] + window_start_ns + 0.5 * bin_width_ns
            )
            return first_bin_center + (start + center_offset) * bin_width_ns

        expected_trigger_time = window_center(first_start if triggered else max_start)
        key = (int(image["event_id"]), int(image["telescope_id"]))
        saved = trigger_by_key.get(key)
        if saved is None:
            failures.append(f"missing trigger row for event/telescope {key}")
            continue
        saved_triggered = bool(saved["triggered"])
        saved_multiplicity = int(saved["n_pixels_above_threshold"])
        saved_trigger_time = float(saved["trigger_time_ns"])
        if saved_triggered != triggered:
            failures.append(f"{key}: triggered {saved_triggered} != {triggered}")
        if saved_multiplicity != int(np.max(multiplicities)):
            failures.append(
                f"{key}: multiplicity {saved_multiplicity} != "
                f"{int(np.max(multiplicities))}"
            )
        if not math.isclose(
            saved_trigger_time, expected_trigger_time, rel_tol=0.0, abs_tol=2.0e-5
        ):
            failures.append(
                f"{key}: trigger time {saved_trigger_time} != "
                f"{expected_trigger_time}"
            )
        checked += 1
    return {
        "available": True,
        "checked": checked,
        "passed": not failures,
        "failure_count": len(failures),
        "failures": failures[:20],
    }


def root_trigger_rows(root_path: Path) -> dict[tuple[int, int], dict[str, Any]]:
    import ROOT  # type: ignore

    root_file = ROOT.TFile.Open(str(root_path))
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"failed to open ROOT file: {root_path}")
    tree = root_file.Get("observations")
    if tree is None:
        raise RuntimeError("ROOT file has no observations tree")
    rows: dict[tuple[int, int], dict[str, Any]] = {}
    for entry in tree:
        key = (int(entry.event_id), int(entry.telescope_id))
        rows[key] = {
            "triggered": bool(entry.triggered),
            "n_pixels_above_threshold": int(entry.n_pixels_above_threshold),
            "trigger_time_ns": float(entry.trigger_time_ns),
            "trigger_first_time_ns": float(entry.trigger_first_time_ns),
            "trigger_max_multiplicity_time_ns": float(
                entry.trigger_max_multiplicity_time_ns
            ),
        }
    root_file.Close()
    return rows


def compare_root_hdf5_triggers(
    root_path: Path, hdf5_rows: np.ndarray
) -> dict[str, Any]:
    root_rows = root_trigger_rows(root_path)
    failures: list[str] = []
    checked = 0
    for row in hdf5_rows:
        key = (int(row["event_id"]), int(row["telescope_id"]))
        root_row = root_rows.get(key)
        if root_row is None:
            failures.append(f"ROOT missing event/telescope {key}")
            continue
        if bool(row["triggered"]) != root_row["triggered"]:
            failures.append(f"{key}: triggered differs")
        if (
            int(row["n_pixels_above_threshold"])
            != root_row["n_pixels_above_threshold"]
        ):
            failures.append(f"{key}: maximum multiplicity differs")
        hdf5_triggered = bool(row["triggered"])
        hdf5_trigger_time = float(row["trigger_time_ns"])
        if hdf5_triggered and not math.isclose(
            hdf5_trigger_time,
            root_row["trigger_time_ns"],
            rel_tol=0.0,
            abs_tol=2.0e-5,
        ):
            failures.append(
                f"{key}: trigger time HDF5={hdf5_trigger_time} "
                f"ROOT={root_row['trigger_time_ns']}"
            )
        if not hdf5_triggered and math.isfinite(root_row["trigger_time_ns"]):
            failures.append(
                f"{key}: untriggered ROOT trigger time is not NaN"
            )
        if not math.isclose(
            float(row["trigger_max_multiplicity_time_ns"]),
            root_row["trigger_max_multiplicity_time_ns"],
            rel_tol=0.0,
            abs_tol=2.0e-5,
        ):
            failures.append(f"{key}: maximum-multiplicity time differs")
        hdf5_first = float(row["trigger_first_time_ns"])
        root_first = root_row["trigger_first_time_ns"]
        if hdf5_triggered:
            if not math.isclose(
                hdf5_first, root_first, rel_tol=0.0, abs_tol=2.0e-5
            ):
                failures.append(f"{key}: first trigger time differs")
        elif math.isfinite(hdf5_first) or math.isfinite(root_first):
            failures.append(f"{key}: untriggered first time is not NaN")
        checked += 1
    extra = sorted(set(root_rows) - {
        (int(row["event_id"]), int(row["telescope_id"])) for row in hdf5_rows
    })
    failures.extend(f"ROOT has extra event/telescope {key}" for key in extra[:20])
    return {
        "checked": checked,
        "passed": not failures,
        "failure_count": len(failures),
        "failures": failures[:20],
    }


def parse_event_totals(log_path: Path) -> dict[str, Any]:
    totals: dict[str, float] = {}
    shower_lines = 0
    key_value = re.compile(r"([A-Za-z_]+)=([-+0-9.eE]+)")
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "event_total" not in line:
            continue
        shower_lines += 1
        for key, text in key_value.findall(line):
            if key in {"output_events", "telescopes"}:
                continue
            totals[key] = totals.get(key, 0.0) + float(text)
    return {"completed_shower_lines": shower_lines, "sums": totals}


def analyze(
    hdf5_path: Path, log_path: Path | None, root_path: Path | None
) -> dict[str, Any]:
    report: dict[str, Any] = {"hdf5_path": str(hdf5_path)}
    with h5py.File(hdf5_path, "r") as handle:
        report["provenance"] = {
            key: scalar(handle.attrs.get(key, ""))
            for key in ("producer_version", "source_path", "source_sha256")
        }
        shapes = dataset_shapes(handle)
        report["dataset_shapes"] = shapes
        report["counts"] = {
            "corsika_event_rows": int(handle["/events/corsika"].shape[0]),
            "corsika_shower_headers": int(handle["/events/corsika_showers"].shape[0]),
            "image_rows": int(handle["/images/index"].shape[0]),
            "camera_pixels": int(handle["/camera/pixels"].shape[0]),
            "telescopes": int(handle["/telescopes/table"].shape[0]),
            "telescope_trigger_rows": int(handle["/trigger/telescope"].shape[0]),
            "array_trigger_rows": int(handle["/trigger/array"].shape[0]),
        }

        telescope_triggered = structured_field(
            handle, "/trigger/telescope", "triggered"
        ).astype(bool)
        array_triggered = structured_field(
            handle, "/trigger/array", "array_triggered"
        ).astype(bool)
        report["trigger"] = {
            "telescope_triggered": int(np.count_nonzero(telescope_triggered)),
            "telescope_total": int(telescope_triggered.size),
            "telescope_fraction": (
                float(np.mean(telescope_triggered)) if telescope_triggered.size else 0.0
            ),
            "array_triggered": int(np.count_nonzero(array_triggered)),
            "array_total": int(array_triggered.size),
            "array_fraction": (
                float(np.mean(array_triggered)) if array_triggered.size else 0.0
            ),
        }

        pe = np.asarray(handle["/images/dense/pe"][()], dtype=np.float64)
        cherenkov = np.asarray(
            handle["/images/dense/cherenkov_pe"][()], dtype=np.float64
        )
        nsb = np.asarray(handle["/images/dense/nsb_pe"][()], dtype=np.float64)
        decomposition = pe - cherenkov - nsb
        report["images"] = {
            "pe": stats(pe),
            "cherenkov_pe": stats(cherenkov),
            "nsb_pe": stats(nsb),
            "per_image_total_pe": stats(np.sum(pe, axis=1, dtype=np.float64)),
            "per_image_cherenkov_pe": stats(
                np.sum(cherenkov, axis=1, dtype=np.float64)
            ),
            "per_image_nsb_pe": stats(np.sum(nsb, axis=1, dtype=np.float64)),
            "decomposition_max_abs_error": (
                float(np.max(np.abs(decomposition))) if decomposition.size else 0.0
            ),
        }

        waveform_report: dict[str, Any] = {"available": False}
        if "/waveforms/samples" in handle:
            waveform = handle["/waveforms/samples"][()]
            waveform_report = {
                "available": True,
                "rows": int(waveform.size),
            }
            if waveform.dtype.names:
                for field in ("pe", "cherenkov_pe", "nsb_pe"):
                    if field in waveform.dtype.names:
                        waveform_report[field] = stats(waveform[field])
        report["waveforms"] = waveform_report
        report["trigger_from_saved_waveform"] = (
            verify_hdf5_trigger_from_waveform(handle)
        )
        if root_path is not None:
            report["root_hdf5_trigger"] = compare_root_hdf5_triggers(
                root_path, handle["/trigger/telescope"][()]
            )

        failures: list[str] = []
        for name, values in (
            ("pe", pe),
            ("cherenkov_pe", cherenkov),
            ("nsb_pe", nsb),
        ):
            if not np.all(np.isfinite(values)):
                failures.append(f"{name} contains non-finite values")
            if np.any(values < 0.0):
                failures.append(f"{name} contains negative values")
        if report["images"]["decomposition_max_abs_error"] > 1.0e-5:
            failures.append("pe != cherenkov_pe + nsb_pe")
        expected_sha = str(report["provenance"]["source_sha256"])
        if len(expected_sha) != 64:
            failures.append("source_sha256 is absent or malformed")
        trigger_check = report["trigger_from_saved_waveform"]
        if trigger_check.get("available") and not trigger_check.get("passed"):
            failures.append("saved waveform does not reproduce HDF5 trigger")
        root_check = report.get("root_hdf5_trigger")
        if root_check is not None and not root_check.get("passed"):
            failures.append("ROOT and HDF5 trigger rows differ")
        report["sanity"] = {"passed": not failures, "failures": failures}

    if log_path is not None:
        report["stdout_event_totals"] = parse_event_totals(log_path)
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("hdf5", type=Path)
    parser.add_argument("--stdout-log", type=Path)
    parser.add_argument("--root", type=Path)
    parser.add_argument("--pretty", action="store_true")
    args = parser.parse_args()

    report = analyze(args.hdf5, args.stdout_log, args.root)
    print(
        json.dumps(
            report,
            ensure_ascii=False,
            indent=2 if args.pretty else None,
            sort_keys=True,
            allow_nan=False,
        )
    )
    return 0 if report["sanity"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
