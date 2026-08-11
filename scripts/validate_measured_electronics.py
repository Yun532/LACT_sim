#!/usr/bin/env python3
"""Validate the measured-electronics event 1909 / telescope 19 outputs.

The script deliberately reads the public ROOT, HDF5 and CSV products instead
of internal C++ objects.  It therefore checks the same files consumed by
downstream analysis.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import uproot


IMAGE_FIELDS = {
    "pe": ("image_pe", "fired_pe"),
    "primary_cherenkov_pe": ("image_primary_cherenkov_pe", "primary_cherenkov_pe"),
    "primary_nsb_pe": ("image_primary_nsb_pe", "primary_nsb_pe"),
    "primary_dark_pe": ("image_primary_dark_pe", "primary_dark_pe"),
    "fired_cherenkov_pe": ("image_fired_cherenkov_pe", "fired_cherenkov_pe"),
    "fired_nsb_pe": ("image_fired_nsb_pe", "fired_nsb_pe"),
    "fired_dark_pe": ("image_fired_dark_pe", "fired_dark_pe"),
    "gap_lost_pe": ("image_gap_lost_pe", "gap_lost_pe"),
    "saturation_lost_pe": ("image_saturation_lost_pe", "saturation_lost_pe"),
}

MODES = (
    "saturation_off_waveform_off",
    "saturation_on_waveform_off",
    "measured_waveform_no_nsb",
    "nsb_only",
    "cherenkov_plus_nsb",
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _dense_root(
    branches: dict[str, np.ndarray], branch: str, pixel_axis: np.ndarray
) -> np.ndarray:
    dense = np.zeros(len(pixel_axis), dtype=np.float64)
    columns = {int(pixel_id): column for column, pixel_id in enumerate(pixel_axis)}
    indices = [columns[int(pixel_id)] for pixel_id in branches["pixel_id"][0]]
    dense[indices] = np.asarray(branches[branch][0], dtype=np.float64)
    return dense


def _dense_csv(frame: pd.DataFrame, column: str, pixel_axis: np.ndarray) -> np.ndarray:
    dense = np.zeros(len(pixel_axis), dtype=np.float64)
    columns = {int(pixel_id): index for index, pixel_id in enumerate(pixel_axis)}
    indices = [columns[int(pixel_id)] for pixel_id in frame["pixel_id"]]
    dense[indices] = frame[column].to_numpy(dtype=np.float64)
    return dense


def _root_sparse(tree: uproot.behaviors.TTree.TTree, value: str) -> np.ndarray:
    arrays = tree.arrays(["pixel_id", "time_bin", value], library="np")
    rows = np.empty(
        len(arrays["pixel_id"][0]),
        dtype=[("pixel_id", "<i4"), ("time_bin", "<i4"), ("value", "<f8")],
    )
    rows["pixel_id"] = arrays["pixel_id"][0]
    rows["time_bin"] = arrays["time_bin"][0]
    rows["value"] = arrays[value][0]
    return np.sort(rows, order=["pixel_id", "time_bin"])


def _table_sparse(pixel_id: np.ndarray, time_bin: np.ndarray, value: np.ndarray) -> np.ndarray:
    rows = np.empty(
        len(pixel_id),
        dtype=[("pixel_id", "<i4"), ("time_bin", "<i4"), ("value", "<f8")],
    )
    rows["pixel_id"] = pixel_id
    rows["time_bin"] = time_bin
    rows["value"] = value
    return np.sort(rows, order=["pixel_id", "time_bin"])


def _compare_sparse(left: np.ndarray, right: np.ndarray, label: str) -> float:
    _require(len(left) == len(right), f"{label}: sparse row count differs")
    _require(np.array_equal(left["pixel_id"], right["pixel_id"]), f"{label}: pixel ids differ")
    _require(np.array_equal(left["time_bin"], right["time_bin"]), f"{label}: time bins differ")
    error = float(np.max(np.abs(left["value"] - right["value"]), initial=0.0))
    _require(np.allclose(left["value"], right["value"], rtol=2e-6, atol=2e-7),
             f"{label}: sample values differ (max abs {error})")
    return error


def validate_mode(path: Path) -> dict[str, object]:
    root_path = path / "lact_events.root"
    h5_path = path / "lact_events.h5"
    csv_path = path / "camera_pixels.csv"
    for required in (root_path, h5_path, csv_path):
        _require(required.exists(), f"missing {required}")

    with uproot.open(root_path) as root, h5py.File(h5_path, "r") as h5:
        observations = root["observations"].arrays(library="np")
        _require(int(observations["event_id"][0]) == 1909, "ROOT event id is not 1909")
        _require(int(observations["telescope_id"][0]) == 19, "ROOT telescope id is not 19")
        n_pixels = int(observations["n_pixels_camera"][0])
        _require(n_pixels == 1664, f"camera has {n_pixels}, expected 1664 pixels")

        event_table = h5["events/table"][:]
        image_index = h5["images/index"][:]
        _require(len(event_table) == len(image_index) == 1, "expected one HDF5 image")
        _require(int(event_table["event_id"][0]) == 1909, "HDF5 event id is not 1909")
        _require(int(image_index["telescope_id"][0]) == 19, "HDF5 telescope id is not 19")

        camera_csv = pd.read_csv(csv_path)
        _require(set(camera_csv["event_id"]) == {1909}, "CSV event id is not 1909")
        _require(set(camera_csv["telescope_id"]) == {19}, "CSV telescope id is not 19")

        field_errors: dict[str, float] = {}
        images: dict[str, np.ndarray] = {}
        pixel_axis = np.asarray(h5["images/dense/pixel_id_axis"][:], dtype=np.int64)
        _require(len(pixel_axis) == n_pixels, "HDF5 pixel axis length differs")
        for h5_name, (root_name, csv_name) in IMAGE_FIELDS.items():
            root_image = _dense_root(observations, root_name, pixel_axis)
            h5_image = np.asarray(h5[f"images/dense/{h5_name}"][0], dtype=np.float64)
            csv_image = _dense_csv(camera_csv, csv_name, pixel_axis)
            root_h5 = float(np.max(np.abs(root_image - h5_image), initial=0.0))
            root_csv = float(np.max(np.abs(root_image - csv_image), initial=0.0))
            _require(np.allclose(root_image, h5_image, rtol=0, atol=1e-6),
                     f"{h5_name}: ROOT/HDF5 mismatch")
            _require(np.allclose(root_image, csv_image, rtol=0, atol=1e-9),
                     f"{h5_name}: ROOT/CSV mismatch")
            field_errors[h5_name] = max(root_h5, root_csv)
            images[h5_name] = root_image

        fired_total = (
            images["fired_cherenkov_pe"]
            + images["fired_nsb_pe"]
            + images["fired_dark_pe"]
        )
        primary_total = (
            images["primary_cherenkov_pe"]
            + images["primary_nsb_pe"]
            + images["primary_dark_pe"]
        )
        _require(np.array_equal(images["pe"], fired_total), "final image is not fired-component sum")
        _require(
            np.array_equal(
                primary_total,
                fired_total + images["gap_lost_pe"] + images["saturation_lost_pe"],
            ),
            "primary p.e. conservation failed",
        )

        result: dict[str, object] = {
            "event_id": 1909,
            "telescope_id": 19,
            "n_pixels": n_pixels,
            "triggered": bool(observations["triggered"][0]),
            "primary_cherenkov_pe": float(images["primary_cherenkov_pe"].sum()),
            "primary_nsb_pe": float(images["primary_nsb_pe"].sum()),
            "fired_cherenkov_pe": float(images["fired_cherenkov_pe"].sum()),
            "fired_nsb_pe": float(images["fired_nsb_pe"].sum()),
            "gap_lost_pe": float(images["gap_lost_pe"].sum()),
            "saturation_lost_pe": float(images["saturation_lost_pe"].sum()),
            "root_hdf5_csv_max_abs_error": max(field_errors.values()),
            "waveform_enabled": "waveforms" in root,
        }

        fired_root = root["fired_pe_hits"].arrays(library="np")
        fired_h5 = h5["electronics/fired_pe_hits"][:]
        fired_csv = pd.read_csv(path / "fired_pe.csv")
        _require(len(fired_h5) == len(fired_csv) == len(fired_root["pixel_id"]),
                 "fired-hit row counts differ")
        for field in ("pixel_id", "channel_id", "microcell_id", "origin"):
            _require(np.array_equal(fired_root[field], fired_h5[field]),
                     f"fired {field}: ROOT/HDF5 mismatch")
            _require(np.array_equal(fired_root[field], fired_csv[field].to_numpy()),
                     f"fired {field}: ROOT/CSV mismatch")
        for field in ("time_ns", "fired_pe", "charge_factor", "time_jitter_ns"):
            _require(np.allclose(fired_root[field], fired_h5[field], rtol=0, atol=1e-12),
                     f"fired {field}: ROOT/HDF5 mismatch")
            _require(np.allclose(fired_root[field], fired_csv[field], rtol=0, atol=1e-12),
                     f"fired {field}: ROOT/CSV mismatch")
        result["fired_hit_rows"] = len(fired_h5)

        if "waveforms" in root:
            root_sparse = _root_sparse(root["waveforms"], "sample_value")
            h5_samples = h5["waveforms/samples"][:]
            h5_sparse = _table_sparse(
                h5_samples["pixel_id"], h5_samples["time_bin"], h5_samples["sample_value"]
            )
            waveform_csv = pd.read_csv(path / "waveforms.csv")
            csv_sparse = _table_sparse(
                waveform_csv["pixel_id"].to_numpy(),
                waveform_csv["time_bin"].to_numpy(),
                waveform_csv["sample_value"].to_numpy(),
            )
            result["waveform_root_hdf5_max_abs_error_mv"] = _compare_sparse(
                root_sparse, h5_sparse, "ROOT/HDF5 waveform"
            )
            result["waveform_root_csv_max_abs_error_mv"] = _compare_sparse(
                root_sparse, csv_sparse, "ROOT/CSV waveform"
            )

            root_cfg = root["waveform_config"].arrays(library="np")
            area = float(root_cfg["single_pe_area_mv_ns"][0])
            unit = str(root_cfg["sample_unit"][0])
            if unit.startswith("b'"):
                unit = unit[2:-1]
            h5_unit_raw = h5["waveforms"].attrs["sample_unit"]
            h5_unit = h5_unit_raw.decode() if isinstance(h5_unit_raw, bytes) else str(h5_unit_raw)
            h5_area = float(h5["waveforms"].attrs["single_pe_area_mv_ns"])
            _require(unit == h5_unit == "mV", f"unexpected waveform units: {unit}, {h5_unit}")
            _require(np.isclose(area, h5_area, rtol=0, atol=1e-12), "SPE area differs")
            root_pulse_t = np.asarray(root_cfg["reference_pulse_time_ns"][0], dtype=np.float64)
            root_pulse_a = np.asarray(root_cfg["reference_pulse_amplitude"][0], dtype=np.float64)
            _require(np.array_equal(root_pulse_t, h5["waveforms/reference_pulse_time_ns"][:]),
                     "reference pulse time differs")
            _require(np.array_equal(root_pulse_a, h5["waveforms/reference_pulse_amplitude"][:]),
                     "reference pulse amplitude differs")

            sample_width_ns = float(np.diff(root_cfg["time_edges_ns"][0])[0])
            waveform_charge_pe = float(root_sparse["value"].sum() * sample_width_ns / area)
            fired_charge_pe = float(
                np.sum(fired_root["fired_pe"] * fired_root["charge_factor"])
            )
            closure = waveform_charge_pe - fired_charge_pe
            result.update(
                {
                    "sample_unit": unit,
                    "single_pe_area_mv_ns": area,
                    "reference_pulse_points": len(root_pulse_t),
                    "waveform_charge_pe": waveform_charge_pe,
                    "fired_charge_factor_sum_pe": fired_charge_pe,
                    "full_waveform_direct_closure_pe": closure,
                    "full_waveform_direct_relative_error": closure / fired_charge_pe,
                }
            )
            # The no-NSB Cherenkov pulse is fully contained by the configured
            # window and must close.  Constant-rate NSB is generated across the
            # full readout window, so pulses born near either boundary are
            # intentionally clipped and their all-hit charge is not a closure
            # reference for the finite waveform.
            if path.name == "measured_waveform_no_nsb":
                _require(abs(closure / fired_charge_pe) < 2e-5,
                         f"direct full-waveform charge closure failed: {closure} p.e.")

    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "base",
        nargs="?",
        type=Path,
        default=Path("validation/measured_electronics"),
        help="directory containing the five validation-mode directories",
    )
    parser.add_argument("--json", type=Path, help="optional JSON report path")
    args = parser.parse_args()

    report = {mode: validate_mode(args.base / mode) for mode in MODES}
    _require(report["nsb_only"]["primary_cherenkov_pe"] == 0.0,
             "NSB-only sample contains Cherenkov p.e.")
    _require(report["measured_waveform_no_nsb"]["primary_nsb_pe"] == 0.0,
             "no-NSB sample contains NSB p.e.")
    _require(
        report["saturation_on_waveform_off"]["fired_cherenkov_pe"]
        == report["measured_waveform_no_nsb"]["fired_cherenkov_pe"],
        "enabling the measured waveform changed fired Cherenkov p.e.",
    )
    _require(
        report["cherenkov_plus_nsb"]["primary_cherenkov_pe"]
        == report["measured_waveform_no_nsb"]["primary_cherenkov_pe"],
        "enabling NSB changed primary Cherenkov p.e.",
    )
    text = json.dumps(report, ensure_ascii=False, indent=2)
    print(text)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
