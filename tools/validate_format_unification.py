#!/usr/bin/env python3
"""Cross-check one canonical electronics result in ROOT, HDF5 and CSV."""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import uproot


def key(event_id: int, telescope_id: int) -> tuple[int, int]:
    return int(event_id), int(telescope_id)


def assert_close(name: str, left, right, atol: float = 2e-5) -> None:
    a = np.asarray(left)
    b = np.asarray(right)
    if a.shape != b.shape or not np.allclose(a, b, rtol=0.0, atol=atol,
                                              equal_nan=True):
        delta = np.nanmax(np.abs(a - b)) if a.shape == b.shape and a.size else np.nan
        raise AssertionError(
            f"{name}: shape/value mismatch {a.shape} vs {b.shape}; max delta={delta}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    directory = args.directory

    pixels_csv = pd.read_csv(directory / "camera_pixels.csv")
    waveform_csv = pd.read_csv(directory / "waveforms.csv")
    trigger_csv = pd.read_csv(directory / "triggers.csv")
    summary_csv = pd.read_csv(directory / "trace_summary.csv")

    with uproot.open(directory / "lact_events.root") as root:
        observations = root["observations"].arrays(library="ak")
        root_waveforms = root["waveforms"].arrays(library="ak")
        root_primary = root["primary_pe_hits"].arrays(library="np")
        root_fired = root["fired_pe_hits"].arrays(library="np")
        waveform_config = root["waveform_config"].arrays(library="np")

    with h5py.File(directory / "lact_events.h5", "r") as h5:
        image_index = h5["images/index"][()]
        pixel_axis = h5["images/dense/pixel_id_axis"][()].astype(int)
        h5_image = {
            name: h5[f"images/dense/{name}"][()]
            for name in (
                "pe", "primary_cherenkov_pe", "primary_nsb_pe",
                "primary_dark_pe",
                "fired_cherenkov_pe", "fired_nsb_pe", "fired_dark_pe",
                "gap_lost_pe", "saturation_lost_pe",
            )
        }
        h5_trigger = h5["trigger/telescope"][()]
        h5_waveform_samples = h5["waveforms/samples"][()]
        h5_waveform_reference = h5["waveforms/reference_time_ns"][()]
        h5_time_centers = h5["waveforms/time_centers_ns"][()]
        h5_primary = h5["electronics/primary_pe_hits"][()]
        h5_fired = h5["electronics/fired_pe_hits"][()]

    root_observation_by_key = {
        key(observations.event_id[i], observations.telescope_id[i]): i
        for i in range(len(observations))
    }
    h5_image_by_key = {
        key(row["event_id"], row["telescope_id"]): int(row["image_index"])
        for row in image_index
    }
    csv_keys = set(zip(trigger_csv.event_id.astype(int),
                       trigger_csv.telescope_id.astype(int)))
    expected_keys = set(root_observation_by_key)
    if expected_keys != set(h5_image_by_key) or expected_keys != csv_keys:
        raise AssertionError("event/telescope key sets differ across formats")

    component_mapping = {
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
    max_image_delta = 0.0
    for event_key in sorted(expected_keys):
        root_row = root_observation_by_key[event_key]
        h5_row = h5_image_by_key[event_key]
        csv_rows = pixels_csv[
            (pixels_csv.event_id == event_key[0]) &
            (pixels_csv.telescope_id == event_key[1])
        ].set_index("pixel_id")
        root_pixel_ids = np.asarray(observations.pixel_id[root_row], dtype=int)

        for h5_name, (root_name, csv_name) in component_mapping.items():
            h5_values = np.asarray(h5_image[h5_name][h5_row], dtype=float)
            csv_values = np.zeros(len(pixel_axis), dtype=float)
            for column, pixel_id in enumerate(pixel_axis):
                if pixel_id in csv_rows.index:
                    csv_values[column] = float(csv_rows.at[pixel_id, csv_name])
            assert_close(f"{event_key} HDF5/CSV {h5_name}", h5_values, csv_values)
            if root_name is not None:
                root_values = np.zeros(len(pixel_axis), dtype=float)
                stored = np.asarray(observations[root_name][root_row], dtype=float)
                for pixel_id, value in zip(root_pixel_ids, stored):
                    root_values[np.searchsorted(pixel_axis, pixel_id)] = value
                assert_close(f"{event_key} ROOT/HDF5 {h5_name}", root_values, h5_values)
                if root_values.size:
                    max_image_delta = max(
                        max_image_delta,
                        float(np.max(np.abs(root_values - h5_values))),
                    )

        root_reference = float(observations.reference_time_ns[root_row])
        h5_reference = float(image_index[h5_row]["reference_time_ns"])
        csv_reference = float(trigger_csv[
            (trigger_csv.event_id == event_key[0]) &
            (trigger_csv.telescope_id == event_key[1])
        ].iloc[0].reference_time_ns)
        assert_close(f"{event_key} reference time",
                     [root_reference, h5_reference],
                     [csv_reference, csv_reference], atol=1e-9)

    def root_waveform_rows() -> pd.DataFrame:
        records = []
        for row in range(len(root_waveforms)):
            reference = float(root_waveforms.reference_time_ns[row])
            for pixel_id, time_bin, value in zip(
                root_waveforms.pixel_id[row],
                root_waveforms.time_bin[row],
                root_waveforms.sample_value[row],
            ):
                records.append((
                    int(root_waveforms.event_id[row]),
                    int(root_waveforms.telescope_id[row]),
                    reference, int(pixel_id), int(time_bin), float(value),
                ))
        return pd.DataFrame(records, columns=(
            "event_id", "telescope_id", "reference_time_ns", "pixel_id",
            "time_bin", "sample_value",
        ))

    root_waveform = root_waveform_rows()
    h5_image_key_by_index = {
        int(row["image_index"]): key(row["event_id"], row["telescope_id"])
        for row in image_index
    }
    h5_waveform_records = []
    for sample in h5_waveform_samples:
        event_key = h5_image_key_by_index[int(sample["image_index"])]
        row = h5_image_by_key[event_key]
        h5_waveform_records.append((
            event_key[0], event_key[1], float(h5_waveform_reference[row]),
            int(sample["pixel_id"]), int(sample["time_bin"]),
            float(sample["sample_value"]),
        ))
    h5_waveform = pd.DataFrame(h5_waveform_records, columns=root_waveform.columns)

    columns = ["event_id", "telescope_id", "pixel_id", "time_bin"]
    root_waveform = root_waveform.sort_values(columns).reset_index(drop=True)
    h5_waveform = h5_waveform.sort_values(columns).reset_index(drop=True)
    waveform_csv_compare = waveform_csv[[
        "event_id", "telescope_id", "reference_time_ns", "pixel_id",
        "time_bin", "sample_value",
    ]].sort_values(columns).reset_index(drop=True)
    if not root_waveform[columns].equals(h5_waveform[columns]) or not \
            root_waveform[columns].equals(waveform_csv_compare[columns]):
        raise AssertionError("waveform sparse coordinates differ across formats")
    assert_close("waveform ROOT/HDF5", root_waveform.sample_value,
                 h5_waveform.sample_value)
    assert_close("waveform ROOT/CSV", root_waveform.sample_value,
                 waveform_csv_compare.sample_value)
    assert_close("waveform references ROOT/HDF5",
                 root_waveform.reference_time_ns,
                 h5_waveform.reference_time_ns, atol=1e-9)
    assert_close("waveform references ROOT/CSV",
                 root_waveform.reference_time_ns,
                 waveform_csv_compare.reference_time_ns, atol=1e-9)

    primary_csv = pd.read_csv(directory / "primary_pe.csv")
    fired_csv = pd.read_csv(directory / "fired_pe.csv")
    if len(root_primary["event_id"]) != len(h5_primary) or \
            len(h5_primary) != len(primary_csv):
        raise AssertionError("primary p.e. sequence lengths differ")
    if len(root_fired["event_id"]) != len(h5_fired) or \
            len(h5_fired) != len(fired_csv):
        raise AssertionError("fired p.e. sequence lengths differ")
    for name in ("event_id", "telescope_id", "pixel_id", "origin"):
        if not np.array_equal(root_primary[name], h5_primary[name]) or not \
                np.array_equal(root_primary[name], primary_csv[name].to_numpy()):
            raise AssertionError(f"primary p.e. {name} differs")
    for name in ("reference_time_ns", "time_ns", "sensor_x_m", "sensor_y_m",
                 "wavelength_nm", "primary_pe"):
        assert_close(f"primary p.e. {name} ROOT/HDF5",
                     root_primary[name], h5_primary[name], atol=1e-9)
        assert_close(f"primary p.e. {name} ROOT/CSV",
                     root_primary[name], primary_csv[name], atol=1e-9)
    for name in ("event_id", "telescope_id", "pixel_id", "channel_id",
                 "microcell_id", "origin"):
        if not np.array_equal(root_fired[name], h5_fired[name]) or not \
                np.array_equal(root_fired[name], fired_csv[name].to_numpy()):
            raise AssertionError(f"fired p.e. {name} differs")
    for name in ("reference_time_ns", "time_ns", "fired_pe"):
        assert_close(f"fired p.e. {name} ROOT/HDF5",
                     root_fired[name], h5_fired[name], atol=1e-9)
        assert_close(f"fired p.e. {name} ROOT/CSV",
                     root_fired[name], fired_csv[name], atol=1e-9)

    h5_trigger_by_key = {
        key(row["event_id"], row["telescope_id"]): row
        for row in h5_trigger
    }
    for event_key in sorted(expected_keys):
        root_row = root_observation_by_key[event_key]
        h5_row = h5_trigger_by_key[event_key]
        csv_row = trigger_csv[
            (trigger_csv.event_id == event_key[0]) &
            (trigger_csv.telescope_id == event_key[1])
        ].iloc[0]
        values = (
            ("triggered", bool(observations.triggered[root_row]),
             bool(h5_row["triggered"]), bool(csv_row.triggered)),
            ("n_pixels_above_threshold",
             int(observations.n_pixels_above_threshold[root_row]),
             int(h5_row["n_pixels_above_threshold"]),
             int(csv_row.n_pixels_above_threshold)),
        )
        for name, root_value, h5_value, csv_value in values:
            if root_value != h5_value or root_value != csv_value:
                raise AssertionError(f"{event_key} {name} differs")
        for name in ("trigger_time_ns", "geometric_delay_ns",
                     "coincidence_time_ns"):
            assert_close(f"{event_key} {name}",
                         [float(observations[name][root_row]), float(h5_row[name])],
                         [float(csv_row[name]), float(csv_row[name])], atol=1e-9)

    if not np.allclose(
        np.asarray(waveform_config["time_centers_ns"][0]),
        h5_time_centers, rtol=0.0, atol=1e-12,
    ):
        raise AssertionError("ROOT/HDF5 waveform time axes differ")

    print(f"events_telescopes={len(expected_keys)}")
    print(f"pixels_csv_rows={len(pixels_csv)}")
    print(f"waveform_sparse_samples={len(root_waveform)}")
    print(f"primary_pe_hits={len(h5_primary)}")
    print(f"fired_pe_hits={len(h5_fired)}")
    print(f"triggered_telescopes={int(trigger_csv.triggered.sum())}")
    print(f"max_root_hdf5_image_delta={max_image_delta:.9g}")
    print("format_unification=PASS")


if __name__ == "__main__":
    main()
