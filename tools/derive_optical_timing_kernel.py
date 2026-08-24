#!/usr/bin/env python3
"""Compress a large optical-ray hit table into a reproducible timing mixture.

The output keeps one Gaussian component per mirror facet.  Component weights
come from the accepted on-axis ray counts; within-facet means and standard
deviations preserve the optical path-time distribution without committing the
full ray table to git.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("hits_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--provenance-json", type=Path, required=True)
    parser.add_argument("--source-config", type=Path)
    parser.add_argument("--mirror-csv", type=Path)
    args = parser.parse_args()

    columns = ["hit_surface", "mirror_id", "time_ns"]
    chunks = []
    for frame in pd.read_csv(args.hits_csv, usecols=columns, chunksize=200_000):
        frame = frame.loc[(frame.hit_surface == 1) & frame.time_ns.notna()]
        grouped = frame.groupby("mirror_id").time_ns.agg(
            count="count", sum_time="sum", sum_time2=lambda x: (x * x).sum(),
            min_time_ns="min", max_time_ns="max")
        chunks.append(grouped)
    reduced = pd.concat(chunks).groupby(level=0).sum()
    # Min/max cannot be summed across chunks; recover them in one bounded pass.
    extrema = pd.concat([
        frame.loc[(frame.hit_surface == 1) & frame.time_ns.notna()]
        .groupby("mirror_id").time_ns.agg(min_time_ns="min", max_time_ns="max")
        for frame in pd.read_csv(args.hits_csv, usecols=columns, chunksize=200_000)
    ]).groupby(level=0).agg({"min_time_ns": "min", "max_time_ns": "max"})
    reduced[["min_time_ns", "max_time_ns"]] = extrema
    reduced["mean_time_ns"] = reduced.sum_time / reduced["count"]
    variance = (
        reduced.sum_time2 - reduced.sum_time**2 / reduced["count"]
    ) / (reduced["count"] - 1).clip(lower=1)
    reduced["std_time_ns"] = variance.clip(lower=0).pow(0.5)
    reduced["weight"] = reduced["count"] / reduced["count"].sum()
    output = reduced.reset_index()[[
        "mirror_id", "count", "weight", "mean_time_ns", "std_time_ns",
        "min_time_ns", "max_time_ns",
    ]]
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output_csv, index=False, float_format="%.12g")

    provenance = {
        "description": "On-axis 400 nm parallel-ray optical arrival-time mixture; one component per facet.",
        "source_hits_csv": str(args.hits_csv.resolve()),
        "source_hits_sha256": sha256(args.hits_csv),
        "source_rows_hit_surface": int(output["count"].sum()),
        "components": int(len(output)),
        "derivation_script": str(Path(__file__).resolve()),
    }
    for label, path in (("source_config", args.source_config),
                        ("mirror_csv", args.mirror_csv)):
        if path is not None:
            provenance[label] = str(path.resolve())
            provenance[f"{label}_sha256"] = sha256(path)
    args.provenance_json.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_json.write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8")


if __name__ == "__main__":
    main()
