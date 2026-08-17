#!/usr/bin/env python3
"""Generate a reproducible external-style photon CSV for array tests."""

import argparse
import csv
import math
import random
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Generate array photon CSV test input.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--telescopes", default="0,1,2")
    parser.add_argument("--events", default="1")
    parser.add_argument("--photons-per-telescope", type=int, default=50000)
    parser.add_argument("--beam-radius-m", type=float, default=4.0)
    parser.add_argument("--source-z-m", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=1229)
    args = parser.parse_args()

    telescope_ids = [int(x.strip()) for x in args.telescopes.split(",") if x.strip()]
    event_ids = [int(x.strip()) for x in args.events.split(",") if x.strip()]
    if not telescope_ids:
        raise SystemExit("no telescope ids requested")
    if not event_ids:
        raise SystemExit("no event ids requested")
    if args.photons_per_telescope <= 0:
        raise SystemExit("--photons-per-telescope must be > 0")

    rng = random.Random(args.seed)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "x_cm",
                "y_cm",
                "z_cm",
                "ux",
                "uy",
                "uz",
                "time_ns",
                "wavelength_nm",
                "weight",
                "multiplicity",
                "event",
                "tel_id",
            ],
        )
        writer.writeheader()
        for event_id in event_ids:
            for tel_id in telescope_ids:
                for _ in range(args.photons_per_telescope):
                    rho = args.beam_radius_m * math.sqrt(rng.random())
                    phi = 2.0 * math.pi * rng.random()
                    x_m = rho * math.cos(phi)
                    y_m = rho * math.sin(phi)
                    writer.writerow({
                        "x_cm": f"{x_m * 100.0:.9g}",
                        "y_cm": f"{y_m * 100.0:.9g}",
                        "z_cm": f"{args.source_z_m * 100.0:.9g}",
                        "ux": "0",
                        "uy": "0",
                        "uz": "-1",
                        "time_ns": "0",
                        "wavelength_nm": "400",
                        "weight": "1",
                        "multiplicity": "1",
                        "event": str(event_id),
                        "tel_id": str(tel_id),
                    })

    total = len(event_ids) * len(telescope_ids) * args.photons_per_telescope
    print(f"wrote {out_path} with {total} photons")


if __name__ == "__main__":
    main()
