#!/usr/bin/env python3
"""生成 LACT2 带日期的完整仰角—镜片标定表。

输入是已经复核的 w30 标定产物；输出对每个仰角锚点、每片镜子保留一行，
只写入运行时真正需要且随仰角或镜片变化的参数。
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source = args.source
    geometry = read_rows(source / "support_deformation_plus_fixed_pointing.csv")
    facets = {
        int(row["id"]): row for row in read_rows(source / "mirror_facets.csv")
    }
    expected_ids = set(range(54))
    if set(facets) != expected_ids:
        raise RuntimeError("expected facet ids 0--53 in the mirror table")
    anchors = sorted({float(row["elevation_deg"]) for row in geometry})
    if len(geometry) != len(anchors) * 54:
        raise RuntimeError("elevation series does not contain 54 facets per anchor")

    fieldnames = [
        "elevation_deg",
        "id",
        "center_x",
        "center_y",
        "center_z",
        "normal_x",
        "normal_y",
        "normal_z",
        "radius_of_curvature",
    ]
    output_rows = []
    for geometry_row in geometry:
        facet_id = int(geometry_row["id"])
        facet = facets[facet_id]
        row = {name: geometry_row.get(name, "") for name in fieldnames}
        row.update(
            {
                "radius_of_curvature": facet["radius_of_curvature"],
            }
        )
        output_rows.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    print(
        f"wrote {len(output_rows)} rows, {len(anchors)} elevations to "
        f"{args.output}"
    )


if __name__ == "__main__":
    main()
