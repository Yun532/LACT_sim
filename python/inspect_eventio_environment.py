#!/usr/bin/env python3
"""Inspect local CORSIKA/EventIO tools and example files."""

import argparse
import shutil
from pathlib import Path


PATTERNS = [
    "*eventio*",
    "*corsika*",
    "*simtel*",
    "*iact*",
    "read_iact",
    "select_iact",
]


def main():
    parser = argparse.ArgumentParser(description="Find local EventIO/CORSIKA tooling.")
    parser.add_argument("--root", default=".")
    parser.add_argument("--max-results", type=int, default=80)
    args = parser.parse_args()

    root = Path(args.root)
    print(f"Search root = {root}")
    for exe in ["read_iact", "select_iact", "extract_simtel", "python3"]:
        found = shutil.which(exe)
        print(f"PATH {exe}: {found or 'not found'}")

    seen = set()
    results = []
    for pattern in PATTERNS:
        for path in root.rglob(pattern):
            if path in seen:
                continue
            seen.add(path)
            results.append(path)
            if len(results) >= args.max_results:
                break
        if len(results) >= args.max_results:
            break

    print(f"Matched files/directories = {len(results)}")
    for path in results:
        kind = "dir" if path.is_dir() else "file"
        print(f"  [{kind}] {path}")


if __name__ == "__main__":
    main()
