#!/usr/bin/env python3
"""Execute a notebook in place with an explicit working directory."""

import argparse
from pathlib import Path

import nbformat
from nbclient import NotebookClient


parser = argparse.ArgumentParser()
parser.add_argument("notebook", type=Path)
parser.add_argument("--cwd", type=Path, default=Path.cwd())
parser.add_argument("--timeout", type=int, default=300)
args = parser.parse_args()

notebook = nbformat.read(args.notebook, as_version=4)
NotebookClient(
    notebook,
    timeout=args.timeout,
    kernel_name="python3",
    resources={"metadata": {"path": str(args.cwd.resolve())}},
).execute()
nbformat.write(notebook, args.notebook)
print(f"Executed {len(notebook.cells)} cells: {args.notebook}")
