#!/usr/bin/env python3
"""在指定工作目录从头执行Notebook；即使单元失败也保存已执行输出供定位。"""

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
try:
    NotebookClient(
        notebook,
        timeout=args.timeout,
        kernel_name="python3",
        resources={"metadata": {"path": str(args.cwd.resolve())}},
        on_cell_executed=lambda **event: print(
            f"已执行单元 {event['cell_index'] + 1}/{len(notebook.cells)}", flush=True),
    ).execute()
finally:
    nbformat.write(notebook, args.notebook)
print(f"Executed {len(notebook.cells)} cells: {args.notebook}")
