#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HESSIO_DIR="$ROOT_DIR/external/hessioxxx/source"
ARCHIVE="$ROOT_DIR/external/hessioxxx/hessioxxx.tar.gz"

# Git 只保存上游压缩包；首次构建时按需展开，避免仓库出现近百个供应商文件。
if [[ ! -f "$HESSIO_DIR/Makefile" ]]; then
    mkdir -p "$HESSIO_DIR"
    tar -xzf "$ARCHIVE" -C "$HESSIO_DIR" --strip-components=1
fi

mkdir -p "$HESSIO_DIR/out" "$HESSIO_DIR/bin" "$HESSIO_DIR/lib"
make -C "$HESSIO_DIR" "$@"
