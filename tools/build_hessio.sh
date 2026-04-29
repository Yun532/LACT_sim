#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HESSIO_DIR="$ROOT_DIR/external/hessioxxx/source"

mkdir -p "$HESSIO_DIR/out" "$HESSIO_DIR/bin" "$HESSIO_DIR/lib"
make -C "$HESSIO_DIR" "$@"
