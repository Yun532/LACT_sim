#!/usr/bin/env sh
set -eu

HESSIO_ROOT="${HESSIO_ROOT:-external/hessioxxx/source}"
OUT="${OUT:-build_eventio/eventio_to_photon_csv}"

mkdir -p "$(dirname "$OUT")"

cc tools/eventio_to_photon_csv.c \
  -I"$HESSIO_ROOT/include" \
  -L"$HESSIO_ROOT/lib" \
  -lhessio \
  -lm \
  -o "$OUT"

printf 'built %s\n' "$OUT"
printf 'run with: DYLD_LIBRARY_PATH=%s/lib %s input.zst output.csv\n' "$HESSIO_ROOT" "$OUT"
