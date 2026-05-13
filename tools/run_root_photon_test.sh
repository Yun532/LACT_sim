#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${LACT_BUILD_DIR:-$ROOT_DIR/build}"

if [[ $# -lt 1 ]]; then
  echo "usage: tools/run_root_photon_test.sh /path/to/test.root" >&2
  exit 2
fi

INPUT_ROOT="$1"
OUT_DIR="$ROOT_DIR/run_logs/root_photon_test"
mkdir -p "$OUT_DIR"

cd "$ROOT_DIR"

cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -DLACT_ENABLE_HESSIO=OFF
cmake --build "$BUILD_DIR" -j "${LACT_BUILD_JOBS:-4}"

python3 python/root_photons_to_photon_csv.py \
  "$INPUT_ROOT" \
  "$OUT_DIR/photons_local.csv" \
  --telescope-az-deg 100.156 \
  --telescope-zenith-deg 25.802 \
  --metadata-output "$OUT_DIR/root_conversion.log" \
  2>&1 | tee "$OUT_DIR/root_conversion.stdout.log"

"$BUILD_DIR/run_optical_sim" configs/root_photon_test/root_photon_camera.cfg \
  2>&1 | tee "$OUT_DIR/run.log"

python3 python/export_trace_hdf5.py \
  --pixel-csv "$OUT_DIR/camera_pixels.csv" \
  --config configs/root_photon_test/root_photon_camera.cfg \
  --camera-csv configs/cameras/new_camera_pixels.csv \
  --storage dense \
  --output "$OUT_DIR/camera_dense.h5"

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_hdf5_camera.py \
  "$OUT_DIR/camera_dense.h5" \
  --event-id 0 \
  --telescope-id 0 \
  --quantity pe \
  --output "$OUT_DIR/camera_pe.png"

echo "ROOT photon camera test completed."
echo "Outputs: $OUT_DIR"
