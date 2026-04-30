#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${LACT_BUILD_DIR:-$ROOT_DIR/build}"
RUN_CORSIKA=1
CORSIKA_FILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --no-corsika)
      RUN_CORSIKA=0
      shift
      ;;
    --corsika-file)
      CORSIKA_FILE="${2:-}"
      if [[ -z "$CORSIKA_FILE" ]]; then
        echo "--corsika-file requires a path." >&2
        exit 2
      fi
      shift 2
      ;;
    -*)
      echo "Unknown option: $1" >&2
      echo "usage: tools/run_official_tests.sh [--no-corsika] [--corsika-file FILE]" >&2
      exit 2
      ;;
    *)
      CORSIKA_FILE="$1"
      shift
      ;;
  esac
done

cd "$ROOT_DIR"
mkdir -p run_logs/official_tests/perfect_parallel
mkdir -p run_logs/official_tests/point_900m
mkdir -p run_logs/official_tests/deformation_scan
mkdir -p run_logs/official_tests/corsika/plots/shower1_array0_whiteboard
mkdir -p run_logs/official_tests/corsika/plots/shower1_array0_camera
mkdir -p run_logs/official_tests/corsika/plots/shower1_array0_nsb_trigger
mkdir -p run_logs/official_tests/collector_angular_response

if [[ "$RUN_CORSIKA" -eq 0 ]]; then
  cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -DLACT_ENABLE_HESSIO=OFF
else
  "$ROOT_DIR/tools/build_hessio.sh"
  cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release \
    -DHESSIO_ROOT="$ROOT_DIR/external/hessioxxx/source"
fi
cmake --build "$BUILD_DIR" -j "${LACT_BUILD_JOBS:-4}"
ctest --test-dir "$BUILD_DIR" --output-on-failure

"$BUILD_DIR/run_optical_sim" configs/official_tests/perfect_parallel_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/perfect_parallel/run.log
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_spot_histogram.py \
  run_logs/official_tests/perfect_parallel/hits.csv \
  --output run_logs/official_tests/perfect_parallel/spot.png \
  --max-bins 520 --dpi 350

"$BUILD_DIR/run_optical_sim" configs/official_tests/perfect_point_900m_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/point_900m/run.log
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_spot_histogram.py \
  run_logs/official_tests/point_900m/hits.csv \
  --output run_logs/official_tests/point_900m/spot.png \
  --max-bins 520 --dpi 350

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/run_elevation_parallel_scan.py \
  --config configs/official_tests/deformation_parallel_whiteboard.cfg \
  --run-binary "$BUILD_DIR/run_optical_sim" \
  --elevations 0,10,20,30,40,50,60,70,80,90 \
  --n-bunches 100000 \
  --output-dir run_logs/official_tests/deformation_scan

"$BUILD_DIR/scan_light_collector_angular_response" \
  --photons-per-angle 2000 \
  --angle-step-deg 1 \
  --max-angle-deg 90 \
  --output run_logs/official_tests/collector_angular_response/collector_angular_response.csv \
  2>&1 | tee run_logs/official_tests/collector_angular_response/run.log
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_collector_angular_response.py \
  run_logs/official_tests/collector_angular_response/collector_angular_response.csv \
  --output run_logs/official_tests/collector_angular_response/collector_angular_response.png \
  --dpi 350

if [[ "$RUN_CORSIKA" -eq 0 ]]; then
  echo "Skipping CORSIKA/EventIO tests (--no-corsika)."
  exit 0
fi

if [[ -z "$CORSIKA_FILE" ]]; then
  echo "A CORSIKA/EventIO file is required for CORSIKA tests." >&2
  echo "usage: tools/run_official_tests.sh --corsika-file /path/to/input.zst" >&2
  exit 2
fi
if [[ ! -x "$BUILD_DIR/run_corsika_trace" ]]; then
  echo "run_corsika_trace was not built. Build external/hessioxxx/source first." >&2
  exit 2
fi

if [[ "$(uname -s)" == "Darwin" ]]; then
  export DYLD_LIBRARY_PATH="$ROOT_DIR/external/hessioxxx/source/lib:${DYLD_LIBRARY_PATH:-}"
else
  export LD_LIBRARY_PATH="$ROOT_DIR/external/hessioxxx/source/lib:${LD_LIBRARY_PATH:-}"
fi

"$BUILD_DIR/run_corsika_trace" configs/official_tests/corsika_whiteboard.cfg "$CORSIKA_FILE" \
  2>&1 | tee run_logs/official_tests/corsika/whiteboard_run.log
"$BUILD_DIR/run_corsika_trace" configs/official_tests/corsika_new_camera.cfg "$CORSIKA_FILE" \
  2>&1 | tee run_logs/official_tests/corsika/camera_run.log
"$BUILD_DIR/run_corsika_trace" configs/official_tests/corsika_nsb_trigger_camera.cfg "$CORSIKA_FILE" \
  2>&1 | tee run_logs/official_tests/corsika/camera_nsb_trigger_run.log

if [[ ! -s run_logs/official_tests/corsika/whiteboard_hits.csv ]]; then
  echo "Expected whiteboard CSV was not created: run_logs/official_tests/corsika/whiteboard_hits.csv" >&2
  echo "Check run_logs/official_tests/corsika/whiteboard_run.log" >&2
  exit 1
fi
if [[ ! -s run_logs/official_tests/corsika/camera_dense.h5 ]]; then
  echo "Expected dense camera HDF5 was not created: run_logs/official_tests/corsika/camera_dense.h5" >&2
  echo "Check run_logs/official_tests/corsika/camera_run.log" >&2
  exit 1
fi
if [[ ! -s run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5 ]]; then
  echo "Expected NSB+trigger HDF5 was not created: run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5" >&2
  echo "Check run_logs/official_tests/corsika/camera_nsb_trigger_run.log" >&2
  exit 1
fi

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_corsika_trace_output.py \
  run_logs/official_tests/corsika/whiteboard_hits.csv \
  --summary-csv run_logs/official_tests/corsika/whiteboard_summary.csv \
  --shower-event-number 1 \
  --array-id 0 \
  --output-dir run_logs/official_tests/corsika/plots/shower1_array0_whiteboard

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --shower-event-number 1 \
  --array-id 0 \
  --telescope-id "${LACT_PLOT_TELESCOPE_ID:-7}" \
  --quantity pe \
  --output run_logs/official_tests/corsika/plots/shower1_array0_camera/telescope_${LACT_PLOT_TELESCOPE_ID:-7}_pe.png

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5 \
  --shower-event-number 1 \
  --array-id 0 \
  --telescope-id "${LACT_PLOT_TELESCOPE_ID:-7}" \
  --quantity pe \
  --output run_logs/official_tests/corsika/plots/shower1_array0_nsb_trigger/telescope_${LACT_PLOT_TELESCOPE_ID:-7}_pe.png

echo "Official tests completed. See run_logs/official_tests/ for outputs and README.md for commands."
