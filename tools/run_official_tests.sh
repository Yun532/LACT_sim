#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${LACT_BUILD_DIR:-$ROOT_DIR/build}"
RUN_CORSIKA=1
CORSIKA_FILE=""
DEFAULT_CORSIKA_FILE="${LACT_DEFAULT_CORSIKA_FILE:-}"

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
mkdir -p run_logs/official_tests/raytrace_structure_parallel
mkdir -p run_logs/official_tests/raytrace_structure_point_30m
mkdir -p run_logs/official_tests/deformation_scan
mkdir -p run_logs/official_tests/corsika/plots
mkdir -p run_logs/official_tests/collector_angular_response
mkdir -p run_logs/official_tests/efficiency_curves
mkdir -p run_logs/official_tests/nsb_spectral

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
  --config configs/official_tests/perfect_parallel_whiteboard.cfg \
  --output run_logs/official_tests/perfect_parallel/spot.png \
  --max-bins 520 --dpi 350

"$BUILD_DIR/run_optical_sim" configs/official_tests/perfect_point_900m_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/point_900m/run.log
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_spot_histogram.py \
  run_logs/official_tests/point_900m/hits.csv \
  --config configs/official_tests/perfect_point_900m_whiteboard.cfg \
  --output run_logs/official_tests/point_900m/spot.png \
  --max-bins 520 --dpi 350

"$BUILD_DIR/run_optical_sim" configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/raytrace_structure_parallel/run.log
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_spot_histogram.py \
  run_logs/official_tests/raytrace_structure_parallel/hits.csv \
  --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  --output run_logs/official_tests/raytrace_structure_parallel/spot.png \
  --max-bins 520 --dpi 350 \
  --title "Parallel beam with 3D obstruction"
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_mirror_hit_map.py \
  run_logs/official_tests/raytrace_structure_parallel/hits.csv \
  --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  --require-surface \
  --overlay-facets \
  --sky-up \
  --output run_logs/official_tests/raytrace_structure_parallel/mirror_hits_with_facet_outlines.png \
  --dpi 350 \
  --title "Parallel beam: mirror hit points with facet outlines"
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_optical_layout_3d.py \
  --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  --show-obstruction \
  --output run_logs/official_tests/raytrace_structure_parallel/layout_3d.png \
  --dpi 350
python3 python/plot_optical_layout_html.py \
  --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  --output run_logs/official_tests/raytrace_structure_parallel/layout_3d.html

"$BUILD_DIR/run_optical_sim" configs/official_tests/point_30m_structure_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/raytrace_structure_point_30m/run.log
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_spot_histogram.py \
  run_logs/official_tests/raytrace_structure_point_30m/hits.csv \
  --config configs/official_tests/point_30m_structure_whiteboard.cfg \
  --output run_logs/official_tests/raytrace_structure_point_30m/spot.png \
  --max-bins 520 --dpi 350 \
  --title "30 m point source with 3D obstruction"
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_mirror_hit_map.py \
  run_logs/official_tests/raytrace_structure_point_30m/hits.csv \
  --config configs/official_tests/point_30m_structure_whiteboard.cfg \
  --require-surface \
  --overlay-facets \
  --sky-up \
  --output run_logs/official_tests/raytrace_structure_point_30m/mirror_hits_with_facet_outlines.png \
  --dpi 350 \
  --title "30 m point source: mirror hit points with facet outlines"

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

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_efficiency_curves.py \
  --output-dir run_logs/official_tests/efficiency_curves \
  --dpi 350 \
  2>&1 | tee run_logs/official_tests/efficiency_curves/run.log

"$BUILD_DIR/compute_nsb_rate" configs/nsb/spectral_rate_check_with_obstruction.cfg \
  2>&1 | tee run_logs/official_tests/nsb_spectral/run.log
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_nsb_spectral_rate.py \
  --effective-area-m2 22.606448 \
  --output run_logs/official_tests/nsb_spectral/nsb_spectral_response.png \
  --diagnostic-csv run_logs/official_tests/nsb_spectral/diagnostic.csv \
  --summary run_logs/official_tests/nsb_spectral/summary.txt \
  2>&1 | tee run_logs/official_tests/nsb_spectral/plot.log

if [[ "$RUN_CORSIKA" -eq 0 ]]; then
  echo "Skipping CORSIKA/EventIO tests (--no-corsika)."
  exit 0
fi

if [[ -z "$CORSIKA_FILE" ]]; then
  if [[ -s "$DEFAULT_CORSIKA_FILE" ]]; then
    CORSIKA_FILE="$DEFAULT_CORSIKA_FILE"
    echo "Using default CORSIKA/EventIO test file: $CORSIKA_FILE"
  else
    echo "A CORSIKA/EventIO file is required for CORSIKA tests." >&2
    echo "usage: tools/run_official_tests.sh --corsika-file /path/to/input.zst" >&2
    echo "or set LACT_DEFAULT_CORSIKA_FILE=/path/to/input.zst" >&2
    exit 2
  fi
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
"$BUILD_DIR/run_corsika_trace" configs/official_tests/corsika_obstruction_nsb_trigger_camera.cfg "$CORSIKA_FILE" \
  2>&1 | tee run_logs/official_tests/corsika/camera_obstruction_nsb_trigger_run.log
"$BUILD_DIR/run_corsika_trace" configs/official_tests/corsika_full_response_camera.cfg "$CORSIKA_FILE" \
  2>&1 | tee run_logs/official_tests/corsika/camera_full_response_run.log

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
if [[ ! -s run_logs/official_tests/corsika/camera_obstruction_nsb_trigger_dense.h5 ]]; then
  echo "Expected obstruction+NSB+trigger HDF5 was not created: run_logs/official_tests/corsika/camera_obstruction_nsb_trigger_dense.h5" >&2
  echo "Check run_logs/official_tests/corsika/camera_obstruction_nsb_trigger_run.log" >&2
  exit 1
fi
if [[ ! -s run_logs/official_tests/corsika/camera_full_response_dense.h5 ]]; then
  echo "Expected full-response HDF5 was not created: run_logs/official_tests/corsika/camera_full_response_dense.h5" >&2
  echo "Check run_logs/official_tests/corsika/camera_full_response_run.log" >&2
  exit 1
fi

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/select_hdf5_event.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5 \
  run_logs/official_tests/corsika/camera_obstruction_nsb_trigger_dense.h5 \
  run_logs/official_tests/corsika/camera_full_response_dense.h5 \
  --quantity pe \
  --min-images 1 \
  --output-env run_logs/official_tests/corsika/plots/selected_event.env \
  --output-summary run_logs/official_tests/corsika/plots/selected_event.txt \
  2>&1 | tee run_logs/official_tests/corsika/plots/select_event.log

# shellcheck disable=SC1091
source run_logs/official_tests/corsika/plots/selected_event.env
PLOT_EVENT_DIR="run_logs/official_tests/corsika/plots/event_${LACT_SELECTED_EVENT_ID}"
LACT_SELECTED_NSB_GIF_TELESCOPE_ID="$(
python3 -c 'import h5py, sys
path, event_id = sys.argv[1], int(sys.argv[2])
with h5py.File(path, "r") as h5:
    rows = h5["images/index"][:]
    rows = rows[rows["event_id"] == event_id]
    if len(rows) == 0:
        raise SystemExit("no image rows for selected event")
    row = max(rows, key=lambda r: float(r["total_pe"]))
    print(int(row["telescope_id"]))' \
  run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5 \
  "$LACT_SELECTED_EVENT_ID"
)"
mkdir -p "$PLOT_EVENT_DIR/whiteboard" \
  "$PLOT_EVENT_DIR/camera" \
  "$PLOT_EVENT_DIR/nsb_trigger" \
  "$PLOT_EVENT_DIR/obstruction_nsb_trigger" \
  "$PLOT_EVENT_DIR/full_response" \
  "$PLOT_EVENT_DIR/layout"

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_corsika_trace_output.py \
  run_logs/official_tests/corsika/whiteboard_hits.csv \
  --summary-csv run_logs/official_tests/corsika/whiteboard_summary.csv \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --config configs/official_tests/corsika_whiteboard.cfg \
  --output-dir "$PLOT_EVENT_DIR/whiteboard"

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --quantity pe \
  --output "$PLOT_EVENT_DIR/camera/all_tel_pe"

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --quantity cherenkov_pe \
  --output "$PLOT_EVENT_DIR/nsb_trigger/all_tel_cherenkov_pe"
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_hdf5_waveform_gif.py \
  run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --quantity cherenkov_pe \
  --output-dir "$PLOT_EVENT_DIR/nsb_trigger/waveform_cherenkov_pe_frames" \
  --gif "$PLOT_EVENT_DIR/nsb_trigger/all_tel_cherenkov_pe_gif" \
  --stride 1

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --quantity nsb_pe \
  --output "$PLOT_EVENT_DIR/nsb_trigger/all_tel_nsb_pe"
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_hdf5_waveform_gif.py \
  run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --telescope-id "$LACT_SELECTED_NSB_GIF_TELESCOPE_ID" \
  --quantity nsb_pe \
  --output-dir "$PLOT_EVENT_DIR/nsb_trigger/waveform_nsb_pe_frames" \
  --gif "$PLOT_EVENT_DIR/nsb_trigger/one_tel_nsb_pe.gif" \
  --stride 1

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --quantity pe \
  --output "$PLOT_EVENT_DIR/nsb_trigger/all_tel_final_pe"
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_hdf5_waveform_gif.py \
  run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --quantity pe \
  --output-dir "$PLOT_EVENT_DIR/nsb_trigger/waveform_final_pe_frames" \
  --gif "$PLOT_EVENT_DIR/nsb_trigger/all_tel_final_pe_gif" \
  --stride 1

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_obstruction_nsb_trigger_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --quantity cherenkov_pe \
  --output "$PLOT_EVENT_DIR/obstruction_nsb_trigger/all_tel_cherenkov_pe"

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_obstruction_nsb_trigger_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --quantity nsb_pe \
  --output "$PLOT_EVENT_DIR/obstruction_nsb_trigger/all_tel_nsb_pe"

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_obstruction_nsb_trigger_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --quantity pe \
  --output "$PLOT_EVENT_DIR/obstruction_nsb_trigger/all_tel_final_pe"

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_full_response_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --quantity pe \
  --output "$PLOT_EVENT_DIR/full_response/all_tel_pe"

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_hdf5_array_layout.py \
  run_logs/official_tests/corsika/camera_full_response_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --quantity pe \
  --log-color \
  --output "$PLOT_EVENT_DIR/layout/core_and_array_pe.png" \
  --dpi 350

echo "Official tests completed. See run_logs/official_tests/ for outputs and README.md for commands."
