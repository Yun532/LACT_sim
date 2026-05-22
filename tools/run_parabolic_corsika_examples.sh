#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${LACT_BUILD_DIR:-${BUILD_DIR:-$ROOT_DIR/build}}"
CORSIKA_FILE="${LACT_DEFAULT_CORSIKA_FILE:-}"
ONLY="both"

usage() {
  cat >&2 <<USAGE
usage: tools/run_parabolic_corsika_examples.sh --corsika-file FILE [--only ideal|full|both]

Runs two standalone 6 m parabolic CORSIKA camera examples:
  ideal : no PDE/efficiency weighting, plots photon_count and 1 ns/bin GIF
  full  : full response with obstruction, plots p.e. and 1 ns/bin GIF
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --corsika-file)
      CORSIKA_FILE="${2:-}"
      if [[ -z "$CORSIKA_FILE" ]]; then
        echo "--corsika-file requires a path." >&2
        exit 2
      fi
      shift 2
      ;;
    --only)
      ONLY="${2:-}"
      if [[ "$ONLY" != "ideal" && "$ONLY" != "full" && "$ONLY" != "both" ]]; then
        echo "--only must be ideal, full, or both." >&2
        exit 2
      fi
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage
      exit 2
      ;;
    *)
      CORSIKA_FILE="$1"
      shift
      ;;
  esac
done

if [[ -z "$CORSIKA_FILE" ]]; then
  echo "A CORSIKA/EventIO .zst file is required." >&2
  usage
  exit 2
fi
if [[ ! -s "$CORSIKA_FILE" ]]; then
  echo "CORSIKA/EventIO file does not exist or is empty: $CORSIKA_FILE" >&2
  exit 2
fi

cd "$ROOT_DIR"

if [[ ! -x "$BUILD_DIR/run_corsika_trace" ]]; then
  echo "run_corsika_trace not found in $BUILD_DIR; running make BUILD_DIR=$BUILD_DIR build"
  make BUILD_DIR="$BUILD_DIR" build
fi

if [[ "$(uname -s)" == "Darwin" ]]; then
  export DYLD_LIBRARY_PATH="$ROOT_DIR/external/hessioxxx/source/lib:${DYLD_LIBRARY_PATH:-}"
else
  export LD_LIBRARY_PATH="$ROOT_DIR/external/hessioxxx/source/lib:${LD_LIBRARY_PATH:-}"
fi

mkdir -p run_logs/parabolic_corsika_examples/ideal/plots
mkdir -p run_logs/parabolic_corsika_examples/full_response_obstruction/plots

select_brightest_telescope() {
  local h5="$1"
  local event_id="$2"
  python3 -c 'import h5py, sys
path, event_id = sys.argv[1], int(sys.argv[2])
with h5py.File(path, "r") as h5:
    rows = h5["images/index"][:]
    rows = rows[rows["event_id"] == event_id]
    if len(rows) == 0:
        raise SystemExit("no image rows for selected event")
    row = max(rows, key=lambda r: float(r["total_pe"]) if "total_pe" in rows.dtype.names else float(r["total_signal"]))
    print(int(row["telescope_id"]))' "$h5" "$event_id"
}

run_ideal() {
  local cfg="configs/parabolic_corsika_examples/ideal_parabolic_camera.cfg"
  local h5="run_logs/parabolic_corsika_examples/ideal/camera_ideal_dense.h5"
  local plot_dir="run_logs/parabolic_corsika_examples/ideal/plots"

  "$BUILD_DIR/run_corsika_trace" "$cfg" "$CORSIKA_FILE" \
    2>&1 | tee run_logs/parabolic_corsika_examples/ideal/run.log

  if [[ ! -s "$h5" ]]; then
    echo "Expected HDF5 was not created: $h5" >&2
    exit 1
  fi

  MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
  python3 python/select_hdf5_event.py "$h5" \
    --quantity photon_count \
    --min-images 1 \
    --output-env "$plot_dir/selected_event.env" \
    --output-summary "$plot_dir/selected_event.txt" \
    2>&1 | tee "$plot_dir/select_event.log"

  # shellcheck disable=SC1091
  source "$plot_dir/selected_event.env"
  local gif_telescope_id
  gif_telescope_id="$(select_brightest_telescope "$h5" "$LACT_SELECTED_EVENT_ID")"
  echo "Ideal GIF telescope_id=${gif_telescope_id} (1-based label $((gif_telescope_id + 1)))"

  MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
  python3 python/plot_hdf5_camera.py "$h5" \
    --event-id "$LACT_SELECTED_EVENT_ID" \
    --quantity photon_count \
    --output "$plot_dir/all_tel_photon_count"

  MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
  python3 python/plot_hdf5_waveform_gif.py "$h5" \
    --event-id "$LACT_SELECTED_EVENT_ID" \
    --telescope-id "$gif_telescope_id" \
    --quantity photon_count \
    --output-dir "$plot_dir/waveform_photon_count_frames" \
    --gif "$plot_dir/brightest_tel_photon_count.gif" \
    --stride 1

  MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
  python3 python/plot_hdf5_time_histogram.py "$h5" \
    --event-id "$LACT_SELECTED_EVENT_ID" \
    --quantity photon_count \
    --output "$plot_dir/time_hist_photon_count.png"
}

run_full() {
  local cfg="configs/parabolic_corsika_examples/full_response_obstruction_parabolic_camera.cfg"
  local h5="run_logs/parabolic_corsika_examples/full_response_obstruction/camera_full_response_obstruction_dense.h5"
  local plot_dir="run_logs/parabolic_corsika_examples/full_response_obstruction/plots"

  "$BUILD_DIR/run_corsika_trace" "$cfg" "$CORSIKA_FILE" \
    2>&1 | tee run_logs/parabolic_corsika_examples/full_response_obstruction/run.log

  if [[ ! -s "$h5" ]]; then
    echo "Expected HDF5 was not created: $h5" >&2
    exit 1
  fi

  MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
  python3 python/select_hdf5_event.py "$h5" \
    --quantity pe \
    --min-images 1 \
    --output-env "$plot_dir/selected_event.env" \
    --output-summary "$plot_dir/selected_event.txt" \
    2>&1 | tee "$plot_dir/select_event.log"

  # shellcheck disable=SC1091
  source "$plot_dir/selected_event.env"
  local gif_telescope_id
  gif_telescope_id="$(select_brightest_telescope "$h5" "$LACT_SELECTED_EVENT_ID")"
  echo "Full-response GIF telescope_id=${gif_telescope_id} (1-based label $((gif_telescope_id + 1)))"

  MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
  python3 python/plot_hdf5_camera.py "$h5" \
    --event-id "$LACT_SELECTED_EVENT_ID" \
    --quantity pe \
    --output "$plot_dir/all_tel_pe"

  MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
  python3 python/plot_hdf5_waveform_gif.py "$h5" \
    --event-id "$LACT_SELECTED_EVENT_ID" \
    --telescope-id "$gif_telescope_id" \
    --quantity pe \
    --output-dir "$plot_dir/waveform_pe_frames" \
    --gif "$plot_dir/brightest_tel_pe.gif" \
    --stride 1

  MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
  python3 python/plot_hdf5_time_histogram.py "$h5" \
    --event-id "$LACT_SELECTED_EVENT_ID" \
    --quantity pe \
    --output "$plot_dir/time_hist_pe.png"
}

case "$ONLY" in
  ideal)
    run_ideal
    ;;
  full)
    run_full
    ;;
  both)
    run_ideal
    run_full
    ;;
esac

echo "Parabolic CORSIKA examples completed. See run_logs/parabolic_corsika_examples/."
