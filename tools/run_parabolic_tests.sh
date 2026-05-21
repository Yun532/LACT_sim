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
      echo "usage: tools/run_parabolic_tests.sh [--no-corsika] [--corsika-file FILE]" >&2
      exit 2
      ;;
    *)
      CORSIKA_FILE="$1"
      shift
      ;;
  esac
done

cd "$ROOT_DIR"

CFG_DIR="run_logs/parabolic_tests/configs"
mkdir -p "$CFG_DIR"
mkdir -p run_logs/parabolic_tests/perfect_parallel
mkdir -p run_logs/parabolic_tests/point_900m
mkdir -p run_logs/parabolic_tests/raytrace_structure_parallel
mkdir -p run_logs/parabolic_tests/raytrace_structure_point_30m
mkdir -p run_logs/parabolic_tests/elevation_scan
mkdir -p run_logs/parabolic_tests/corsika

write_common_whiteboard_cfg() {
  local path="$1"
  local source_cfg="$2"
  local output_cfg="$3"
  local output_csv="$4"
  cat > "$path" <<CFG
telescope.config=../../../configs/official_tests/telescope_1229_minimal.cfg
mirror.config=../../../configs/mirrors/mirror_6m_parabolic.cfg
source.config=../../../configs/sources/${source_cfg}
output.config=../../../configs/outputs/${output_cfg}
propagation.speed_of_light_m_per_ns=0.299792458
output.csv=${output_csv}
CFG
}

write_common_whiteboard_cfg \
  "$CFG_DIR/perfect_parallel_whiteboard.cfg" \
  "parallel_1M_on_axis.cfg" \
  "whiteboard_f8.cfg" \
  "run_logs/parabolic_tests/perfect_parallel/hits.csv"

write_common_whiteboard_cfg \
  "$CFG_DIR/perfect_point_900m_whiteboard.cfg" \
  "point_900m_on_axis.cfg" \
  "whiteboard_point_900m_focus.cfg" \
  "run_logs/parabolic_tests/point_900m/hits.csv"

write_common_whiteboard_cfg \
  "$CFG_DIR/perfect_parallel_raytrace_structure_whiteboard.cfg" \
  "parallel_1M_on_axis.cfg" \
  "whiteboard_f8.cfg" \
  "run_logs/parabolic_tests/raytrace_structure_parallel/hits.csv"
cat >> "$CFG_DIR/perfect_parallel_raytrace_structure_whiteboard.cfg" <<CFG
obstruction.config=../../../configs/obstructions/raytrace_final_structure.cfg
CFG

write_common_whiteboard_cfg \
  "$CFG_DIR/point_30m_structure_whiteboard.cfg" \
  "point_30m_from_whiteboard_on_axis.cfg" \
  "whiteboard_f8.cfg" \
  "run_logs/parabolic_tests/raytrace_structure_point_30m/hits.csv"
cat >> "$CFG_DIR/point_30m_structure_whiteboard.cfg" <<CFG
obstruction.config=../../../configs/obstructions/raytrace_final_structure.cfg
CFG

cat > "$CFG_DIR/elevation_parallel_whiteboard.cfg" <<CFG
telescope.config=../../../configs/official_tests/telescope_1229_minimal.cfg
mirror.config=../../../configs/mirrors/mirror_6m_parabolic.cfg
source.config=../../../configs/sources/parallel_1M_on_axis.cfg
output.config=../../../configs/outputs/whiteboard_f8.cfg
propagation.speed_of_light_m_per_ns=0.299792458
output.csv=run_logs/parabolic_tests/elevation_scan/hits.csv
CFG

if [[ "$RUN_CORSIKA" -eq 0 ]]; then
  cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -DLACT_ENABLE_HESSIO=OFF
else
  "$ROOT_DIR/tools/build_hessio.sh"
  cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release \
    -DHESSIO_ROOT="$ROOT_DIR/external/hessioxxx/source"
fi
cmake --build "$BUILD_DIR" -j "${LACT_BUILD_JOBS:-4}"
ctest --test-dir "$BUILD_DIR" --output-on-failure

"$BUILD_DIR/run_optical_sim" "$CFG_DIR/perfect_parallel_whiteboard.cfg" \
  2>&1 | tee run_logs/parabolic_tests/perfect_parallel/run.log
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_spot_histogram.py \
  run_logs/parabolic_tests/perfect_parallel/hits.csv \
  --config "$CFG_DIR/perfect_parallel_whiteboard.cfg" \
  --output run_logs/parabolic_tests/perfect_parallel/spot.png \
  --max-bins 520 --dpi 350 \
  --title "Parallel beam with 6 m parabolic mirror"

"$BUILD_DIR/run_optical_sim" "$CFG_DIR/perfect_point_900m_whiteboard.cfg" \
  2>&1 | tee run_logs/parabolic_tests/point_900m/run.log
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_spot_histogram.py \
  run_logs/parabolic_tests/point_900m/hits.csv \
  --config "$CFG_DIR/perfect_point_900m_whiteboard.cfg" \
  --output run_logs/parabolic_tests/point_900m/spot.png \
  --max-bins 520 --dpi 350 \
  --title "900 m point source with 6 m parabolic mirror"

"$BUILD_DIR/run_optical_sim" "$CFG_DIR/perfect_parallel_raytrace_structure_whiteboard.cfg" \
  2>&1 | tee run_logs/parabolic_tests/raytrace_structure_parallel/run.log
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_spot_histogram.py \
  run_logs/parabolic_tests/raytrace_structure_parallel/hits.csv \
  --config "$CFG_DIR/perfect_parallel_raytrace_structure_whiteboard.cfg" \
  --output run_logs/parabolic_tests/raytrace_structure_parallel/spot.png \
  --max-bins 520 --dpi 350 \
  --title "Parallel beam with 6 m parabolic mirror and 3D obstruction"
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_mirror_hit_map.py \
  run_logs/parabolic_tests/raytrace_structure_parallel/hits.csv \
  --config "$CFG_DIR/perfect_parallel_raytrace_structure_whiteboard.cfg" \
  --require-surface \
  --overlay-facets \
  --sky-up \
  --output run_logs/parabolic_tests/raytrace_structure_parallel/mirror_hits_with_facet_outlines.png \
  --dpi 350 \
  --title "6 m parabolic mirror: hit points with aperture outline"
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_optical_layout_3d.py \
  --config "$CFG_DIR/perfect_parallel_raytrace_structure_whiteboard.cfg" \
  --show-obstruction \
  --output run_logs/parabolic_tests/raytrace_structure_parallel/layout_3d.png \
  --dpi 350
python3 python/plot_optical_layout_html.py \
  --config "$CFG_DIR/perfect_parallel_raytrace_structure_whiteboard.cfg" \
  --output run_logs/parabolic_tests/raytrace_structure_parallel/layout_3d.html

"$BUILD_DIR/run_optical_sim" "$CFG_DIR/point_30m_structure_whiteboard.cfg" \
  2>&1 | tee run_logs/parabolic_tests/raytrace_structure_point_30m/run.log
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_spot_histogram.py \
  run_logs/parabolic_tests/raytrace_structure_point_30m/hits.csv \
  --config "$CFG_DIR/point_30m_structure_whiteboard.cfg" \
  --output run_logs/parabolic_tests/raytrace_structure_point_30m/spot.png \
  --max-bins 520 --dpi 350 \
  --title "30 m point source with 6 m parabolic mirror and 3D obstruction"
MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/plot_mirror_hit_map.py \
  run_logs/parabolic_tests/raytrace_structure_point_30m/hits.csv \
  --config "$CFG_DIR/point_30m_structure_whiteboard.cfg" \
  --require-surface \
  --overlay-facets \
  --sky-up \
  --output run_logs/parabolic_tests/raytrace_structure_point_30m/mirror_hits_with_facet_outlines.png \
  --dpi 350 \
  --title "30 m point source: 6 m parabolic mirror hit points"

MPLBACKEND=Agg MPLCONFIGDIR="${MPLCONFIGDIR:-/tmp}" \
python3 python/run_elevation_parallel_scan.py \
  --config "$CFG_DIR/elevation_parallel_whiteboard.cfg" \
  --run-binary "$BUILD_DIR/run_optical_sim" \
  --elevations 0,10,20,30,40,50,60,70,80,90 \
  --n-bunches 100000 \
  --output-dir run_logs/parabolic_tests/elevation_scan

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
    echo "usage: tools/run_parabolic_tests.sh --corsika-file /path/to/input.zst" >&2
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

write_corsika_cfg() {
  local path="$1"
  local output_path="$2"
  local nsb_cfg="$3"
  local trigger_cfg="$4"
  cat > "$path" <<CFG
telescope.config=../../../configs/official_tests/telescope_1229_minimal.cfg
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70
mirror.config=../../../configs/mirrors/mirror_6m_parabolic.cfg
output.config=../../../configs/outputs/focal_plane_f8.cfg
camera.config=../../../configs/cameras/new_camera.cfg
sipm.config=../../../configs/sipm/ideal_sipm.cfg
atmosphere.config=../../../configs/atmosphere/ideal.cfg
nsb.config=../../../configs/nsb/${nsb_cfg}
trigger.config=../../../configs/trigger/${trigger_cfg}
source.mode=EventIO
source.eventio_path=
source.event_id_mode=event_array100
source.eventio_coordinate_frame=corsika_iact
source.max_shower_events=1
propagation.speed_of_light_m_per_ns=0.299792458
output.format=hdf5
output.hdf5_path=${output_path}
output.hdf5_storage=dense
output.hdf5_write_components=true
output.save_only_triggered=false
output.write_pixel_time_stats=true
waveform.enabled=true
waveform.source=pe
waveform.time_reference=image_mean
waveform.time_bin_width_ns=1
waveform.time_window_start_ns=-20
waveform.time_window_end_ns=120
CFG
}

write_corsika_cfg "$CFG_DIR/corsika_new_camera.cfg" \
  "run_logs/parabolic_tests/corsika/camera_dense.h5" \
  "ideal.cfg" "disabled.cfg"
write_corsika_cfg "$CFG_DIR/corsika_nsb_trigger_camera.cfg" \
  "run_logs/parabolic_tests/corsika/camera_nsb_trigger_dense.h5" \
  "example_constant_rate.cfg" "example_simple_multiplicity.cfg"
cat >> "$CFG_DIR/corsika_nsb_trigger_camera.cfg" <<CFG
output.save_only_triggered=true
trigger.pixel_threshold_pe=10
CFG

"$BUILD_DIR/run_corsika_trace" "$CFG_DIR/corsika_new_camera.cfg" "$CORSIKA_FILE" \
  2>&1 | tee run_logs/parabolic_tests/corsika/camera_run.log
"$BUILD_DIR/run_corsika_trace" "$CFG_DIR/corsika_nsb_trigger_camera.cfg" "$CORSIKA_FILE" \
  2>&1 | tee run_logs/parabolic_tests/corsika/camera_nsb_trigger_run.log

if [[ ! -s run_logs/parabolic_tests/corsika/camera_dense.h5 ]]; then
  echo "Expected dense camera HDF5 was not created: run_logs/parabolic_tests/corsika/camera_dense.h5" >&2
  exit 1
fi
if [[ ! -s run_logs/parabolic_tests/corsika/camera_nsb_trigger_dense.h5 ]]; then
  echo "Expected NSB+trigger HDF5 was not created: run_logs/parabolic_tests/corsika/camera_nsb_trigger_dense.h5" >&2
  exit 1
fi
