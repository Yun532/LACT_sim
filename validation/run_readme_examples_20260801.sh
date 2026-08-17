#!/usr/bin/env bash

set -u -o pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 /path/to/corsika.zst" >&2
    exit 2
fi

corsika_file=$1
build_dir=${BUILD_DIR:-build_official_electronics_enabled_20260801}
run_root=run_logs/readme_examples_20260801
figure_dir=$run_root/figures
status_file=$run_root/status.tsv
mkdir -p "$run_root/logs" "$figure_dir"
printf 'case\texit_code\n' > "$status_file"

run_case() {
    local name=$1
    shift
    echo "[RUN] $name"
    "$@" > "$run_root/logs/$name.log" 2>&1
    local code=$?
    printf '%s\t%s\n' "$name" "$code" >> "$status_file"
    if [[ $code -eq 0 ]]; then
        echo "[PASS] $name"
    else
        echo "[FAIL] $name (exit $code; see $run_root/logs/$name.log)"
    fi
}

# README: synthetic optical-source entries.
run_case perfect_parallel \
    "$build_dir/run_optical_sim" \
    configs/official_tests/perfect_parallel_whiteboard.cfg
run_case perfect_parallel_plot \
    python3 python/plot_spot_histogram.py \
    run_logs/official_tests/perfect_parallel/hits.csv \
    --config configs/official_tests/perfect_parallel_whiteboard.cfg \
    --output "$figure_dir/perfect_parallel_spot.png"

run_case point_900m \
    "$build_dir/run_optical_sim" \
    configs/official_tests/perfect_point_900m_whiteboard.cfg
run_case point_900m_plot \
    python3 python/plot_spot_histogram.py \
    run_logs/official_tests/point_900m/hits.csv \
    --config configs/official_tests/perfect_point_900m_whiteboard.cfg \
    --output "$figure_dir/point_900m_spot.png"

run_case obstruction_parallel \
    "$build_dir/run_optical_sim" \
    configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg
run_case obstruction_parallel_plot \
    python3 python/plot_spot_histogram.py \
    run_logs/official_tests/raytrace_structure_parallel/hits.csv \
    --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
    --output "$figure_dir/obstruction_parallel_spot.png"
run_case obstruction_layout_plot \
    python3 python/plot_optical_layout_3d.py \
    --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
    --show-obstruction \
    --output "$figure_dir/obstruction_layout_3d.png"

run_case deformation_parallel \
    "$build_dir/run_optical_sim" \
    configs/official_tests/deformation_parallel_whiteboard.cfg
run_case deformation_parallel_plot \
    python3 python/plot_spot_histogram.py \
    run_logs/official_tests/deformation_scan/hits.csv \
    --config configs/official_tests/deformation_parallel_whiteboard.cfg \
    --output "$figure_dir/deformation_parallel_spot.png"

run_case lact2_measured_parallel \
    "$build_dir/run_optical_sim" configs/examples/lact2_measured_parallel.cfg
run_case lact2_measured_parallel_plot \
    python3 python/plot_spot_histogram.py \
    run_logs/lact2_measured_20260622_parallel/hits.csv \
    --config configs/examples/lact2_measured_parallel.cfg \
    --output "$figure_dir/lact2_measured_parallel_spot.png"

run_case photon_csv_local \
    "$build_dir/run_optical_sim" configs/examples/photon_csv_local_whiteboard.cfg
run_case photon_csv_local_plot \
    python3 python/plot_spot_histogram.py \
    run_logs/examples/photon_csv_local_whiteboard/hits.csv \
    --config configs/examples/photon_csv_local_whiteboard.cfg \
    --output "$figure_dir/photon_csv_local_spot.png"

# README: minimal event-1909 Photon CSV checks.
run_case photon_csv_minimal_optics \
    "$build_dir/run_optical_sim" configs/examples/photon_csv_minimal_optics.cfg
run_case photon_csv_minimal_optics_plot \
    python3 python/plot_minimal_photon_csv_outputs.py \
    --mode optics \
    --hits run_logs/examples/photon_csv_minimal/whiteboard_hits.csv \
    --photon-pixels run_logs/examples/photon_csv_minimal/camera_photon_counts.csv \
    --camera configs/cameras/new_camera_pixels.csv \
    --output-dir "$figure_dir/photon_csv_minimal"

run_case photon_csv_full_camera_root \
    "$build_dir/run_corsika_trace" configs/examples/photon_csv_full_camera_root.cfg
run_case photon_csv_full_camera_root_plot \
    python3 python/plot_photon_csv_root_pylast.py \
    run_logs/examples/photon_csv_full_camera/lact_events.root \
    --event-id 1909 \
    --output "$figure_dir/photon_csv_event1909_camera_pe.png"

# README: EventIO/CORSIKA entries.
run_case corsika_whiteboard \
    "$build_dir/run_corsika_trace" \
    configs/official_tests/corsika_whiteboard.cfg "$corsika_file"
run_case corsika_whiteboard_plot \
    python3 python/plot_corsika_trace_output.py \
    run_logs/official_tests/corsika/whiteboard_hits.csv \
    --summary-csv run_logs/official_tests/corsika/whiteboard_summary.csv \
    --event-id 107 \
    --output-dir "$figure_dir/corsika_whiteboard"

run_case corsika_ideal_camera \
    "$build_dir/run_corsika_trace" \
    configs/official_tests/corsika_new_camera.cfg "$corsika_file"
run_case corsika_ideal_camera_plot \
    python3 python/plot_hdf5_camera.py \
    run_logs/official_tests/corsika/camera_dense.h5 \
    --image-index 0 --quantity pe \
    --output "$figure_dir/corsika_ideal_camera_first.png"

run_case corsika_nsb_trigger_camera \
    "$build_dir/run_corsika_trace" \
    configs/official_tests/corsika_nsb_trigger_camera.cfg "$corsika_file"
run_case corsika_nsb_trigger_camera_plot \
    python3 python/plot_hdf5_camera.py \
    run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5 \
    --image-index 0 --quantity pe \
    --output "$figure_dir/corsika_nsb_trigger_camera_first.png"

run_case corsika_full_response_camera \
    "$build_dir/run_corsika_trace" \
    configs/official_tests/corsika_full_response_camera.cfg "$corsika_file"
run_case corsika_full_response_camera_plot \
    python3 python/plot_hdf5_camera.py \
    run_logs/official_tests/corsika/camera_full_response_dense.h5 \
    --image-index 0 --quantity pe \
    --output "$figure_dir/corsika_full_response_camera_first.png"

run_case lactroot_quicklook \
    "$build_dir/run_corsika_trace" \
    configs/examples/corsika_lact_root_quicklook.cfg "$corsika_file"
run_case lactroot_quicklook_plot \
    python3 python/plot_lact_root_output.py \
    run_logs/lact_root_quicklook/lact_events.root \
    --outdir "$figure_dir/lactroot_quicklook"

if awk -F '\t' 'NR > 1 && $2 != 0 {failed=1} END {exit failed}' "$status_file"; then
    echo "All README examples and plot commands passed."
    exit 0
fi

echo "One or more README examples failed; inspect $status_file." >&2
exit 1
