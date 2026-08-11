# LACT_sim User Guide

[中文手册](user_guide_zh.md) | [Back to English README](../README_EN.md)

This guide covers the complete user workflow from installation and compilation to standard runs and plots. Detailed physics conventions and file schemas are linked at the end.

## Install and build

```bash
git clone https://github.com/Yun532/LACT_sim.git
cd LACT_sim

# Ubuntu/Debian
sudo apt-get update
sudo apt-get install cmake g++ make zlib1g-dev libhdf5-dev python3 python3-pip

python3 -m pip install numpy pandas matplotlib h5py

make
make test
```

The normal build first compiles the vendored `external/hessioxxx` and then builds LACT_sim under `build/`. The principal executables are:

```text
build/run_optical_sim
build/run_corsika_trace
build/compute_nsb_rate
```

For synthetic sources without EventIO support:

```bash
make no-hessio
ctest --test-dir build --output-on-failure
```

To keep EventIO/HDF5 but disable ROOT:

```bash
make no-root
```

For an installed ROOT 6.24 or newer:

```bash
source /path/to/root/bin/thisroot.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DHESSIO_ROOT="$PWD/external/hessioxxx/source" \
  -DLACT_ENABLE_ROOT=ON \
  -DROOT_DIR="$(root-config --cmakedir)"
cmake --build build -j4
```

On a headless server:

```bash
export MPLBACKEND=Agg
export MPLCONFIGDIR=/tmp
```

## Override cfg values on the command line

`run_optical_sim`, `run_corsika_trace`, `compute_nsb_rate`, and
`run_camera_electronics` accept repeatable sim_telarray-style `-C key=value`
overrides:

```bash
build/run_corsika_trace configs/examples/photon_csv_full_camera_root.cfg \
  -C nsb.enabled=true \
  -C nsb.seed=20260811 \
  -C electronics.single_pe.enabled=false \
  -C waveform.enabled=false \
  -C output.lact_root_path=run_logs/batch/run_20260811.root
```

Overrides are applied after all component cfg files are expanded, so they have
the highest precedence and the last repeated value wins. Compact `-Ckey=value`,
`--set key=value`, and `--set=key=value` forms are also accepted.

## Official-test entry point

Without CORSIKA:

```bash
tools/run_official_tests.sh --no-corsika \
  2>&1 | tee run_logs/official_tests/run_no_corsika.log
```

With CORSIKA/EventIO:

```bash
tools/run_official_tests.sh --corsika-file /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/run_with_corsika.log
```

The script configures and builds the code, runs CTest, executes the standard cases, and writes plots and data under `run_logs/official_tests/`. See the [official-test reference](official_tests.md) for every case and validation criterion.

## Synthetic-source tests

### Ideal parallel light

```bash
mkdir -p run_logs/official_tests/perfect_parallel
build/run_optical_sim configs/official_tests/perfect_parallel_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/perfect_parallel/run.log

MPLBACKEND=Agg python3 python/plot_spot_histogram.py \
  run_logs/official_tests/perfect_parallel/hits.csv \
  --config configs/official_tests/perfect_parallel_whiteboard.cfg \
  --output run_logs/official_tests/perfect_parallel/spot.png \
  --max-bins 520 --dpi 350
```

This cfg combines the complete 54-facet mirror, one million on-axis parallel photons, and the 8 m whiteboard.

### 900 m point source

```bash
mkdir -p run_logs/official_tests/point_900m
build/run_optical_sim configs/official_tests/perfect_point_900m_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/point_900m/run.log

MPLBACKEND=Agg python3 python/plot_spot_histogram.py \
  run_logs/official_tests/point_900m/hits.csv \
  --config configs/official_tests/perfect_point_900m_whiteboard.cfg \
  --output run_logs/official_tests/point_900m/spot.png
```

### Parallel light with structural obstruction

```bash
mkdir -p run_logs/official_tests/raytrace_structure_parallel
build/run_optical_sim \
  configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg

MPLBACKEND=Agg python3 python/plot_optical_layout_3d.py \
  --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  --show-obstruction \
  --output run_logs/official_tests/raytrace_structure_parallel/layout_3d.png

python3 python/plot_optical_layout_html.py \
  --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  --output run_logs/official_tests/raytrace_structure_parallel/layout_3d.html
```

The obstruction is enabled by:

```ini
obstruction.config=../obstructions/raytrace_final_structure.cfg
```

The one-command official script also runs the elevation-dependent deformation scan. See [test 5](official_tests.md#5-支架形变仰角扫描) for its standalone controls.

## Photon CSV

The minimum input is:

```csv
x_m,y_m,z_m,dir_x,dir_y,dir_z
0.0,0.0,1.0,0.0,0.0,-1.0
```

Run the supplied example:

```bash
build/run_optical_sim configs/examples/photon_csv_local_whiteboard.cfg
```

After copying the cfg, users normally edit:

```ini
source.csv_path=/path/to/photons.csv
source.coordinate_frame=telescope_local
output.csv=run_logs/my_photons/hits.csv
```

To save the photon position and direction actually passed into tracing:

```ini
output.whiteboard_input_photon=true
```

Plot these directions on the 3D telescope layout:

```bash
python3 python/plot_optical_layout_html.py \
  --config configs/examples/photon_csv_local_whiteboard.cfg \
  --input-photon-csv run_logs/examples/photon_csv_local_whiteboard/hits.csv \
  --output run_logs/examples/photon_csv_local_whiteboard/layout_with_photons.html
```

See the [Photon CSV format](photon_csv_format.md) for optional columns, input frames, and normalization tools.

## CORSIKA whiteboard

Use this before enabling the camera to check EventIO input, pointing, and focal-plane spots:

```bash
mkdir -p run_logs/official_tests/corsika/whiteboard_plots
build/run_corsika_trace \
  configs/official_tests/corsika_whiteboard.cfg \
  /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/whiteboard_run.log
```

Main outputs:

```text
run_logs/official_tests/corsika/whiteboard_hits.csv
run_logs/official_tests/corsika/whiteboard_summary.csv
```

Plot a known output event id:

```bash
MPLBACKEND=Agg python3 python/plot_corsika_trace_output.py \
  run_logs/official_tests/corsika/whiteboard_hits.csv \
  --summary-csv run_logs/official_tests/corsika/whiteboard_summary.csv \
  --event-id OUTPUT_EVENT_ID \
  --output-dir run_logs/official_tests/corsika/whiteboard_plots
```

The one-command script selects a common event automatically and writes `run_logs/official_tests/corsika/plots/selected_event.env`.

## CORSIKA camera and full response

### Ideal camera HDF5

```bash
build/run_corsika_trace \
  configs/official_tests/corsika_new_camera.cfg \
  /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_run.log
```

This cfg uses the standard mirror, real camera pixels, and light collector while keeping atmosphere, NSB, trigger, and PDE losses disabled. It writes:

```text
run_logs/official_tests/corsika/camera_dense.h5
```

Plot the first image:

```bash
MPLBACKEND=Agg python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --image-index 0 --quantity pe \
  --output run_logs/official_tests/corsika/camera_first.png
```

### NSB and trigger

```bash
build/run_corsika_trace \
  configs/official_tests/corsika_nsb_trigger_camera.cfg \
  /path/to/input.zst
```

The cfg adds measured optical efficiencies, SiPM PDE, spectral NSB, trigger tables, and waveform output.

### Full response

```bash
build/run_corsika_trace \
  configs/official_tests/corsika_full_response_camera.cfg \
  /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_full_response_run.log
```

This cfg enables mirror deformation, random errors, reflectivity, filter transmission, MODTRAN atmosphere, SiPM PDE, structural obstruction, and trigger. It writes:

```text
run_logs/official_tests/corsika/camera_full_response_dense.h5
```

Camera plot:

```bash
MPLBACKEND=Agg python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_full_response_dense.h5 \
  --image-index 0 --quantity pe \
  --output run_logs/official_tests/corsika/full_response_first.png
```

For all telescopes in one event, use `python/select_hdf5_event.py` or the selected-event file produced by the official script.

### ROOT/pylast

```bash
build/run_corsika_trace \
  configs/examples/lactroot_only.cfg \
  /path/to/input.zst
```

Default output:

```text
run_logs/lactroot_only/lact_events.root
```

Quick-look plots:

```bash
python3 python/plot_lact_root_output.py \
  run_logs/lactroot_only/lact_events.root \
  --outdir run_logs/lactroot_only/root_quicklook
```

The cfg uses `timeseries_pe`, processes all showers in the input file, and stores
only telescope events that pass the trigger. See the
[ROOT server check](server_root_output_check_zh.md), [pylast data levels](pylast_event_data_levels_zh.md), and `notebooks/lact_root_to_pylast_visualize.ipynb`.

## Common cfg fields

```ini
# Telescope pointing
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70

# Input
source.mode=EventIO
source.eventio_path=
source.max_shower_events=10

# Output
output.format=hdf5
output.hdf5_path=run_logs/my_run/output.h5
output.save_only_triggered=false
```

The EventIO file is normally passed as the second command-line argument, so `source.eventio_path` can remain empty. Use `run_optical_sim` for synthetic/Photon CSV sources and `run_corsika_trace` for EventIO.

## Troubleshooting

- HDF5 is not found: install the HDF5 C development package, remove `build/`, and rebuild.
- ROOT is not found: activate the ROOT environment before CMake, or use `make no-root`.
- Plotting fails on a server: set `MPLBACKEND=Agg` and `MPLCONFIGDIR=/tmp`.
- CORSIKA output is empty: start with `corsika_whiteboard.cfg`, disable triggered-only output, and reduce filters.
- The event id is unknown: use the official script or `python/select_hdf5_event.py`.

## Detailed references

- [All official tests](official_tests.md)
- [Photon CSV format](photon_csv_format.md)
- [Coordinate conventions](coordinate_systems.md)
- [CORSIKA/EventIO adapter](corsika_eventio_adapter.md)
- [HDF5 output format](hdf5_output_format.md)
- [Camera timing and waveforms](camera_timing_waveform_zh.md)
- [NSB spectral model](nsb_spectral_model_zh.md)
- [Structural obstruction](obstruction_primitives_csv_zh.md)
