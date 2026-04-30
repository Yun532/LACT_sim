# LACT_sim

LACT_sim is a C++ optical ray-tracing program for the LACT 1229 mirror layout.
It supports:

- perfect optical tests with synthetic parallel light and finite-distance point sources
- elevation-dependent structural deformation of mirror facets
- direct CORSIKA/EventIO photon-bunch input through the vendored hessioxxx library
- whiteboard output-plane images
- the current `new_camera` square-pixel camera, light collector, SiPM, and p.e. output

Generated build folders, logs, CSV outputs, and plots are intentionally not part
of the source package. They can be regenerated with the commands below.

## Repository Layout

```text
apps/        C++ executables and tests
configs/     mirror, source, output, camera, electronics, and official-test configs
docs/        coordinate-system and EventIO notes
include/     public C++ headers
python/      plotting and conversion utilities
src/         C++ implementation
tools/       build and reproduction scripts
external/    vendored hessioxxx source and metadata
```

## Server Build

From the repository root:

```bash
cd /path/to/LACT_sim
```

For optical tests that do not need CORSIKA/EventIO:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLACT_ENABLE_HESSIO=OFF
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

For CORSIKA/EventIO input and dense HDF5 camera output, install the HDF5 C
development package on the server first.  Common examples:

```bash
# Ubuntu/Debian
sudo apt-get install cmake g++ make zlib1g-dev libhdf5-dev

# CentOS/RHEL-like systems
sudo yum install cmake gcc-c++ make zlib-devel hdf5-devel
```

Then build hessioxxx and LACT_sim:

```bash
tools/build_hessio.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DHESSIO_ROOT="$PWD/external/hessioxxx/source"
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

Native HDF5 camera output is enabled automatically when the HDF5 C library is
available. Check the CMake output for a line like `Found HDF5`. To force a
CSV-only build, add:

```bash
-DLACT_ENABLE_HDF5=OFF
```

If you build hessioxxx manually, make sure these folders exist first:

```bash
mkdir -p external/hessioxxx/source/out \
         external/hessioxxx/source/bin \
         external/hessioxxx/source/lib
make -C external/hessioxxx/source
```

Before running EventIO executables, set the runtime library path.

Linux:

```bash
export LD_LIBRARY_PATH="$PWD/external/hessioxxx/source/lib:${LD_LIBRARY_PATH:-}"
```

macOS:

```bash
export DYLD_LIBRARY_PATH="$PWD/external/hessioxxx/source/lib:${DYLD_LIBRARY_PATH:-}"
```

For server plotting without a display, use:

```bash
export MPLBACKEND=Agg
export MPLCONFIGDIR=/tmp
```

## One-Command Reproduction

The helper script rebuilds the C++ code, runs tests, runs the official optical
cases, and plots the standard figures.

Without CORSIKA:

```bash
tools/run_official_tests.sh --no-corsika
```

With CORSIKA/EventIO:

```bash
tools/run_official_tests.sh --corsika-file /path/to/input.zst
```

Outputs are written under:

```text
run_logs/official_tests/
```

The sections below show the same tests as explicit commands.

## Configuring CORSIKA Pointing

For CORSIKA/EventIO runs, the input file gives photon positions and directions
in the CORSIKA IACT horizontal frame. LACT_sim then rotates those photons into
the telescope-local optical frame using the telescope pointing in the cfg:

```ini
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70
source.eventio_coordinate_frame=corsika_iact
```

Use the CORSIKA run direction to set these values. If the CORSIKA file name or
input card says:

```text
zenith = Z deg
azimuth = A deg
```

then use:

```ini
telescope.pointing_el_deg = 90 - Z
telescope.pointing_az_deg = A
```

For example:

```text
zenith_20 azimuth_0
```

becomes:

```ini
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70
```

The azimuth convention is the hessio/sim_telarray-compatible one documented in
`docs/coordinate_systems.md`: `0 deg` points along array `+X` / North, and
`90 deg` points toward East, represented as array `-Y` in the CORSIKA frame.

For the smallest normal CORSIKA camera run, copy this template and edit only
pointing plus input/output paths:

```text
configs/templates/minimal_corsika_camera.cfg
```

For a more verbose example with all important sections visible, copy:

```text
configs/examples/corsika_new_user_full.cfg
```

Example:

```bash
cp configs/examples/corsika_new_user_full.cfg my_corsika_run.cfg
```

Then edit:

```ini
telescope.pointing_az_deg=...
telescope.pointing_el_deg=...
output.format=hdf5
output.hdf5_path=run_logs/my_corsika_run/corsika_trace.h5
```

Optional atmosphere transmission is configured as an extra factor after the
CORSIKA photons arrive at the telescope plane:

```ini
atmosphere.transmission=none
atmosphere.transmission=0.92
atmosphere.transmission=configs/atmosphere/my_transmission.csv
```

The SiPM/electronics block is currently an interface layer. It converts collected
photons into integrated p.e. through `electronics.pe_conversion`; detailed
waveform, saturation, crosstalk, afterpulse, and trigger electronics are reserved
for later replacement with the real implementation.

Run:

```bash
build/run_corsika_trace my_corsika_run.cfg /path/to/input.zst
```

For a quick server smoke test on a large CORSIKA file, add one of:

```ini
source.max_shower_events=1
source.filter_shower_event_id=327666
```

`source.max_shower_events=N` streams only the first `N` matching shower events.
`source.filter_shower_event_id=ID` keeps only one original CORSIKA shower event.

## Test 1: Perfect Optics Whiteboard Spots

This test uses the ideal 1229 mirror geometry, no structural deformation, and a
whiteboard output plane at the 8 m focal plane.

### Parallel Light

Run ray tracing:

```bash
mkdir -p run_logs/official_tests/perfect_parallel
build/run_optical_sim configs/official_tests/perfect_parallel_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/perfect_parallel/run.log
```

Plot the whiteboard spot:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_spot_histogram.py \
  run_logs/official_tests/perfect_parallel/hits.csv \
  --output run_logs/official_tests/perfect_parallel/spot.png \
  --max-bins 520 \
  --dpi 350 \
  --title "Perfect optics: on-axis parallel light"
```

Main outputs:

```text
run_logs/official_tests/perfect_parallel/hits.csv
run_logs/official_tests/perfect_parallel/run.log
run_logs/official_tests/perfect_parallel/spot.png
```

### 900 m Point Source

Run ray tracing:

```bash
mkdir -p run_logs/official_tests/point_900m
build/run_optical_sim configs/official_tests/perfect_point_900m_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/point_900m/run.log
```

Plot the whiteboard spot:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_spot_histogram.py \
  run_logs/official_tests/point_900m/hits.csv \
  --output run_logs/official_tests/point_900m/spot.png \
  --max-bins 520 \
  --dpi 350 \
  --title "Perfect optics: 900 m point source"
```

Main outputs:

```text
run_logs/official_tests/point_900m/hits.csv
run_logs/official_tests/point_900m/run.log
run_logs/official_tests/point_900m/spot.png
```

## Test 2: Structural Deformation Scan

This test enables the elevation-dependent mirror series:

```text
configs/mirror_1229_elevation_series.csv
configs/errors/structural_deformation_1229.cfg
configs/official_tests/deformation_parallel_whiteboard.cfg
```

The scan interpolates the mirror state at the requested elevation and runs a
parallel-light whiteboard trace for each elevation.

Run the 0 to 90 degree scan:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/run_elevation_parallel_scan.py \
  --config configs/official_tests/deformation_parallel_whiteboard.cfg \
  --run-binary build/run_optical_sim \
  --elevations 0,10,20,30,40,50,60,70,80,90 \
  --n-bunches 100000 \
  --output-dir run_logs/official_tests/deformation_scan \
  --dpi 350
```

Main outputs:

```text
run_logs/official_tests/deformation_scan/elevation_scan_metrics.csv
run_logs/official_tests/deformation_scan/elevation_spot_grid.png
run_logs/official_tests/deformation_scan/elevation_metrics.png
run_logs/official_tests/deformation_scan/hits_el_*.csv
run_logs/official_tests/deformation_scan/run_el_*.log
```

Plot the 3D mirror and whiteboard layout for selected elevations:

```bash
mkdir -p run_logs/official_tests/deformation_scan/layouts

for el in 0 20 45 60 90; do
  MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_optical_layout_3d.py \
    --config configs/official_tests/deformation_parallel_whiteboard.cfg \
    --elevation-deg "$el" \
    --output "run_logs/official_tests/deformation_scan/layouts/layout_el_${el}.png" \
    --dpi 350 \
    --ray-stride 2 \
    --view 32,-58
done
```

For a top-view mirror check:

```bash
for el in 0 20 45 60 90; do
  MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_mirror_top_view.py \
    --config configs/official_tests/deformation_parallel_whiteboard.cfg \
    --elevation-deg "$el" \
    --frame local \
    --output "run_logs/official_tests/deformation_scan/layouts/top_el_${el}.png" \
    --dpi 350
done
```

## Test 3: CORSIKA/EventIO to Whiteboard and Camera Images

This test requires the hessio-enabled build. The input file is passed on the
command line; the configs keep `source.eventio_path` empty.

Recommended CORSIKA coordinate setting:

```ini
source.eventio_coordinate_frame=corsika_iact
```

This means photon positions are telescope-relative CORSIKA IACT coordinates,
then rotated into the LACT local optical frame using the telescope pointing
configuration. See `docs/coordinate_systems.md` for the full convention.

### Whiteboard Output

Run the CORSIKA whiteboard trace. This official whiteboard config intentionally
uses CSV because it saves per-photon hit positions for spot inspection:

```bash
mkdir -p run_logs/official_tests/corsika
build/run_corsika_trace \
  configs/official_tests/corsika_whiteboard.cfg \
  /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/whiteboard_run.log
```

Main outputs:

```text
run_logs/official_tests/corsika/whiteboard_hits.csv
run_logs/official_tests/corsika/whiteboard_summary.csv
run_logs/official_tests/corsika/whiteboard_run.log
```

The required setting is already in
`configs/official_tests/corsika_whiteboard.cfg`:

```ini
output.format=csv
output.hits_csv=run_logs/official_tests/corsika/whiteboard_hits.csv
output.summary_csv=run_logs/official_tests/corsika/whiteboard_summary.csv
```

Plot by CORSIKA shower-event order and array id:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_corsika_trace_output.py \
  run_logs/official_tests/corsika/whiteboard_hits.csv \
  --summary-csv run_logs/official_tests/corsika/whiteboard_summary.csv \
  --shower-event-number 1 \
  --array-id 0 \
  --output-dir run_logs/official_tests/corsika/plots/whiteboard_shower1_array0 \
  --dpi 350
```

Add `--telescope-id N` if you only want one telescope. If omitted, the plotting
script draws all telescopes with matching rows.

If you already know the exact output event id, you can also select it directly:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_corsika_trace_output.py \
  run_logs/official_tests/corsika/whiteboard_hits.csv \
  --event-id OUTPUT_EVENT_ID \
  --output-dir run_logs/official_tests/corsika/plots/whiteboard_event_id \
  --dpi 350
```

### Dense HDF5 Pixel Camera Output

Run the CORSIKA camera trace. This is the recommended production path: the
output file contains the static camera/mirror/telescope metadata once, plus one
full 1616-pixel camera image per event/telescope stream.

```bash
mkdir -p run_logs/official_tests/corsika
build/run_corsika_trace \
  configs/official_tests/corsika_new_camera.cfg \
  /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_run.log
```

Main outputs:

```text
run_logs/official_tests/corsika/camera_dense.h5
run_logs/official_tests/corsika/camera_run.log
```

The required setting is already in
`configs/official_tests/corsika_new_camera.cfg`:

```ini
output.format=hdf5
output.hdf5_path=run_logs/official_tests/corsika/camera_dense.h5
output.hdf5_storage=dense
```

The dense HDF5 camera image is already after pixel entrance, light collector,
SiPM, and p.e. conversion when these components are enabled in the config.
Inside the file:

```text
/images/dense/pe            [n_images, 1616]
/images/dense/signal        [n_images, 1616]
/images/dense/photon_count  [n_images, 1616]
/camera/pixels              camera geometry
/telescopes/table           telescope metadata
/events/table               event_index, event_id
```

Plot one telescope by CORSIKA shower-event order:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --shower-event-number 1 \
  --array-id 0 \
  --telescope-id 7 \
  --quantity pe \
  --output run_logs/official_tests/corsika/camera_shower1_array0_tel7_pe.png \
  --dpi 350
```

If you already know the output event id, select it directly:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --event-id OUTPUT_EVENT_ID \
  --telescope-id 7 \
  --quantity pe \
  --output run_logs/official_tests/corsika/camera_event_id_tel7_pe.png \
  --dpi 350
```

For CORSIKA text debugging, set `output.format=csv` or `output.format=both` and
provide `output.pixel_csv`/`output.summary_csv`.

For synthetic sources, `run_optical_sim` also supports compact camera output:

```ini
output.mode=pixel
output.pixel_csv=run_logs/my_camera/camera_pixels.csv
```

Use `output.mode=hits` for a full per-photon whiteboard/camera hit table, or
`output.mode=both` to write both files.

## HDF5 Packed Output

CSV remains useful for debugging. For CORSIKA camera runs, the C++ program can
write camera images and static geometry/config directly into one HDF5 file:

```ini
output.format=hdf5
output.hdf5_path=run_logs/my_corsika_run/corsika_trace.h5
output.hdf5_storage=dense
```

Use `output.format=csv` for text-only debug output, or `output.format=both` to
write HDF5 plus `output.pixel_csv`/`output.summary_csv` for regression checks.
For production camera output, `dense` is recommended: each event/telescope image
is stored as a full 1616-pixel vector, which is the natural format once NSB and
electronics noise are added. `sparse` is still available for compact no-noise
debugging, and `both` stores both views.

The Python exporter is still available for converting older CSV outputs:

```bash
python3 python/export_trace_hdf5.py \
  --pixel-csv run_logs/half_mirror_50m/camera_pixels.csv \
  --config configs/experiments/half_mirror_point_50m_new_camera.cfg \
  --output run_logs/half_mirror_50m/camera_sparse.h5 \
  --storage sparse
```

Use `--storage sparse` before NSB/noise, when only hit pixels need to be stored.
Use `--storage dense` after NSB/noise, when every pixel has a value. Use
`--storage both` when you want both views.

Plot directly from HDF5, without a separate camera pixel CSV:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --shower-event-number 1 \
  --array-id 0 \
  --telescope-id 7 \
  --quantity pe \
  --output run_logs/official_tests/corsika/camera_from_h5.png
```

See `docs/hdf5_output_format.md` for the file layout.

## Configs Used by the Official Tests

```text
configs/official_tests/perfect_parallel_whiteboard.cfg
configs/official_tests/perfect_point_900m_whiteboard.cfg
configs/official_tests/deformation_parallel_whiteboard.cfg
configs/official_tests/corsika_whiteboard.cfg
configs/official_tests/corsika_new_camera.cfg
```

Core shared components:

```text
configs/mirrors/mirror_1229_imported.cfg
configs/mirror_1229_facets.csv
configs/mirror_1229_elevation_series.csv
configs/outputs/whiteboard_f8.cfg
configs/cameras/new_camera.cfg
configs/cameras/new_camera_pixels.csv
configs/sipm/new_camera_sipm.cfg
configs/electronics/new_camera_pe.cfg
```

Reusable starting templates:

```text
configs/templates/perfect_whiteboard.cfg
configs/templates/real_camera_pixel.cfg
configs/templates/corsika_camera.cfg
```

## Coordinate Notes

Read this before interpreting off-axis images:

```text
docs/coordinate_systems.md
```

Important short version:

- synthetic `source.beam_theta_deg/source.beam_phi_deg` describe photon
  propagation direction, not the sky-source position
- `phi=0 deg` starts at local `+x`; `phi=90 deg` points to local `+y`
- a sky source in `+y` corresponds to photons traveling with a `-y` component
- the reflecting optical system forms an inverted image on the focal plane
- camera/whiteboard image coordinates are defined by
  `output.plane_u_axis` and `output.plane_v_axis`
