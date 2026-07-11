# LACT_sim

[**English**](README_EN.md) | [中文](README.md)

LACT_sim is a C++ optical ray-tracing program for the LACT 1229 mirror layout.
It supports:

- perfect optical tests with synthetic parallel light and finite-distance point sources
- elevation-dependent structural deformation of mirror facets
- direct CORSIKA/EventIO photon-bunch input through the vendored hessioxxx library
- whiteboard output-plane images
- the current `new_camera` square-pixel camera, light collector, SiPM, and p.e. output

Generated build folders, logs, CSV outputs, and plots are intentionally not part
of the source package. They can be regenerated with the commands below.

## Chinese Overview

For a compact Chinese guide to the program flow, coordinate conventions, core
functions, output files, and current physics boundaries, see:

```text
docs/program_overview_zh.md
```

The default repository landing page is now the Chinese guide:
[`README.md`](README.md).

## Photon Input Coordinate Frames

`run_optical_sim` always converts source photons into the existing
telescope-local optical ray-tracing frame before tracing. The input setting
does not redefine the mirror, camera, output-plane, or plotting axes.

```ini
source.coordinate_frame=telescope_local
```

Supported values are:

- `telescope_local`: positions and directions are already in the selected
  telescope's local optical frame; this is the recommended Photon CSV default.
- `corsika_nwu_relative`: CORSIKA NWU vectors with positions already relative
  to the selected telescope.
- `corsika_nwu_global`: absolute CORSIKA NWU array positions; the configured
  telescope position is subtracted exactly once.
- `lact_generic_global`: the legacy LACT global frame whose azimuth is measured
  from global `+x` toward `+y`.

For whiteboard hit output, the post-conversion photon position and direction
actually passed into ray tracing can be appended directly to the CSV:

```ini
output.whiteboard_input_photon=true
```

This adds `input_local_x_m`, `input_local_y_m`, `input_local_z_m`, and
`input_local_dir_x/y/z`. These are captured from the ray-tracing input and are
not transformed again for output. See
[`docs/photon_csv_format.md`](docs/photon_csv_format.md) and
[`docs/coordinate_systems.md`](docs/coordinate_systems.md).

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

For a normal server build with CORSIKA/EventIO support:

```bash
make
make test
```

This builds the vendored hessioxxx library first, configures CMake in `build/`,
builds LACT_sim, and runs the CTest suite. The generated executables are linked
with a build-tree runtime path pointing at:

```text
external/hessioxxx/source/lib
```

so in the normal in-place build you do not need to re-export `LD_LIBRARY_PATH`
in every new terminal. If you move the compiled `build/` directory to another
location, rebuild from the repository root.

For optical tests that do not need CORSIKA/EventIO:

```bash
make no-hessio
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

The `make` command above is equivalent to:

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

Optional pylast-oriented ROOT output is enabled when ROOT is available at CMake
configure time. ROOT is not vendored under `external/`; use a system, module, or
conda ROOT installation. LACT_sim requires ROOT 6.24 or newer for this output.
Common setup patterns:

```bash
# Cluster module, if available
module load root
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLACT_ENABLE_ROOT=ON

# Explicit ROOT install
source /path/to/root/bin/thisroot.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DLACT_ENABLE_ROOT=ON \
  -DROOT_DIR="$(root-config --cmakedir)"

# If your root-config has no --cmakedir option
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DLACT_ENABLE_ROOT=ON \
  -DCMAKE_PREFIX_PATH="$(root-config --prefix)"

# Conda/mamba ROOT plus quick-look plotting tools
mamba create -n lact-root -c conda-forge root cmake compilers hdf5 zlib uproot matplotlib numpy
mamba activate lact-root
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DLACT_ENABLE_ROOT=ON \
  -DROOT_DIR="$(root-config --cmakedir)"
```

If ROOT is not found, the build continues with ROOT output disabled. To force a
build without ROOT probing:

```bash
make no-root
# or:
cmake -S . -B build -DLACT_ENABLE_ROOT=OFF
```

The ROOT output is an additional `run_corsika_trace` product for pylast adapter
workflows:

```ini
output.lact_root_enabled=true
output.lact_root_path=run_logs/my_corsika_run/lact_events.root
output.lact_profile=image_pe       # or timeseries_pe
output.lact_root_write_components=false
```

`image_pe` stores integrated p.e. camera images. `timeseries_pe` additionally
stores sparse p.e. time-series waveforms and per-pixel `time_peak_ns`. Both are
intended to map to pylast readout layers; pylast may derive DL0 from R1
waveforms. ROOT observations also include `image_cherenkov_pe` for pylast
`simulation.true_image`. Set `output.lact_root_write_components=true` to write
optional diagnostic components such as `image_nsb_pe`.

For server setup and quick-look validation plots, see
`docs/server_root_output_check_zh.md`. The plotting helper is:

```bash
python/plot_lact_root_output.py run_logs/my_corsika_run/lact_events.root \
  --outdir run_logs/my_corsika_run/root_quicklook
```

A runnable ROOT quick-look config is available at
`configs/examples/corsika_lact_root_quicklook.cfg`:

```bash
./build/run_corsika_trace \
  configs/examples/corsika_lact_root_quicklook.cfg \
  /path/to/input.simtel.zst
```

If you build hessioxxx manually, make sure these folders exist first:

```bash
mkdir -p external/hessioxxx/source/out \
         external/hessioxxx/source/bin \
         external/hessioxxx/source/lib
make -C external/hessioxxx/source
```

If a system strips or ignores executable RPATHs, set the runtime library path
manually before running EventIO executables.

Linux:

```bash
export LD_LIBRARY_PATH="$PWD/external/hessioxxx/source/lib:${LD_LIBRARY_PATH:-}"
```

macOS:

```bash
export DYLD_LIBRARY_PATH="$PWD/external/hessioxxx/source/lib:${DYLD_LIBRARY_PATH:-}"
```

To rebuild from scratch:

```bash
make clean      # remove the LACT_sim CMake build directory
make            # configure and build again
```

To also remove hessioxxx build products:

```bash
make distclean
make
```

For server plotting without a display, use:

```bash
export MPLBACKEND=Agg
export MPLCONFIGDIR=/tmp
```

## One-Command Reproduction

The helper script rebuilds the C++ code, runs tests, runs the official optical
cases, and plots the standard figures.

For a cfg-by-cfg explanation of every official test, including independent run
commands and expected outputs, see
[`docs/official_tests.md`](docs/official_tests.md).

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
the telescope-local optical frame using the telescope pointing in the cfg.
Use `run_corsika_trace` for these runs; `run_optical_sim` is for synthetic and
PhotonCsv optical debugging and rejects `source.mode=EventIO`.

```ini
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70
source.coordinate_frame=corsika_nwu_relative
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

The generic coordinate convention used by `run_optical_sim` is unchanged:
`buildTelescopeFrame()` measures azimuth from global `+x` toward `+y`.
`source.coordinate_frame` only controls how input rows are converted into the
existing runtime ray-tracing frame; it does not redefine the mirror, camera,
output, or plotting coordinates.

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

Optional atmosphere transmission can still be configured as a simple
wavelength-only extra factor:

```ini
atmosphere.transmission=none
atmosphere.transmission=0.92
atmosphere.transmission=configs/atmosphere/my_transmission.csv
```

For height- and wavelength-dependent MODTRAN total optical depth, use:

```ini
atmosphere.config=configs/atmosphere/modtran_4400_desert.cfg
```

This applies `T=exp(-tau_total)` before telescope optics. EventIO photon
bunches use their emission height when present; PhotonCsv input can provide an
optional `emission_altitude_km` column.

The SiPM block converts collected photons into integrated p.e. through
`sipm.pde`, which can be omitted, set to a constant, or set to a wavelength
table. Detailed waveform, saturation, crosstalk, afterpulse, and trigger
electronics are reserved for later replacement with the real implementation.

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

To measure where a long CORSIKA run spends time, enable the lightweight
profiler:

```ini
profile.enabled=true
```

The final log will include a `[Profile]` block with coarse timings for EventIO
streaming, coordinate transforms, optical tracing, obstruction checks, camera
response, accumulation, and HDF5 writing. `eventio_stream_s` is wall time for
streaming plus callback processing; the sub-stage timings inside the callback
are not exclusive.

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

## Test 2: Imported 3D Obstruction

This test uses the external simplified 3D obstruction model converted into:

```text
configs/obstructions/raytrace_final_structure_primitives.csv
```

It checks both incoming and reflected-ray obstruction. The official smoke test
keeps four plots: a parallel-light whiteboard spot, a parallel-light mirror-hit
distribution with projected mirror-facet outlines, a 30 m point-source
whiteboard spot, and the point-source mirror-hit distribution with projected
mirror-facet outlines.

Parallel-light obstruction run:

```bash
mkdir -p run_logs/official_tests/raytrace_structure_parallel
build/run_optical_sim configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/raytrace_structure_parallel/run.log

MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_spot_histogram.py \
  run_logs/official_tests/raytrace_structure_parallel/hits.csv \
  --output run_logs/official_tests/raytrace_structure_parallel/spot.png \
  --max-bins 520 \
  --dpi 350 \
  --title "Parallel beam with 3D obstruction"

MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_mirror_hit_map.py \
  run_logs/official_tests/raytrace_structure_parallel/hits.csv \
  --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  --require-surface \
  --overlay-facets \
  --output run_logs/official_tests/raytrace_structure_parallel/mirror_hits_with_facet_outlines.png \
  --dpi 350 \
  --title "Parallel beam: mirror hit points with facet outlines"
```

30 m point-source obstruction run:

```bash
mkdir -p run_logs/official_tests/raytrace_structure_point_30m
build/run_optical_sim configs/official_tests/point_30m_structure_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/raytrace_structure_point_30m/run.log

MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_spot_histogram.py \
  run_logs/official_tests/raytrace_structure_point_30m/hits.csv \
  --output run_logs/official_tests/raytrace_structure_point_30m/spot.png \
  --max-bins 520 \
  --dpi 350 \
  --title "30 m point source with 3D obstruction"
```

Mirror-hit distribution with projected facet outlines:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_mirror_hit_map.py \
  run_logs/official_tests/raytrace_structure_point_30m/hits.csv \
  --config configs/official_tests/point_30m_structure_whiteboard.cfg \
  --require-surface \
  --overlay-facets \
  --output run_logs/official_tests/raytrace_structure_point_30m/mirror_hits_with_facet_outlines.png \
  --dpi 350 \
  --title "30 m point source: mirror hit points with facet outlines"
```

Main outputs:

```text
run_logs/official_tests/raytrace_structure_parallel/spot.png
run_logs/official_tests/raytrace_structure_parallel/mirror_hits_with_facet_outlines.png
run_logs/official_tests/raytrace_structure_parallel/layout_3d.png
run_logs/official_tests/raytrace_structure_parallel/layout_3d.html
run_logs/official_tests/raytrace_structure_point_30m/spot.png
run_logs/official_tests/raytrace_structure_point_30m/mirror_hits_with_facet_outlines.png
```

The obstruction run log also reports true before/after obstruction statistics
and equivalent collection areas, for example:

```text
hit_output_before_obstruction
hit_output_plane
output_transmission_after_obstruction
output_loss_fraction_from_obstruction
source_sampling_area_m2
output_collecting_area_before_obstruction_m2
output_collecting_area_after_obstruction_m2
output_collecting_area_loss_from_obstruction_m2
```

## Test 3: Structural Deformation Scan

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

## Test 4: CORSIKA/EventIO to Whiteboard and Camera Images

This test requires the hessio-enabled build. The input file is passed on the
command line; the configs keep `source.eventio_path` empty.

Recommended CORSIKA coordinate setting:

```ini
source.coordinate_frame=corsika_nwu_relative
source.missing_wavelength_model=cherenkov
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

Plot by CORSIKA shower-event order and array id. Here `--array-id` is the
CORSIKA `CSCAT` / `MC_TELOFF` array-use index, not the telescope ID. With
`source.event_id_mode=event_array100`, output `event_id` is encoded as:

```text
event_id = shower_event * 100 + array_id
```

For a file with `CSCAT 10 ...`, the valid array IDs are normally `0..9`.
The official script now auto-selects one bright event common to all camera HDF5
outputs and writes it to:

```text
run_logs/official_tests/corsika/plots/selected_event.env
```

After sourcing that file, all CORSIKA plots can use the same `event_id`:

```bash
source run_logs/official_tests/corsika/plots/selected_event.env
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_corsika_trace_output.py \
  run_logs/official_tests/corsika/whiteboard_hits.csv \
  --summary-csv run_logs/official_tests/corsika/whiteboard_summary.csv \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --output-dir "run_logs/official_tests/corsika/plots/event_${LACT_SELECTED_EVENT_ID}/whiteboard" \
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
/metadata/nsb               NSB settings
/metadata/trigger           trigger settings
/trigger/telescope          telescope trigger rows
/trigger/array              array trigger rows
```

Plot all telescope camera images for the selected CORSIKA event. Omitting
`--telescope-id` writes one camera PNG for every telescope image available in
the selected event:

```bash
source run_logs/official_tests/corsika/plots/selected_event.env
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --quantity pe \
  --output "run_logs/official_tests/corsika/plots/event_${LACT_SELECTED_EVENT_ID}/camera/all_tel_pe" \
  --dpi 350
```

Add `--telescope-id N` only if you want one telescope:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --event-id "$LACT_SELECTED_EVENT_ID" \
  --telescope-id 3 \
  --quantity pe \
  --output "run_logs/official_tests/corsika/plots/event_${LACT_SELECTED_EVENT_ID}/camera/tel3_pe.png" \
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

You can also select by the original CORSIKA shower event id, which is often the
number printed in the EventIO log:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --shower-event-id 468898 \
  --array-id 2 \
  --telescope-id 3 \
  --quantity pe \
  --output run_logs/official_tests/corsika/camera_shower468898_array2_tel3_pe.png \
  --dpi 350
```

Plot the full telescope array layout, with North up and East right:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_array_layout.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --output run_logs/official_tests/corsika/array_layout.png \
  --dpi 350
```

Plot one event with total p.e. per telescope shown by a logarithmic colorbar.
For CORSIKA HDF5 files written by the current `run_corsika_trace`, the script
also marks the core position and draws the horizontal arrival-direction arrow
from `/events/corsika`. The telescope labels are one-based (`T1`, `T2`, ...),
while the stored telescope IDs remain zero-based:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_array_layout.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --shower-event-id 468898 \
  --array-id 2 \
  --quantity pe \
  --output run_logs/official_tests/corsika/array_event46889802_pe.png \
  --dpi 350
```

For the same event, `--event-id 46889802` is equivalent to
`--shower-event-id 468898 --array-id 2` because the default CORSIKA setting uses
`event_id = shower_event * 100 + array_id`. So `46889800` is not a different
shower from `468898`; it is shower `468898` with array/core-offset stream `0`.
This is why HDF5 event IDs can look like the log shower ID with two extra
digits appended.

To inspect several shower cores near the array center, use:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_array_layout.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --show-nearby-cores 6 \
  --array-only-limits \
  --output run_logs/official_tests/corsika/array_nearby_cores.png \
  --dpi 350
```

For CORSIKA text debugging, set `output.format=csv` or `output.format=both` and
provide `output.pixel_csv`/`output.summary_csv`.

To run the official NSB + simple-trigger smoke test, use:

```bash
build/run_corsika_trace \
  configs/official_tests/corsika_nsb_trigger_camera.cfg \
  /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_nsb_trigger_run.log
```

This writes `camera_nsb_trigger_dense.h5`, where `/images/dense/pe` includes
per-pixel integrated Cherenkov plus Poisson NSB p.e. This official smoke test enables
`output.save_only_triggered=true`, so `/images/index` and `/images/dense/*`
contain only telescope images that pass the simple camera trigger. The full
trigger decision table is still kept in `/trigger/telescope` and
`/trigger/array`. With `output.hdf5_write_components=true`, the file also
contains `/images/dense/cherenkov_pe` and `/images/dense/nsb_pe`.
The proxy waveform uses `waveform.time_reference=image_first`, so compact GIF
frames are saved relative to each telescope image first Cherenkov photon time
T0; the subtracted absolute reference is stored in
`/waveforms/reference_time_ns`.
For the benchmark file used in development,
`lact_prod1_corsika_particle_gamma_energy_1000.0_10000.0_zenith_20.0_azimuth_0.0_run_2418_event_468898.zst`,
the official script plots all available/triggered telescope images for
`shower-event-number=1, array-id=2`. The NSB+trigger file saves only telescope
images that pass the trigger, so the number of PNGs may be smaller than the
full array.

To run the official full-response CORSIKA smoke test, use:

```bash
build/run_corsika_trace \
  configs/official_tests/corsika_full_response_camera.cfg \
  /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_full_response_run.log
```

This writes `camera_full_response_dense.h5`. It enables the current optical
error paths, the wavelength-dependent mirror reflectivity, filter transmission,
SiPM PDE curve, MODTRAN atmosphere absorption, and the simple multiplicity
trigger. SiPM PDE is configured only once in
`configs/sipm/new_camera_sipm.cfg` as `sipm.pde`; the electronics module
is still only a placeholder for future waveform/electronics effects. This test uses
`output.save_only_triggered=true`, so `/images/dense/*`
contains only telescope images that pass the trigger. The complete trigger
decision tables remain in `/trigger/telescope` and `/trigger/array`.

The official one-command script also plots the triggered camera images and the
array/core layout for the same automatically selected event:

```text
run_logs/official_tests/corsika/plots/event_<event_id>/full_response/all_tel_pe/
run_logs/official_tests/corsika/plots/event_<event_id>/layout/core_and_array_pe.png
```

The random optical error values in
`configs/errors/full_response_1229.cfg` are reproducible smoke-test values, not
a calibrated alignment model.

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
  --quantity pe \
  --output run_logs/official_tests/corsika/camera_from_h5.png
```

Omit `--telescope-id` to plot all telescope images for that event. Add
`--telescope-id N` when you want one camera only. In these HDF5 plotting
commands, `--array-id` has the same meaning: it selects the CORSIKA
array-use/core-offset stream for the chosen original shower event. It is
separate from `--telescope-id`.

See `docs/hdf5_output_format.md` for the file layout.

## Efficiency Curve Checks

Mirror reflectivity, filter transmission, SiPM PDE, and optional atmosphere
transmission are wavelength-dependent multiplicative factors. To inspect the
input tables and the total response curve:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_efficiency_curves.py \
  --output-dir run_logs/official_tests/efficiency_curves \
  2>&1 | tee run_logs/official_tests/efficiency_curves/run.log
```

This writes one plot per component, `total_efficiency.png`, sampled values in
`efficiency_curve_samples.csv`, and a short report. See
`docs/efficiency_validation.md` for the exact interpolation and duplicate-point
rules.

## NSB Spectral Rate Check

The official NSB smoke tests now use a SkyCalc LoNS spectrum instead of a
hand-written constant rate. To inspect the standalone rate calculation:

```bash
mkdir -p run_logs/official_tests/nsb_spectral

build/compute_nsb_rate configs/nsb/spectral_rate_check_with_obstruction.cfg \
  2>&1 | tee run_logs/official_tests/nsb_spectral/run.log

MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_nsb_spectral_rate.py \
  --effective-area-m2 22.606448 \
  --output run_logs/official_tests/nsb_spectral/nsb_spectral_response.png \
  --diagnostic-csv run_logs/official_tests/nsb_spectral/diagnostic.csv \
  --summary run_logs/official_tests/nsb_spectral/summary.txt \
  2>&1 | tee run_logs/official_tests/nsb_spectral/plot.log
```

The default dark-sky, obstruction-aware result is about
`0.074375315 pe/ns/pixel`. See `docs/nsb_spectral_model_zh.md` for the
LoNS unit conversion, fixed effective areas, and validation plot details.

## Configs Used by the Official Tests

```text
configs/official_tests/perfect_parallel_whiteboard.cfg
configs/official_tests/perfect_point_900m_whiteboard.cfg
configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg
configs/official_tests/point_30m_structure_whiteboard.cfg
configs/official_tests/deformation_parallel_whiteboard.cfg
configs/official_tests/corsika_whiteboard.cfg
configs/official_tests/corsika_new_camera.cfg
configs/official_tests/corsika_nsb_trigger_camera.cfg
configs/official_tests/corsika_obstruction_nsb_trigger_camera.cfg
configs/official_tests/corsika_full_response_camera.cfg
```

Core shared components:

```text
configs/mirrors/mirror_1229_imported.cfg
configs/mirror_1229_facets.csv
configs/mirror_1229_elevation_series.csv
configs/outputs/whiteboard_f8.cfg
configs/outputs/focal_plane_f8.cfg
configs/cameras/new_camera.cfg
configs/cameras/new_camera_pixels.csv
configs/sipm/ideal_sipm.cfg
configs/sipm/new_camera_sipm.cfg
configs/electronics/ideal_pe.cfg
configs/efficiency/curves_all.cfg
configs/errors/full_response_1229.cfg
configs/atmosphere/modtran_4400_desert.cfg
configs/nsb/ideal.cfg
configs/trigger/disabled.cfg
```

Use `configs/sipm/ideal_sipm.cfg` for pure optical/camera tests where the SiPM
aperture is present but PDE losses are disabled. Use
`configs/sipm/new_camera_sipm.cfg` for full-response runs with the measured SiPM
PDE curve.

Reusable starting templates:

```text
configs/templates/perfect_whiteboard.cfg
configs/templates/real_camera_pixel.cfg
configs/templates/corsika_camera.cfg
configs/templates/minimal_corsika_camera.cfg
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
- standard 2D plotters use a sky-up display convention when enough pointing
  metadata or `--config` is provided: display `+y` is the projection of global
  `+z/up` onto the plotted plane. Use `--raw-camera-xy` in
  `plot_hdf5_camera.py` when you need stored camera x/y without display
  rotation.
