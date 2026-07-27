# LACT_sim

[中文](README.md) | [**English**](README_EN.md)

LACT_sim is a C++ optical ray tracer using the default LACT telescope mirror layout. It supports:

- parallel beams and finite-distance point sources;
- custom Photon CSV input;
- CORSIKA/EventIO photon input;
- mirror deformation, optical errors, and structural obstruction;
- camera, light collector, SiPM, NSB, and trigger response;
- whiteboard CSV, camera HDF5, and LACT ROOT/pylast output.

This page contains only the main user workflow and cfg settings. See the [complete user guide](docs/user_guide_en.md) and [official-test reference](docs/official_tests.md) for full commands and outputs.

## Install and build

```bash
git clone https://github.com/Yun532/LACT_sim.git
cd LACT_sim

# Ubuntu/Debian
sudo apt-get install cmake g++ make zlib1g-dev libhdf5-dev

# Basic plotting environment
python3 -m pip install numpy pandas matplotlib h5py

make
make test
```

The main executables are:

```text
build/run_optical_sim      parallel, point, and Photon CSV sources
build/run_corsika_trace    CORSIKA/EventIO, camera, HDF5, and ROOT
```

For synthetic sources without EventIO support:

```bash
make no-hessio
```

See the [complete build guide](docs/user_guide_en.md#install-and-build) for ROOT/pylast setup and troubleshooting.

## Run the official tests

Without CORSIKA:

```bash
tools/run_official_tests.sh --no-corsika
```

With a CORSIKA/EventIO file:

```bash
tools/run_official_tests.sh --corsika-file /path/to/input.zst
```

The script builds the code, runs CTest, traces the standard cases, and creates plots under:

```text
run_logs/official_tests/
```

## Parallel light

```bash
mkdir -p run_logs/official_tests/perfect_parallel
build/run_optical_sim configs/official_tests/perfect_parallel_whiteboard.cfg
```

Main cfg entries:

```ini
telescope.config=telescope_1229_minimal.cfg
mirror.config=../mirrors/mirror_1229_imported.cfg
source.config=../sources/parallel_1M_on_axis.cfg
output.config=../outputs/whiteboard_f8.cfg
output.csv=run_logs/official_tests/perfect_parallel/hits.csv
```

Plot the focal-plane spot:

```bash
MPLBACKEND=Agg python3 python/plot_spot_histogram.py \
  run_logs/official_tests/perfect_parallel/hits.csv \
  --config configs/official_tests/perfect_parallel_whiteboard.cfg \
  --output run_logs/official_tests/perfect_parallel/spot.png
```

Use `configs/official_tests/perfect_point_900m_whiteboard.cfg` for the finite-distance point-source example. Obstruction, deformation, and 3D-layout commands are in the [synthetic-source workflow](docs/user_guide_en.md#synthetic-source-tests).

## Photon CSV

```bash
build/run_optical_sim configs/examples/photon_csv_local_whiteboard.cfg
```

Typical user settings:

```ini
source.mode=PhotonCsv
source.csv_path=/path/to/photons.csv
source.coordinate_frame=telescope_local
output.mode=hits
output.csv=run_logs/my_photons/hits.csv
```

Required columns:

```csv
x_m,y_m,z_m,dir_x,dir_y,dir_z
3.9014025878906251,0.16619310379028321,0,-0.33789920806884766,-0.014721360988914967,-0.94106716376520094
```

See
[`configs/sources/photon_csv_six_column_example.csv`](configs/sources/photon_csv_six_column_example.csv)
for a small reference file. The plotting input,
[`configs/sources/event1909_tel19_minimal_photons.csv`](configs/sources/event1909_tel19_minimal_photons.csv),
contains event 1909, telescope 19 in the same six-column format.
Its `z_m=0` retains the raw CORSIKA 2D bunch reference plane. Both example
cfg files set `source.eventio_2d=true`, so the programs apply the same
telescope-local `-16 m` input-plane default and signed ray intersection as
direct EventIO input.

```bash
# Pure optics: whiteboard hits plus per-pixel photon counts; two diagnostic plots
build/run_optical_sim configs/examples/photon_csv_minimal_optics.cfg
python3 python/plot_minimal_photon_csv_outputs.py \
  --mode optics \
  --hits run_logs/examples/photon_csv_minimal/whiteboard_hits.csv \
  --photon-pixels run_logs/examples/photon_csv_minimal/camera_photon_counts.csv \
  --camera configs/cameras/new_camera_pixels.csv \
  --output-dir run_logs/examples/photon_csv_minimal/plots

# Full camera chain: write ROOT, then plot one p.e. camera image with pyLAST
build/run_corsika_trace configs/examples/photon_csv_full_camera_root.cfg
python3 python/plot_photon_csv_root_pylast.py \
  run_logs/examples/photon_csv_full_camera/lact_events.root \
  --event-id 1909 --telescope-id 19 \
  --output run_logs/examples/photon_csv_full_camera/camera_pe.png
```

The full example shares the normal `run_corsika_trace` camera processing and
ROOT writer with EventIO input. When the CSV omits `multiplicity`, the cfg
fallback is used; this example sets it to one photon per row. The cfg also
supplies a fixed `400 nm` wavelength and disables NSB, trigger, and waveform
generation. The two optical diagnostics retain the LACT_sim focal-plane
orientation; only the ROOT/pyLAST plot uses the sky-view conversion
(`pix_x=-x_m`, `pix_y=-y_m`). See the
[minimal Photon CSV example](docs/minimal_photon_csv.md) for details.

To save the photon position and direction actually passed into ray tracing:

```ini
output.whiteboard_input_photon=true
```

See the [Photon CSV format](docs/photon_csv_format.md) and [user guide](docs/user_guide_en.md#photon-csv) for optional columns, other input frames, and 3D HTML plotting.

## CORSIKA whiteboard

```bash
mkdir -p run_logs/official_tests/corsika
build/run_corsika_trace \
  configs/official_tests/corsika_whiteboard.cfg \
  /path/to/input.zst
```

Main cfg entries:

```ini
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70
source.mode=EventIO
source.coordinate_frame=corsika_nwu_relative
output.format=csv
output.hits_csv=run_logs/official_tests/corsika/whiteboard_hits.csv
```

See the [CORSIKA whiteboard workflow](docs/user_guide_en.md#corsika-whiteboard) for event selection and plotting.

## CORSIKA camera and full response

Ideal camera HDF5:

```bash
build/run_corsika_trace \
  configs/official_tests/corsika_new_camera.cfg \
  /path/to/input.zst
```

Full optical response:

```bash
build/run_corsika_trace \
  configs/official_tests/corsika_full_response_camera.cfg \
  /path/to/input.zst
```

Main cfg entries:

```ini
telescope.pointing_az_deg=0
telescope.pointing_el_deg=70
camera.config=../cameras/new_camera.cfg
source.mode=EventIO
source.coordinate_frame=corsika_nwu_relative
output.format=hdf5
output.hdf5_storage=dense
```

Plot the first camera image:

```bash
MPLBACKEND=Agg python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --image-index 0 \
  --quantity pe \
  --output run_logs/official_tests/corsika/camera_first.png
```

See the [complete CORSIKA workflow](docs/user_guide_en.md#corsika-camera-and-full-response) for NSB, trigger, obstruction, array plots, waveforms, and ROOT/pylast output.

## LACT ROOT / pylast

Write LACT ROOT only:

```bash
build/run_corsika_trace \
  configs/examples/lactroot_only.cfg \
  /path/to/input.zst
```

Default output:

```text
run_logs/lactroot_only/lact_events.root
```

The cfg stores integrated p.e. images, CORSIKA truth, telescope/camera metadata,
and sparse `timeseries_pe` waveforms. By default it processes all showers in the
input file and stores only telescope events that pass the trigger. Create
quick-look plots with:

```bash
python3 python/plot_lact_root_output.py \
  run_logs/lactroot_only/lact_events.root \
  --outdir run_logs/lactroot_only/root_quicklook
```

See the [ROOT and pylast check guide](docs/server_root_output_check_zh.md) for
pylast reading, ROOT trees, and the notebook workflow.

## Common cfg files

| Purpose | cfg |
|---|---|
| Ideal parallel light | `configs/official_tests/perfect_parallel_whiteboard.cfg` |
| 900 m point source | `configs/official_tests/perfect_point_900m_whiteboard.cfg` |
| Parallel light with obstruction | `configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg` |
| Photon CSV | `configs/examples/photon_csv_local_whiteboard.cfg` |
| CORSIKA whiteboard | `configs/official_tests/corsika_whiteboard.cfg` |
| CORSIKA ideal camera | `configs/official_tests/corsika_new_camera.cfg` |
| CORSIKA with NSB and trigger | `configs/official_tests/corsika_nsb_trigger_camera.cfg` |
| CORSIKA full response | `configs/official_tests/corsika_full_response_camera.cfg` |
| pylast ROOT-only | `configs/examples/lactroot_only.cfg` |

## Documentation

- [Documentation index by topic](docs/README.md)
- [Complete user guide](docs/user_guide_en.md)
- [All official tests](docs/official_tests.md)
- [Photon CSV format](docs/photon_csv_format.md)
- [Coordinate conventions](docs/coordinate_systems.md)
- [HDF5 output format](docs/hdf5_output_format.md)
- [Camera timing and waveforms](docs/camera_timing_waveform_zh.md)
- [NSB spectral model](docs/nsb_spectral_model_zh.md)
- [ROOT/pylast data levels](docs/pylast_event_data_levels_zh.md)
