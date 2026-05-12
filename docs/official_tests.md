# Official Test Suite

This document explains the reproducible tests run by
`tools/run_official_tests.sh`. Each official cfg under `configs/official_tests/`
is meant to be small and explicit: it lists only the non-default modules needed
for that test.

Run non-CORSIKA tests:

```bash
tools/run_official_tests.sh --no-corsika 2>&1 | tee run_logs/official_tests/run_no_corsika.log
```

Run the full suite:

```bash
tools/run_official_tests.sh --corsika-file /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/run_all.log
```

## Shared Config Blocks

- `telescope.config`: telescope metadata and pointing.
- `mirror.config`: mirror facet layout.
- `source.config`: synthetic source for non-CORSIKA tests.
- `output.config`: whiteboard or focal-plane geometry.
- `camera.config`: pixel layout and collector. `new_camera.cfg` enables the
  Bezier square-cone collector with `true_reflect`.
- `sipm.config`: SiPM size and optional PDE curve.
- `atmosphere.config`: extra transmission after CORSIKA; ideal is off.
- `nsb.config`: night-sky background. Ideal is disabled.
- `trigger.config`: simple multiplicity trigger. Disabled unless specified.
- `obstruction.config`: imported 3D obstruction primitives.
- `error.config`: structural deformation and random optical errors.

`electronics.config=ideal_pe` is intentionally not used by official cfgs now.
Electronics is still a placeholder; SiPM PDE is configured through `sipm.pde`.

## Optical Whiteboard Tests

### Perfect Parallel Whiteboard

Cfg: `configs/official_tests/perfect_parallel_whiteboard.cfg`

Content:

- `telescope.config=telescope_1229_minimal.cfg`: vertical 1229 metadata.
- `mirror.config=../mirrors/mirror_1229_imported.cfg`: all 54 ideal facets.
- `source.config=../sources/parallel_1M_on_axis.cfg`: 1,000,000 on-axis
  parallel photons sampled over a 4 m radius disk.
- `output.config=../outputs/whiteboard_f8.cfg`: virtual whiteboard at `z=-8 m`.
- `output.csv=.../perfect_parallel/hits.csv`: photon hit table.

Run:

```bash
mkdir -p run_logs/official_tests/perfect_parallel
build/run_optical_sim configs/official_tests/perfect_parallel_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/perfect_parallel/run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_spot_histogram.py \
  run_logs/official_tests/perfect_parallel/hits.csv \
  --output run_logs/official_tests/perfect_parallel/spot.png \
  --max-bins 520 --dpi 350
```

Expected: no obstruction and no errors; `hit_output_plane` equals
`hit_output_before_obstruction`.

### Perfect 900 m Point Source

Cfg: `configs/official_tests/perfect_point_900m_whiteboard.cfg`

Content:

- Standard telescope and 54 ideal facets.
- `source.config=../sources/point_900m_on_axis.cfg`: point source at local
  `z=900 m`.
- `output.config=../outputs/whiteboard_point_900m_focus.cfg`: finite-distance
  focal plane at about `8.07175 m`.

Run:

```bash
mkdir -p run_logs/official_tests/point_900m
build/run_optical_sim configs/official_tests/perfect_point_900m_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/point_900m/run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_spot_histogram.py \
  run_logs/official_tests/point_900m/hits.csv \
  --output run_logs/official_tests/point_900m/spot.png \
  --max-bins 520 --dpi 350
```

Expected: finite-distance spot differs from the parallel PSF.

## Imported 3D Obstruction Tests

### Parallel Whiteboard With Obstruction

Cfg: `configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg`

Content:

- Same ideal parallel source as the perfect whiteboard test.
- `obstruction.config=../obstructions/raytrace_final_structure.cfg`: imported
  camera/support primitives are enabled.
- Normal physics mode: blocked photons are discarded.

Run:

```bash
mkdir -p run_logs/official_tests/raytrace_structure_parallel
build/run_optical_sim configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/raytrace_structure_parallel/run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_spot_histogram.py \
  run_logs/official_tests/raytrace_structure_parallel/hits.csv \
  --output run_logs/official_tests/raytrace_structure_parallel/spot.png \
  --max-bins 520 --dpi 350 \
  --title "Parallel beam with 3D obstruction"
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_optical_layout_3d.py \
  --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  --show-obstruction \
  --output run_logs/official_tests/raytrace_structure_parallel/layout_3d.png \
  --dpi 350
python3 python/plot_optical_layout_html.py \
  --config configs/official_tests/perfect_parallel_raytrace_structure_whiteboard.cfg \
  --output run_logs/official_tests/raytrace_structure_parallel/layout_3d.html
```

Expected: log reports blocked counts, transmission after obstruction, and
before/after equivalent collection areas.

### 30 m Point Source With Obstruction

Cfg: `configs/official_tests/point_30m_structure_whiteboard.cfg`

Content:

- `source.config=../sources/point_30m_from_whiteboard_on_axis.cfg`: source is
  30 m in front of the standard `z=-8 m` whiteboard, so local source `z=22 m`.
- Standard whiteboard remains fixed at `z=-8 m`.
- Obstruction primitives are enabled.

Run:

```bash
mkdir -p run_logs/official_tests/raytrace_structure_point_30m
build/run_optical_sim configs/official_tests/point_30m_structure_whiteboard.cfg \
  2>&1 | tee run_logs/official_tests/raytrace_structure_point_30m/run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_spot_histogram.py \
  run_logs/official_tests/raytrace_structure_point_30m/hits.csv \
  --output run_logs/official_tests/raytrace_structure_point_30m/spot.png \
  --max-bins 520 --dpi 350 \
  --title "30 m point source with 3D obstruction"
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_mirror_hit_map.py \
  run_logs/official_tests/raytrace_structure_point_30m/hits.csv \
  --config configs/official_tests/point_30m_structure_whiteboard.cfg \
  --require-surface --overlay-facets \
  --output run_logs/official_tests/raytrace_structure_point_30m/mirror_hits_with_facet_outlines.png \
  --dpi 350 \
  --title "30 m point source: mirror hit points with facet outlines"
```

Expected: spot and mirror-hit maps show finite-distance obstruction projection.

## Structural Deformation Scan

Cfg: `configs/official_tests/deformation_parallel_whiteboard.cfg`

Content:

- Standard ideal parallel setup.
- `error.config=../errors/structural_deformation_1229.cfg`: uses the
  elevation-dependent mirror series.
- The scan runner overrides elevation and photon count.

Run:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/run_elevation_parallel_scan.py \
  --config configs/official_tests/deformation_parallel_whiteboard.cfg \
  --run-binary build/run_optical_sim \
  --elevations 0,10,20,30,40,50,60,70,80,90 \
  --n-bunches 100000 \
  --output-dir run_logs/official_tests/deformation_scan
```

Expected: one PSF and one mirror-layout diagnostic per elevation.

## Light Collector Angular Response

This test uses a dedicated binary rather than an assembly cfg.

Model:

- One 2.44 cm square pixel.
- Bezier square-cone collector.
- `true_reflect` material.
- 1.30 cm square SiPM.

Run:

```bash
mkdir -p run_logs/official_tests/collector_angular_response
build/scan_light_collector_angular_response \
  --photons-per-angle 2000 \
  --angle-step-deg 1 \
  --max-angle-deg 90 \
  --output run_logs/official_tests/collector_angular_response/collector_angular_response.csv \
  2>&1 | tee run_logs/official_tests/collector_angular_response/run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_collector_angular_response.py \
  run_logs/official_tests/collector_angular_response/collector_angular_response.csv \
  --output run_logs/official_tests/collector_angular_response/collector_angular_response.png \
  --dpi 350
```

Expected: CSV and plot show geometric acceptance, weighted acceptance,
collector photon weight, and mean reflection count versus incidence angle.

## CORSIKA/EventIO Tests

All CORSIKA tests read the input file from the command line:

```bash
build/run_corsika_trace CONFIG.cfg /path/to/input.zst
```

All use:

- `source.mode=EventIO`
- `source.event_id_mode=event_array100`, so `event_id = shower_event * 100 + array_id`
- `source.eventio_coordinate_frame=corsika_iact`

### Whiteboard Debug

Cfg: `configs/official_tests/corsika_whiteboard.cfg`

Content:

- CORSIKA photons traced to a whiteboard.
- `output.format=csv`: saves `whiteboard_hits.csv` and `whiteboard_summary.csv`.
- No camera, collector, NSB, trigger, efficiency curves, errors, or obstruction.

Run:

```bash
build/run_corsika_trace configs/official_tests/corsika_whiteboard.cfg /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/whiteboard_run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_corsika_trace_output.py \
  run_logs/official_tests/corsika/whiteboard_hits.csv \
  --summary-csv run_logs/official_tests/corsika/whiteboard_summary.csv \
  --shower-event-number 1 --array-id 0 \
  --output-dir run_logs/official_tests/corsika/plots/shower1_array0_whiteboard
```

### Perfect Camera Dense HDF5

Cfg: `configs/official_tests/corsika_new_camera.cfg`

Content:

- `telescope.pointing_az_deg=0`, `telescope.pointing_el_deg=70`.
- Standard ideal mirrors and ideal atmosphere.
- `camera.config=../cameras/new_camera.cfg`: real pixel map plus Bezier
  collector.
- `sipm.config=../sipm/ideal_sipm.cfg`: SiPM size only; PDE off.
- NSB and trigger disabled.
- Dense HDF5 stores every pixel.

Run:

```bash
build/run_corsika_trace configs/official_tests/corsika_new_camera.cfg /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --shower-event-number 1 --array-id 2 --quantity pe \
  --output run_logs/official_tests/corsika/plots/shower1_array2_camera/all_tel_pe.png
```

### Camera + NSB + Trigger

Cfg: `configs/official_tests/corsika_nsb_trigger_camera.cfg`

Content:

- Same camera chain as the perfect camera test.
- `nsb.config=../nsb/example_constant_rate.cfg`: Poisson NSB with
  `0.05 p.e./ns/pixel * 16 ns = 0.8 p.e./pixel` mean.
- `trigger.config=../trigger/example_simple_multiplicity.cfg` and
  `trigger.pixel_threshold_pe=10`.
- `output.hdf5_write_components=true`: writes Cherenkov, NSB, and final p.e.
- `output.save_only_triggered=true`: keeps triggered telescope images.

Run:

```bash
build/run_corsika_trace configs/official_tests/corsika_nsb_trigger_camera.cfg /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_nsb_trigger_run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_nsb_trigger_dense.h5 \
  --shower-event-number 1 --array-id 2 --quantity pe \
  --output run_logs/official_tests/corsika/plots/shower1_array2_nsb_trigger/all_tel_final_pe.png
```

### Camera + Obstruction + NSB + Trigger

Cfg: `configs/official_tests/corsika_obstruction_nsb_trigger_camera.cfg`

Content:

- Same as the NSB+trigger camera test.
- Adds `obstruction.config=../obstructions/raytrace_final_structure.cfg`.
- Blocked photons are discarded before camera pixelization.
- This is the smoke test for the current end-to-end chain with structure
  shadowing.

Run:

```bash
build/run_corsika_trace configs/official_tests/corsika_obstruction_nsb_trigger_camera.cfg /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_obstruction_nsb_trigger_run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_obstruction_nsb_trigger_dense.h5 \
  --shower-event-number 1 --array-id 2 --quantity pe \
  --output run_logs/official_tests/corsika/plots/shower1_array2_obstruction_nsb_trigger/all_tel_final_pe.png
```

Expected: the log includes obstruction statistics and HDF5 includes triggered
final camera images.

### Full Response Smoke Test

Cfg: `configs/official_tests/corsika_full_response_camera.cfg`

Content:

- Structural deformation and small optical errors through
  `error.config=../errors/full_response_1229.cfg`.
- Mirror reflectivity and filter transmission through
  `efficiency.config=../efficiency/curves_all.cfg`.
- SiPM PDE through `sipm.config=../sipm/new_camera_sipm.cfg`.
- NSB is disabled to isolate optical/efficiency response plus trigger.
- Trigger threshold is set to 10 p.e.

Run:

```bash
build/run_corsika_trace configs/official_tests/corsika_full_response_camera.cfg /path/to/input.zst \
  2>&1 | tee run_logs/official_tests/corsika/camera_full_response_run.log
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/official_tests/corsika/camera_full_response_dense.h5 \
  --shower-event-number 1 --array-id 2 --quantity pe \
  --output run_logs/official_tests/corsika/plots/shower1_array2_full_response/all_tel_pe.png
```

Expected: images differ from the ideal camera test because deformation, random
errors, mirror/filter curves, and SiPM PDE are active.

## Diagnostic Cfgs Not Run By Default

- `configs/official_tests/inner_two_rings_parallel_whiteboard.cfg`: inner 18
  facets, no obstruction.
- `configs/official_tests/outer_ring_parallel_structure_mark_only_whiteboard.cfg`:
  outer 18 facets, obstruction marked but photons still propagated.
- `configs/official_tests/perfect_parallel_structure_mark_only_whiteboard.cfg`:
  all 54 facets, obstruction marked but photons still propagated.

For mark-only runs, use:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_obstruction_marked_hits.py \
  run_logs/official_tests/raytrace_structure_parallel_mark_obstructed/hits.csv \
  --space mirror \
  --config configs/official_tests/perfect_parallel_structure_mark_only_whiteboard.cfg \
  --overlay-facets \
  --output run_logs/official_tests/raytrace_structure_parallel_mark_obstructed/mirror_hits_marked_obstructed.png
```
