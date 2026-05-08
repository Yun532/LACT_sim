# LACT_sim HDF5 Output Format

CSV remains the primary debugging format. HDF5 is the recommended packed format
for batch camera-image analysis and is the default target for formal CORSIKA
camera runs.

Native C++ output from `run_corsika_trace`:

```ini
output.format=hdf5
output.hdf5_path=run_logs/my_corsika_run/corsika_trace.h5
output.hdf5_storage=dense
```

Use `output.format=both` to write HDF5 and CSV together for validation.  Use
`output.format=csv` on machines where HDF5 is not installed.

Create an HDF5 file from an existing pixel CSV:

```bash
python3 python/export_trace_hdf5.py \
  --pixel-csv run_logs/half_mirror_50m/camera_pixels.csv \
  --config configs/experiments/half_mirror_point_50m_new_camera.cfg \
  --output run_logs/half_mirror_50m/camera_sparse.h5 \
  --storage sparse
```

Supported storage modes:

```text
sparse  only non-zero / hit pixels are stored per image
dense   every image stores a full [n_pixels] camera vector
both    stores both sparse and dense views
```

The native C++ writer supports `sparse`, `dense`, and `both` for camera images.
For production CORSIKA camera output, use `dense`: each image is a full
event/telescope vector with one value per camera pixel, currently 1616 pixels
for `new_camera`.

## Layout

```text
/config
  main_config_text
  attrs: expanded key/value configuration
  /components
    camera_text
    mirror_text
    output_text
    source_text
    ...

/metadata/electronics
  attrs: placeholder model

/metadata/sipm
  attrs: size_m, pde

/metadata/efficiency
  attrs: mirror/filter/atmosphere/funnel settings

/metadata/nsb
  attrs: enabled, model, rate_pe_per_ns_per_pixel, window_ns, seed

/metadata/trigger
  attrs: enabled, simple multiplicity thresholds

/metadata/coordinates
  attrs: array and pointing coordinate conventions

/telescopes/table
  telescope_id
  x_m, y_m, z_m                       # compatibility aliases
  array_x_north_m, array_y_west_m, array_z_up_m
  radius_m
  pointing_az_deg, pointing_el_deg
  focal_length_m

/camera/pixels
  pixel_id
  x_m, y_m
  size_m
  shape_code

/mirrors/facets
  mirror_id
  center_x/y/z_m
  normal_x/y/z
  radius_of_curvature_m
  size1_m, size2_m
  aperture_rotation_rad
  shape_code

/events/table
  event_index
  event_id

/events/corsika
  event_id
  shower_event_id
  array_id
  primary_type
  energy_gev
  theta_deg
  phi_deg
  azimuth_north_to_east_deg
  core_x_north_m
  core_y_west_m
  array_rotation_deg

For CORSIKA files with multiple IACT array offsets (`CSCAT`), these
`/events/corsika` core coordinates are the effective core for the specific
`event_id = shower_event * 100 + array_id`. EventIO `MC_TELOFF` stores the
detector-array offset with respect to the shower core, so LACT_sim writes the
core in the input/telescope array frame as:

```text
core_x_north_m = -MC_TELOFF.xoff
core_y_west_m  = -MC_TELOFF.yoff
```

The stored coordinates are still CORSIKA IACT/NWU coordinates, not ENU.

/events/corsika_showers
  shower_event_id
  primary_type
  energy_gev
  theta_deg
  phi_deg
  azimuth_north_to_east_deg
  core_x_north_m
  core_y_west_m
  array_rotation_deg

`/events/corsika_showers` keeps the raw shower-header core coordinates once per
shower before selecting an array-offset stream.

/images/index
  image_index
  event_id
  telescope_id
  start
  count
  total_photons
  total_pe
  total_signal
  time_mean_ns
  time_rms_ns

/images/sparse/pixels
  pixel_id
  photon_count
  pe
  signal
  time_mean_ns
  time_rms_ns

/images/dense/signal
/images/dense/pe
/images/dense/photon_count
/images/dense/pixel_id_axis

# optional when output.hdf5_write_components=true
/images/dense/cherenkov_pe
/images/dense/nsb_pe

/trigger/telescope
  event_id
  telescope_id
  triggered
  n_pixels_above_threshold
  total_pe
  trigger_time_ns

/trigger/array
  event_id
  array_triggered
  n_triggered_telescopes
```

`/camera`, `/mirrors`, `/telescopes`, and `/config` are static for the file and
are stored once. Event/telescope images only reference them by `event_id` and
`telescope_id`. The formal output intentionally uses one event identifier.
When `source.event_id_mode=event_array100`, that `event_id` is already the
combined CORSIKA shower/array stream id used by LACT_sim.

## Coordinate Metadata

For CORSIKA/EventIO outputs, telescope positions use the CORSIKA IACT horizontal
array frame:

```text
array_x_north_m : CORSIKA magnetic North positive
array_y_west_m  : West positive
array_z_up_m    : Up positive
```

This is a magnetic-North-West-Up (NWU) frame. LACT_sim keeps CORSIKA/EventIO
provenance tables in this frame and does not apply a magnetic-declination
rotation to true geographic north. Plotting tools that show conventional East-North axes
convert only for display with:

```text
plot_east_m  = -array_y_west_m
plot_north_m =  array_x_north_m
```

Pointing angles use azimuth measured from North toward East and elevation above
the horizon:

```text
pointing_az_deg = 0   -> North / +array_x
pointing_az_deg = 90  -> East / -array_y
pointing_el_deg = 90 - zenith_deg
```

The legacy `x_m,y_m,z_m` columns in `/telescopes/table` are kept as aliases for
older scripts. New analysis should prefer the explicit
`array_x_north_m,array_y_west_m,array_z_up_m` fields and read
`/metadata/coordinates` for the convention.

`radius_m` is the EventIO telescope/detector radius. Some CORSIKA producers
write `TELESCOPE x y 0 400` in the input card while hessio's `MC_TELPOS` block
reports `array_z_up_m=4 m` and `radius_m=4 m`; the HDF5 table preserves the
EventIO values actually used by the file.

Camera image coordinates are the output-plane `u/v` axes. The camera pixel table
stores those coordinates as `x_m,y_m`, matching the plot axes in
`python/plot_hdf5_camera.py`.

`/events/table` remains the minimal formal event index. `/events/corsika` is an
optional provenance table for CORSIKA/EventIO runs, useful for plotting array
core positions and arrival directions without rereading the EventIO file.

## Sparse Images

Sparse storage stores only pixels with signal. `/images/index` tells you where
one image lives inside `/images/sparse/pixels`:

```python
row = h5["images/index"][i]
start = row["start"]
count = row["count"]
pixels = h5["images/sparse/pixels"][start:start + count]
```

This is efficient before NSB/noise is added, when most camera pixels are zero.
It is mainly useful for debugging and cross-checking.

## Dense Images

Dense storage stores a full image vector for every event/telescope:

```python
signal = h5["images/dense/signal"][image_index, :]
pixel_ids = h5["images/dense/pixel_id_axis"][:]
```

This is the recommended production format. It is preferable after
NSB/background/electronics noise is added, when every pixel has a baseline
value, and it also avoids needing an external pixel CSV when reading images.

`/images/dense/pe` is the final integrated p.e. image. If NSB is enabled, it
already includes the Poisson NSB contribution. Set
`output.hdf5_write_components=true` to also save the Cherenkov-only and NSB-only
components for debugging.

## Trigger Tables

The first trigger implementation is a simple multiplicity model. Telescope rows
record the number of pixels above `trigger.pixel_threshold_pe`; array rows
record how many telescopes triggered for the same `event_id`.

When `output.save_only_triggered=true`, only triggered telescope images are
written under `/images`. The `/trigger` tables remain complete and should be
used to inspect non-triggered telescope decisions.

## Plotting

The HDF5 plotter reads camera geometry from the same file:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/half_mirror_50m/camera_sparse.h5 \
  --event-id 0 \
  --telescope-id 0 \
  --output run_logs/half_mirror_50m/hdf5_camera.png
```

For CORSIKA/EventIO HDF5 files, camera and array plots support the same event
selection forms:

```text
--event-id 46889802
--shower-event-id 468898 --array-id 2
--shower-event-number 1 --array-id 2
```

With the recommended `source.event_id_mode=event_array100`, these identify the
same output event because:

```text
event_id = shower_event * 100 + array_id
```

The original CORSIKA shower id may therefore look like the HDF5 event id with
the last two digits removed. Omit `--telescope-id` in
`python/plot_hdf5_camera.py` to write one camera image for every telescope in
the selected event.

Plot the telescope distribution only:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_array_layout.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --output run_logs/official_tests/corsika/array_layout.png
```

Plot total p.e. per telescope for one event. The figure uses East on the x-axis
and North on the y-axis, while reading telescope coordinates from the
CORSIKA-compatible HDF5 fields:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_array_layout.py \
  run_logs/official_tests/corsika/camera_dense.h5 \
  --event-id 46889802 \
  --quantity pe \
  --output run_logs/official_tests/corsika/array_event46889802_pe.png
```

For p.e. maps, the colorbar is logarithmic by default. Add `--linear-color` for
a linear p.e. colorbar during debugging.
