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

/telescopes/table
  telescope_id
  x_m, y_m, z_m
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
```

`/camera`, `/mirrors`, `/telescopes`, and `/config` are static for the file and
are stored once. Event/telescope images only reference them by `event_id` and
`telescope_id`. The formal output intentionally uses one event identifier.
When `source.event_id_mode=event_array100`, that `event_id` is already the
combined CORSIKA shower/array stream id used by LACT_sim.

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

## Plotting

The HDF5 plotter reads camera geometry from the same file:

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/tmp python3 python/plot_hdf5_camera.py \
  run_logs/half_mirror_50m/camera_sparse.h5 \
  --event-id 0 \
  --telescope-id 0 \
  --output run_logs/half_mirror_50m/hdf5_camera.png
```
