# Minimal photon CSV

The simple `PhotonCsv` input is intended for optical tests and camera-image
generation without requiring CORSIKA event metadata.

## Required columns

```csv
x_m,y_m,z_m,dir_x,dir_y,dir_z
```

Each row represents one photon ray:

- `x_m,y_m,z_m`: photon position in metres.
- `dir_x,dir_y,dir_z`: unit propagation direction.
- The coordinate frame is selected by `source.coordinate_frame`.
- `telescope_local` is the simplest choice for hand-written optical tests.
- A CSV stripped from CORSIKA should remain in `corsika_nwu_relative`; the
  simulator then uses its existing CORSIKA coordinate adapter instead of a
  second conversion implemented in a helper script.

No bunch column is required. In this mode one CSV row contributes one photon.

## Optional quantities

Wavelength and time do not have to be present in the CSV:

- For a pure optical test, use ideal/disabled efficiency curves. Wavelength is
  then irrelevant.
- For a monochromatic camera test, set one value such as
  `source.wavelength_nm=400` in the configuration.
- If waveform simulation is not requested, photon arrival time is irrelevant
  and may be omitted.

The loader still accepts optional columns for advanced replay, but they are not
part of the minimal user interface.

## Included examples

Pure optics, whiteboard hits, and photon-count camera:

```text
configs/examples/photon_csv_minimal_optics.cfg
```

Expected-photoelectron camera image with wavelength-dependent detector
efficiency, but without waveform generation:

```text
configs/examples/photon_csv_minimal_pe_camera.cfg
```

The full expected-p.e. figure is rendered through
`pylast.visualize.plot_camera_image`, using the same angular axes and camera
orientation as the normal CORSIKA/ROOT pyLAST display. The whiteboard remains a
focal-plane optical diagnostic rather than a pyLAST camera plot.

At the LACT ROOT/pyLAST boundary, pyLAST converts focal-plane hit coordinates
to source-offset camera coordinates with
`pix_x=-x_m` and `pix_y=-y_m`. The minimal-CSV plotting helper applies the same
conversion before calling `plot_camera_image`; pixel-id order is unchanged.

Run them from the `LACT_sim` directory:

```bash
./build/run_optical_sim configs/examples/photon_csv_minimal_optics.cfg
./build/run_optical_sim configs/examples/photon_csv_minimal_pe_camera.cfg

# Pure optics: whiteboard plus raw photon-count camera, no display reversal.
python python/plot_minimal_photon_csv_outputs.py \
  --mode optics \
  --hits run_logs/examples/photon_csv_minimal/whiteboard_hits.csv \
  --photon-pixels run_logs/examples/photon_csv_minimal/camera_photon_counts.csv \
  --camera configs/cameras/new_camera_pixels.csv \
  --output-dir run_logs/examples/photon_csv_minimal/plots

# Full camera: one expected-PE image in the pyLAST sky-view convention.
python python/plot_minimal_photon_csv_outputs.py \
  --mode camera \
  --pe-pixels run_logs/examples/photon_csv_minimal/camera_expected_pe.csv \
  --camera configs/cameras/new_camera_pixels.csv \
  --output-dir run_logs/examples/photon_csv_minimal/plots
```

## Making a minimal CSV from CORSIKA photons

The helper

```text
python/corsika_photon_csv_to_minimal.py
```

selects one telescope and removes all columns except position and direction.
It deliberately keeps the original CORSIKA North-West-Up coordinates. This
avoids duplicating the coordinate conversion already implemented and tested in
the simulator.

Example:

```bash
python python/corsika_photon_csv_to_minimal.py \
  event1909_photons.csv \
  configs/sources/event1909_tel19_minimal_photons.csv \
  --telescope-id 19
```

For this CORSIKA 2-D example, `z_m=0` is already a valid upstream reference
plane and does not need to be reconstructed. If another producer supplies
points downstream of the mirror, that producer should move each point along
its own ray to an upstream plane before writing the minimal CSV.

The supplied event-1909 example deliberately treats every CORSIKA bunch record
as one representative photon row. It is therefore suitable for checking image
shape and coordinate transformations, but its total photon count is not the
original multiplicity-weighted CORSIKA photon count.
