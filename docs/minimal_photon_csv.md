# Minimal Photon CSV

The minimal user input is a CSV with six columns:

```csv
x_m,y_m,z_m,dir_x,dir_y,dir_z
3.9014025878906251,0.16619310379028321,0,-0.33789920806884766,-0.014721360988914967,-0.94106716376520094
```

A small file is included at
[`configs/sources/photon_csv_six_column_example.csv`](../configs/sources/photon_csv_six_column_example.csv).
When `multiplicity` is absent, `source.multiplicity` supplies the value; the
included examples set it to 1, so each row represents one photon.

- `x_m,y_m,z_m` are the photon position in metres.
- `dir_x,dir_y,dir_z` are its propagation direction.
- `source.coordinate_frame` defines the axes. Use `telescope_local` for
  hand-written optical tests and `corsika_nwu_relative` for the supplied
  CORSIKA-derived example.

Wavelength and time are optional. A pure optical test can use ideal efficiency
files, while a monochromatic camera test can set, for example,
`source.wavelength_nm=400`. Time is unnecessary when waveform processing is
disabled.

## Included event-1909 example

The plotting input is
[`configs/sources/event1909_tel19_minimal_photons.csv`](../configs/sources/event1909_tel19_minimal_photons.csv).
It contains only the six required columns.

Pure optics produces a whiteboard image and a raw photon-count camera image:

```bash
build/run_optical_sim configs/examples/photon_csv_minimal_optics.cfg
python3 python/plot_minimal_photon_csv_outputs.py \
  --mode optics \
  --hits run_logs/examples/photon_csv_minimal/whiteboard_hits.csv \
  --photon-pixels run_logs/examples/photon_csv_minimal/camera_photon_counts.csv \
  --camera configs/cameras/new_camera_pixels.csv \
  --output-dir run_logs/examples/photon_csv_minimal/plots
```

The full camera example uses the normal `run_corsika_trace` event/camera path,
writes LACT ROOT, and is read back through `pylast.io.LactEventSource`:

```bash
build/run_corsika_trace configs/examples/photon_csv_full_camera_root.cfg
python3 python/plot_photon_csv_root_pylast.py \
  run_logs/examples/photon_csv_full_camera/lact_events.root \
  --event-id 1909 --telescope-id 19 \
  --output run_logs/examples/photon_csv_full_camera/camera_pe.png
```

This cfg assigns every row a fixed wavelength of 400 nm, uses discrete p.e.
sampling, and disables NSB, trigger, and waveform processing. Those detector
features can still be enabled with the same keys used by an EventIO run.
EventIO metadata is optional for PhotonCsv; unavailable shower-truth fields in
the ROOT file remain at their documented default or non-finite values.

The optical plots retain the native focal-plane orientation. At the ROOT/pyLAST
boundary, `LactEventSource` applies the normal source-offset display convention:

```text
pyLAST pix_x = -LACT_sim camera_x
pyLAST pix_y = -LACT_sim camera_y
```

## Creating the six-column input

`python/corsika_photon_csv_to_minimal.py` selects one telescope and removes
optional columns:

```bash
python3 python/corsika_photon_csv_to_minimal.py \
  event1909_photons.csv \
  configs/sources/event1909_tel19_minimal_photons.csv \
  --telescope-id 19
```

The helper retains CORSIKA North-West-Up coordinates, so use
`source.coordinate_frame=corsika_nwu_relative`. The supplied file treats every
CORSIKA bunch record as one photon row; it is intended for image-shape and
coordinate checks, not for preserving the original bunch-weighted intensity.
