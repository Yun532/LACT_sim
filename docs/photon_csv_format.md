# LACT_sim PhotonCsv Format

`PhotonCsv` is the common photon-table format used before the optical ray
tracer. External sources such as CORSIKA/EventIO adapters should convert their
output into this table first.

## Required Columns

```text
x_m,y_m,z_m,dir_x,dir_y,dir_z
```

- `x_m,y_m,z_m`: photon position in meters.
- `dir_x,dir_y,dir_z`: photon direction vector. It does not have to be exactly
  normalized; the C++ loader normalizes it before tracing.

## Optional Columns

```text
time_ns,wavelength_nm,weight,multiplicity,event_id,telescope_id,emission_altitude_km
```

If an optional column is absent, the source config supplies the default value:

```ini
time_ns=0
wavelength_nm=400
photon_weight=1
multiplicity=1
event_id=0
telescope_id=0
```

`weight * multiplicity` is applied to the photon before optical propagation.
`emission_altitude_km` is used by MODTRAN tau atmosphere tables; omit it when
atmosphere modeling is disabled or when the atmosphere config supplies an
explicit default emission altitude.

## Coordinate Convention

The source config controls the coordinate frame:

```ini
local_telescope_frame=true
```

means the photon rows are interpreted in each telescope's local optical frame
and then transformed by that telescope's position and pointing.

```ini
local_telescope_frame=false
```

means the photon rows are already expressed in the generic LACT global frame
used by `buildTelescopeFrame()`: azimuth is measured from global `+x` toward
global `+y`. This is not the CORSIKA NWU frame.

Raw CORSIKA/EventIO photon bunches should be read with `run_corsika_trace`.
Do not feed raw EventIO-derived rows into `run_optical_sim` as generic global
PhotonCsv unless they have first been explicitly converted into this LACT frame.

## Array Distribution

For multi-telescope input, include `telescope_id` in the CSV. During an array
run, `python/run_array_sim.py` writes:

```ini
source.filter_telescope_id=<current telescope_id>
```

so each telescope traces only its own rows.

For event-level slicing, `run_array_sim.py` can also write:

```ini
source.filter_event_id=<current event_id>
```

Run multiple events with:

```bash
python3 python/run_array_sim.py ... --event-ids 1,2,3
```

Each event is written into its own `eventNNN/` subdirectory, and the top-level
`array_run_summary.csv` combines all event/telescope rows.

## Normalizing External Files

Use `python/normalize_photon_csv.py` when an external file has different column
names or units:

```bash
python3 python/normalize_photon_csv.py external.csv photons.csv \
  --map x_m=x_cm,y_m=y_cm,z_m=z_cm,dir_x=ux,dir_y=uy,dir_z=uz,telescope_id=tel_id,event_id=event \
  --defaults wavelength_nm=400,weight=1,multiplicity=1,time_ns=0 \
  --scale-position 0.01 \
  --normalize-direction \
  --fail-on-nonfinite
```

The output can then be used with:

```ini
mode=PhotonCsv
csv_path=...
local_telescope_frame=true
```
