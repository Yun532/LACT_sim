# CORSIKA ROOT Photon-Bunch Adapter

The file:

```text
/path/to/photon_E500_th0_run000001_corsika.root
```

contains a ROOT `bunch` tree that is already close to LACT_sim `PhotonCsv`.

Important branches:

```text
bunch_x, bunch_y   arrival position at telescope plane [m]
cx, cy             direction cosines
time               arrival time [ns]
p_height           emission height-like CORSIKA quantity, not used by raytrace yet
lambda             wavelength [nm], 0 when unspecified
nbunch             photons represented by this bunch
itel               telescope id inside the ROOT conversion
runid              event/array identifier from the ROOT conversion
```

Convert it with:

```bash
python3 python/root_bunch_to_photon_csv.py \
  /path/to/photon_E500_th0_run000001_corsika.root \
  run_logs/root_bunch_adapter/photon_E500_th0_photons.csv
```

By default:

```text
x_m = bunch_x
y_m = bunch_y
z_m = 0
dir_x = cx
dir_y = cy
dir_z = -sqrt(1 - cx^2 - cy^2)
multiplicity = nbunch
event_id = runid
telescope_id = itel
```

If an efficiency model needs actual wavelengths, note that this ROOT file has
`lambda=0`; use `--replace-zero-wavelength-nm 400` only for controlled tests,
not as a physical claim.

The validation run is recorded in:

```text
run_logs/root_bunch_adapter/README.md
```
