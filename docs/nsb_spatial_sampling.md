# NSB spatial sampling

## Scope

The NSB response is a pre-electronics photoelectron-arrival model. It produces
integer Poisson p.e. counts in camera pixels and time bins. It does not yet
model SiPM pulses, dark counts, crosstalk, afterpulses, gain fluctuations, or
ADC electronics.

The scalar NSB rate still comes from either `constant_rate` or `spectral_flux`.
The spatial model controls how that scalar is distributed over the camera:

```ini
nsb.spatial_model=uniform
```

uses the same mean rate in every pixel. Fluctuations remain statistically
independent between disjoint pixel/time cells.

## Per-pixel rate scale

Use a strict CSV rate map when pixels or telescopes need different means:

```ini
nsb.spatial_model=pixel_scale
nsb.pixel_scale_csv=path/to/nsb_pixel_scale.csv
```

The required columns are:

```csv
telescope_id,pixel_id,relative_scale
-1,0,1.00
-1,1,0.97
3,1,1.12
```

`telescope_id=-1` is the fallback for every telescope. An exact telescope row
overrides the fallback. Every camera pixel must resolve to exactly one scale;
missing rows, duplicate keys, negative/non-finite scales, unknown columns, or
malformed values stop the run with an error. A zero scale disables NSB for that
pixel. The final rate is

```text
pixel_rate = base_rate_pe_per_ns_per_pixel * relative_scale
```

This format deliberately represents a rate map, not a sky-coordinate model.
Star catalogues, pointing-dependent sky maps, and off-axis acceptance belong in
a later spatial-model implementation and do not require changes to the sampler
or output interfaces.

## Canonical realization

For one `(seed, event_id, telescope_id, pixel-rate map, time grid)` tuple, the
sampler returns one deterministic sparse realization plus its integrated image.
The total count follows a Poisson distribution with

```text
expected_total_pe = sum(pixel_rate) * bin_width_ns * n_bins
```

Poisson splitting then gives the requested per-pixel means and zero covariance
between disjoint cells. Repeating the same tuple reproduces the same result;
changing the event or telescope changes the random stream.

ROOT and HDF5 use this same realization for the integrated image, p.e. sequence,
and trigger decision. With `waveform.source=pe`, the serialized image is exactly
the sum over the configured waveform time bins, including both Cherenkov and
NSB p.e. The HDF5 dense and sparse image representations are also equivalent.

## Output metadata

HDF5 `/metadata/nsb` records `spatial_model`, `pixel_scale_csv`, and the parsed
row count in addition to the resolved scalar rate. The run summary prints the
same values so a rate-map run can be audited without opening the CSV manually.
