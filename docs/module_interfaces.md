# LACT_sim Module Interfaces

This document records the intended boundaries for modules that are not yet full
detector implementations.

## Event Identity

Formal output uses one event key:

```text
event_id
```

For CORSIKA/EventIO runs with `source.event_id_mode=event_array100`, LACT_sim
sets:

```text
event_id = shower_event * 100 + array_id
```

The HDF5 output stores only `event_id` in `/events/table` and `/images/index`.
Intermediate logs may still print decoded shower/array information for human
debugging, but downstream analysis should use `event_id`.

For example, original CORSIKA shower `468898` with `array_id=0` is stored as
`event_id=46889800`; with `array_id=2` it is stored as `event_id=46889802`.
Those trailing two digits are the array/core-offset stream, not an extra shower
number.

## Atmosphere

CORSIKA photons are interpreted as photons already at the telescope plane. The
atmosphere module is therefore an optional extra transmission factor, not a
second full shower propagation.

Supported first-stage configuration:

```ini
atmosphere.transmission=none
atmosphere.transmission=0.92
atmosphere.transmission=configs/atmosphere/transmission.csv
```

CSV tables use:

```text
wavelength_nm,transmission
```

For formal production, this CSV should come from the same atmospheric profile
used in the CORSIKA run, or from a clearly documented site measurement/model.
Until that source is fixed, the recommended benchmark setting is
`atmosphere.transmission=none` to avoid adding uncalibrated physics.

Atmospheric refractive-index timing is reserved but disabled by default:

```ini
atmosphere.refractive_index_model=none
atmosphere.speed_of_light_correction=false
```

Current optical-path timing uses:

```ini
propagation.speed_of_light_m_per_ns=0.299792458
```

## SiPM And Electronics

The current SiPM/electronics path is an interface layer:

```text
collector photon -> pe conversion -> integrated pixel image
```

The active setting is:

```ini
electronics.pe_conversion=none
electronics.pe_conversion=0.35
electronics.pe_conversion=configs/efficiency/sipm_pde.csv
```

Detailed waveform, saturation, crosstalk, afterpulse, dark count, gain
fluctuation, and trigger electronics should be implemented behind this boundary
later. The HDF5 format already keeps integrated `pe`, `signal`, `photon_count`,
and timing summaries per event/telescope image.

For each collected Cherenkov photon or photon bunch, the integrated weight is
handled multiplicatively:

```text
p.e. contribution =
  CORSIKA bunch multiplicity
  * optical photon weight
  * mirror/filter/atmosphere efficiency factors
  * light-collector acceptance and reflection weight
  * electronics.pe_conversion
```

`electronics.pe_conversion` can be disabled, a constant, or a wavelength table.
The current output is an integrated p.e. image. It is not yet a waveform and it
does not model SiPM saturation, crosstalk, afterpulse, dark count, or gain
fluctuations.

## NSB

The first NSB implementation is a constant-rate Poisson model:

```ini
nsb.enabled=false
nsb.model=constant_rate
nsb.rate_pe_per_ns_per_pixel=0.05
nsb.window_ns=16
nsb.seed=12345
```

It applies only to dense pixel-camera HDF5 output. The final
`/images/dense/pe` includes the NSB p.e. contribution. Use
`output.hdf5_write_components=true` to additionally write `cherenkov_pe` and
`nsb_pe`.

## Trigger

The first trigger implementation is intentionally simple:

```ini
trigger.enabled=false
trigger.pixel_threshold_pe=5
trigger.camera_multiplicity=3
trigger.array_multiplicity=2
trigger.coincidence_window_ns=50
```

Camera trigger counts pixels above threshold in one telescope image. Array
trigger counts triggered telescopes for the same `event_id`. Neighbor topology,
waveform thresholding, and hardware-specific trigger logic are reserved for the
future electronics implementation.
