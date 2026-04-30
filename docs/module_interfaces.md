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
