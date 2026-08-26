# Photon response modes

`run_corsika_trace` supports two detector-response modes:

- `response.mode=expectation` traces one ray per EventIO bunch and
  carries the bunch multiplicity and efficiencies as a continuous weight.
- `response.mode=stochastic_pe` (default) expands the bunch into represented photons,
  samples a wavelength and random stream for each one, and writes discrete
  photoelectrons. A fractional final photon is handled as a Bernoulli trial.

Set `response.seed` to make stochastic output reproducible. Stochastic mode
requires the source photon weight to be one; bunch multiplicity remains the
only represented-photon count.

For time-resolved PhotonCsv backgrounds, repeated photons from one weighted
row share the same input geometry and time. Use independently sampled rows
with `weight=1,multiplicity=1` when waveform timing, microcell occupancy, or
trigger correlations matter. Weighted rows remain appropriate for
`expectation` images over arbitrarily long exposures.

For speed, wavelength-only factors (atmosphere, mirror/filter curves and PDE)
are sampled before ray tracing. This is statistically equivalent for the
camera p.e. result because those factors do not depend on the ray geometry.
Facet reflectivity scale, funnel acceptance and collector loss remain after
geometry. Consequently, stochastic-mode trace counters such as `hit_mirror`
describe the thinned rays that were actually traced, while `input_photons`
still describes the original bunch multiplicity.
