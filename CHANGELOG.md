# Changelog

## [1.2.0] - 2026-05-08

### Added
- `luminous_flux`: total flux in lm/klm by trapezoidal integration (CIE 190 method).
  Mean over all C-planes per γ angle, exact solid angle formula `2π|cos(γ_k) − cos(γ_{k+1})|`.
  Returns `None` if matrix or γ-grid is empty.
- `lorl_computed`: LORL in % from the intensity matrix (`flux / 10`). Returns `None` if
  `luminous_flux` returns `None`.
- `dff_computed`: downward flux fraction in %, with exact handling of zones straddling γ=90°.
  Returns `None` if total flux is zero.
- `check_photometric_consistency`: compares matrix-computed LORL and DFF against header
  declared values. Returns a dict with computed values, header values, signed deltas, and
  absolute flux (`flux_lm_abs=None` when `num_lamps=-1` or `lamp_flux=0`).

## [1.1.0] - 2026-04-22

### Added
- `resample`: resamples a `Ldt` object to a target angular resolution (`c_step`, `g_step`).
  Sequential linear interpolation — γ pass per C-plane, then C pass with circular continuity
  at 360°. Returns `None` for invalid inputs (step ≤ 0, fewer than 2 planes/angles, or target
  resolution finer than source). Default target: 15°×5°.

## [1.0.0] - 2026-04-07

### Added
- `half_angle`: computes the half-angle at half maximum (HAHM) per C-plane,
  using CubicSpline interpolation and Brent root-finding.
