# Changelog

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
