# eulumdat-analysis — Flux integration over an arbitrary gamma window

## Purpose

`luminous_flux_range` computes the luminous flux emitted within an arbitrary
angular cone defined by two gamma boundaries. Use it to quantify the fraction
of flux in a beam lobe, a UGR zone, or any photometric window — without the
limitation of the fixed 0°–180° integral used by `luminous_flux`.

---

## Basic usage

```python
from pyldt import LdtReader
from ldt_analysis import luminous_flux, luminous_flux_range, dff_computed

ldt = LdtReader.read("data/input/sample_14.ldt")

# Total flux (full sphere)
phi_total = luminous_flux(ldt)                        # e.g. 943.7 lm/klm

# Flux in the 0°–60° cone (main beam of a downlight)
phi_cone = luminous_flux_range(ldt, 0.0, 60.0)       # e.g. 756.2 lm/klm

# Fraction of total flux within the cone
fraction = phi_cone / phi_total * 100                 # e.g. 80.1 %

# Lower hemisphere (equivalent to dff_computed but in lm/klm)
phi_down = luminous_flux_range(ldt, 0.0, 90.0)

# Arbitrary intermediate zone
phi_zone = luminous_flux_range(ldt, 30.0, 75.0)
```

---

## Parameters

| Parameter | Type    | Description |
|-----------|---------|-------------|
| `ldt`     | `Ldt`   | EULUMDAT data object with a full expanded intensity matrix |
| `g_min`   | `float` | Lower gamma boundary in degrees. Must be in [0.0, 180.0] |
| `g_max`   | `float` | Upper gamma boundary in degrees. Must be in [0.0, 180.0] and strictly greater than `g_min` |

---

## Return value

| Case | Return value |
|------|-------------|
| Window overlaps the gamma grid | `float` — flux in lm/klm over [g_min, g_max] |
| Window does not overlap the gamma grid | `0.0` — not an error; the file simply does not cover that range |
| `g_min >= g_max` | `None` |
| `g_min < 0` or `g_max > 180` | `None` |
| Intensity matrix or gamma grid empty | `None` |

---

## Edge cases

- **Boundary not on a grid angle** — the intensity at the boundary is estimated
  by linear interpolation between the two adjacent grid angles. The trapezoidal
  zone is then computed only over the partial interval.
- **`luminous_flux_range(ldt, 0.0, 180.0)`** is always equal to `luminous_flux(ldt)`
  (identity by construction).
- **Additivity** — `flux_range(0, 45) + flux_range(45, 90)` ≈ `flux_range(0, 90)`,
  with a small second-order interpolation error when 45° is not a grid angle.
- **`dff_computed`** is internally implemented as
  `luminous_flux_range(ldt, 0.0, 90.0) / luminous_flux(ldt) * 100`.

---

## Notes

The algorithm is strictly equivalent to the CIE 190 trapezoidal method used in
`luminous_flux`, restricted to the requested window. All intensities are in cd/klm
(the EULUMDAT normalisation), so the result is in lm/klm. To convert to absolute
lumens, multiply by `num_lamps[0] * lamp_flux[0] / 1000`.
