# Example 02 — Resampling angular resolution with `resample`

This example shows how to resample a `Ldt` object to a coarser angular resolution using `eulumdat-analysis`.

---

## Setup

```bash
pip install eulumdat-analysis
```

---

## Basic usage — raw measurement to 15°×5°

A TechnoTeam RiGO801 goniometer typically measures at 2.5°×1°.
`resample` converts to the standard 15°×5° grid used for archiving and simulation:

```python
from pyldt import LdtReader
from ldt_analysis import resample

ldt = LdtReader.read("TB238123-C-1-1_Plvk_raw.ldt")  # 144 C-planes × 181 γ-angles

ldt_15x5 = resample(ldt)  # default: c_step=15.0, g_step=5.0
# → 24 C-planes × 37 γ-angles

print(ldt_15x5.header.mc)        # 24
print(ldt_15x5.header.ng)        # 37
print(ldt_15x5.header.c_angles)  # [0.0, 15.0, 30.0, ..., 345.0]
print(ldt_15x5.header.g_angles)  # [0.0, 5.0, 10.0, ..., 180.0]
```

The source `ldt` object is never modified.

---

## Custom target resolution

```python
ldt_10x2 = resample(ldt, c_step=10.0, g_step=2.0)
# → 36 C-planes × 91 γ-angles
```

---

## Guard: target finer than source returns None

`resample` refuses to invent angular data points. If the requested resolution
is finer than the source on any dimension, it returns `None`:

```python
ldt_src = LdtReader.read("coarse_15x5.ldt")  # 15°×5° source

result = resample(ldt_src, c_step=5.0, g_step=2.0)
if result is None:
    print("Target resolution is finer than source — resampling refused")
```

---

## Return values

| Case | Return value |
|------|-------------|
| Success | `Ldt` — new object at target resolution |
| `c_step <= 0` or `g_step <= 0` | `None` |
| Source has fewer than 2 C-planes or γ-angles | `None` |
| Target finer than source (`c_step < source_c_step` or `g_step < source_g_step`) | `None` |

---

## What is preserved

All non-angular header fields are copied unchanged:

- `luminaire_name`, `company`, `luminaire_number`
- `lamp_flux`, `lamp_cct`, `lamp_cri`, `lamp_watt`
- `isym`, `lorl`, `dff`, `conv_factor`
- Physical dimensions (`length`, `width`, `height`, etc.)

Only `mc`, `ng`, `dc`, `dg`, `c_angles`, `g_angles` are updated.

---

## Algorithm note

Linear interpolation (`numpy.interp`) in two sequential passes:

1. **γ pass** — for each source C-plane, interpolate onto the target γ grid
2. **C pass** — for each target γ angle, interpolate across C with circular continuity at 360°

Linear interpolation is used instead of cubic to avoid overshoot between physically measured intensity values — a requirement for ISO 17025 traceability.
