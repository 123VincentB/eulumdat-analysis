# Example 01 — Basic usage of `half_angle`

This example shows how to compute the **half-angle at half maximum** for a parsed EULUMDAT file using `eulumdat-analysis`.

---

## Setup

```bash
pip install eulumdat-analysis
```

---

## Read an LDT file and compute half-angles

```python
from pyldt import LdtReader
from ldt_analysis import half_angle

# Parse the LDT file
ldt = LdtReader.read("luminaire.ldt")

# Compute half-angle for the four main C-planes
result = half_angle(ldt, [0.0, 90.0, 180.0, 270.0])

print(result)
# {0.0: 35.4, 90.0: 36.1, 180.0: 35.8, 270.0: 36.0}
```

The function returns the **absolute gamma angle** (in degrees, from nadir) where the intensity first drops to 50 % of the C-plane maximum, searching from `γ_max` toward 90°.

---

## Compute for all C-planes in the file

```python
from pyldt import LdtReader
from ldt_analysis import half_angle

ldt = LdtReader.read("luminaire.ldt")

# Use every C-plane declared in the file
result = half_angle(ldt, ldt.header.c_angles)

for c_plane, angle in result.items():
    if angle is not None:
        print(f"C={c_plane:6.1f}°  →  half-angle = {angle:.2f}°")
    else:
        print(f"C={c_plane:6.1f}°  →  undefined (multi-peak or I_max = 0)")
```

---

## Return value

| Case | Return value |
|------|-------------|
| Normal beam | `float` — crossing angle in degrees |
| C-plane not found in file (±0.01°) | `None` |
| `I_max = 0` (dark or inactive plane) | `None` |
| Intensity never drops to half-max within [γ_max, 90°] | `None` |
| Multi-peak distribution (secondary peak prominence > 5 % of I_max) | `None` |

---

## Note on ISYM = 1 (full rotational symmetry)

For luminaires with full rotational symmetry (`ldt.header.isym == 1`), only one
C-plane is stored in the file. `half_angle` automatically applies that single
C-plane to any requested angle:

```python
ldt = LdtReader.read("rotationally_symmetric.ldt")  # ISYM=1, one C-plane

result = half_angle(ldt, [0.0, 45.0, 90.0, 180.0])
# All four planes return the same value
# {0.0: 35.2, 45.0: 35.2, 90.0: 35.2, 180.0: 35.2}
```

---

## FWHM from two complementary C-planes

`half_angle` is a **unilateral** measurement (one side only). To compute the
Full Width at Half Maximum (FWHM) of a beam using two complementary C-planes:

```python
ha_c0   = result[0.0]    # crossing angle on the C=0° side
ha_c180 = result[180.0]  # crossing angle on the C=180° side

if ha_c0 is not None and ha_c180 is not None:
    fwhm = ha_c0 + ha_c180
    print(f"FWHM (C0/C180 plane) = {fwhm:.1f}°")
```

For a rotationally symmetric luminaire (ISYM=1):

```python
ha = result[0.0]
if ha is not None:
    fwhm = 2 * ha
    print(f"FWHM = {fwhm:.1f}°")
```
