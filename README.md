# eulumdat-analysis

[![PyPI](https://img.shields.io/pypi/v/eulumdat-analysis)](https://pypi.org/project/eulumdat-analysis/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/eulumdat-analysis)](https://pypi.org/project/eulumdat-analysis/)
[![License: MIT](https://img.shields.io/github/license/123VincentB/eulumdat-analysis)](https://github.com/123VincentB/eulumdat-analysis/blob/main/LICENSE)

Practical photometric calculations on EULUMDAT data.

Part of the `eulumdat-*` Python ecosystem built on [`eulumdat-py`](https://github.com/123VincentB/eulumdat-py).  
Developed in an ISO 17025 accredited photometry laboratory.

---

## Features

- `half_angle` — half-angle at half maximum (HAHM) per C-plane, with CubicSpline interpolation
- `resample` — resample a `Ldt` to a coarser angular resolution (linear interpolation, ISO 17025 traceable)
- `luminous_flux` — total flux in lm/klm by trapezoidal integration (CIE 190 method)
- `lorl_computed` — LORL in % computed from the intensity matrix
- `dff_computed` — downward flux fraction in % computed from the intensity matrix
- `check_photometric_consistency` — compare matrix-computed LORL and DFF against header declared values
- Handles all ISYM symmetry types including full rotational symmetry (ISYM=1)
- Automatically rejects multi-peak distributions (secondary peak prominence > 5 % of I_max)
- Returns `None` for undefined cases — never raises unhandled exceptions

---

## Installation

```bash
pip install eulumdat-analysis
```

> For development:
>
> ```bash
> git clone https://github.com/123VincentB/eulumdat-analysis.git
> cd eulumdat-analysis
> pip install -e ".[dev]"
> ```

---

## Quick start

```python
from pyldt import LdtReader
from ldt_analysis import half_angle, resample, luminous_flux, lorl_computed, dff_computed, check_photometric_consistency

ldt = LdtReader.read("luminaire.ldt")

# Half-angle at half maximum per C-plane
result = half_angle(ldt, [0.0, 90.0, 180.0, 270.0])
print(result)
# {0.0: 35.4, 90.0: 36.1, 180.0: 35.8, 270.0: 36.0}

# Resample a raw 2.5°×1° measurement to standard 15°×5°
ldt_raw = LdtReader.read("luminaire_raw.ldt")   # 144 C-planes × 181 γ-angles
ldt_15x5 = resample(ldt_raw)                    # → 24 C-planes × 37 γ-angles

# Flux and photometric consistency check (CIE 190)
flux = luminous_flux(ldt)                        # e.g. 865.3 lm/klm
lorl = lorl_computed(ldt)                        # e.g. 86.53 %
dff  = dff_computed(ldt)                         # e.g. 99.8 %
report = check_photometric_consistency(ldt)
# {"lorl_header": 87.0, "lorl_computed": 86.53, "lorl_delta": -0.47, ...}
```

---

## Examples

| File | Description |
|------|-------------|
| [`examples/01_basic_usage.md`](https://github.com/123VincentB/eulumdat-analysis/blob/main/examples/01_basic_usage.md) | `half_angle` — basic usage, return values, FWHM |
| [`examples/02_resample.md`](https://github.com/123VincentB/eulumdat-analysis/blob/main/examples/02_resample.md) | `resample` — resolution resampling, guards, preserved fields |

---

## Project structure

```
eulumdat-analysis/
├── src/
│   └── ldt_analysis/
│       ├── __init__.py
│       ├── half_angle.py
│       ├── resample.py
│       └── flux.py
├── examples/
│   ├── 01_basic_usage.md
│   └── 02_resample.md
├── tests/
│   ├── test_half_angle.py
│   ├── test_resample.py
│   └── test_flux.py
├── CHANGELOG.md
├── LICENSE
└── README.md
```

---

## Dependencies

- [`eulumdat-py`](https://pypi.org/project/eulumdat-py/) ≥ 1.0.0 — EULUMDAT parser
- [`scipy`](https://scipy.org/) ≥ 1.7 — CubicSpline interpolation and Brent root-finding

---

## `half_angle` — return values

| Case | Return value |
|------|-------------|
| Normal beam | `float` — crossing angle in degrees |
| C-plane not found in file (±0.01°) | `None` |
| `I_max = 0` (dark or inactive plane) | `None` |
| Intensity never drops to half-max within [γ_max, 90°] | `None` |
| Multi-peak distribution | `None` |

---

## `resample` — return values

| Case | Return value |
|------|-------------|
| Success | `Ldt` — new object at target resolution |
| `c_step <= 0` or `g_step <= 0` | `None` |
| Source has fewer than 2 C-planes or γ-angles | `None` |
| Target finer than source | `None` |

---

## eulumdat-* ecosystem

> **New to the ecosystem?** [eulumdat-quickstart](https://github.com/123VincentB/eulumdat-quickstart) — a step-by-step guide covering all 8 packages with working examples.

| Package | Description |
|---|---|
| [eulumdat-py](https://pypi.org/project/eulumdat-py/) | Read / write EULUMDAT files |
| [eulumdat-symmetry](https://pypi.org/project/eulumdat-symmetry/) | Symmetrise and detect ISYM |
| [eulumdat-plot](https://pypi.org/project/eulumdat-plot/) | Polar intensity diagram (SVG/PNG) |
| [eulumdat-luminance](https://pypi.org/project/eulumdat-luminance/) | Luminance table and polar diagram |
| [eulumdat-ugr](https://pypi.org/project/eulumdat-ugr/) | UGR catalogue (CIE 117/190) |
| **`eulumdat-analysis`** | **Beam half-angle, flux integration, photometric consistency — this package** |
| [eulumdat-report](https://pypi.org/project/eulumdat-report/) | Full photometric datasheet (HTML/PDF) |
| [eulumdat-ies](https://pypi.org/project/eulumdat-ies/) | LDT ↔ IES LM-63-2002 conversion |

---

## License

MIT — see [LICENSE](https://github.com/123VincentB/eulumdat-analysis/blob/main/LICENSE).

---

## Context

This package was developed as a practical tool in the context of ISO 17025 accredited photometric testing. It is shared as open-source in the hope that it will be useful to others working with EULUMDAT files in Python.
