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
from ldt_analysis import half_angle

ldt = LdtReader.read("luminaire.ldt")

result = half_angle(ldt, [0.0, 90.0, 180.0, 270.0])
print(result)
# {0.0: 35.4, 90.0: 36.1, 180.0: 35.8, 270.0: 36.0}
```

The function returns the **absolute gamma angle** (degrees from nadir) where intensity drops to 50 % of the C-plane maximum, searching from `γ_max` toward 90°.

---

## Examples

| File | Description |
|------|-------------|
| [`examples/01_basic_usage.md`](https://github.com/123VincentB/eulumdat-analysis/blob/main/examples/01_basic_usage.md) | Basic usage, return values, FWHM |

---

## Project structure

```
eulumdat-analysis/
├── src/
│   └── ldt_analysis/
│       ├── __init__.py
│       └── half_angle.py
├── examples/
│   └── 01_basic_usage.md
├── tests/
│   └── test_half_angle.py
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

## License

MIT — see [LICENSE](https://github.com/123VincentB/eulumdat-analysis/blob/main/LICENSE).

---

## Context

This package was developed as a practical tool in the context of ISO 17025 accredited photometric testing. It is shared as open-source in the hope that it will be useful to others working with EULUMDAT files in Python.
