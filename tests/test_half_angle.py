import math
import pytest
from pyldt import Ldt, LdtHeader
from ldt_analysis import half_angle


def make_ldt(c_angles, g_angles, intensities):
    """Build a minimal Ldt object for testing."""
    header = LdtHeader(
        c_angles=list(c_angles),
        g_angles=list(g_angles),
        mc=len(c_angles),
        ng=len(g_angles),
    )
    return Ldt(header=header, intensities=intensities)


def cosine_intensities(c_angles, g_angles, scale=1000.0):
    """Lambertian distribution: I(gamma) = scale * cos(gamma), restricted to [0,90]."""
    result = []
    for _ in c_angles:
        row = [scale * math.cos(math.radians(g)) if g <= 90.0 else 0.0 for g in g_angles]
        result.append(row)
    return result


# ---------------------------------------------------------------------------
# 1. Standard Lambertian (max at gamma=0°) — analytical result: half-angle = 60°
# ---------------------------------------------------------------------------
def test_lambertian_half_angle():
    g_angles = [float(g) for g in range(0, 185, 5)]  # 0 to 180 step 5
    c_angles = [0.0]
    intensities = cosine_intensities(c_angles, g_angles)
    ldt = make_ldt(c_angles, g_angles, intensities)

    result = half_angle(ldt, [0.0])
    assert result[0.0] is not None
    assert abs(result[0.0] - 60.0) < 0.5  # within 0.5° of analytical


# ---------------------------------------------------------------------------
# 2. Off-axis beam (max at gamma != 0°) — half-angle measured from gamma_max
# ---------------------------------------------------------------------------
def test_offaxis_half_angle():
    # Gaussian centered at 30°: I(gamma) = exp(-((gamma-30)/15)^2)
    # Analytical crossing toward 90°: gamma_root = 30 + 15*sqrt(ln2) ≈ 42.49°
    # half_angle returns the absolute crossing angle, not the distance from gamma_max
    g_angles = [float(g) for g in range(0, 91, 5)]
    c_angles = [0.0]
    intensities = [
        [1000.0 * math.exp(-((g - 30) / 15) ** 2) for g in g_angles]
    ]
    ldt = make_ldt(c_angles, g_angles, intensities)

    result = half_angle(ldt, [0.0])
    assert result[0.0] is not None
    expected = 30.0 + 15.0 * math.sqrt(math.log(2))  # ≈ 42.49°
    assert abs(result[0.0] - expected) < 1.0  # within 1° given 5° sampling


# ---------------------------------------------------------------------------
# 3. Beam that never drops to half-maximum within 90° → None
# ---------------------------------------------------------------------------
def test_beam_never_reaches_half_maximum():
    # Flat intensity of 800 cd/klm — I_half = 400, never crossed
    g_angles = [float(g) for g in range(0, 91, 5)]
    c_angles = [0.0]
    intensities = [[800.0] * len(g_angles)]
    ldt = make_ldt(c_angles, g_angles, intensities)

    result = half_angle(ldt, [0.0])
    assert result[0.0] is None


# ---------------------------------------------------------------------------
# 4. Dark C-plane (I_max = 0) → None
# ---------------------------------------------------------------------------
def test_dark_plane():
    g_angles = [float(g) for g in range(0, 91, 5)]
    c_angles = [0.0]
    intensities = [[0.0] * len(g_angles)]
    ldt = make_ldt(c_angles, g_angles, intensities)

    result = half_angle(ldt, [0.0])
    assert result[0.0] is None


# ---------------------------------------------------------------------------
# 5. Invalid C-plane (not in c_angles) → None
# ---------------------------------------------------------------------------
def test_invalid_c_plane():
    g_angles = [float(g) for g in range(0, 91, 5)]
    c_angles = [0.0]
    intensities = cosine_intensities(c_angles, g_angles)
    ldt = make_ldt(c_angles, g_angles, intensities)

    result = half_angle(ldt, [45.0])
    assert result[45.0] is None


# ---------------------------------------------------------------------------
# 6. Coarse angular resolution (15° steps) — interpolation should give ~60°
# ---------------------------------------------------------------------------
def test_coarse_resolution_lambertian():
    g_angles = [float(g) for g in range(0, 181, 15)]  # 0, 15, 30, ..., 180
    c_angles = [0.0]
    intensities = cosine_intensities(c_angles, g_angles)
    ldt = make_ldt(c_angles, g_angles, intensities)

    result = half_angle(ldt, [0.0])
    assert result[0.0] is not None
    assert abs(result[0.0] - 60.0) < 2.0  # looser tolerance for coarse grid


# ---------------------------------------------------------------------------
# 7. Multiple C-planes returned in a single call
# ---------------------------------------------------------------------------
def test_multiple_c_planes():
    g_angles = [float(g) for g in range(0, 91, 5)]
    c_angles = [0.0, 90.0, 180.0, 270.0]
    intensities = cosine_intensities(c_angles, g_angles)
    ldt = make_ldt(c_angles, g_angles, intensities)

    result = half_angle(ldt, [0.0, 90.0, 999.0])
    assert result[0.0] is not None
    assert result[90.0] is not None
    assert result[999.0] is None
    assert abs(result[0.0] - 60.0) < 0.5
    assert abs(result[90.0] - 60.0) < 0.5
