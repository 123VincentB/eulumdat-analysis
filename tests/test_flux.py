import math
import pathlib
import pytest
from pyldt import Ldt, LdtHeader, LdtReader
from ldt_analysis import (
    luminous_flux, luminous_flux_range,
    lorl_computed, dff_computed, check_photometric_consistency,
)

DATA_DIR = pathlib.Path(__file__).parent.parent / "data" / "input"
SAMPLE_FILES = sorted(DATA_DIR.glob("*.ldt"))


def make_ldt(c_angles, g_angles, intensities, *, lorl=80.0, dff=100.0,
             num_lamps=1, lamp_flux=3500.0):
    header = LdtHeader(
        c_angles=list(c_angles),
        g_angles=list(g_angles),
        mc=len(c_angles),
        ng=len(g_angles),
        lorl=lorl,
        dff=dff,
        n_sets=1,
        num_lamps=[num_lamps],
        lamp_flux=[lamp_flux],
    )
    return Ldt(header=header, intensities=intensities)


# ---------------------------------------------------------------------------
# Smoke — all 4 functions run on every real LDT file without exception
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", SAMPLE_FILES, ids=[p.name for p in SAMPLE_FILES])
def test_smoke(path):
    ldt = LdtReader.read(path)
    luminous_flux(ldt)
    lorl_computed(ldt)
    dff_computed(ldt)
    check_photometric_consistency(ldt)


# ---------------------------------------------------------------------------
# luminous_flux
# ---------------------------------------------------------------------------

def test_isotropic_sphere():
    """Uniform sphere (I = 1000/4π cd/klm): flux must equal 1000 lm/klm exactly."""
    i_val = 1000.0 / (4.0 * math.pi)
    g_angles = [float(g) for g in range(0, 185, 5)]
    c_angles = [0.0]
    intensities = [[i_val] * len(g_angles)]
    ldt = make_ldt(c_angles, g_angles, intensities)
    phi = luminous_flux(ldt)
    assert phi is not None
    assert abs(phi - 1000.0) < 1e-9


def test_luminous_flux_empty_g_angles():
    ldt = make_ldt([0.0], [], [])
    assert luminous_flux(ldt) is None


def test_luminous_flux_empty_intensities():
    ldt = make_ldt([], [0.0, 45.0, 90.0], [])
    assert luminous_flux(ldt) is None


def test_luminous_flux_isym1_equivalent():
    """1 C-plane vs 24 identical C-planes must yield the same flux (ISYM=1 transparency)."""
    g_angles = [float(g) for g in range(0, 95, 5)]
    i_row = [1000.0 * math.cos(math.radians(g)) for g in g_angles]
    ldt_1 = make_ldt([0.0], g_angles, [list(i_row)])
    ldt_24 = make_ldt(
        [float(c) for c in range(0, 360, 15)],
        g_angles,
        [list(i_row) for _ in range(24)],
    )
    phi_1 = luminous_flux(ldt_1)
    phi_24 = luminous_flux(ldt_24)
    assert phi_1 is not None and phi_24 is not None
    assert abs(phi_1 - phi_24) < 1e-9


# ---------------------------------------------------------------------------
# lorl_computed — bounds on real files
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", SAMPLE_FILES, ids=[p.name for p in SAMPLE_FILES])
def test_lorl_in_bounds(path):
    """LORL must be non-negative. Values slightly above 100% are acceptable in real
    files due to measurement uncertainty and rounded lamp flux declarations."""
    ldt = LdtReader.read(path)
    lorl = lorl_computed(ldt)
    if lorl is not None:
        assert 0.0 <= lorl <= 110.0


# ---------------------------------------------------------------------------
# dff_computed — bounds and correctness
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", SAMPLE_FILES, ids=[p.name for p in SAMPLE_FILES])
def test_dff_in_bounds(path):
    ldt = LdtReader.read(path)
    dff = dff_computed(ldt)
    if dff is not None:
        assert 0.0 <= dff <= 100.0


def test_dff_downward_only():
    """g_angles stops at 90°: all flux is downward, DFF must be 100%."""
    g_angles = [float(g) for g in range(0, 95, 5)]  # 0° to 90°
    c_angles = [0.0, 90.0]
    intensities = [[500.0] * len(g_angles) for _ in range(2)]
    ldt = make_ldt(c_angles, g_angles, intensities, dff=100.0)
    assert abs(dff_computed(ldt) - 100.0) < 1e-9


def test_dff_symmetric_luminaire():
    """Uniform sphere 0°–180°: DFF must equal 50%."""
    g_angles = [float(g) for g in range(0, 185, 5)]
    c_angles = [0.0]
    intensities = [[100.0] * len(g_angles)]
    ldt = make_ldt(c_angles, g_angles, intensities)
    assert abs(dff_computed(ldt) - 50.0) < 1e-9


def test_dff_straddling_zone():
    """Zone straddling 90° exactly (g=[0,85,95,180], uniform I): DFF must be 50%."""
    g_angles = [0.0, 85.0, 95.0, 180.0]
    c_angles = [0.0]
    intensities = [[100.0, 100.0, 100.0, 100.0]]
    ldt = make_ldt(c_angles, g_angles, intensities)
    result = dff_computed(ldt)
    assert result is not None
    assert abs(result - 50.0) < 1e-9


def test_dff_zero_flux():
    """Zero intensities: DFF must be None (division by zero avoided)."""
    g_angles = [0.0, 45.0, 90.0, 135.0, 180.0]
    c_angles = [0.0]
    intensities = [[0.0] * len(g_angles)]
    ldt = make_ldt(c_angles, g_angles, intensities)
    assert dff_computed(ldt) is None


# ---------------------------------------------------------------------------
# check_photometric_consistency
# ---------------------------------------------------------------------------

def test_consistency_dict_keys():
    g_angles = [float(g) for g in range(0, 185, 5)]
    c_angles = [0.0, 90.0, 180.0, 270.0]
    intensities = [[100.0] * len(g_angles) for _ in range(4)]
    ldt = make_ldt(c_angles, g_angles, intensities, lorl=50.0, dff=50.0)
    result = check_photometric_consistency(ldt)
    assert result is not None
    for key in ("lorl_header", "lorl_computed", "lorl_delta",
                "dff_header", "dff_computed", "dff_delta",
                "flux_lm_klm", "flux_lm_abs"):
        assert key in result


def test_consistency_delta_sign():
    """lorl_delta = lorl_computed - lorl_header (signed)."""
    g_angles = [float(g) for g in range(0, 185, 5)]
    c_angles = [0.0]
    intensities = [[1000.0 / (4.0 * math.pi)] * len(g_angles)]
    ldt = make_ldt(c_angles, g_angles, intensities, lorl=90.0)
    result = check_photometric_consistency(ldt)
    assert result is not None
    assert abs(result["lorl_delta"] - (result["lorl_computed"] - result["lorl_header"])) < 1e-12


def test_consistency_num_lamps_minus1():
    """num_lamps=-1 (unspecified): dict is returned but flux_lm_abs is None."""
    g_angles = [float(g) for g in range(0, 185, 5)]
    c_angles = [0.0]
    intensities = [[100.0] * len(g_angles)]
    ldt = make_ldt(c_angles, g_angles, intensities, num_lamps=-1)
    result = check_photometric_consistency(ldt)
    assert result is not None
    assert result["flux_lm_abs"] is None


def test_consistency_lamp_flux_zero():
    """lamp_flux=0.0: flux_lm_abs is None."""
    g_angles = [float(g) for g in range(0, 185, 5)]
    c_angles = [0.0]
    intensities = [[100.0] * len(g_angles)]
    ldt = make_ldt(c_angles, g_angles, intensities, lamp_flux=0.0)
    result = check_photometric_consistency(ldt)
    assert result is not None
    assert result["flux_lm_abs"] is None


def test_consistency_flux_lm_abs_value():
    """flux_lm_abs = flux_lm_klm * num_lamps * lamp_flux / 1000."""
    g_angles = [float(g) for g in range(0, 185, 5)]
    c_angles = [0.0]
    intensities = [[100.0] * len(g_angles)]
    ldt = make_ldt(c_angles, g_angles, intensities, num_lamps=2, lamp_flux=4000.0)
    result = check_photometric_consistency(ldt)
    assert result is not None
    assert result["flux_lm_abs"] is not None
    expected = result["flux_lm_klm"] * 2 * 4000.0 / 1000.0
    assert abs(result["flux_lm_abs"] - expected) < 1e-9


@pytest.mark.parametrize("path", SAMPLE_FILES, ids=[p.name for p in SAMPLE_FILES])
def test_lorl_delta_within_tolerance(path):
    """For real files with declared LORL > 0, |lorl_delta| < 5 pp."""
    ldt = LdtReader.read(path)
    result = check_photometric_consistency(ldt)
    if result is None or result["lorl_header"] <= 0:
        return
    assert abs(result["lorl_delta"]) < 5.0


# ---------------------------------------------------------------------------
# luminous_flux_range — invalid parameters
# ---------------------------------------------------------------------------

def test_flux_range_invalid_reversed():
    ldt = make_ldt([0.0], [0.0, 45.0, 90.0], [[100.0, 100.0, 100.0]])
    assert luminous_flux_range(ldt, 90.0, 0.0) is None


def test_flux_range_invalid_equal():
    ldt = make_ldt([0.0], [0.0, 45.0, 90.0], [[100.0, 100.0, 100.0]])
    assert luminous_flux_range(ldt, 60.0, 60.0) is None


def test_flux_range_invalid_g_min_negative():
    ldt = make_ldt([0.0], [0.0, 45.0, 90.0], [[100.0, 100.0, 100.0]])
    assert luminous_flux_range(ldt, -5.0, 60.0) is None


def test_flux_range_invalid_g_max_over_180():
    ldt = make_ldt([0.0], [0.0, 45.0, 90.0], [[100.0, 100.0, 100.0]])
    assert luminous_flux_range(ldt, 0.0, 200.0) is None


# ---------------------------------------------------------------------------
# luminous_flux_range — window outside gamma grid
# ---------------------------------------------------------------------------

def test_flux_range_window_above_grid():
    """g_max <= first g_angle: no overlap, must return 0.0."""
    ldt = make_ldt([0.0], [100.0, 110.0, 120.0], [[100.0, 100.0, 100.0]])
    assert luminous_flux_range(ldt, 0.0, 90.0) == 0.0


def test_flux_range_window_below_grid():
    """g_min >= last g_angle: no overlap, must return 0.0."""
    ldt = make_ldt([0.0], [0.0, 10.0, 20.0], [[100.0, 100.0, 100.0]])
    assert luminous_flux_range(ldt, 30.0, 80.0) == 0.0


# ---------------------------------------------------------------------------
# luminous_flux_range — coherence with luminous_flux on real files
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", SAMPLE_FILES, ids=[p.name for p in SAMPLE_FILES])
def test_flux_range_full_equals_total(path):
    """luminous_flux_range(ldt, 0, 180) must equal luminous_flux(ldt) to 1e-6."""
    ldt = LdtReader.read(path)
    phi_total = luminous_flux(ldt)
    phi_range = luminous_flux_range(ldt, 0.0, 180.0)
    assert phi_total is not None and phi_range is not None
    assert abs(phi_range - phi_total) < 1e-6


# ---------------------------------------------------------------------------
# luminous_flux_range — coherence with dff_computed on real files
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", SAMPLE_FILES, ids=[p.name for p in SAMPLE_FILES])
def test_flux_range_dff_coherence(path):
    """luminous_flux_range(0, 90) / luminous_flux * 100 must equal dff_computed to 1e-6."""
    ldt = LdtReader.read(path)
    phi_total = luminous_flux(ldt)
    phi_down = luminous_flux_range(ldt, 0.0, 90.0)
    dff = dff_computed(ldt)
    if dff is None:
        return
    assert phi_total is not None and phi_down is not None
    assert abs(phi_down / phi_total * 100.0 - dff) < 1e-6


# ---------------------------------------------------------------------------
# luminous_flux_range — monotonicity on real files
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", SAMPLE_FILES, ids=[p.name for p in SAMPLE_FILES])
def test_flux_range_monotone(path):
    """flux(0-45) <= flux(0-90) <= flux(0-180) on every real file."""
    ldt = LdtReader.read(path)
    phi_45 = luminous_flux_range(ldt, 0.0, 45.0)
    phi_90 = luminous_flux_range(ldt, 0.0, 90.0)
    phi_180 = luminous_flux_range(ldt, 0.0, 180.0)
    assert phi_45 is not None and phi_90 is not None and phi_180 is not None
    assert phi_45 <= phi_90 + 1e-12
    assert phi_90 <= phi_180 + 1e-12


# ---------------------------------------------------------------------------
# luminous_flux_range — additivity
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", SAMPLE_FILES, ids=[p.name for p in SAMPLE_FILES])
def test_flux_range_additivity(path):
    """flux(0-45) + flux(45-90) ≈ flux(0-90) to 1e-4 (interpolation at boundary)."""
    ldt = LdtReader.read(path)
    phi_a = luminous_flux_range(ldt, 0.0, 45.0)
    phi_b = luminous_flux_range(ldt, 45.0, 90.0)
    phi_ab = luminous_flux_range(ldt, 0.0, 90.0)
    assert phi_a is not None and phi_b is not None and phi_ab is not None
    assert abs(phi_a + phi_b - phi_ab) < 1e-4


# ---------------------------------------------------------------------------
# dff_computed — non-regression: new implementation must match old algorithm
# ---------------------------------------------------------------------------

def _dff_old(ldt: Ldt) -> float | None:
    """Old dff_computed algorithm, preserved verbatim for non-regression."""
    h = ldt.header
    mc = len(ldt.intensities)
    ng = len(h.g_angles)
    if mc == 0 or ng == 0:
        return None
    i_mean = [sum(ldt.intensities[c][g] for c in range(mc)) / mc for g in range(ng)]
    g_rad = [math.radians(g) for g in h.g_angles]
    g90 = math.pi / 2.0
    phi_total = 0.0
    phi_down = 0.0
    for k in range(len(g_rad) - 1):
        i_zone = (i_mean[k] + i_mean[k + 1]) / 2.0
        omega = 2.0 * math.pi * abs(math.cos(g_rad[k]) - math.cos(g_rad[k + 1]))
        phi_total += i_zone * omega
        gk = h.g_angles[k]
        gk1 = h.g_angles[k + 1]
        if gk1 <= 90.0:
            phi_down += i_zone * omega
        elif gk < 90.0:
            omega_down = 2.0 * math.pi * abs(math.cos(g_rad[k]) - math.cos(g90))
            phi_down += i_zone * omega_down
    if phi_total == 0.0:
        return None
    return phi_down / phi_total * 100.0


@pytest.mark.parametrize("path", SAMPLE_FILES, ids=[p.name for p in SAMPLE_FILES])
def test_dff_refactoring_noop(path):
    """New dff_computed (via luminous_flux_range) must match the old algorithm to 1e-6 %."""
    ldt = LdtReader.read(path)
    old = _dff_old(ldt)
    new = dff_computed(ldt)
    if old is None:
        assert new is None
    else:
        assert new is not None
        assert abs(new - old) < 1e-6
