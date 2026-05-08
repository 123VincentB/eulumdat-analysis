import math
import pathlib
import pytest
from pyldt import Ldt, LdtHeader, LdtReader
from ldt_analysis import luminous_flux, lorl_computed, dff_computed, check_photometric_consistency

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
