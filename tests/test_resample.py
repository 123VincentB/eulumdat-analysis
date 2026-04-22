"""
Tests unitaires pour ldt_analysis.resample.

Stratégie : objets Ldt synthétiques (pas de fichier fixture),
en cohérence avec test_half_angle.py.

Tests couverts :
    - Forme de sortie (mc, ng, c_angles, g_angles)
    - Conservation des champs header non angulaires
    - Source non modifiée
    - Conservation des valeurs sur une matrice uniforme
    - Conservation aux angles communs source/cible (±tolérance interpolation)
    - Continuité circulaire C (raccordement 345°→0°)
    - Guards : step <= 0, résolution plus fine, source < 2 plans
    - Conservation approchée du flux intégré
"""
import dataclasses

import numpy as np
import pytest
from pyldt import Ldt, LdtHeader

from ldt_analysis import resample


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_ldt(
    c_step: float = 2.5,
    g_step: float = 1.0,
    intensity_value: float = 100.0,
    mc: int | None = None,
    ng: int | None = None,
) -> Ldt:
    """
    Ldt synthétique à distribution uniforme.
    mc et ng sont déduits de c_step/g_step si non fournis.
    """
    if mc is None:
        mc = round(360.0 / c_step)      # ex: 2.5° → 144
    if ng is None:
        ng = round(180.0 / g_step) + 1  # ex: 1.0° → 181

    c_angles = [i * c_step for i in range(mc)]
    g_angles = [i * g_step for i in range(ng)]

    header = LdtHeader(
        company        = "TestCo",
        luminaire_name = "Test luminaire",
        luminaire_number = "TEST-001",
        ityp  = 2,
        isym  = 0,
        mc    = mc,
        dc    = c_step,
        ng    = ng,
        dg    = g_step,
        lorl  = 85.0,
        lamp_flux = [5000.0],
        c_angles  = c_angles,
        g_angles  = g_angles,
    )
    intensities = [[intensity_value] * ng for _ in range(mc)]
    return Ldt(header=header, intensities=intensities)


def make_ldt_gradient(c_step: float = 2.5, g_step: float = 1.0) -> Ldt:
    """Ldt avec intensité variant linéairement sur γ (0 à 100 cd/klm)."""
    mc = round(360.0 / c_step)
    ng = round(180.0 / g_step) + 1
    c_angles = [i * c_step for i in range(mc)]
    g_angles = [i * g_step for i in range(ng)]
    header = LdtHeader(
        mc=mc, dc=c_step, ng=ng, dg=g_step,
        c_angles=c_angles, g_angles=g_angles,
    )
    row = [100.0 * (1.0 - g / 180.0) for g in g_angles]
    intensities = [row[:] for _ in range(mc)]
    return Ldt(header=header, intensities=intensities)


# ---------------------------------------------------------------------------
# Forme de sortie
# ---------------------------------------------------------------------------

class TestOutputShape:
    def test_mc_24_for_c_step_15(self):
        r = resample(make_ldt(), c_step=15.0, g_step=5.0)
        assert r is not None
        assert r.header.mc == 24
        assert len(r.header.c_angles) == 24
        assert len(r.intensities) == 24

    def test_ng_37_for_g_step_5(self):
        r = resample(make_ldt(), c_step=15.0, g_step=5.0)
        assert r is not None
        assert r.header.ng == 37
        assert len(r.header.g_angles) == 37
        assert all(len(row) == 37 for row in r.intensities)

    def test_c_angles_values(self):
        r = resample(make_ldt(), c_step=15.0, g_step=5.0)
        expected = list(np.arange(0.0, 360.0, 15.0))
        assert r.header.c_angles == pytest.approx(expected)

    def test_g_angles_values(self):
        r = resample(make_ldt(), c_step=15.0, g_step=5.0)
        expected = list(np.arange(0.0, 185.0, 5.0))
        assert r.header.g_angles == pytest.approx(expected)

    def test_dc_dg_updated_in_header(self):
        r = resample(make_ldt(), c_step=15.0, g_step=5.0)
        assert r.header.dc == pytest.approx(15.0)
        assert r.header.dg == pytest.approx(5.0)

    def test_custom_step_10x2(self):
        r = resample(make_ldt(), c_step=10.0, g_step=2.0)
        # c: 0,10,...,350 → 36 ; g: 0,2,...,180 → 91
        assert r.header.mc == 36
        assert r.header.ng == 91


# ---------------------------------------------------------------------------
# Conservation des champs header non angulaires
# ---------------------------------------------------------------------------

class TestHeaderPreservation:
    def test_luminaire_name_preserved(self):
        ldt = make_ldt()
        r = resample(ldt, c_step=15.0, g_step=5.0)
        assert r.header.luminaire_name == ldt.header.luminaire_name

    def test_company_preserved(self):
        ldt = make_ldt()
        r = resample(ldt, c_step=15.0, g_step=5.0)
        assert r.header.company == ldt.header.company

    def test_lamp_flux_preserved(self):
        ldt = make_ldt()
        r = resample(ldt, c_step=15.0, g_step=5.0)
        assert r.header.lamp_flux == ldt.header.lamp_flux

    def test_isym_preserved(self):
        ldt = make_ldt()
        r = resample(ldt, c_step=15.0, g_step=5.0)
        assert r.header.isym == ldt.header.isym

    def test_lorl_preserved(self):
        ldt = make_ldt()
        r = resample(ldt, c_step=15.0, g_step=5.0)
        assert r.header.lorl == pytest.approx(ldt.header.lorl)


# ---------------------------------------------------------------------------
# Source non modifiée
# ---------------------------------------------------------------------------

class TestSourceNotModified:
    def test_c_angles_not_mutated(self):
        ldt = make_ldt()
        c_orig = ldt.header.c_angles[:]
        resample(ldt, c_step=15.0, g_step=5.0)
        assert ldt.header.c_angles == c_orig

    def test_g_angles_not_mutated(self):
        ldt = make_ldt()
        g_orig = ldt.header.g_angles[:]
        resample(ldt, c_step=15.0, g_step=5.0)
        assert ldt.header.g_angles == g_orig

    def test_intensities_not_mutated(self):
        ldt = make_ldt(intensity_value=77.0)
        resample(ldt, c_step=15.0, g_step=5.0)
        assert ldt.intensities[0][0] == pytest.approx(77.0)

    def test_mc_not_mutated(self):
        ldt = make_ldt()
        mc_orig = ldt.header.mc
        resample(ldt, c_step=15.0, g_step=5.0)
        assert ldt.header.mc == mc_orig


# ---------------------------------------------------------------------------
# Valeurs interpolées
# ---------------------------------------------------------------------------

class TestInterpolatedValues:
    def test_uniform_matrix_preserved(self):
        r = resample(make_ldt(intensity_value=123.0), c_step=15.0, g_step=5.0)
        assert all(
            pytest.approx(v, abs=0.01) == 123.0
            for row in r.intensities for v in row
        )

    def test_values_within_source_range(self):
        ldt = make_ldt_gradient()
        r = resample(ldt, c_step=15.0, g_step=5.0)
        src_min = min(v for row in ldt.intensities for v in row)
        src_max = max(v for row in ldt.intensities for v in row)
        for row in r.intensities:
            for v in row:
                assert src_min - 0.01 <= v <= src_max + 0.01

    def test_common_grid_point_preserved(self):
        """Aux angles communs source/cible, la valeur doit être conservée (±0.5%)."""
        ldt = make_ldt_gradient()
        r = resample(ldt, c_step=15.0, g_step=5.0)
        # G=0° est commun (source 0° → valeur 100.0)
        assert r.intensities[0][0] == pytest.approx(100.0, rel=0.005)
        # G=180° est commun (source 180° → valeur 0.0)
        assert r.intensities[0][-1] == pytest.approx(0.0, abs=0.5)

    def test_g90_value_correct(self):
        """G=90° sur matrice gradient : valeur attendue = 50.0 (mi-chemin 0–180)."""
        ldt = make_ldt_gradient()
        r = resample(ldt, c_step=15.0, g_step=5.0)
        idx_g90 = r.header.g_angles.index(90.0)
        assert r.intensities[0][idx_g90] == pytest.approx(50.0, rel=0.005)


# ---------------------------------------------------------------------------
# Cohérence aux frontières C
# ---------------------------------------------------------------------------

class TestCBoundary:
    def test_c0_and_c345_values_on_asymmetric_distribution(self):
        """
        Source asymétrique : C=0°=100, C=357.5°=200, autres=150.
        Après rééchantillonnage, C=0° et C=345° cibles doivent avoir
        des valeurs distinctes reflétant la distribution source.
        """
        mc = 144
        ng = 181
        c_step_src, g_step_src = 2.5, 1.0
        c_angles = [i * c_step_src for i in range(mc)]
        g_angles = [i * g_step_src for i in range(ng)]
        intensities = [[150.0] * ng for _ in range(mc)]
        intensities[0]  = [100.0] * ng   # C=0°
        intensities[-1] = [200.0] * ng   # C=357.5°
        header = LdtHeader(mc=mc, dc=c_step_src, ng=ng, dg=g_step_src,
                           c_angles=c_angles, g_angles=g_angles)
        ldt = Ldt(header=header, intensities=intensities)

        r = resample(ldt, c_step=15.0, g_step=5.0)
        assert r is not None

        # C=0° doit rester à 100 (point source exact)
        assert r.intensities[0][0] == pytest.approx(100.0, abs=0.5)
        # C=345° est entre C=345° source (150) et C=357.5° (200) — valeur interpolée
        idx_c345 = r.header.c_angles.index(345.0)
        # 345° source = exactement 150.0 (pas d'interpolation entre plans)
        assert r.intensities[idx_c345][0] == pytest.approx(150.0, abs=0.5)


# ---------------------------------------------------------------------------
# Guards — valeurs invalides
# ---------------------------------------------------------------------------

class TestGuards:
    def test_none_on_c_step_zero(self):
        assert resample(make_ldt(), c_step=0.0, g_step=5.0) is None

    def test_none_on_g_step_zero(self):
        assert resample(make_ldt(), c_step=15.0, g_step=0.0) is None

    def test_none_on_negative_c_step(self):
        assert resample(make_ldt(), c_step=-15.0, g_step=5.0) is None

    def test_none_on_negative_g_step(self):
        assert resample(make_ldt(), c_step=15.0, g_step=-5.0) is None

    def test_none_on_finer_c_step(self):
        # Source à 2.5°, cible à 1.0° → refusé
        assert resample(make_ldt(c_step=2.5), c_step=1.0, g_step=5.0) is None

    def test_none_on_finer_g_step(self):
        # Source à 1.0°, cible à 0.5° → refusé
        assert resample(make_ldt(g_step=1.0), c_step=15.0, g_step=0.5) is None

    def test_equal_step_is_valid(self):
        # Même résolution que la source → retourne un Ldt valide (copie rééchantillonnée)
        ldt = make_ldt(c_step=15.0, g_step=5.0)
        r = resample(ldt, c_step=15.0, g_step=5.0)
        assert r is not None

    def test_none_on_single_c_plane(self):
        header = LdtHeader(mc=1, dc=0.0, ng=37, dg=5.0,
                           c_angles=[0.0], g_angles=list(np.arange(0, 185, 5)))
        ldt = Ldt(header=header, intensities=[[100.0] * 37])
        assert resample(ldt) is None

    def test_none_on_single_g_angle(self):
        header = LdtHeader(mc=24, dc=15.0, ng=1, dg=0.0,
                           c_angles=list(np.arange(0, 360, 15)),
                           g_angles=[0.0])
        ldt = Ldt(header=header, intensities=[[100.0]] * 24)
        assert resample(ldt) is None


# ---------------------------------------------------------------------------
# Conservation approchée du flux intégré
# ---------------------------------------------------------------------------

class TestFluxConservation:
    @staticmethod
    def _zonal_flux(ldt: Ldt) -> float:
        """Flux intégré approché par somme de Riemann."""
        g = np.array(ldt.header.g_angles) * np.pi / 180.0
        c = np.array(ldt.header.c_angles) * np.pi / 180.0
        dg = np.full(len(g), ldt.header.dg * np.pi / 180.0)
        dc = np.full(len(c), ldt.header.dc * np.pi / 180.0)
        I  = np.array(ldt.intensities)
        return float(np.sum(I * np.sin(g)[None, :] * dg[None, :] * dc[:, None]))

    def test_flux_conserved_within_5_percent(self):
        """Tolérance ±5% — limitation connue de l'interpolation linéaire."""
        ldt = make_ldt_gradient()
        r = resample(ldt, c_step=15.0, g_step=5.0)
        flux_src = self._zonal_flux(ldt)
        flux_dst = self._zonal_flux(r)
        assert abs(flux_dst - flux_src) / flux_src < 0.05

    def test_flux_conserved_uniform(self):
        """Sur matrice uniforme, flux conservé à ±1%."""
        ldt = make_ldt(intensity_value=200.0)
        r = resample(ldt, c_step=15.0, g_step=5.0)
        flux_src = self._zonal_flux(ldt)
        flux_dst = self._zonal_flux(r)
        assert abs(flux_dst - flux_src) / flux_src < 0.01
