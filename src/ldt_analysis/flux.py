import math
from pyldt import Ldt


def _mean_intensities(ldt: Ldt) -> list[float] | None:
    mc = len(ldt.intensities)
    ng = len(ldt.header.g_angles)
    if mc == 0 or ng == 0:
        return None
    return [
        sum(ldt.intensities[c][g] for c in range(mc)) / mc
        for g in range(ng)
    ]


def luminous_flux(ldt: Ldt) -> float | None:
    """Return total luminous flux in lm/klm by trapezoidal integration over the full matrix.

    The result is normalised to the lamp flux used in the file (cd/klm convention).
    Multiply by num_lamps * lamp_flux / 1000 to obtain absolute lumens.
    Returns None if the intensity matrix or gamma grid is empty.
    """
    h = ldt.header
    i_mean = _mean_intensities(ldt)
    if i_mean is None:
        return None
    g_rad = [math.radians(g) for g in h.g_angles]
    phi = 0.0
    for k in range(len(g_rad) - 1):
        i_zone = (i_mean[k] + i_mean[k + 1]) / 2.0
        omega = 2.0 * math.pi * abs(math.cos(g_rad[k]) - math.cos(g_rad[k + 1]))
        phi += i_zone * omega
    return phi


def lorl_computed(ldt: Ldt) -> float | None:
    """Return light output ratio computed from the intensity matrix, in %.

    Since intensities are in cd/klm, the flux integral is in lm/klm.
    1000 lm/klm corresponds to LORL = 100%, so LORL [%] = flux [lm/klm] / 10.
    Returns None if luminous_flux() returns None.
    """
    phi = luminous_flux(ldt)
    if phi is None:
        return None
    return phi / 10.0


def dff_computed(ldt: Ldt) -> float | None:
    """Return downward flux fraction computed from the intensity matrix, in %.

    Downward hemisphere: gamma in [0°, 90°] (nadir to horizontal).
    Zones straddling 90° are split exactly at the boundary using the exact
    solid angle formula: omega_down = 2*pi * |cos(gamma_k) - cos(90°)|.
    Returns None if the total flux is zero or the matrix is empty.
    """
    h = ldt.header
    i_mean = _mean_intensities(ldt)
    if i_mean is None:
        return None
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


def check_photometric_consistency(ldt: Ldt) -> dict | None:
    """Compare matrix-computed LORL and DFF against the header declared values.

    Returns a dict with computed values, header values, and signed deltas
    (positive means computed > declared). Returns None if the intensity matrix
    is empty or any base computation fails.

    flux_lm_abs is None when num_lamps[0] is -1 (unspecified) or lamp_flux[0] is 0.
    """
    h = ldt.header
    phi = luminous_flux(ldt)
    if phi is None:
        return None
    lorl_c = lorl_computed(ldt)
    dff_c = dff_computed(ldt)
    if lorl_c is None or dff_c is None:
        return None

    flux_lm_abs = None
    if h.num_lamps and h.lamp_flux:
        n = h.num_lamps[0]
        f = h.lamp_flux[0]
        if n != -1 and n > 0 and f > 0.0:
            flux_lm_abs = phi * n * f / 1000.0

    return {
        "lorl_header": h.lorl,
        "lorl_computed": lorl_c,
        "lorl_delta": lorl_c - h.lorl,
        "dff_header": h.dff,
        "dff_computed": dff_c,
        "dff_delta": dff_c - h.dff,
        "flux_lm_klm": phi,
        "flux_lm_abs": flux_lm_abs,
    }
