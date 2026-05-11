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
    """Return total luminous flux in lm/klm by trapezoidal integration (CIE 190 method).

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


def luminous_flux_range(ldt: Ldt, g_min: float, g_max: float) -> float | None:
    """Return luminous flux integrated over the gamma window [g_min, g_max], in lm/klm.

    Uses the same trapezoidal method as luminous_flux() (CIE 190).
    Zones partially overlapping the window boundaries are split by linear
    interpolation of the mean intensity at the boundary angle.

    Parameters
    ----------
    ldt : Ldt
        EULUMDAT data object (full expanded matrix expected).
    g_min : float
        Lower gamma boundary in degrees. Must be in [0, 180].
    g_max : float
        Upper gamma boundary in degrees. Must be in [0, 180] and > g_min.

    Returns
    -------
    float
        Flux in lm/klm over [g_min, g_max].
    float (0.0)
        If the window does not overlap the file's gamma grid at all.
    None
        If g_min >= g_max, parameters are out of [0, 180], or the matrix is empty.
    """
    if g_min < 0.0 or g_max > 180.0 or g_min >= g_max:
        return None
    h = ldt.header
    i_mean = _mean_intensities(ldt)
    if i_mean is None:
        return None
    g_angles = h.g_angles
    ng = len(g_angles)
    if g_max <= g_angles[0] or g_min >= g_angles[-1]:
        return 0.0
    phi = 0.0
    for k in range(ng - 1):
        gk = g_angles[k]
        gk1 = g_angles[k + 1]
        if gk1 <= g_min or gk >= g_max:
            continue
        lo = max(gk, g_min)
        hi = min(gk1, g_max)
        span = gk1 - gk
        i_lo = i_mean[k] if lo == gk else i_mean[k] + (i_mean[k + 1] - i_mean[k]) * (lo - gk) / span
        i_hi = i_mean[k + 1] if hi == gk1 else i_mean[k] + (i_mean[k + 1] - i_mean[k]) * (hi - gk) / span
        i_zone = (i_lo + i_hi) / 2.0
        omega = 2.0 * math.pi * abs(math.cos(math.radians(lo)) - math.cos(math.radians(hi)))
        phi += i_zone * omega
    return phi


def dff_computed(ldt: Ldt) -> float | None:
    """Return downward flux fraction computed from the intensity matrix, in %.

    Downward hemisphere: gamma in [0°, 90°] (nadir to horizontal).
    Returns None if the total flux is zero or the matrix is empty.
    """
    phi_total = luminous_flux(ldt)
    if phi_total is None or phi_total == 0.0:
        return None
    phi_down = luminous_flux_range(ldt, 0.0, 90.0)
    if phi_down is None:
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
