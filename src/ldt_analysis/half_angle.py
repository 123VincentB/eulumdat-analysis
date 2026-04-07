from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from pyldt import Ldt


def half_angle(ldt: Ldt, c_planes: list[float]) -> dict[float, float | None]:
    """
    Compute the half-angle at half maximum (HAHM) for one or more C-planes.

    For each requested C-plane, the half-angle is the absolute gamma angle
    (measured from nadir) where intensity first drops to 50% of the C-plane
    maximum, searching from gamma_max toward 90°. The measurement is
    unilateral: from gamma_max toward 90° only.

    Parameters
    ----------
    ldt : Ldt
        Parsed EULUMDAT object (from pyldt.LdtReader).
    c_planes : list[float]
        C-plane angles in degrees to process. Each must match a value in
        ldt.header.c_angles within ±0.01°. For ISYM=1 (full rotational
        symmetry), any requested C-plane is resolved to the single available
        C-plane.

    Returns
    -------
    dict[float, float | None]
        Maps each requested C-plane angle to its half-angle in degrees,
        or None if undefined.
    """
    return {c: _half_angle_single(ldt, c) for c in c_planes}


def _half_angle_single(ldt: Ldt, c_plane: float) -> float | None:
    try:
        # --- resolve C-plane index ---
        c_index = None
        for i, c in enumerate(ldt.header.c_angles):
            if abs(c - c_plane) <= 0.01:
                c_index = i
                break

        # ISYM=1: full rotational symmetry — all C-planes are identical.
        # pyldt keeps only 1 C-plane in this case; use it for any request.
        if c_index is None and ldt.header.isym == 1 and len(ldt.header.c_angles) == 1:
            c_index = 0

        if c_index is None:
            return None

        # --- restrict gamma to [0°, 90°] ---
        g_all = ldt.header.g_angles
        i_all = ldt.intensities[c_index]

        g_vals = []
        i_vals = []
        for j, g in enumerate(g_all):
            if g <= 90.0:
                g_vals.append(g)
                i_vals.append(i_all[j])

        if len(g_vals) < 2:
            return None

        # --- reject multi-peak distributions ---
        # A secondary peak is real if its prominence toward the global maximum
        # exceeds 5% of I_max. Tiny fluctuations (noise, rounding) are ignored.
        i_max_val = max(i_vals)
        g_max_idx_full = i_vals.index(i_max_val)
        prominence_threshold = 0.05 * i_max_val

        n = len(i_vals)
        significant_peaks = 1  # the global maximum always counts

        def _local_max(k):
            if k == 0:
                return n > 1 and i_vals[0] > i_vals[1]
            if k == n - 1:
                return n > 1 and i_vals[-1] > i_vals[-2]
            return i_vals[k] > i_vals[k - 1] and i_vals[k] > i_vals[k + 1]

        for k in range(n):
            if k == g_max_idx_full:
                continue
            if not _local_max(k):
                continue
            # prominence = peak value minus deepest valley on path to global max
            if k < g_max_idx_full:
                valley = min(i_vals[k:g_max_idx_full + 1])
            else:
                valley = min(i_vals[g_max_idx_full:k + 1])
            if i_vals[k] - valley > prominence_threshold:
                significant_peaks += 1

        if significant_peaks > 1:
            return None

        # --- find gamma_max (smallest gamma in case of ties) ---
        i_max = max(i_vals)
        if i_max == 0:
            return None
        g_max_idx = i_vals.index(i_max)
        g_max = g_vals[g_max_idx]

        i_half = i_max / 2.0

        # --- build search domain: from gamma_max to 90° ---
        search_g = g_vals[g_max_idx:]
        search_i = i_vals[g_max_idx:]

        if len(search_g) < 2:
            return None

        # --- find first bracket where intensity crosses i_half ---
        bracket_left = None
        for k in range(len(search_i) - 1):
            if search_i[k] >= i_half and search_i[k + 1] <= i_half:
                bracket_left = k
                break
        if bracket_left is None:
            return None

        # --- fit CubicSpline on the search domain ---
        cs = CubicSpline(search_g, search_i)

        # --- find root of cs(gamma) - i_half using brentq ---
        g_a = search_g[bracket_left]
        g_b = search_g[bracket_left + 1]
        g_root = brentq(lambda g: cs(g) - i_half, g_a, g_b)

        if g_root > 0 and (g_root == g_root):  # positive and not NaN
            return float(g_root)
        return None

    except Exception:
        return None
