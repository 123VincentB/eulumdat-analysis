"""
resample.py — Rééchantillonnage d'un objet Ldt vers une résolution angulaire cible.

Algorithme : interpolation linéaire numpy en deux passes séquentielles.
    Passe 1 — sur γ (pour chaque plan C source).
    Passe 2 — sur C (pour chaque γ cible), avec continuité circulaire period=360°.

Choix métrologique : l'interpolation linéaire est conservative (pas d'overshoot)
et traçable lors d'un audit ISO 17025. L'interpolation cubique peut introduire
des dépassements non physiques entre deux valeurs d'intensité lumineuse.
"""
from __future__ import annotations

import dataclasses

import numpy as np
from pyldt import Ldt


def resample(
    ldt: Ldt,
    c_step: float = 15.0,
    g_step: float = 5.0,
) -> Ldt | None:
    """
    Rééchantillonne un objet ``Ldt`` vers la résolution angulaire cible.

    L'objet source n'est jamais modifié. Tous les champs du header (nom,
    flux, ISYM, dimensions, etc.) sont conservés ; seuls ``mc``, ``ng``,
    ``dc``, ``dg``, ``c_angles`` et ``g_angles`` sont mis à jour.

    Paramètres
    ----------
    ldt :
        Objet ``Ldt`` source (produit par ``pyldt.LdtReader`` ou construit
        en mémoire). Doit exposer la matrice complète d'intensités.
    c_step :
        Pas des plans C cibles en degrés. Défaut : 15.0°.
    g_step :
        Pas des angles γ cibles en degrés. Défaut : 5.0°.

    Retourne
    --------
    ``Ldt``
        Nouvel objet avec la résolution cible.
    ``None`` si :
        - ``c_step <= 0`` ou ``g_step <= 0``
        - ``ldt`` contient moins de 2 plans C ou moins de 2 angles γ
        - la résolution cible est plus fine que la source sur au moins
          une dimension (on n'invente pas de points)

    Exemples
    --------
    ::

        from pyldt import LdtReader
        from ldt_analysis import resample

        ldt = LdtReader.read("luminaire_raw.ldt")   # 2.5° × 1°
        ldt_15x5 = resample(ldt)                    # → 24 plans C × 37 γ
        if ldt_15x5 is None:
            raise ValueError("résolution source plus fine que la cible")
    """
    # ── Guards ────────────────────────────────────────────────────────────────
    if c_step <= 0 or g_step <= 0:
        return None

    c_src = ldt.header.c_angles
    g_src = ldt.header.g_angles

    if len(c_src) < 2 or len(g_src) < 2:
        return None

    src_c_step = c_src[1] - c_src[0]
    src_g_step = g_src[1] - g_src[0]

    if c_step < src_c_step or g_step < src_g_step:
        return None

    # ── Axes cibles ───────────────────────────────────────────────────────────
    c_new = list(np.arange(0.0, 360.0, c_step))
    g_new = list(np.arange(0.0, 180.0 + g_step, g_step))

    c_src_arr = np.asarray(c_src, dtype=float)
    g_src_arr = np.asarray(g_src, dtype=float)
    c_new_arr = np.asarray(c_new, dtype=float)
    g_new_arr = np.asarray(g_new, dtype=float)

    intensities_src = np.asarray(ldt.intensities, dtype=float)
    mc_src = len(c_src)
    mc_new = len(c_new_arr)
    ng_new = len(g_new_arr)

    # ── Passe 1 : rééchantillonnage γ pour chaque plan C source ──────────────
    temp = np.zeros((mc_src, ng_new))
    for i in range(mc_src):
        temp[i] = np.interp(g_new_arr, g_src_arr, intensities_src[i])

    # ── Passe 2 : rééchantillonnage C avec continuité circulaire ─────────────
    # period=360 : numpy gère le raccordement C=345°→C=0° sans extension manuelle
    result = np.zeros((mc_new, ng_new))
    for j in range(ng_new):
        result[:, j] = np.interp(c_new_arr, c_src_arr, temp[:, j], period=360.0)

    # ── Construction du nouvel objet Ldt ─────────────────────────────────────
    new_header = dataclasses.replace(
        ldt.header,
        mc=mc_new,
        ng=ng_new,
        dc=c_step,
        dg=g_step,
        c_angles=c_new,
        g_angles=g_new,
    )
    new_intensities = [list(np.round(row, 2)) for row in result]

    return Ldt(header=new_header, intensities=new_intensities)
