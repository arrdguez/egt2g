"""Referencia analítica (cuasiestática) del modelo de dos niveles.

No hace falta minimización numérica: para un centro a distancia r y temperatura
T, la fracción de HS es la solución exacta del ensemble canónico,
n_HS = 1 / (1 + exp(ΔE/(k_B T))).
"""

from __future__ import annotations

from .config import Config
from .model import crossing_distance, delta_o, energy_gap, hs_fraction, ligand_field


def equilibrium_scan(steps: list[dict], cfg: Config) -> list[dict]:
    """Evalúa la fracción analítica de HS para cada paso del barrido."""
    rows = []
    for step in steps:
        if cfg.protocol == "distance":
            r, T = step["r_A"], cfg.T_fixed
        else:
            r, T = cfg.r_fixed, step["T_K"]
        rows.append(
            {
                **step,
                "r_A": r,
                "T_K": T,
                "Dq_cm1": float(ligand_field(r, cfg)),
                "delta_o_cm1": float(delta_o(r, cfg)),
                "gap_cm1": float(energy_gap(r, cfg)),
                "hs_fraction": float(hs_fraction(r, T, cfg)),
            }
        )
    return rows


def crossing_info(cfg: Config) -> dict:
    """Resumen del cruce HS/LS (útil para reportes y figuras)."""
    r_c = crossing_distance(cfg)
    if r_c is None:
        return {"r_c_A": None, "Dq_c_cm1": None, "delta_o_c_cm1": None}
    return {
        "r_c_A": r_c,
        "Dq_c_cm1": float(ligand_field(r_c, cfg)),
        "delta_o_c_cm1": float(delta_o(r_c, cfg)),
    }
