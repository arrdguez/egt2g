"""Genera los pasos del barrido (control externo).

Dos controles físicos del mismo fenómeno:
  - "distance":    variar la distancia Fe–N → variar Dq(r) → cruce HS/LS
                   (el "fenómeno de raíz" del código `prueba`).
  - "temperature": variar la temperatura a distancia fija → población térmica
                   (cruce suave, no cooperativo, de un centro aislado).
"""

from __future__ import annotations

import numpy as np

from .config import Config


def _descending(start: float, end: float, increment: float) -> list[float]:
    """Serie estrictamente decreciente desde `start` hasta `end` (inclusive).

    Se redondea a 10 decimales para eliminar la deriva de punto flotante de
    `np.arange` con paso no entero (p. ej. 2.3800000000000003 → 2.38).
    """
    values = list(
        np.round(
            np.arange(start, end - 1e-9 * abs(increment), -increment, dtype=float), 10
        )
    )
    if not values or not np.isclose(values[-1], end, atol=1e-6):
        values.append(float(end))
    return values


def build_steps(cfg: Config) -> list[dict]:
    """Devuelve la lista de pasos según el protocolo configurado."""
    if cfg.protocol == "distance":
        radii = _descending(cfg.R_start, cfg.R_end, cfg.delta_R)
        return [{"step": i, "r_A": r} for i, r in enumerate(radii)]
    temperatures = _descending(cfg.T_start, cfg.T_end, cfg.delta_T)
    return [{"step": i, "T_K": t} for i, t in enumerate(temperatures)]
