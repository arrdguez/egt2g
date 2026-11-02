"""Externally imposed geometric contraction; nominal T is NOT thermodynamic T."""

from __future__ import annotations

import numpy as np

from .config import Config


def _descending(start: float, end: float, increment: float) -> list[float]:
    values = list(np.arange(start, end, -increment, dtype=float))
    if not values or not np.isclose(values[-1], end):
        values.append(float(end))
    return values


def cooling_steps(cfg: Config) -> list[dict]:
    if cfg.protocol == "nominal":
        values = _descending(cfg.nominal_start, cfg.nominal_end, cfg.delta_T)
        return [{"step": i, "nominal_T_K": t,
                 "R_red_A": cfg.R_0K + (cfg.R_300K - cfg.R_0K) * t / 300.0}
                for i, t in enumerate(values)]
    values = _descending(cfg.R_start, cfg.R_end, cfg.delta_R)
    return [{"step": i, "nominal_T_K": float("nan"), "R_red_A": r}
            for i, r in enumerate(values)]
