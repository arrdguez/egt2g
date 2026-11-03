"""Conduce el Monte Carlo a temperatura finita sobre los pasos del barrido.

Para cada paso del protocolo se termaliza (pasos descartados) y se mide la
fracción de HS con Metropolis (espín propuesto y aceptado con Boltzmann),
promediando sobre réplicas independientes.
"""

from __future__ import annotations

import numpy as np

from .config import Config
from .metropolis import sample_hs_fraction


def run_mc(steps: list[dict], cfg: Config) -> list[dict]:
    """Devuelve una fila por paso con la fracción de HS promediada entre réplicas."""
    rows = []
    for step_index, step in enumerate(steps):
        if cfg.protocol == "distance":
            r, T = step["r_A"], cfg.T_fixed
        else:
            r, T = cfg.r_fixed, step["T_K"]

        # Cada paso usa un SeedSequence propio, para que los pasos no compartan
        # la misma secuencia aleatoria (reproducible: misma semilla → mismo resultado).
        streams = np.random.SeedSequence([cfg.seed, step_index]).spawn(cfg.replicas)

        fractions = []
        for stream in streams:
            rng = np.random.default_rng(stream)
            sample_hs_fraction(r, T, cfg, cfg.equilibration, rng)  # termalización
            fractions.append(sample_hs_fraction(r, T, cfg, cfg.mc_steps, rng))

        rows.append({**step, "hs_fraction": float(np.mean(fractions))})
    return rows
