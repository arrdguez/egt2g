"""Metropolis a temperatura finita para el grado de libertad de espín.

Corrige el defecto principal del prototipo que hizo Gemini: allí solo se
aceptaban movimientos "cuesta abajo" (0 K). Aquí el espín se PROPONE y se acepta
con la probabilidad de Boltzmann min(1, exp(-ΔE/(k_B T))), recuperando el
ensemble térmico del código original `docs/prueba/main.cpp`.
"""

from __future__ import annotations

import numpy as np

from .config import Config
from .constants import KB_CM1
from .model import energy_gap


def sample_hs_fraction(
    r: float,
    T: float,
    cfg: Config,
    steps: int,
    rng: np.random.Generator,
    start_hs: bool = True,
) -> float:
    """Muestrea el espín de un centro por Metropolis y devuelve la fracción de HS.

    El sistema es de dos niveles: σ = 1 (HS) con energía E_HS y σ = 0 (LS) con
    E_LS. Al voltear el espín, ΔE = E_nuevo - E_viejo. Se acepta si ΔE ≤ 0 o si
    rand < exp(-ΔE/(k_B T)). En T = 0 solo se aceptan los movimientos que bajan
    la energía (función escalón, sin ensanchamiento térmico).
    """
    gap = float(energy_gap(r, cfg))  # E_HS - E_LS (constante a r fija)
    beta = np.inf if T <= 0 else 1.0 / (KB_CM1 * T)

    sigma = 1 if start_hs else 0
    accumulated = 0.0
    for _ in range(steps):
        # ΔE al pasar σ → 1-σ: desde HS es E_LS-E_HS = -gap; desde LS es +gap.
        delta = -gap if sigma == 1 else gap
        if delta <= 0 or rng.random() < np.exp(-delta * beta):
            sigma = 1 - sigma
        accumulated += sigma
    return accumulated / steps
