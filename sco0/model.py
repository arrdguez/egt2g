"""Modelo de campo de ligandos: Dq(r) provoca el cruce HS↔LS.

El espíritu del modelo (heredado de `docs/prueba/main.cpp`):

    variar la distancia Fe–N (r)  ⟶  cambia el campo de ligandos Dq(r)  ⟶  cruza HS/LS.

La temperatura NO entra en estas energías; solo modifica la población térmica
a través del muestreo de Metropolis (ver `metropolis.py`).

Energías (linealización del diagrama de Tanabe–Sugano para d⁶):

    E_HS(r) = -4·Dq(r)          (t₂g⁴ e_g²)
    E_LS(r) = -24·Dq(r) + P     (t₂g⁶, con energía de apareamiento P)

El campo de ligandos decae con la distancia (modelo de cargas puntuales,
10·Dq ∝ r⁻⁵): Dq(r) = Dq0·(r0/r)⁵.

Unidades: energía en cm⁻¹, distancia en Å.
"""

from __future__ import annotations

import numpy as np

from .config import Config
from .constants import KB_CM1


def ligand_field(r: float | np.ndarray, cfg: Config) -> float | np.ndarray:
    """Campo de ligandos Dq(r) = Dq0·(r0/r)⁵ (cm⁻¹)."""
    r = np.asarray(r, dtype=float)
    return cfg.Dq0 * (cfg.r0 / r) ** 5


def delta_o(r: float | np.ndarray, cfg: Config) -> float | np.ndarray:
    """Desdoblamiento octaédrico Δo(r) = 10·Dq(r) (cm⁻¹)."""
    return 10.0 * ligand_field(r, cfg)


def electronic_energy(
    r: float | np.ndarray, hs: bool, cfg: Config
) -> float | np.ndarray:
    """Energía electrónica de un centro: -4·Dq (HS) o -24·Dq + P (LS)."""
    dq = ligand_field(r, cfg)
    return -4.0 * dq if hs else -24.0 * dq + cfg.P


def energy_gap(r: float | np.ndarray, cfg: Config) -> float | np.ndarray:
    """ΔE = E_HS - E_LS = 20·Dq(r) - P. Positivo ⇒ LS es más estable."""
    return 20.0 * ligand_field(r, cfg) - cfg.P


def electronic_derivative(
    r: float | np.ndarray, hs: bool, cfg: Config
) -> float | np.ndarray:
    """dE/dr. Como d(Dq)/dr = -5·Dq/r ⇒ dE_HS/dr = 20·Dq/r y dE_LS/dr = 120·Dq/r."""
    r = np.asarray(r, dtype=float)
    dq = ligand_field(r, cfg)
    return (20.0 * dq / r) if hs else (120.0 * dq / r)


def crossing_distance(cfg: Config) -> float | None:
    """Distancia de cruce r_c donde E_HS = E_LS ⇒ Dq = P/20.

    r_c = r0·(20·Dq0/P)^(1/5). Devuelve None si cae fuera de [r_min, r_max].
    """
    if cfg.P <= 0:
        return None
    r_c = cfg.r0 * (20.0 * cfg.Dq0 / cfg.P) ** 0.2
    return float(r_c) if cfg.r_min <= r_c <= cfg.r_max else None


def hs_fraction(r: float | np.ndarray, T: float, cfg: Config) -> float | np.ndarray:
    """Fracción térmica de HS (solución analítica del sistema de dos niveles).

    n_HS = 1 / (1 + exp(ΔE/(k_B T))). En T = 0 es la función escalón
    (HS si ΔE ≤ 0, LS si ΔE > 0).
    """
    if T <= 0:
        return np.where(energy_gap(r, cfg) <= 0, 1.0, 0.0)
    return 1.0 / (1.0 + np.exp(energy_gap(r, cfg) / (KB_CM1 * T)))
