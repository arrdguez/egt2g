"""Configuración del modelo de campo de ligandos (energías en cm⁻¹, distancias en Å).

El único grado de libertad microscópico es la distancia Fe–N `r`. El campo de
ligandos Dq(r) es lo que provoca, en su raíz, la transición HS↔LS. La
temperatura actúa como control externo de los barridos (no entra en el
Hamiltoniano: solo modifica la población térmica vía Metropolis).
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, fields
from pathlib import Path


@dataclass(frozen=True)
class Config:
    # --- Modelo de campo de ligandos (un solo centro Fe) ---
    Dq0: float = 1500.0  # campo de ligandos en r0 (cm⁻¹); Δo = 10·Dq0
    r0: float = 2.0  # distancia Fe–N de referencia (Å)
    P: float = 20000.0  # energía de apareamiento (cm⁻¹)

    # --- Cota de la distancia Fe–N ---
    r_min: float = 1.80  # Å
    r_max: float = 2.40  # Å

    # --- Protocolo de barrido ---
    protocol: str = "distance"  # "distance" | "temperature"

    # Barrido de distancia (a T_fixed constante)
    R_start: float = 2.40  # Å (campo débil → HS)
    R_end: float = 1.80  # Å (campo fuerte → LS)
    delta_R: float = 0.005  # Å (paso)
    T_fixed: float = 300.0  # K, temperatura del barrido de distancia

    # Barrido de temperatura (a r_fixed constante)
    r_fixed: float = 2.17  # Å, distancia fijada para el barrido térmico
    T_start: float = 400.0  # K
    T_end: float = 0.0  # K
    delta_T: float = 4.0  # K (paso)

    # --- Monte Carlo ---
    equilibration: int = 2000  # pasos de termalización (descartados)
    mc_steps: int = 20000  # pasos de medición por réplica
    replicas: int = 4  # réplicas independientes (promediadas)
    seed: int = 2026  # semilla reproducible

    def __post_init__(self) -> None:
        if self.protocol not in ("distance", "temperature"):
            raise ValueError("protocol debe ser 'distance' o 'temperature'")
        for key in ("Dq0", "r0", "P", "delta_R", "delta_T"):
            if not math.isfinite(getattr(self, key)) or getattr(self, key) <= 0:
                raise ValueError(f"{key} debe ser finito y positivo")
        for key in ("T_fixed", "T_start", "T_end", "R_start", "R_end"):
            if not math.isfinite(getattr(self, key)) or getattr(self, key) < 0:
                raise ValueError(f"{key} debe ser finito y no negativo")
        if not 0 < self.r_min < self.r_max:
            raise ValueError("Se requiere 0 < r_min < r_max")
        for key in ("r_fixed", "R_start", "R_end"):
            if not self.r_min <= getattr(self, key) <= self.r_max:
                raise ValueError(f"{key} debe estar dentro de [r_min, r_max]")
        if self.equilibration < 1 or self.mc_steps < 1 or self.replicas < 1:
            raise ValueError("equilibration, mc_steps y replicas deben ser >= 1")
        if self.seed < 0:
            raise ValueError("seed debe ser no negativo")


def load_config(path: str | Path) -> Config:
    """Carga la configuración desde un JSON y rechaza claves desconocidas."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("La configuración JSON debe ser un objeto")
    allowed = {field.name for field in fields(Config)}
    unknown = set(data) - allowed
    if unknown:
        raise ValueError(f"Claves desconocidas: {', '.join(sorted(unknown))}")
    return Config(**data)
