"""Configuration of the illustrative, zero-temperature model."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, fields
from pathlib import Path


@dataclass(frozen=True)
class Config:
    N_Fe: int = 1
    E_scale: float = 15.0  # eV
    P_LS: float = 4.7  # eV; illustrative, not an experimental Racah parameter
    r_ref: float = 2.2  # Å
    ligand_strength: float = 1.0
    geometry: str = "radial"  # radial | six_radial
    initial_branch: str = "HS"  # starting GEOMETRY near HS/LS minimum, not a spin constraint
    K_red: float = 20.0  # eV/Å²
    K_shape: float = 10.0  # eV/Å²; radial shape penalty in six_radial mode
    K_bridge: float = 0.0  # eV/Å²
    L_bridge: float = 5.0  # Å; effective length, not an atomistic geometry
    R_bridge: float = 9.2  # Å
    R_300K: float = 2.30  # Å; independent from HS equilibrium distance
    R_0K: float = 1.96  # Å; independent from LS equilibrium distance
    nominal_start: float = 300.0
    nominal_end: float = 0.0
    delta_T: float = 5.0
    protocol: str = "nominal"  # nominal | distance
    R_start: float = 2.30
    R_end: float = 1.96
    delta_R: float = 0.005  # Å
    delta_r: float = 0.025  # Å, local proposals
    sweeps_fast: int = 20
    sweeps_slow: int = 400
    replicas: int = 4
    seed: int = 2026
    r_min: float = 1.65  # Å
    r_max: float = 2.55  # Å

    def __post_init__(self) -> None:
        if self.N_Fe not in (1, 2):
            raise ValueError("N_Fe must be 1 or 2")
        if self.geometry not in ("radial", "six_radial"):
            raise ValueError("geometry must be 'radial' or 'six_radial'")
        if self.initial_branch not in ("HS", "LS"):
            raise ValueError("initial_branch must be 'HS' or 'LS'")
        for key in ("E_scale", "r_ref", "ligand_strength", "delta_T", "delta_R", "delta_r"):
            if not math.isfinite(getattr(self, key)) or getattr(self, key) <= 0:
                raise ValueError(f"{key} must be finite and positive")
        for key in ("K_red", "K_bridge", "K_shape"):
            if not math.isfinite(getattr(self, key)) or getattr(self, key) < 0:
                raise ValueError(f"{key} must be finite and nonnegative")
        for key in ("P_LS", "L_bridge", "R_bridge", "R_300K", "R_0K", "nominal_start", "nominal_end", "R_start", "R_end", "r_min", "r_max"):
            if not math.isfinite(getattr(self, key)):
                raise ValueError(f"{key} must be finite")
        if not 0 < self.r_min < self.r_max:
            raise ValueError("Require 0 < r_min < r_max")
        if not self.r_min <= self.R_0K <= self.r_max or not self.r_min <= self.R_300K <= self.r_max:
            raise ValueError("R_0K and R_300K must lie within the distance bounds")
        if self.protocol not in ("nominal", "distance"):
            raise ValueError("protocol must be 'nominal' or 'distance'")
        if self.protocol == "nominal" and self.nominal_start <= self.nominal_end:
            raise ValueError("nominal_start must exceed nominal_end")
        if self.protocol == "distance" and not (self.r_min <= self.R_end < self.R_start <= self.r_max):
            raise ValueError("Require r_min <= R_end < R_start <= r_max")
        if self.replicas < 1 or self.sweeps_fast < 1 or self.sweeps_slow < 1:
            raise ValueError("replicas and sweeps must be positive integers")
        if self.seed < 0:
            raise ValueError("seed must be nonnegative")


def load_config(path: str | Path) -> Config:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("JSON configuration must be an object")
    allowed = {field.name for field in fields(Config)}
    unknown = set(data) - allowed
    if unknown:
        raise ValueError(f"Unknown configuration keys: {', '.join(sorted(unknown))}")
    return Config(**data)
