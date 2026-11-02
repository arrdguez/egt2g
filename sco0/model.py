"""Illustrative electronic energy surfaces (eV) and radial forces (eV/Å)."""

from __future__ import annotations

import numpy as np

from .config import Config


def electronic_energy(r: float | np.ndarray, hs: bool, cfg: Config):
    x = np.asarray(r) / cfg.r_ref
    rho = 2 if hs else 12
    energy = cfg.E_scale * (1 / (9 * x**9) - 0.02 * cfg.ligand_strength * rho / x**5 - 1 / x)
    return energy if hs else energy + cfg.P_LS


def electronic_derivative(r: float | np.ndarray, hs: bool, cfg: Config):
    x = np.asarray(r) / cfg.r_ref
    rho = 2 if hs else 12
    return cfg.E_scale / cfg.r_ref * (-1 / x**10 + 0.1 * cfg.ligand_strength * rho / x**6 + 1 / x**2)


def crossing_distance(cfg: Config) -> float | None:
    """Analytic intersection of the unrestrained surfaces, if within bounds."""
    if cfg.P_LS <= 0:
        return None
    r = cfg.r_ref * (0.2 * cfg.E_scale * cfg.ligand_strength / cfg.P_LS) ** 0.2
    return r if cfg.r_min <= r <= cfg.r_max else None


def ligand_delta(r: float | np.ndarray, cfg: Config) -> float:
    """Effective octahedral Δ_o (eV); mean r^-5, NOT (mean r)^-5."""
    x = np.asarray(r, dtype=float) / cfg.r_ref
    return 0.1 * cfg.E_scale * cfg.ligand_strength * float(np.mean(x ** -5))


def ligand_positions(distances: np.ndarray, cfg: Config) -> np.ndarray:
    """Cartesian N coordinates (N_Fe, 6, 3) relative to Fe at the origin."""
    if cfg.geometry != "six_radial":
        raise ValueError("Explicit ligand positions require six_radial geometry")
    radii = _radii(distances, cfg)
    axes = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0],
                     [0, 0, 1], [0, 0, -1]], dtype=float)
    return radii[:, :, None] * axes[None, :, :]


def _radii(distances: np.ndarray, cfg: Config) -> np.ndarray:
    r = np.asarray(distances, dtype=float)
    expected = (cfg.N_Fe,) if cfg.geometry == "radial" else (cfg.N_Fe, 6)
    if r.shape != expected:
        raise ValueError(f"Expected distances with shape {expected}")
    if not np.all(np.isfinite(r)) or np.any(r < cfg.r_min) or np.any(r > cfg.r_max):
        raise ValueError("Fe–N distance outside configured bounds or not finite")
    return r if cfg.geometry == "six_radial" else r[:, None]


def ground_state(distances: np.ndarray, R_red: float, cfg: Config) -> tuple[float, np.ndarray]:
    """Electrons instantly occupy the lower fixed-geometry surface; spin is not sampled."""
    r = _radii(distances, cfg)
    hs = np.array([np.mean(electronic_energy(site, True, cfg)) <=
                   np.mean(electronic_energy(site, False, cfg)) for site in r], dtype=bool)
    return total_energy(distances, hs, R_red, cfg), hs


def total_energy(distances: np.ndarray, spins_hs: np.ndarray, R_red: float, cfg: Config) -> float:
    r = _radii(distances, cfg)
    spins = np.asarray(spins_hs, dtype=bool)
    if spins.shape != (cfg.N_Fe,):
        raise ValueError("Expected one spin for each Fe")
    energy = sum(float(np.mean(electronic_energy(ri, bool(si), cfg))) for ri, si in zip(r, spins))
    energy += 0.5 * cfg.K_red * float(np.sum(np.mean((r - R_red) ** 2, axis=1)))
    if cfg.geometry == "six_radial":
        energy += 0.5 * cfg.K_shape * float(np.sum((r - r.mean(axis=1, keepdims=True)) ** 2))
    if cfg.N_Fe == 2:
        energy += 0.5 * cfg.K_bridge * float((r.mean(axis=1).sum() + cfg.L_bridge - cfg.R_bridge) ** 2)
    return energy


def gradient(distances: np.ndarray, spins_hs: np.ndarray, R_red: float, cfg: Config) -> np.ndarray:
    """Derivative of total energy with respect to the Fe–N distances; force = -gradient."""
    r = _radii(distances, cfg)
    spins = np.asarray(spins_hs, dtype=bool)
    if spins.shape != (cfg.N_Fe,):
        raise ValueError("Expected one spin for each Fe")
    count = r.shape[1]
    grad = np.stack([electronic_derivative(ri, bool(si), cfg) for ri, si in zip(r, spins)]) / count
    grad += cfg.K_red * (r - R_red) / count
    if cfg.geometry == "six_radial":
        grad += cfg.K_shape * (r - r.mean(axis=1, keepdims=True))
    if cfg.N_Fe == 2:
        grad += cfg.K_bridge * (r.mean(axis=1).sum() + cfg.L_bridge - cfg.R_bridge) / count
    return grad if cfg.geometry == "six_radial" else grad[:, 0]
