"""Quasistatic numerical reference: compare minima of every spin configuration."""

from __future__ import annotations

from itertools import product

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from .config import Config
from .model import gradient, ligand_delta, ligand_positions, total_energy


def _six_minimum(red: float, spins: np.ndarray, cfg: Config) -> tuple[float, np.ndarray]:
    shape = (cfg.N_Fe, 6)
    midpoint = (cfg.r_min + cfg.r_max) / 2
    starts = [np.full(shape, value) for value in (cfg.r_min + 0.1, midpoint, cfg.r_max - 0.1)]
    distorted = np.full(shape, midpoint)
    distorted[:, :4] += 0.05
    distorted[:, 4:] -= 0.05
    starts.append(distorted)
    candidates = []
    for start in starts:
        objective = lambda flat: total_energy(flat.reshape(shape), spins, red, cfg)
        jacobian = lambda flat: gradient(flat.reshape(shape), spins, red, cfg).ravel()
        result = minimize(objective, start.ravel(), jac=jacobian,
                          bounds=[(cfg.r_min, cfg.r_max)] * start.size, method="L-BFGS-B",
                          options={"ftol": 1e-12, "maxiter": 2000})
        if result.success:
            candidates.append((float(result.fun), result.x.reshape(shape)))
    if not candidates:
        raise RuntimeError("No fixed-spin minimization converged for the six-ligand geometry")
    return min(candidates, key=lambda item: item[0])


def surfaces_at(R_red: float, cfg: Config) -> list[dict]:
    results = []
    for combination in product((False, True), repeat=cfg.N_Fe):
        spins = np.array(combination, dtype=bool)
        if cfg.geometry == "six_radial":
            energy, distances = _six_minimum(R_red, spins, cfg)
        elif cfg.N_Fe == 1:
            objective = lambda r: total_energy(np.array([r]), spins, R_red, cfg)
            optimum = minimize_scalar(objective, bounds=(cfg.r_min, cfg.r_max), method="bounded",
                                      options={"xatol": 1e-11})
            candidates = [(objective(cfg.r_min), cfg.r_min), (objective(cfg.r_max), cfg.r_max),
                          (optimum.fun, optimum.x)]
            energy, r = min(candidates)
            distances = np.array([r])
        else:
            starts = [np.array(start) for start in product(
                (cfg.r_min, (cfg.r_min + cfg.r_max) / 2, cfg.r_max), repeat=2)]
            candidates = []
            for start in starts:
                result = minimize(total_energy, start, args=(spins, R_red, cfg),
                                  jac=gradient, bounds=[(cfg.r_min, cfg.r_max)] * 2,
                                  method="L-BFGS-B", options={"ftol": 1e-13, "gtol": 1e-9})
                if not result.success and np.max(np.abs(gradient(result.x, spins, R_red, cfg))) > 1e-4:
                    raise RuntimeError(f"Equilibrium minimization failed: {result.message}")
                candidates.append((total_energy(result.x, spins, R_red, cfg), result.x))
            energy, distances = min(candidates, key=lambda item: item[0])
        row = {"spin_pattern": "-".join("HS" if s else "LS" for s in spins),
               "energy_eV": float(energy), "hs_fraction": float(np.mean(spins))}
        positions = ligand_positions(distances, cfg) if cfg.geometry == "six_radial" else None
        for i, (r, s) in enumerate(zip(distances, spins), 1):
            row[f"r{i}_A"] = float(np.mean(r))
            row[f"spin{i}_HS"] = int(s)
            row[f"delta{i}_eV"] = ligand_delta(r, cfg)
            if cfg.geometry == "six_radial":
                for j, radius in enumerate(r, 1):
                    row[f"fe{i}_n{j}_A"] = float(radius)
                    for axis, coordinate in zip("xyz", positions[i - 1, j - 1]):
                        row[f"fe{i}_n{j}_{axis}_A"] = float(coordinate)
        results.append(row)
    return results


def equilibrium_scan(steps: list[dict], cfg: Config) -> tuple[list[dict], list[dict]]:
    all_surfaces = []
    ground_states = []
    for step in steps:
        candidates = [{**step, **row} for row in surfaces_at(step["R_red_A"], cfg)]
        all_surfaces.extend(candidates)
        ground_states.append(min(candidates, key=lambda row: row["energy_eV"]))
    return ground_states, all_surfaces
