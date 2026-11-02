"""Geometry-only stochastic descent on the instantaneous ground-state envelope."""

from __future__ import annotations

import csv
from contextlib import nullcontext
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from .config import Config
from .model import electronic_energy, ground_state, ligand_delta, ligand_positions


def initial_state(cfg: Config) -> tuple[np.ndarray, np.ndarray]:
    minimum = minimize_scalar(lambda r: electronic_energy(r, cfg.initial_branch == "HS", cfg),
                              bounds=(cfg.r_min, cfg.r_max), method="bounded")
    shape = (cfg.N_Fe,) if cfg.geometry == "radial" else (cfg.N_Fe, 6)
    distances = np.full(shape, minimum.x)
    return distances, ground_state(distances, cfg.R_300K, cfg)[1]


def attempt(distances: np.ndarray, spins: np.ndarray, red: float,
            cfg: Config, rng: np.random.Generator, details: dict | None = None) -> bool:
    """Trial geometry is temporary; the spin is recomputed, never proposed."""
    index = int(rng.integers(cfg.N_Fe))
    ligand = int(rng.integers(6)) if cfg.geometry == "six_radial" else 0
    new_distances = distances.copy()
    old_r = float(distances[index] if cfg.geometry == "radial" else distances[index, ligand])
    proposed_r = old_r + float(rng.uniform(-cfg.delta_r, cfg.delta_r))
    old_energy, current_spins = ground_state(distances, red, cfg)
    spins[:] = current_spins
    trial_energy = float("nan")
    trial_spin = int(spins[index])
    accepted = False
    if cfg.r_min <= proposed_r <= cfg.r_max:
        if cfg.geometry == "radial":
            new_distances[index] = proposed_r
        else:
            new_distances[index, ligand] = proposed_r
        trial_energy, new_spins = ground_state(new_distances, red, cfg)
        trial_spin = int(new_spins[index])
        accepted = trial_energy <= old_energy
    if accepted:
        distances[:] = new_distances
        spins[:] = new_spins
    if details is not None:
        details.update({"center": index + 1, "ligand": ligand + 1,
                        "r_before_A": old_r, "r_proposed_A": proposed_r,
                        "r_after_A": proposed_r if accepted else old_r,
                        "spin_before_HS": int(current_spins[index]),
                        "spin_proposed_HS": trial_spin if cfg.r_min <= proposed_r <= cfg.r_max else "",
                        "spin_after_HS": int(spins[index]), "energy_before_eV": old_energy,
                        "energy_proposed_eV": trial_energy,
                        "energy_after_eV": trial_energy if accepted else old_energy,
                        "accepted": int(accepted)})
    return accepted


def run_mc(steps: list[dict], cfg: Config, attempts_csv: Path | None = None) -> tuple[list[dict], list[dict]]:
    rows = []
    convergence = []
    fields = ["speed", "replica", "direction", "step", "nominal_T_K", "R_red_A", "sweep",
              "attempt_in_sweep", "center", "ligand", "r_before_A", "r_proposed_A", "r_after_A",
              "spin_before_HS", "spin_proposed_HS", "spin_after_HS", "energy_before_eV",
              "energy_proposed_eV", "energy_after_eV", "accepted"]
    if attempts_csv is not None:
        attempts_csv.parent.mkdir(parents=True, exist_ok=True)
    stream = open(attempts_csv, "w", newline="", encoding="utf-8") if attempts_csv else nullcontext(None)
    with stream as log:
        writer = csv.DictWriter(log, fieldnames=fields) if log is not None else None
        if writer is not None:
            writer.writeheader()
        streams = iter(np.random.SeedSequence(cfg.seed).spawn(2 * cfg.replicas))
        for speed, sweeps in (("fast", cfg.sweeps_fast), ("slow", cfg.sweeps_slow)):
            for replica in range(cfg.replicas):
                rng = np.random.default_rng(next(streams))
                distances, spins = initial_state(cfg)
                for direction, path in (("cooling", steps), ("heating", list(reversed(steps)))):
                    for position, step in enumerate(path):
                        # Re-evaluate the electronic state when the external constraint changes.
                        _, spins[:] = ground_state(distances, step["R_red_A"], cfg)
                        accepted = 0
                        for sweep in range(sweeps):
                            for trial in range(cfg.N_Fe * (6 if cfg.geometry == "six_radial" else 1)):
                                details = {} if writer is not None else None
                                accepted += attempt(distances, spins, step["R_red_A"], cfg, rng, details)
                                if writer is not None:
                                    writer.writerow({"speed": speed, "replica": replica,
                                                     "direction": direction, "step": position,
                                                     "nominal_T_K": step["nominal_T_K"],
                                                     "R_red_A": step["R_red_A"], "sweep": sweep + 1,
                                                     "attempt_in_sweep": trial + 1, **details})
                            if sweep == 0 or (sweep + 1) % max(1, sweeps // 10) == 0 or sweep == sweeps - 1:
                                convergence.append({"speed": speed, "replica": replica,
                                                    "direction": direction, "step": position,
                                                    "sweep": sweep + 1, "energy_eV": ground_state(
                                                        distances, step["R_red_A"], cfg)[0]})
                        energy, spins[:] = ground_state(distances, step["R_red_A"], cfg)
                        row = {**step, "step": position, "speed": speed, "replica": replica,
                               "direction": direction, "sweeps": sweeps,
                               "attempts": sweeps * cfg.N_Fe * (6 if cfg.geometry == "six_radial" else 1),
                               "accepted": accepted, "energy_eV": energy,
                               "hs_fraction": float(np.mean(spins))}
                        positions = ligand_positions(distances, cfg) if cfg.geometry == "six_radial" else None
                        for i in range(cfg.N_Fe):
                            site = distances[i]
                            row[f"r{i + 1}_A"] = float(np.mean(site))
                            row[f"spin{i + 1}_HS"] = int(spins[i])
                            row[f"delta{i + 1}_eV"] = ligand_delta(site, cfg)
                            if cfg.geometry == "six_radial":
                                for j, r in enumerate(site, 1):
                                    row[f"fe{i + 1}_n{j}_A"] = float(r)
                                    for axis, coordinate in zip("xyz", positions[i, j - 1]):
                                        row[f"fe{i + 1}_n{j}_{axis}_A"] = float(coordinate)
                        rows.append(row)
    return rows, convergence
