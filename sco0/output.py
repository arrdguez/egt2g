"""CSV and PNG reports of zero-temperature model calculations."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from .config import Config
from .model import crossing_distance, electronic_energy, ligand_delta


def export(out: Path, cfg: Config, equilibrium: list[dict], surfaces: list[dict],
           samples: list[dict], convergence: list[dict]) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / "config_used.json").write_text(json.dumps(asdict(cfg), indent=2), encoding="utf-8")
    eq = pd.DataFrame(equilibrium)
    eq.to_csv(out / "equilibrium.csv", index=False)
    pd.DataFrame(surfaces).to_csv(out / "spin_surfaces.csv", index=False)
    df = pd.DataFrame(samples)
    df.to_csv(out / "trajectories.csv", index=False)
    pd.DataFrame(convergence).to_csv(out / "convergence.csv", index=False)
    numeric = ["energy_eV", "hs_fraction", "accepted", "attempts"] + [
        value for i in range(1, cfg.N_Fe + 1)
        for value in (f"r{i}_A", f"spin{i}_HS", f"delta{i}_eV")]
    if cfg.geometry == "six_radial":
        numeric.extend(f"fe{i}_n{j}_A" for i in range(1, cfg.N_Fe + 1) for j in range(1, 7))
    averaged = df.groupby(["speed", "direction", "step", "R_red_A", "nominal_T_K"],
                          dropna=False, as_index=False)[numeric].mean()
    averaged.to_csv(out / "averages.csv", index=False)

    radii = np.linspace(cfg.r_min, cfg.r_max, 600)
    fig, axes = plt.subplots(1, 3, figsize=(16, 4))
    for hs, name in ((True, "HS"), (False, "LS")):
        energy = electronic_energy(radii, hs, cfg)
        axes[0].plot(radii, energy, label=name)
        minimum = minimize_scalar(lambda r: electronic_energy(r, hs, cfg),
                                  bounds=(cfg.r_min, cfg.r_max), method="bounded")
        axes[0].scatter([minimum.x], [minimum.fun], s=28)
    axes[0].set(xlabel="Fe–N (Å)", ylabel="E electrónica (eV)", title="Superficies sin red")
    axes[0].legend()
    axes[1].plot(radii, electronic_energy(radii, False, cfg) - electronic_energy(radii, True, cfg))
    axes[1].axhline(0, color="black", linewidth=0.8)
    crossing = crossing_distance(cfg)
    if crossing is not None:
        axes[0].axvline(crossing, linestyle=":", color="gray")
        axes[1].axvline(crossing, linestyle=":", color="gray", label=f"Cruce {crossing:.3f} Å")
        axes[1].legend()
    axes[1].set(xlabel="Fe–N (Å)", ylabel="E LS − E HS (eV)", title="Competencia electrónica")
    axes[2].plot(radii, [ligand_delta(r, cfg) for r in radii])
    axes[2].set(xlabel="Fe–N (Å)", ylabel="Δ_o efectivo (eV)", title="Campo ligando ∝ r⁻⁵")
    _save(fig, out / "surfaces.png")

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for i in range(1, cfg.N_Fe + 1):
        axes[0, 0].plot(eq["R_red_A"], eq[f"r{i}_A"], color="black", linewidth=2,
                        label=f"equilibrio Fe{i}" if cfg.N_Fe == 2 else "equilibrio")
        axes[0, 1].plot(eq["R_red_A"], eq[f"spin{i}_HS"], color="black", linewidth=1)
    for (speed, direction), group in averaged.groupby(["speed", "direction"]):
        group = group.sort_values("step")
        label = f"{speed} {direction}"
        style = "-" if direction == "cooling" else "--"
        for i in range(1, cfg.N_Fe + 1):
            suffix = f" Fe{i}" if cfg.N_Fe == 2 else ""
            axes[0, 0].plot(group["R_red_A"], group[f"r{i}_A"], style, label=label + suffix)
            axes[0, 1].plot(group["R_red_A"], group[f"spin{i}_HS"], style, label=label + suffix)
        axes[1, 0].plot(group["R_red_A"], group["energy_eV"], style, label=label)
        if cfg.protocol == "nominal":
            axes[1, 1].plot(group["nominal_T_K"], group["hs_fraction"], style, label=label)
        else:
            axes[1, 1].plot(group["R_red_A"], group["hs_fraction"], style, label=label)
    axes[0, 0].set(xlabel="R_red (Å)", ylabel="r_Fe–N (Å)", title="Contracción: enlace frente a red")
    axes[0, 1].set(xlabel="R_red (Å)", ylabel="HS = 1; LS = 0 (media MC)", title="Estado de cada Fe")
    axes[1, 0].set(xlabel="R_red (Å)", ylabel="E total (eV)", title="Energía de la trayectoria")
    axes[1, 1].set(xlabel="Temperatura nominal (K)" if cfg.protocol == "nominal" else "R_red (Å)",
                   ylabel="Fracción HS", title="Barridos inversos y velocidades")
    for ax in axes.flat:
        ax.legend(fontsize=7)
    _save(fig, out / "trajectories.png")

    if cfg.protocol == "nominal":
        fig, ax = plt.subplots(figsize=(8, 5))
        for (speed, direction), group in averaged.groupby(["speed", "direction"]):
            for i in range(1, cfg.N_Fe + 1):
                ax.plot(group["nominal_T_K"], group[f"r{i}_A"],
                        label=f"{speed} {direction} Fe{i}")
        ax.set(xlabel="Temperatura nominal (K), no termodinámica", ylabel="Fe–N (Å)")
        ax.legend(fontsize=8)
        _save(fig, out / "nominal_distance.png")

    if cfg.geometry == "six_radial":
        fig, axes = plt.subplots(cfg.N_Fe, 1, figsize=(9, 4 * cfg.N_Fe), squeeze=False)
        for i, ax in enumerate(axes[:, 0], 1):
            for j in range(1, 7):
                subset = averaged[(averaged["speed"] == "slow") &
                                  (averaged["direction"] == "cooling")]
                ax.plot(subset["R_red_A"], subset[f"fe{i}_n{j}_A"],
                        label=f"N{j} {'ecuatorial' if j <= 4 else 'axial'}")
            ax.set(xlabel="R_red (Å)", ylabel="Distancia Fe–N (Å)",
                   title=f"Fe{i}: seis enlaces, enfriamiento lento (media entre réplicas)")
            ax.legend(fontsize=8, ncol=2)
        _save(fig, out / "ligand_distances.png")

    fig, ax = plt.subplots(figsize=(8, 5))
    example = pd.DataFrame(convergence)
    example = example[(example["speed"] == "slow") & (example["replica"] == 0) &
                      (example["direction"] == "cooling") & (example["step"] == 0)]
    ax.plot(example["sweep"], example["energy_eV"], marker=".")
    ax.set(xlabel="Barridos MC en el primer paso", ylabel="E total (eV)",
           title="Convergencia: ejemplo (todos los pasos en convergence.csv)")
    _save(fig, out / "convergence.png")


def _save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
