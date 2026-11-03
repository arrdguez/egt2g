"""CSV y PNG del modelo de campo de ligandos (numpy + matplotlib).

Genera:
  - config_used.json   configuración efectiva
  - equilibrium.csv    referencia analítica (fracción HS exacta)
  - monte_carlo.csv    fracción HS muestreada por Metropolis
  - surfaces.png       E_HS(r), E_LS(r), Dq(r) y Δo(r) con el cruce
  - crossover.png      fracción HS vs r (o vs T): analítico frente a Monte Carlo
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .config import Config
from .model import crossing_distance, delta_o, electronic_energy, ligand_field


def export(
    out: Path, cfg: Config, equilibrium: list[dict], samples: list[dict]
) -> None:
    """Escribe la configuración, los CSV y las figuras PNG."""
    out.mkdir(parents=True, exist_ok=True)
    (out / "config_used.json").write_text(
        json.dumps(asdict(cfg), indent=2), encoding="utf-8"
    )

    _write_csv(out / "equilibrium.csv", equilibrium)
    _write_csv(out / "monte_carlo.csv", samples)

    _plot_surfaces(out, cfg)
    _plot_scan(out, cfg, equilibrium, samples)


def _write_csv(path: Path, rows: list[dict]) -> None:
    """Escribe una lista de diccionarios como CSV (sin pandas)."""
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _plot_surfaces(out: Path, cfg: Config) -> None:
    """Figura 1: superficies electrónicas y campo de ligandos en función de r."""
    radii = np.linspace(cfg.r_min, cfg.r_max, 600)
    r_c = crossing_distance(cfg)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4))

    axes[0].plot(radii, electronic_energy(radii, True, cfg), label="HS (t₂g⁴e_g²)")
    axes[0].plot(radii, electronic_energy(radii, False, cfg), label="LS (t₂g⁶)")
    if r_c is not None:
        axes[0].axvline(r_c, color="gray", linestyle=":", label=f"r_c = {r_c:.3f} Å")
    axes[0].set(xlabel="Fe–N (Å)", ylabel="E (cm⁻¹)", title="Superficies electrónicas")
    axes[0].legend()

    axes[1].plot(radii, ligand_field(radii, cfg))
    if r_c is not None:
        axes[1].axvline(r_c, color="gray", linestyle=":")
    axes[1].set(
        xlabel="Fe–N (Å)", ylabel="Dq (cm⁻¹)", title="Campo de ligandos Dq(r) ∝ r⁻⁵"
    )

    axes[2].plot(radii, delta_o(radii, cfg))
    if r_c is not None:
        axes[2].axvline(r_c, color="gray", linestyle=":")
    axes[2].set(
        xlabel="Fe–N (Å)", ylabel="Δo = 10·Dq (cm⁻¹)", title="Desdoblamiento octaédrico"
    )

    fig.tight_layout()
    fig.savefig(out / "surfaces.png", dpi=160)
    plt.close(fig)


def _plot_scan(
    out: Path, cfg: Config, equilibrium: list[dict], samples: list[dict]
) -> None:
    """Figura 2: fracción HS del barrido (referencia analítica frente a Monte Carlo)."""
    fig, ax = plt.subplots(figsize=(8, 5))
    if cfg.protocol == "distance":
        x_eq = [row["r_A"] for row in equilibrium]
        x_mc = [row["r_A"] for row in samples]
        ax.set(
            xlabel="Distancia Fe–N (Å)",
            title="Fracción HS vs distancia (el campo de ligandos provoca el cruce)",
        )
    else:
        x_eq = [row["T_K"] for row in equilibrium]
        x_mc = [row["T_K"] for row in samples]
        ax.set(
            xlabel="Temperatura (K)",
            title="Fracción HS vs temperatura (centro aislado, cruce suave)",
        )

    ax.plot(
        x_eq,
        [row["hs_fraction"] for row in equilibrium],
        "-",
        color="black",
        linewidth=2,
        label="Analítico (2 niveles)",
    )
    ax.plot(
        x_mc,
        [row["hs_fraction"] for row in samples],
        "o",
        markersize=3,
        color="crimson",
        label="Monte Carlo (T finita)",
    )
    ax.set(ylabel="Fracción HS", ylim=(-0.05, 1.05))
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "crossover.png", dpi=160)
    plt.close(fig)
