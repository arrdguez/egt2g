"""Ejecución: python -m sco0 --config examples/one_fe.json --out results/one_fe"""

import argparse
from pathlib import Path

from .config import load_config
from .equilibrium import equilibrium_scan
from .monte_carlo import run_mc
from .output import export
from .protocol import build_steps


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Modelo de campo de ligandos para la transición de espín Fe(II) HS/LS"
    )
    parser.add_argument(
        "--config", default="examples/one_fe.json", help="Archivo JSON de configuración"
    )
    parser.add_argument(
        "--out", default="results/one_fe", help="Directorio de salida (CSV y PNG)"
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    steps = build_steps(cfg)
    equilibrium = equilibrium_scan(steps, cfg)
    samples = run_mc(steps, cfg)
    export(Path(args.out), cfg, equilibrium, samples)

    print(
        f"Completados {len(steps)} pasos (protocolo '{cfg.protocol}'), "
        f"{cfg.replicas} réplicas."
    )
    print(f"Resultados en {Path(args.out).resolve()}")


if __name__ == "__main__":
    main()
