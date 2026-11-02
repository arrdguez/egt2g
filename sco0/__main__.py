"""Run using: python -m sco0 --config examples/one_fe.json --out results/one_fe"""

import argparse
from pathlib import Path

from .config import load_config
from .equilibrium import equilibrium_scan
from .monte_carlo import run_mc
from .output import export
from .protocol import cooling_steps


def main() -> None:
    parser = argparse.ArgumentParser(description="Fe(II) electronic/structural spin crossover at 0 K")
    parser.add_argument("--config", default="examples/one_fe.json", help="JSON configuration file")
    parser.add_argument("--out", default="results/one_fe", help="Output directory (CSV and PNG)")
    args = parser.parse_args()
    cfg = load_config(args.config)
    steps = cooling_steps(cfg)
    ground, surfaces = equilibrium_scan(steps, cfg)
    samples, convergence = run_mc(steps, cfg, Path(args.out) / "attempts.csv")
    export(Path(args.out), cfg, ground, surfaces, samples, convergence)
    print(f"Completed {len(steps)} imposed-distance steps, {cfg.replicas} replicas per speed.")
    print(f"Results saved to {Path(args.out).resolve()}")


if __name__ == "__main__":
    main()
