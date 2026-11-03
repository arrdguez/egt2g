"""Pruebas del modelo de campo de ligandos. Ejecutar: python -m unittest discover -s tests -v"""

import unittest
from dataclasses import replace

import numpy as np

from sco0.config import Config
from sco0.equilibrium import equilibrium_scan, crossing_info
from sco0.metropolis import sample_hs_fraction
from sco0.model import (
    crossing_distance,
    delta_o,
    electronic_derivative,
    electronic_energy,
    energy_gap,
    hs_fraction,
    ligand_field,
)
from sco0.monte_carlo import run_mc
from sco0.protocol import build_steps


class ModelTests(unittest.TestCase):
    def setUp(self):
        self.cfg = Config()
        # Distancia de cruce de la configuración por defecto.
        self.r_c = crossing_distance(self.cfg)

    # --- Campo de ligandos y energías ---

    def test_ligand_field_decays_with_distance(self):
        self.assertAlmostEqual(
            ligand_field(self.cfg.r0, self.cfg), self.cfg.Dq0, places=9
        )
        self.assertGreater(ligand_field(1.9, self.cfg), ligand_field(2.1, self.cfg))

    def test_electronic_energy_formulas(self):
        for r in (1.9, 2.1, 2.3):
            self.assertAlmostEqual(
                electronic_energy(r, True, self.cfg),
                -4.0 * ligand_field(r, self.cfg),
                places=9,
            )
            self.assertAlmostEqual(
                electronic_energy(r, False, self.cfg),
                -24.0 * ligand_field(r, self.cfg) + self.cfg.P,
                places=9,
            )

    def test_energy_gap_and_crossing(self):
        # En r_c el gap debe anularse y ambas energías coincidir.
        self.assertAlmostEqual(energy_gap(self.r_c, self.cfg), 0.0, places=6)
        self.assertAlmostEqual(
            electronic_energy(self.r_c, True, self.cfg),
            electronic_energy(self.r_c, False, self.cfg),
            places=6,
        )
        # Por encima de r_c domina HS (gap < 0); por debajo, LS (gap > 0).
        self.assertLess(energy_gap(self.r_c + 0.1, self.cfg), 0.0)
        self.assertGreater(energy_gap(self.r_c - 0.1, self.cfg), 0.0)

    def test_derivative_matches_finite_difference(self):
        r, h = 2.07, 1e-5
        for hs in (True, False):
            numeric = (
                electronic_energy(r + h, hs, self.cfg)
                - electronic_energy(r - h, hs, self.cfg)
            ) / (2 * h)
            self.assertAlmostEqual(
                electronic_derivative(r, hs, self.cfg), numeric, delta=1e-3
            )

    def test_delta_o_is_ten_times_dq(self):
        self.assertAlmostEqual(
            delta_o(2.0, self.cfg), 10.0 * ligand_field(2.0, self.cfg), places=9
        )

    # --- Fracción térmica analítica ---

    def test_hs_fraction_limits(self):
        # Lejos por encima de r_c (campo débil) → casi todo HS.
        self.assertGreater(hs_fraction(self.r_c + 0.2, 300.0, self.cfg), 0.99)
        # Lejos por debajo (campo fuerte) → casi todo LS.
        self.assertLess(hs_fraction(self.r_c - 0.2, 300.0, self.cfg), 0.01)
        # En el cruce exacto la fracción es 1/2 a cualquier T > 0.
        self.assertAlmostEqual(hs_fraction(self.r_c, 300.0, self.cfg), 0.5, places=9)

    def test_hs_fraction_zero_temperature_is_a_step(self):
        # En T = 0 la fracción es 1 (HS) si gap <= 0 y 0 (LS) si gap > 0.
        self.assertEqual(hs_fraction(self.r_c + 0.05, 0.0, self.cfg), 1.0)
        self.assertEqual(hs_fraction(self.r_c - 0.05, 0.0, self.cfg), 0.0)

    # --- Metropolis ---

    def test_metropolis_matches_analytic(self):
        # La media de Metropolis debe acercarse a la fracción analítica exacta.
        rng = np.random.default_rng(7)
        for r in (self.r_c - 0.004, self.r_c, self.r_c + 0.004):
            expected = hs_fraction(r, 300.0, self.cfg)
            measured = sample_hs_fraction(r, 300.0, self.cfg, 40000, rng)
            self.assertAlmostEqual(measured, expected, delta=0.04)

    def test_metropolis_zero_temperature_is_ground_state(self):
        rng = np.random.default_rng(3)
        # gap < 0 ⇒ HS es el suelo ⇒ la fracción debe ser 1 tras termalizar.
        self.assertEqual(
            sample_hs_fraction(self.r_c + 0.1, 0.0, self.cfg, 5000, rng), 1.0
        )
        # gap > 0 ⇒ LS es el suelo ⇒ la fracción debe ser 0.
        self.assertEqual(
            sample_hs_fraction(self.r_c - 0.1, 0.0, self.cfg, 5000, rng), 0.0
        )

    # --- Protocolos y reproducción ---

    def test_distance_protocol_steps(self):
        cfg = replace(
            self.cfg, protocol="distance", R_start=2.4, R_end=2.38, delta_R=0.01
        )
        steps = build_steps(cfg)
        self.assertEqual([s["r_A"] for s in steps], [2.4, 2.39, 2.38])

    def test_temperature_protocol_steps(self):
        cfg = replace(
            self.cfg, protocol="temperature", T_start=10.0, T_end=0.0, delta_T=4.0
        )
        steps = build_steps(cfg)
        self.assertEqual([s["T_K"] for s in steps], [10.0, 6.0, 2.0, 0.0])

    def test_equilibrium_scan_uses_the_right_control(self):
        cfg = replace(
            self.cfg, protocol="distance", R_start=2.4, R_end=2.38, delta_R=0.01
        )
        rows = equilibrium_scan(build_steps(cfg), cfg)
        self.assertEqual(len(rows), 3)
        self.assertEqual(rows[0]["T_K"], cfg.T_fixed)
        self.assertIn("Dq_cm1", rows[0])

    def test_run_mc_is_reproducible(self):
        cfg = replace(
            self.cfg,
            protocol="distance",
            R_start=2.2,
            R_end=2.15,
            delta_R=0.05,
            replicas=2,
            mc_steps=2000,
        )
        steps = build_steps(cfg)
        self.assertEqual(run_mc(steps, cfg), run_mc(steps, cfg))

    def test_crossing_info(self):
        info = crossing_info(self.cfg)
        self.assertAlmostEqual(info["r_c_A"], self.r_c, places=9)
        self.assertAlmostEqual(info["delta_o_c_cm1"], 10.0 * info["Dq_c_cm1"], places=6)


if __name__ == "__main__":
    unittest.main()
