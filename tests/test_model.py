"""Run from the project root: python -m unittest discover -s tests -v"""

import unittest
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

import numpy as np
from scipy.optimize import minimize_scalar

from sco0.config import Config
from sco0.equilibrium import surfaces_at
from sco0.model import (crossing_distance, electronic_derivative, electronic_energy,
                        gradient, ground_state, ligand_delta, ligand_positions, total_energy)
from sco0.monte_carlo import attempt, initial_state, run_mc
from sco0.protocol import cooling_steps


class ModelTests(unittest.TestCase):
    def setUp(self):
        self.cfg = Config()

    def test_jorgensen_minima_and_crossing(self):
        minima = []
        for hs in (True, False):
            optimum = minimize_scalar(lambda r: electronic_energy(r, hs, self.cfg),
                                      bounds=(self.cfg.r_min, self.cfg.r_max), method="bounded")
            minima.append(optimum.x / self.cfg.r_ref)
            self.assertAlmostEqual(electronic_derivative(optimum.x, hs, self.cfg), 0, places=3)
        self.assertAlmostEqual(minima[0], 0.975, delta=0.005)
        self.assertAlmostEqual(minima[1], 0.867, delta=0.005)
        self.assertLess(minima[1], minima[0])
        crossing = crossing_distance(self.cfg)
        self.assertIsNotNone(crossing)
        self.assertAlmostEqual(electronic_energy(crossing, True, self.cfg),
                               electronic_energy(crossing, False, self.cfg), places=9)

    def test_forces_match_finite_difference(self):
        r = 2.07
        h = 1e-5
        for hs in (True, False):
            numeric = (electronic_energy(r + h, hs, self.cfg) -
                       electronic_energy(r - h, hs, self.cfg)) / (2 * h)
            self.assertAlmostEqual(electronic_derivative(r, hs, self.cfg), numeric, places=6)
        cfg = replace(self.cfg, N_Fe=2, K_bridge=8)
        distances = np.array([2.07, 2.15])
        spins = np.array([True, False])
        for i in range(2):
            plus = distances.copy()
            minus = distances.copy()
            plus[i] += h
            minus[i] -= h
            numeric = (total_energy(plus, spins, 2.15, cfg) -
                       total_energy(minus, spins, 2.15, cfg)) / (2 * h)
            self.assertAlmostEqual(gradient(distances, spins, 2.15, cfg)[i], numeric, places=6)

    def test_independent_centers_without_bridge(self):
        cfg = replace(self.cfg, N_Fe=2, K_bridge=0)
        distances = np.array([2.11, 1.99])
        spins = np.array([True, False])
        expected = sum(total_energy(np.array([r]), np.array([s]), 2.10, self.cfg)
                       for r, s in zip(distances, spins))
        self.assertAlmostEqual(total_energy(distances, spins, 2.10, cfg), expected)
        candidates = surfaces_at(2.10, cfg)
        self.assertEqual(len(candidates), 4)
        self.assertAlmostEqual(next(row["energy_eV"] for row in candidates
                                    if row["spin_pattern"] == "HS-LS"),
                               next(row["energy_eV"] for row in candidates
                                    if row["spin_pattern"] == "LS-HS"), places=7)

    def test_equilibrium_changes_spin_for_imposed_contraction(self):
        self.assertEqual(min(surfaces_at(self.cfg.R_300K, self.cfg),
                             key=lambda row: row["energy_eV"])["spin_pattern"], "HS")
        self.assertEqual(min(surfaces_at(self.cfg.R_0K, self.cfg),
                             key=lambda row: row["energy_eV"])["spin_pattern"], "LS")

    def test_bridge_competes_with_first_center_contraction(self):
        cfg = replace(self.cfg, N_Fe=2, K_bridge=3)
        options = {row["spin_pattern"]: row for row in surfaces_at(2.0, cfg)}
        self.assertLess(options["LS-HS"]["r1_A"], options["HS-HS"]["r1_A"])
        self.assertGreater(options["LS-HS"]["r2_A"], options["HS-HS"]["r2_A"])

    def test_spin_is_determined_at_fixed_geometry(self):
        crossover = crossing_distance(self.cfg)
        self.assertTrue(ground_state(np.array([crossover + 0.1]), 2.1, self.cfg)[1][0])
        self.assertFalse(ground_state(np.array([crossover - 0.1]), 2.1, self.cfg)[1][0])
        for r in (1.92, 2.06, 2.20):
            delta = ligand_delta(r, self.cfg)
            difference = electronic_energy(r, False, self.cfg) - electronic_energy(r, True, self.cfg)
            self.assertAlmostEqual(difference, self.cfg.P_LS - 2 * delta)

    def test_downhill_geometry_can_switch_spin_without_spin_proposal(self):
        cfg = replace(self.cfg, P_LS=3.9, delta_r=0.06)
        distances = np.array([2.10])
        _, spins = ground_state(distances, 1.96, cfg)
        self.assertTrue(spins[0])

        class FixedProposal:
            @staticmethod
            def integers(_):
                return 0

            @staticmethod
            def uniform(_lower, _upper):
                return -0.03

        details = {}
        self.assertTrue(attempt(distances, spins, 1.96, cfg, FixedProposal(), details))
        self.assertFalse(spins[0])
        self.assertEqual(details["spin_before_HS"], 1)
        self.assertEqual(details["spin_after_HS"], 0)
        self.assertLess(details["energy_after_eV"], details["energy_before_eV"])

    def test_six_ligands_recover_isotropic_energy_and_force(self):
        cfg = replace(self.cfg, geometry="six_radial")
        for hs in (False, True):
            radii = np.full((1, 6), 2.04)
            self.assertAlmostEqual(total_energy(radii, np.array([hs]), 2.1, cfg),
                                   total_energy(np.array([2.04]), np.array([hs]), 2.1, self.cfg))
            self.assertAlmostEqual(gradient(radii, np.array([hs]), 2.1, cfg).sum(),
                                   gradient(np.array([2.04]), np.array([hs]), 2.1, self.cfg)[0])
        r = np.array([[2.05] * 6])
        r[0, 0] = 1.99
        self.assertNotAlmostEqual(ligand_delta(r[0], cfg), ligand_delta(r.mean(), cfg))
        self.assertEqual(len(surfaces_at(2.05, cfg)), 2)

    def test_six_ligand_force_and_starting_geometries(self):
        cfg = replace(self.cfg, geometry="six_radial", N_Fe=2, K_bridge=3)
        r = np.array([[2.07, 2.08, 2.06, 2.10, 2.03, 2.12],
                      [2.16, 2.13, 2.17, 2.14, 2.18, 2.11]])
        spins = np.array([False, True])
        xyz = ligand_positions(r, cfg)
        np.testing.assert_array_equal(xyz[0, 0], [r[0, 0], 0, 0])
        np.testing.assert_array_equal(xyz[0, 1], [-r[0, 1], 0, 0])
        np.testing.assert_array_equal(xyz[0, 4], [0, 0, r[0, 4]])
        np.testing.assert_array_equal(xyz[0, 5], [0, 0, -r[0, 5]])
        h = 1e-5
        for i, j in ((0, 0), (0, 5), (1, 2)):
            plus, minus = r.copy(), r.copy()
            plus[i, j] += h
            minus[i, j] -= h
            numerical = (total_energy(plus, spins, 2.08, cfg) -
                         total_energy(minus, spins, 2.08, cfg)) / (2 * h)
            self.assertAlmostEqual(gradient(r, spins, 2.08, cfg)[i, j], numerical, places=6)
        r_hs, _ = initial_state(self.cfg)
        ls_config = replace(self.cfg, initial_branch="LS")
        r_ls, spins_ls = initial_state(ls_config)
        self.assertLess(r_ls[0], r_hs[0])
        np.testing.assert_array_equal(spins_ls, ground_state(r_ls, ls_config.R_300K, ls_config)[1])

    def test_six_ligands_moves_only_one_radius(self):
        cfg = replace(self.cfg, geometry="six_radial")
        distances = np.full((1, 6), 2.18)
        _, spins = ground_state(distances, 2.2, cfg)
        rng = np.random.default_rng(42)
        for _ in range(100):
            before = distances.copy()
            old_energy, _ = ground_state(distances, 2.2, cfg)
            details = {}
            accepted = attempt(distances, spins, 2.2, cfg, rng, details)
            self.assertLessEqual(np.count_nonzero(before != distances), 1)
            self.assertEqual(int(accepted), details["accepted"])
            self.assertLessEqual(ground_state(distances, 2.2, cfg)[0], old_energy + 1e-12)
            np.testing.assert_array_equal(spins, ground_state(distances, 2.2, cfg)[1])
            if not accepted:
                np.testing.assert_array_equal(before, distances)

    def test_mc_accepts_only_downhill_and_keeps_rejected_state(self):
        cfg = replace(self.cfg, r_min=1.65, r_max=2.55)
        rng = np.random.default_rng(42)
        distances = np.array([2.12])
        spins = np.array([True])
        for _ in range(1000):
            before = total_energy(distances, spins, 2.07, cfg)
            old_r, old_spin = distances.copy(), spins.copy()
            accepted = attempt(distances, spins, 2.07, cfg, rng)
            self.assertLessEqual(total_energy(distances, spins, 2.07, cfg), before + 1e-12)
            if not accepted:
                np.testing.assert_array_equal(old_r, distances)
                np.testing.assert_array_equal(old_spin, spins)

    def test_protocol_and_reproducibility(self):
        cfg = replace(self.cfg, delta_T=150, sweeps_fast=2, sweeps_slow=3, replicas=1)
        steps = cooling_steps(cfg)
        self.assertEqual([x["nominal_T_K"] for x in steps], [300.0, 150.0, 0.0])
        self.assertEqual(run_mc(steps, cfg), run_mc(steps, cfg))

    def test_distance_protocol_and_two_six_ligand_centers(self):
        cfg = replace(self.cfg, N_Fe=2, geometry="six_radial", K_bridge=0,
                      protocol="distance", R_start=2.2, R_end=2.1, delta_R=0.1,
                      sweeps_fast=1, sweeps_slow=2, replicas=1)
        steps = cooling_steps(cfg)
        self.assertEqual(len(steps), 2)
        self.assertTrue(np.isnan(steps[0]["nominal_T_K"]))
        r = np.array([[2.2] * 6, [1.9] * 6])
        _, spins = ground_state(r, steps[0]["R_red_A"], cfg)
        np.testing.assert_array_equal(spins, [True, False])
        rows, _ = run_mc(steps, cfg)
        self.assertEqual(rows[0]["attempts"], 6 * cfg.N_Fe)
        self.assertIn("fe2_n6_z_A", rows[0])

    def test_each_logged_proposal_preserves_rejection_and_is_downhill(self):
        cfg = replace(self.cfg, delta_T=150, sweeps_fast=4, sweeps_slow=6, replicas=1)
        with TemporaryDirectory() as folder:
            path = Path(folder) / "attempts.csv"
            run_mc(cooling_steps(cfg), cfg, path)
            data = pd.read_csv(path)
            self.assertEqual(len(data), 2 * 3 * (4 + 6))
            rejected = data[data.accepted == 0]
            np.testing.assert_allclose(rejected.r_before_A, rejected.r_after_A)
            np.testing.assert_array_equal(rejected.spin_before_HS, rejected.spin_after_HS)
            np.testing.assert_allclose(rejected.energy_before_eV, rejected.energy_after_eV)
            accepted = data[data.accepted == 1]
            self.assertTrue((accepted.energy_after_eV <= accepted.energy_before_eV).all())


if __name__ == "__main__":
    unittest.main()
