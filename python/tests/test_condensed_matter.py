"""Analytical and independent-reference checks for condensed-matter tools."""

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from numpy.testing import assert_allclose
from atomforge.science.dynamics_analysis import (
    mean_square_displacement, diffusion_coefficient, velocity_autocorrelation, vibrational_spectrum)
from atomforge.science.local_structure import (
    local_strain, centrosymmetry, bond_order, wigner_seitz, static_structure_factor)
from atomforge.science.electronic_properties import band_gap, effective_mass, work_function
from atomforge.science.thermomechanics import (
    equation_of_state, elastic_tensor, phonon_dos, harmonic_thermodynamics)
from atomforge.science.advanced_simulation import migration_path, nvt_dynamics, npt_dynamics, reciprocal_path
from atomforge.science.simulation import from_ase
from ase.calculators.calculator import Calculator, all_changes


class CurvedDoubleWell(Calculator):
    implemented_properties = ["energy", "forces"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        x, y, z = self.atoms.positions.T
        residual = y - .2*(1-x*x)
        self.results = {"energy": float(np.sum((x*x-1)**2 + 5*residual**2 + z*z)),
                        "forces": np.column_stack((-4*x*(x*x-1)-4*x*residual, -10*residual, -2*z))}


class FreeParticles(Calculator):
    implemented_properties = ["energy", "forces", "stress"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.results = {"energy": 0., "forces": np.zeros_like(self.atoms.positions), "stress": np.zeros(6)}


class DynamicsTests(unittest.TestCase):
    def test_ballistic_msd_and_drift(self):
        velocity = np.array([[.1, .2, .3], [-.2, .1, 0]])
        times = np.arange(12) * 2.
        positions = times[:, None, None] * velocity
        result = mean_square_displacement(positions, 2)
        expected = times ** 2 * np.mean(np.sum(velocity ** 2, axis=1))
        assert_allclose(result["msd_A2"], expected, atol=1e-14)
        translation = np.tile(positions[:, :1], (1, 3, 1))
        assert_allclose(mean_square_displacement(translation, 2, remove_drift=True)["msd_A2"], 0, atol=1e-28)

    def test_wrapping_skew_and_partial_pbc(self):
        cell = np.array([[5., 0, 0], [4., 3., 0], [0, 0, 6.]])
        fraction = np.arange(20)[:, None, None] * np.array([[[.13, .07, .01]]])
        unwrapped = fraction @ cell
        wrapped = (fraction % 1) @ cell
        a = mean_square_displacement(unwrapped, 1)
        b = mean_square_displacement(wrapped, 1, wrapped=True, cell=cell)
        assert_allclose(a["msd_A2"], b["msd_A2"], atol=1e-12)
        with self.assertRaises(ValueError):
            mean_square_displacement(wrapped, 1, wrapped=True)

    def test_diffusion_units_and_fit(self):
        time = np.arange(101.) * 10
        result = diffusion_coefficient(time, 6 * .003 * time + .4, fit_range_fs=(100, 800))
        self.assertAlmostEqual(result["D_A2_per_fs"], .003)
        self.assertAlmostEqual(result["D_m2_per_s"], 3e-8)
        self.assertAlmostEqual(result["intercept_A2"], .4)
        self.assertAlmostEqual(result["r_squared"], 1)
        with self.assertRaises(ValueError):
            diffusion_coefficient(time, 30 - time * .01, fit_range_fs=(0, 500))

    def test_brownian_diffusion(self):
        rng = np.random.default_rng(824)
        diffusion = .004
        positions = np.cumsum(rng.normal(size=(400, 300, 3)) * np.sqrt(2 * diffusion), axis=0)
        msd = mean_square_displacement(positions, 1, max_lag=60)
        fit = diffusion_coefficient(msd["lag_fs"], msd["msd_A2"], fit_range_fs=(10, 60))
        self.assertLess(abs(fit["D_A2_per_fs"] / diffusion - 1), .06)

    def test_fft_vacf_matches_direct_correlation(self):
        values = np.random.default_rng(10).normal(size=(35, 4, 3))
        result = velocity_autocorrelation(values, .5)
        direct = [np.mean(np.sum(values[:len(values)-lag] * values[lag:], axis=2))
                  for lag in range(len(values))]
        assert_allclose(result["vacf"], direct, atol=2e-14)
        self.assertAlmostEqual(velocity_autocorrelation(values, 1, normalize=True)["vacf"][0], 1)
        with self.assertRaises(ValueError):
            velocity_autocorrelation(np.zeros_like(values), 1, normalize=True)

    def test_vibrational_frequency_and_normalization(self):
        time = np.arange(1000.)
        velocity = np.zeros((1000, 2, 3))
        velocity[:, :, 0] = np.cos(2 * np.pi * .02 * time)[:, None]
        result = vibrational_spectrum(velocity, 1, masses=[1, 12])
        peak = result["frequency_THz"][np.argmax(result["density_per_THz"])]
        self.assertAlmostEqual(peak, 20)
        self.assertAlmostEqual(np.trapezoid(result["density_per_THz"], result["frequency_THz"]), 1)
        with self.assertRaises(ValueError):
            vibrational_spectrum(np.zeros_like(velocity), 1)


class LocalStructureTests(unittest.TestCase):
    def test_affine_strain_and_rotation(self):
        from ase.build import bulk
        atoms = bulk("Cu", "fcc", a=3.6, cubic=True).repeat((2, 2, 2))
        deformation = np.array([[1.03, .08, 0], [0, .98, .01], [.02, 0, 1.04]])
        expected = (deformation.T @ deformation - np.eye(3)) / 2
        result = local_strain(atoms.positions, atoms.positions @ deformation.T, 2.8,
                              reference_cell=atoms.cell, current_cell=atoms.cell.array @ deformation.T,
                              pbc=(True, True, True))
        self.assertTrue(result["valid"].all())
        assert_allclose(result["green_lagrange_strain"], np.tile(expected, (len(atoms), 1, 1)), atol=1e-14)
        assert_allclose(result["d2min_A2"], 0, atol=1e-25)
        angle = .3
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        rotated = local_strain(atoms.positions, atoms.positions @ rotation.T, 2.8,
                               reference_cell=atoms.cell, current_cell=atoms.cell.array @ rotation.T,
                               pbc=(True, True, True))
        assert_allclose(rotated["green_lagrange_strain"], 0, atol=1e-14)
        perturbed = atoms.positions.copy(); perturbed[0, 0] += .15
        nonaffine = local_strain(atoms.positions, perturbed, 2.8, reference_cell=atoms.cell, pbc=(True, True, True))
        self.assertGreater(nonaffine["d2min_A2"].max(), .01)

    def test_rank_deficient_strain(self):
        planar = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]]
        result = local_strain(planar, planar, 2)
        self.assertFalse(result["valid"].any())

    def test_csp_fcc_and_distortion(self):
        from ase.build import bulk
        atoms = bulk("Cu", "fcc", a=3.6, cubic=True).repeat((2, 2, 2))
        result = centrosymmetry(atoms.positions, 2.8, cell=atoms.cell, pbc=(True, True, True))
        assert_allclose(result["centrosymmetry_A2"], 0, atol=1e-25)
        atoms.positions[0, 0] += .1
        result = centrosymmetry(atoms.positions, 2.8, cell=atoms.cell, pbc=(True, True, True))
        self.assertGreater(result["centrosymmetry_A2"][0], .01)

    def test_fcc_bond_order_and_rotation_invariance(self):
        from ase.build import bulk
        atoms = bulk("Cu", "fcc", a=3.6)  # single-site primitive cell tests periodic self images
        result = bond_order(atoms.positions, 2.8, cell=atoms.cell, pbc=(True, True, True))
        assert_allclose(result["order"]["q4"], .1909406539564933, atol=1e-12)
        assert_allclose(result["order"]["q6"], .5745242597140697, atol=1e-12)
        atoms.rotate(37, (1, 2, 3), rotate_cell=True)
        rotated = bond_order(atoms.positions, 2.8, cell=atoms.cell, pbc=(True, True, True))
        assert_allclose(rotated["order"]["q6"], result["order"]["q6"], atol=1e-12)

    def test_wigner_seitz_frenkel_and_box_mapping(self):
        sites = np.array([[0., 0, 0], [2., 0, 0], [4., 0, 0]])
        moved = sites.copy(); moved[1] = [.2, 0, 0]
        result = wigner_seitz(sites, moved)
        self.assertEqual(result["vacancies"], 1)
        self.assertEqual(result["interstitial_excess"], 1)
        assert_allclose(result["occupancy"], [2, 0, 1])
        self.assertEqual(int(np.sum(result["occupancy"])), len(moved))
        result = wigner_seitz(sites, sites * 1.2, cell=np.eye(3) * 6,
                              current_cell=np.eye(3) * 7.2, pbc=(True, True, True))
        assert_allclose(result["occupancy"], 1)
        self.assertTrue(wigner_seitz(sites, [[1, 0, 0]])["ambiguous"][0])

    def test_static_structure_factor_extinction(self):
        positions = [[0, 0, 0], [.5, .5, .5]]
        q = 2 * np.pi * np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]])
        result = static_structure_factor(positions, q)
        assert_allclose(result["S_q"], [2, 0, 2], atol=1e-14)
        translated = np.asarray(positions) + [7.1, -3.4, 2.6]
        assert_allclose(static_structure_factor(translated, q)["S_q"], result["S_q"], atol=1e-14)


class ElectronicPropertyTests(unittest.TestCase):
    def test_direct_indirect_and_metal(self):
        result = band_gap([[-1, 2], [-.2, 3]], 0)
        self.assertAlmostEqual(result["gap_eV"], 2.2)
        self.assertAlmostEqual(result["direct_gap_eV"], 3)
        self.assertFalse(result["metal_on_sampled_mesh"])
        self.assertTrue(band_gap([[-1, .1], [.1, 2]], 0)["metal_on_sampled_mesh"])
        self.assertTrue(band_gap([[-1, 0, 2]], 0)["metal_on_sampled_mesh"])
        spin = [[[-1, 2], [-.2, 3]], [[-2, 1], [-1, 1.5]]]
        self.assertAlmostEqual(band_gap(spin, 0)["gap_eV"], 1.2)

    def test_effective_mass_full_tensor(self):
        from scipy.spatial.transform import Rotation
        rotation = Rotation.from_rotvec([.3, -.1, .2]).as_matrix()
        mass = rotation @ np.diag([.2, .5, 1.4]) @ rotation.T
        hessian = 7.619964231073853 * np.linalg.inv(mass)
        k = np.random.default_rng(44).normal(size=(60, 3)) * .02
        energy = 2.4 + k @ [.01, -.02, .005] + np.einsum("ni,ij,nj->n", k, hessian, k)/2
        result = effective_mass(k, energy, center_inv_A=[0, 0, 0])
        assert_allclose(result["mass_tensor_m_e"], mass, atol=1e-12)
        assert_allclose(result["gradient_eV_A"], [.01, -.02, .005], atol=1e-12)
        self.assertLess(result["rms_fit_error_eV"], 1e-12)
        with self.assertRaises(ValueError):
            effective_mass(np.column_stack((np.arange(10), np.zeros((10, 2)))), np.arange(10), center_inv_A=[0, 0, 0])

    def test_work_function_and_sloped_vacuum(self):
        x = np.linspace(0, 20, 201)
        potential = np.where(x < 10, 2., 6.)
        result = work_function(x, potential, 1.5, vacuum_range_A=[12, 18])
        self.assertAlmostEqual(result["work_function_eV"], 4.5)
        self.assertTrue(result["flat_vacuum"])
        result = work_function(x, potential + .1*x, 1.5, vacuum_range_A=[12, 18])
        self.assertFalse(result["flat_vacuum"])
        self.assertAlmostEqual(result["slope_eV_per_A"], .1)


class ThermomechanicsTests(unittest.TestCase):
    def test_birch_murnaghan_recovery(self):
        # Independent explicit third-order energy formula.
        volume = np.linspace(14, 18, 13)
        v0, e0, b0, derivative = 16.2, -4.1, .7, 4.3
        eta = (v0/volume)**(2/3)-1
        energy = e0 + 9*v0*b0/16 * (derivative*eta**3 + (6-4*(eta+1))*eta**2)
        result = equation_of_state(volume, energy)
        self.assertAlmostEqual(result["volume_A3"], v0, places=6)
        self.assertAlmostEqual(result["bulk_modulus_GPa"], b0*160.2176634, places=5)
        self.assertAlmostEqual(result["bulk_derivative"], derivative, places=5)
        with self.assertRaises(ValueError):
            equation_of_state(volume, volume)

    def test_isotropic_elasticity_and_singular_inputs(self):
        bulk, shear = 120., 50.
        lame = bulk - 2*shear/3
        tensor = np.zeros((6, 6)); tensor[:3, :3] = lame
        np.fill_diagonal(tensor[:3, :3], lame + 2*shear)
        np.fill_diagonal(tensor[3:, 3:], shear)
        strain = np.vstack((np.zeros((1, 6)), np.eye(6)*.005, -np.eye(6)*.005))
        prestress = np.array([.1, -.2, .05, 0, .03, 0])
        result = elastic_tensor(strain, prestress + strain @ tensor.T)
        assert_allclose(result["stiffness_GPa"], tensor, atol=1e-11)
        self.assertAlmostEqual(result["moduli"]["bulk_hill_GPa"], bulk)
        self.assertAlmostEqual(result["moduli"]["shear_hill_GPa"], shear)
        self.assertAlmostEqual(result["moduli"]["universal_anisotropy"], 0)
        self.assertTrue(result["stable_zero_prestress"])
        with self.assertRaises(ValueError):
            elastic_tensor(np.zeros((10, 6)), np.zeros((10, 6)))
        tensor[5, 5] = -1
        result = elastic_tensor(strain, strain @ tensor.T)
        self.assertFalse(result["stable_zero_prestress"])
        self.assertIsNone(result["moduli"])

    def test_phonon_dos_integral_and_weights(self):
        grid = np.linspace(-.05, .1, 3001)
        result = phonon_dos([[.01, .02, .03], [-.01, .02, .03]], grid, weights=[3, 1])
        self.assertAlmostEqual(result["enclosed_modes"], 3, places=10)
        self.assertAlmostEqual(result["imaginary_weight"], .25)

    def test_harmonic_limits_and_derivatives(self):
        kb = 8.617333262145e-5
        result = harmonic_thermodynamics([[.01, .02, .03]], [0, 299.9, 300, 300.1, 1e7])
        self.assertAlmostEqual(result["zero_point_energy_eV"], .03)
        self.assertEqual(result["heat_capacity_eV_per_K"][0], 0)
        self.assertAlmostEqual(result["heat_capacity_eV_per_K"][-1], 3*kb, places=11)
        df = (result["free_energy_eV"][3] - result["free_energy_eV"][1])/.2
        self.assertAlmostEqual(-df, result["entropy_eV_per_K"][2], places=10)
        du = (result["internal_energy_eV"][3] - result["internal_energy_eV"][1])/.2
        self.assertAlmostEqual(du, result["heat_capacity_eV_per_K"][2], places=10)
        zero = harmonic_thermodynamics([[0, .01, .02]], [300])
        self.assertEqual(zero["omitted_zero_mode_weight"], 1)
        with self.assertRaises(ValueError):
            harmonic_thermodynamics([[-.001, .01]], [300])


class AdvancedSimulationTests(unittest.TestCase):
    def test_neb_curved_double_well(self):
        from ase import Atoms
        initial = from_ase(Atoms("H", positions=[[-1, 0, 0]]))
        final = from_ase(Atoms("H", positions=[[1, 0, 0]]))
        result = migration_path(initial, final, CurvedDoubleWell, images=7, fmax=.002, steps=300)
        self.assertTrue(result["converged"])
        self.assertAlmostEqual(result["forward_barrier_eV"], 1, places=5)
        self.assertAlmostEqual(result["images"][3].atoms[0].y, .2, places=3)
        assert_allclose(result["endpoint_max_forces_eV_per_A"], 0)
        calculator = CurvedDoubleWell()
        with self.assertRaises(ValueError):
            migration_path(initial, final, lambda: calculator)

    def test_nvt_temperature_and_reproducibility(self):
        from ase import Atoms
        structure = from_ase(Atoms("Ar100", positions=np.random.default_rng(1).random((100, 3))*10))
        result = nvt_dynamics(structure, FreeParticles(), steps=400, thermostat_fs=10, seed=27, sample_interval=10)
        self.assertLess(abs(result["temperature_K"][10:].mean()/300 - 1), .1)
        repeat = nvt_dynamics(structure, FreeParticles(), steps=10, thermostat_fs=10, seed=27, sample_interval=10)
        assert_allclose(result["velocities_A_per_fs"][:2], repeat["velocities_A_per_fs"])
        self.assertEqual(len(result["frames"]), 41)
        self.assertEqual(result["time_fs"][-1], 400)

    def test_npt_stress_and_cell_evolution(self):
        from ase.build import bulk
        from ase.calculators.emt import EMT
        structure = from_ase(bulk("Cu", "fcc", a=3.6, cubic=True).repeat((2, 2, 2)))
        result = npt_dynamics(structure, EMT(), steps=30, timestep_fs=.5, temperature_K=100,
                              pressure_GPa=5, thermostat_fs=20, barostat_fs=100, sample_interval=5)
        self.assertTrue(np.isfinite(result["pressure_GPa"]).all())
        self.assertTrue(np.isfinite(result["temperature_K"]).all())
        self.assertLess(result["volume_A3"][-1], result["volume_A3"][0])
        self.assertEqual(len(result["frames"]), 7)

    def test_symmetry_path_basis_consistency(self):
        from ase.build import bulk
        structure = from_ase(bulk("Si", "diamond", a=5.43, cubic=True))
        result = reciprocal_path(structure, spacing_inv_A=.1)
        self.assertEqual(result["spacegroup_number"], 227)
        self.assertEqual(len(result["primitive_structure"]), 2)
        assert_allclose(result["kpoints_fractional"] @ result["reciprocal_lattice_inv_A"],
                        result["kpoints_inv_A"], atol=1e-12)
        self.assertIn("GAMMA", result["labels"])
        self.assertTrue(len(result["segments"]) > 1)


if __name__ == "__main__":
    unittest.main()
