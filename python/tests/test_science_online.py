"""Offline regression checks using pinned upstream scientific example files.

Sources, exact revisions, SHA-256 hashes and upstream licenses accompany the
fixtures. Derived velocities use a stipulated 1 fs interval: these checks test
numerics and units, not a claimed physical diffusion coefficient for that run.
"""

from pathlib import Path
import hashlib
import json
import sys
import tempfile
import unittest
import warnings

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from numpy.testing import assert_allclose
from atomforge.science import *
from atomforge.science.workflows import run_tool, _reference

DATA = Path(__file__).parent / "data" / "science_online"
PMG = DATA / "pymatgen-test-files"
PHONOPY = DATA / "phonopy"


class PublishedDataTests(unittest.TestCase):
    def test_fixture_integrity(self):
        for item in json.loads((DATA / "manifest.json").read_text()):
            with self.subTest(file=item["file"]):
                data = (DATA / item["file"]).read_bytes()
                self.assertEqual(hashlib.sha256(data).hexdigest(), item["sha256"])
                self.assertEqual(len(data), item["bytes"])

    def test_compressed_eigenval_against_pymatgen(self):
        from pymatgen.io.vasp.outputs import Eigenval
        for path in PMG.glob("EIGENVAL*"):
            with self.subTest(file=path.name):
                expected = Eigenval(path)
                parsed = read_bands(path)
                for actual, reference in zip(parsed["energies_eV"].values(), expected.eigenvalues.values()):
                    assert_allclose(actual, reference[:, :, 0], atol=1e-12)
                gap, cbm, vbm, direct = expected.eigenvalue_band_properties
                result = run_tool("band-gap", {"energies_eV": {"file": str(path.resolve())},
                                                "fermi_eV": (cbm + vbm) / 2})
                assert_allclose([result["gap_eV"], result["vbm_eV"], result["cbm_eV"]],
                                [gap, vbm, cbm], atol=1e-10)
                self.assertEqual(abs(result["direct_gap_eV"] - gap) < 1e-8, direct)

    def test_real_band_path_insufficient_for_mass_tensor(self):
        # Real one-k-point and planar data cannot identify a 3D curvature tensor.
        for name in ("EIGENVAL.ispin2.gz", "EIGENVAL_separate_spins.gz"):
            bands = read_bands(PMG / name)
            with self.subTest(file=name), self.assertRaises(ValueError):
                effective_mass(bands["kpoints"], next(iter(bands["energies_eV"].values()))[:, 0],
                               center_inv_A=[0, 0, 0])

    def test_compressed_bands_through_desktop_entry(self):
        import subprocess
        from pymatgen.io.vasp.outputs import Eigenval
        path = PMG / "EIGENVAL.gz"
        gap, cbm, vbm, _ = Eigenval(path).eigenvalue_band_properties
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); request = root / "request.json"
            request.write_text(json.dumps({"energies_eV": {"file": str(path.resolve())}, "fermi_eV": (cbm+vbm)/2}))
            entry = Path(__file__).resolve().parents[1] / "atomforge/science/_desktop_entry.py"
            process = subprocess.run([sys.executable, str(entry), "band-gap", "--input", str(request),
                "--output", str(root/"result.json"), "--report", str(root/"summary.txt")],
                capture_output=True, text=True, timeout=30)
            self.assertEqual(process.returncode, 0, process.stderr)
            result = json.loads((root/"result.json").read_text())["result"]
            self.assertAlmostEqual(result["gap_eV"], gap, places=10)
            self.assertIn("gap_eV", (root/"summary.txt").read_text())

    def test_xdatcar_msd_diffusion_vacf_spectrum(self):
        from pymatgen.io.vasp.outputs import Xdatcar
        from scipy.signal import periodogram
        frames = Xdatcar(PMG / "XDATCAR_traj").structures
        fractional = np.array([frame.frac_coords for frame in frames])
        delta = np.diff(fractional, axis=0)
        # This example is orthogonal; fractional rounding is an independent MIC.
        cell = frames[0].lattice.matrix
        assert_allclose(cell - np.diag(np.diag(cell)), 0, atol=1e-10)
        unwrapped = np.concatenate([fractional[:1], fractional[:1] + np.cumsum(delta - np.rint(delta), axis=0)]) @ cell
        result = run_tool("msd", {"positions": {"file": "XDATCAR_traj"}, "timestep_fs": 1,
            "cell": {"file": "XDATCAR_traj"}, "wrapped": True, "max_lag": 10}, base_directory=PMG)
        expected = [0.] + [np.mean(np.sum((unwrapped[i:] - unwrapped[:-i])**2, axis=2)) for i in range(1, 11)]
        assert_allclose(result["msd_A2"], expected, rtol=1e-10, atol=1e-12)
        fit = diffusion_coefficient(result["lag_fs"], result["msd_A2"], fit_range_fs=[1, 10])
        slope = np.polyfit(np.arange(1, 11), expected[1:], 1)[0]
        self.assertAlmostEqual(fit["D_m2_per_s"], slope / 6 * 1e-5, places=15)
        velocities = np.diff(unwrapped, axis=0)  # A/fs at assumed dt = 1 fs
        correlation = velocity_autocorrelation(velocities, 1, max_lag=10)["vacf"]
        brute = [np.mean(np.sum(velocities[:len(velocities)-i] * velocities[i:], axis=2)) for i in range(11)]
        assert_allclose(correlation, brute, rtol=1e-11)
        spectrum = vibrational_spectrum(velocities, 1)
        frequency, power = periodogram(velocities, fs=1000, window=np.hanning(len(velocities)),
                                        detrend="constant", axis=0)
        power = power.sum(axis=(1, 2)); power /= np.trapezoid(power, frequency)
        assert_allclose(spectrum["frequency_THz"], frequency)
        assert_allclose(spectrum["density_per_THz"], power, rtol=1e-10, atol=1e-12)

    def test_published_copper_local_structure_and_diffraction(self):
        from ase.io import read
        atoms = read(PHONOPY / "POSCAR-05", format="vasp")
        xyz, cell = atoms.positions, atoms.cell.array
        options = {"cell": cell, "pbc": [True]*3}
        cutoff = np.linalg.norm(cell[0]) * .8
        assert_allclose(centrosymmetry(xyz, cutoff, **options)["centrosymmetry_A2"], 0, atol=1e-20)
        order = bond_order(xyz, cutoff, **options)
        assert_allclose(order["order"]["q4"], np.sqrt(7/192), atol=1e-12)
        assert_allclose(order["order"]["q6"], np.sqrt(169/512), atol=1e-12)
        deformation = np.array([[1.02, .03, 0], [.01, .98, .02], [0, .01, 1.04]])
        strain = local_strain(xyz, xyz @ deformation.T, cutoff, reference_cell=cell,
                             current_cell=cell @ deformation.T, pbc=[True]*3)
        assert_allclose(strain["deformation_gradient"], np.tile(deformation, (len(atoms), 1, 1)), atol=1e-12)
        assert_allclose(strain["d2min_A2"], 0, atol=1e-20)
        defects = wigner_seitz(xyz, np.vstack([xyz[1:], xyz[1] + [.01, 0, 0]]), **options)
        self.assertEqual((defects["vacancies"], defects["interstitial_excess"]), (1, 1))
        q = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 1], [2, 0, 0]]) @ (2*np.pi*np.linalg.inv(cell).T)
        assert_allclose(static_structure_factor(xyz, q)["S_q"], [4, 0, 4, 4], atol=1e-12)

    def test_published_energy_volume_curve(self):
        from ase.eos import EquationOfState
        table = np.loadtxt(PHONOPY / "e-v.dat")
        actual = run_tool("equation-of-state", {"volumes_A3": {"file": "e-v.dat", "column": 0},
            "energies_eV": {"file": "e-v.dat", "column": 1}}, base_directory=PHONOPY)
        volume, energy, bulk = EquationOfState(*table.T, eos="birchmurnaghan").fit()
        assert_allclose([actual["volume_A3"], actual["energy_eV"], actual["bulk_modulus_GPa"]],
                        [volume, energy, bulk*160.2176634], rtol=1e-7)
        self.assertLess(actual["rms_error_eV"], .0001)

    def test_real_stress_strain_against_independent_elastic_fit(self):
        from pymatgen.analysis.elasticity.strain import Strain
        from pymatgen.analysis.elasticity.stress import Stress
        from pymatgen.analysis.elasticity.elastic import ElasticTensor
        data = json.loads((PMG / "Sn_def_stress.json").read_text())
        tensors = np.array([Strain.from_deformation(d) for d in data["deformations"]])
        # Upstream VASP stress uses kbar, positive compression: convert to GPa.
        stress = -.1 * np.array(data["stresses"])
        stress = (stress + stress.transpose(0, 2, 1))/2
        strains = np.array([Strain(x).voigt for x in tensors])
        stresses = np.array([Stress(x).voigt for x in stress])
        actual = elastic_tensor(strains, stresses)
        # Center both to eliminate prestress before an independent pseudoinverse.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            reference = ElasticTensor.from_pseudoinverse(tensors-tensors.mean(axis=0), stress-stress.mean(axis=0))
        matrix = (reference.voigt + reference.voigt.T)/2
        assert_allclose(actual["stiffness_GPa"], matrix, atol=1e-10)
        self.assertTrue(actual["stable_zero_prestress"])
        symmetric = ElasticTensor.from_voigt(matrix)
        self.assertAlmostEqual(actual["moduli"]["bulk_hill_GPa"], symmetric.k_vrh, places=9)
        self.assertAlmostEqual(actual["moduli"]["shear_hill_GPa"], symmetric.g_vrh, places=9)

    def test_published_phonon_mesh_dos_and_thermodynamics(self):
        from ruamel.yaml import YAML
        from scipy.constants import physical_constants
        from scipy.stats import norm
        mesh = YAML(typ="safe").load((PHONOPY / "mesh.yaml").read_text())
        h = physical_constants["Planck constant in eV/Hz"][0]
        energy = np.array([[b["frequency"] for b in q["band"]] for q in mesh["phonon"]])*h*1e12
        weights = np.array([q["weight"] for q in mesh["phonon"]]); weights = weights/weights.sum()
        grid = np.linspace(-.01, .06, 1501)
        dos = phonon_dos(energy, grid, sigma_eV=.0007, weights=weights)
        expected = sum(w * norm.pdf(grid[:, None], loc=row, scale=.0007).sum(axis=1) for w, row in zip(weights, energy))
        assert_allclose(dos["dos_per_eV"], expected, rtol=1e-12, atol=1e-12)
        self.assertAlmostEqual(dos["enclosed_modes"], 3, places=7)
        thermal = harmonic_thermodynamics(energy, [0, 100, 300, 1000], weights=weights)
        # Independent partition function from explicitly summing oscillator levels.
        kb = physical_constants["Boltzmann constant in eV/K"][0]
        for i, temperature in enumerate([100, 300, 1000], 1):
            levels = energy[..., None] * (np.arange(5000) + .5)
            probabilities = np.exp(-(levels-levels[..., :1])/(kb*temperature))
            partition = probabilities.sum(axis=2)
            probabilities /= partition[..., None]
            mean = (probabilities*levels).sum(axis=2)
            variance = (probabilities*(levels-mean[..., None])**2).sum(axis=2)
            free = np.sum(weights[:, None]*(energy/2-kb*temperature*np.log(partition)))
            assert_allclose(thermal["free_energy_eV"][i], free, rtol=1e-10)
            assert_allclose(thermal["internal_energy_eV"][i], np.sum(weights[:, None]*mean), rtol=1e-10)
            assert_allclose(thermal["heat_capacity_eV_per_K"][i], np.sum(weights[:, None]*variance)/(kb*temperature**2), rtol=1e-10)

    def test_real_potential_profile_rejects_vacuum_assumption(self):
        from pymatgen.io.vasp.outputs import Locpot
        potential = Locpot.from_file(PMG / "LOCPOT.gz")
        profile = np.asarray(potential.get_average_along_axis(2))
        distance = np.arange(len(profile))*potential.structure.lattice.c/len(profile)
        # This is bulk, not a slab. Test arithmetic and non-flat warning; do not
        # claim a physically meaningful work function from a nonexistent vacuum.
        result = work_function(distance, profile, 0, vacuum_range_A=distance[[0, 7]])
        self.assertAlmostEqual(result["work_function_eV"], profile[:8].mean(), places=12)
        self.assertFalse(result["flat_vacuum"])

    def test_real_structures_reciprocal_and_dynamics(self):
        from ase.io import read
        from ase.calculators.emt import EMT
        from ase import units
        copper = from_ase(read(PHONOPY / "POSCAR-05", format="vasp").repeat((2, 2, 2)))
        silicon = run_tool("reciprocal-path", {"structure": {"file": str((PHONOPY / "POSCAR-unitcell").resolve())}})
        self.assertEqual(silicon["spacegroup_number"], 227)
        assert_allclose(silicon["kpoints_fractional"] @ silicon["reciprocal_lattice_inv_A"], silicon["kpoints_inv_A"], atol=1e-12)
        for function in (nvt_dynamics, npt_dynamics):
            with self.subTest(ensemble=function.__name__):
                result = function(copper, EMT(), steps=50, timestep_fs=.5, sample_interval=5, temperature_K=100)
                self.assertEqual(len(result["frames"]), 11)
                for frame, velocity, temperature in zip(result["frames"], result["velocities_A_per_fs"], result["temperature_K"]):
                    atoms = to_ase(frame); atoms.set_velocities(velocity/units.fs)
                    self.assertAlmostEqual(atoms.get_temperature(), temperature, places=9)
                self.assertTrue(np.isfinite(result["total_energy_eV"]).all())
                if function is nvt_dynamics:
                    assert_allclose(result["volume_A3"], result["volume_A3"][0])

    def test_real_neb_endpoints_and_nonconvergence(self):
        result = run_tool("neb", {"initial": {"file": "Cu_slab_init.cif"},
            "final": {"file": "Cu_slab_fin.cif"}, "images": 5, "steps": 1, "fmax": 1e-8,
            "calculator_factory": {"module": "ase.calculators.emt", "attribute": "EMT"}}, base_directory=PMG)
        self.assertFalse(result["converged"])
        self.assertEqual(result["steps"], 1)
        self.assertTrue(np.isfinite(result["energies_eV"]).all())
        self.assertGreater(max(result["endpoint_max_forces_eV_per_A"]), .01)
        self.assertAlmostEqual(result["forward_barrier_eV"]-result["reverse_barrier_eV"], result["reaction_energy_eV"])


class FileReferenceTests(unittest.TestCase):
    def test_boolean_strings_are_not_silently_true(self):
        velocities = np.ones((5, 2, 3))
        for function, kwargs in ((mean_square_displacement, {"wrapped": "false"}),
                                 (velocity_autocorrelation, {"normalize": "false"}),
                                 (vibrational_spectrum, {"remove_mean": "false"})):
            with self.subTest(tool=function.__name__), self.assertRaisesRegex(ValueError, "boolean"):
                function(velocities, 1, **kwargs)

    def test_all_vacancies(self):
        result = wigner_seitz([[0, 0, 0], [1, 0, 0]], [])
        assert_allclose(result["occupancy"], [0, 0])
        self.assertEqual(result["vacancies"], 2)
        self.assertEqual(result["interstitial_excess"], 0)
        self.assertEqual(result["site_index"].dtype.kind, "i")

    def test_nonperiodic_boxes_and_unsupported_partial_periodicity(self):
        from ase import Atoms
        atoms = Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]], cell=np.eye(3)*5, pbc=False)
        self.assertIsNone(from_ase(atoms).cell)
        atoms.pbc = [True, True, False]
        with self.assertRaisesRegex(ValueError, "periodic"):
            from_ase(atoms)

    def test_trajectory_timing_must_match_analysis(self):
        from ase.io import write
        from ase import Atoms
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            frames = [Atoms("H", positions=[[0, 0, i]]) for i in range(4)]
            for frame, time in zip(frames, [0, 10, 20, 25]): frame.info["time_fs"] = time
            write(root / "frames.extxyz", frames)
            request = {"positions": {"file": "frames.extxyz"}, "timestep_fs": 10}
            with self.assertRaisesRegex(ValueError, "uniform"):
                run_tool("msd", request, base_directory=root)
            frames[-1].info["time_fs"] = 30
            write(root / "frames.extxyz", frames)
            assert_allclose(run_tool("msd", request, base_directory=root)["msd_A2"], [0, 1, 4, 9])

    def test_cli_prevents_output_collisions(self):
        from atomforge.science.__main__ import main
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); source = root / "request.json"; report = root / "report.txt"
            request = '{"energies_eV":[[-1,2]],"fermi_eV":0}'
            source.write_text(request); report.write_text("keep report")
            common = ["band-gap", "--input", str(source), "--output", str(root / "result.json")]
            self.assertEqual(main(common + ["--report", str(report)]), 1)
            self.assertEqual(report.read_text(), "keep report")
            self.assertEqual(main(common + ["--report", str(source), "--overwrite"]), 1)
            self.assertEqual(source.read_text(), request)
            self.assertFalse((root / "result.json").exists())

    def test_column_selection_all_table_formats(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); table = np.arange(12).reshape(4, 3)
            np.save(root / "data.npy", table)
            np.savetxt(root / "data.csv", table, delimiter=",")
            (root / "data.json").write_text(json.dumps({"samples": table.tolist()}))
            for name in ("data.npy", "data.csv", "data.json"):
                options = {"file": name, "column": 1}
                if name.endswith("json"): options["field"] = "samples"
                assert_allclose(_reference(options, "energies_eV", root), table[:, 1])
                for bad in (-1, 1.5, True, 10):
                    with self.subTest(file=name, column=bad), self.assertRaises(ValueError):
                        _reference(dict(options, column=bad), "energies_eV", root)

    def test_velocity_and_cell_file_defaults(self):
        from ase.io import write
        from ase import Atoms, units
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            frames = [Atoms("Cu2", positions=[[0, 0, 0], [2, 0, 0]], cell=np.eye(3)*6, pbc=True) for _ in range(5)]
            for i, frame in enumerate(frames): frame.set_velocities(np.ones((2, 3))*(i+1))
            write(root / "frames.extxyz", frames)
            reference = {"file": "frames.extxyz"}
            assert_allclose(_reference(reference, "velocities", root), np.array([f.get_velocities() for f in frames])*units.fs)
            assert_allclose(_reference(reference, "cell", root), frames[0].cell.array)
            assert_allclose(_reference(reference, "masses", root), frames[0].get_masses())

    def test_numpy_integer_parameters(self):
        result = bond_order([[0, 0, 0], [1, 0, 0]], 2, degrees=np.array([4, 6]))
        assert_allclose(result["order"]["q4"], 1)

    def test_npt_ideal_gas_equilibrium(self):
        from ase import Atoms, units
        from ase.calculators.calculator import Calculator, all_changes
        class IdealGas(Calculator):
            implemented_properties = ["energy", "forces", "stress"]
            def calculate(self, atoms=None, properties=None, system_changes=all_changes):
                super().calculate(atoms, properties, system_changes)
                self.results = {"energy": 0., "forces": np.zeros((len(atoms), 3)), "stress": np.zeros(6)}
        atoms = Atoms("Ar100", positions=np.random.default_rng(10).random((100, 3))*8,
                      cell=np.eye(3)*8, pbc=True)
        result = npt_dynamics(from_ase(atoms), IdealGas(), steps=4000, timestep_fs=1,
            temperature_K=300, pressure_GPa=1, thermostat_fs=20, barostat_fs=100, sample_interval=10, seed=23)
        # NPT ideal-gas distribution has <V>=(N+1)kBT/P. Discard equilibration.
        volume = (len(atoms)+1)*units.kB*300/units.GPa
        self.assertLess(abs(np.mean(result["volume_A3"][100:])/volume-1), .05)
        self.assertLess(abs(np.mean(result["temperature_K"][100:])/300-1), .05)
        self.assertLess(abs(np.mean(result["pressure_GPa"][100:])-1), .05)


if __name__ == "__main__":
    unittest.main()
