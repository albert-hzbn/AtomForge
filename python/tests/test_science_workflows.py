"""End-to-end JSON/CLI/desktop-entry tests for the twenty scientific tools."""

from pathlib import Path
import json
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from atomforge.science.workflows import TOOLS, json_result, run_tool, tool_catalog


def requests():
    from ase.build import bulk
    from ase.eos import birchmurnaghan
    from atomforge.science.simulation import from_ase
    crystal = bulk("Cu", "fcc", a=3.6, cubic=True)
    positions = crystal.positions
    cell = crystal.cell.array
    structure = json_result(from_ase(crystal))
    time = np.arange(40.)
    trajectory = np.zeros((40, 2, 3)); trajectory[:, :, 0] = time[:, None]*.01
    velocity = np.zeros_like(trajectory); velocity[:, :, 0] = np.cos(time[:, None]*.4)
    k = np.random.default_rng(44).normal(size=(30, 3))*.03
    volumes = np.linspace(14, 18, 9)
    strains = np.vstack((np.eye(6)*.01, -np.eye(6)*.01))
    calculator = {"module": "ase.calculators.emt", "attribute": "EMT"}
    initial = {"symbols": ["H", "H"], "positions": [[0, 0, 0], [1, 0, 0]]}
    final = {"symbols": ["H", "H"], "positions": [[0, 0, 0], [1.3, 0, 0]]}
    return {
        "msd": {"positions": trajectory, "timestep_fs": 1},
        "diffusion": {"lag_fs": time, "msd_A2": time*.06, "fit_range_fs": [10, 30]},
        "vacf": {"velocities": velocity, "timestep_fs": 1},
        "vibrational-spectrum": {"velocities": velocity, "timestep_fs": 1},
        "local-strain": {"reference": positions, "current": positions*1.01, "cutoff_A": 2.8,
                         "reference_cell": cell, "current_cell": cell*1.01, "pbc": [True]*3},
        "centrosymmetry": {"positions": positions, "cutoff_A": 2.8, "cell": cell, "pbc": [True]*3},
        "bond-order": {"positions": positions, "cutoff_A": 2.8, "cell": cell, "pbc": [True]*3},
        "wigner-seitz": {"reference_sites": positions, "positions": positions, "cell": cell, "pbc": [True]*3},
        "structure-factor": {"positions": positions, "q_vectors": [[0, 0, 0], [1, 0, 0]]},
        "band-gap": {"energies_eV": [[-1, 2], [-.2, 3]], "fermi_eV": 0},
        "effective-mass": {"kpoints_inv_A": k, "energies_eV": np.sum(k*k, axis=1)*3.8099821155369265,
                           "center_inv_A": [0, 0, 0]},
        "work-function": {"distance_A": time, "potential_eV": np.full(40, 6.), "fermi_eV": 1.5, "vacuum_range_A": [20, 30]},
        "equation-of-state": {"volumes_A3": volumes, "energies_eV": birchmurnaghan(volumes, -4, .7, 4, 16)},
        "elastic-tensor": {"strains": strains, "stresses_GPa": strains*100},
        "phonon-dos": {"energies_eV": [[.01, .02, .03]], "energy_grid_eV": np.linspace(-.01, .05, 101)},
        "harmonic-thermodynamics": {"energies_eV": [[.01, .02, .03]], "temperatures_K": [0, 300, 1000]},
        "neb": {"initial": initial, "final": final, "calculator_factory": {"module": "ase.calculators.morse", "attribute": "MorsePotential"},
                "images": 5, "steps": 3},
        "nvt": {"structure": structure, "calculator": calculator, "steps": 3, "sample_interval": 1},
        "npt": {"structure": structure, "calculator": calculator, "steps": 3, "sample_interval": 1},
        "reciprocal-path": {"structure": structure, "spacing_inv_A": .2},
    }


class WorkflowTests(unittest.TestCase):
    def test_all_twenty_desktop_cli_requests(self):
        entry = Path(__file__).resolve().parents[1] / "atomforge/science/_desktop_entry.py"
        cases = requests()
        self.assertEqual(set(cases), set(TOOLS))
        self.assertEqual(len(cases), 20)
        with tempfile.TemporaryDirectory(prefix="atomforge science spaces ") as directory:
            root = Path(directory)
            for name, parameters in cases.items():
                with self.subTest(tool=name):
                    source = root / (name + " request.json")
                    output = root / (name + " output.json")
                    report = root / (name + " report.txt")
                    structures = root / (name + " frames.extxyz")
                    source.write_text(json.dumps(json_result(parameters)), encoding="utf-8")
                    process = subprocess.run([sys.executable, str(entry), name, "--input", str(source),
                                              "--output", str(output), "--report", str(report), "--structures", str(structures)],
                                             capture_output=True, text=True, timeout=45)
                    self.assertEqual(process.returncode, 0, process.stderr)
                    payload = json.loads(output.read_text(encoding="utf-8"))
                    self.assertEqual(payload["tool"], name)
                    self.assertTrue(payload["result"])
                    self.assertIn(TOOLS[name][0], report.read_text(encoding="utf-8"))
                    if name in ("nvt", "npt", "neb", "reciprocal-path"):
                        from ase.io import read
                        frames = read(structures, index=":")
                        expected = {"nvt": 4, "npt": 4, "neb": 5, "reciprocal-path": 1}[name]
                        self.assertEqual(len(frames), expected)
                        if name in ("nvt", "npt"):
                            self.assertTrue(frames[-1].has("momenta"))
            gap = json.loads((root / "band-gap output.json").read_text())["result"]
            self.assertAlmostEqual(gap["gap_eV"], 2.2)
            diffusion = json.loads((root / "diffusion output.json").read_text())["result"]
            self.assertAlmostEqual(diffusion["D_m2_per_s"], 1e-7)

    def test_data_references_and_invalid_parameters(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            np.savetxt(root / "bands.csv", [[-1, 2], [-.2, 3]], delimiter=",")
            result = run_tool("band-gap", {"energies_eV": {"file": "bands.csv"}, "fermi_eV": 0}, base_directory=root)
            self.assertAlmostEqual(result["gap_eV"], 2.2)
            (root / "previous.json").write_text('{"result": {"msd_A2": [0, 0.06, 0.12, 0.18]}}')
            result = run_tool("diffusion", {"lag_fs": [0, 1, 2, 3],
                "msd_A2": {"file": "previous.json", "field": "result.msd_A2"}, "fit_range_fs": [0, 3]}, base_directory=root)
            self.assertAlmostEqual(result["D_A2_per_fs"], .01)
            (root / "structure.json").write_text(json.dumps(requests()["nvt"]["structure"]))
            result = run_tool("reciprocal-path", {"structure": {"file": "structure.json"}}, base_directory=root)
            self.assertEqual(result["spacegroup_number"], 225)
        with self.assertRaises(TypeError):
            run_tool("band-gap", {"energies_eV": [[-1, 1]], "fermi_eV": 0, "unknown": 3})
        with self.assertRaises(ValueError):
            run_tool("arbitrary-python", {})
        self.assertEqual([item["priority"] for item in tool_catalog()], list(range(1, 21)))

    def test_failed_run_preserves_result(self):
        from atomforge.science.__main__ import main
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "request.json"
            output = Path(directory) / "result.json"
            output.write_text("prior result")
            source.write_text('{"energies_eV": [[1, 2]], "fermi_eV": 0}')
            status = main(["band-gap", "--input", str(source), "--output", str(output), "--overwrite"])
            self.assertEqual(status, 1)
            self.assertEqual(output.read_text(), "prior result")


if __name__ == "__main__":
    unittest.main()
