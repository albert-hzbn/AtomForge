"""Per-atom simulation data (magnetic moments, charges, velocities), layer
stacking and database downloads (with the network mocked)."""
import io
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

if sys.platform == "win32" and os.path.isdir(r"C:\msys64\ucrt64\bin"):
    os.add_dll_directory(r"C:\msys64\ucrt64\bin")  # runtime of a development build
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, os.path.dirname(__file__))
from _optional import missing
import atomforge as af
from atomforge.databases import DatabaseError, fetch_cod, fetch_materials_project
from atomforge.pipeline import Pipeline, PipelineError

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None


def run(text, structure=None):
    return Pipeline.parse(text).run(structure or af.Structure())


def iron(repeat=3):
    s = af.Structure()
    a = 2.87
    for i in range(repeat):
        for j in range(repeat):
            for k in range(repeat):
                s.add_atom("Fe", i * a, j * a, k * a)
                s.add_atom("Fe", (i + .5) * a, (j + .5) * a, (k + .5) * a)
    s.cell = [[a * repeat, 0, 0], [0, a * repeat, 0], [0, 0, a * repeat]]
    return s


@unittest.skipIf(np is None, "numpy is required")
class MomentsChargesVelocities(unittest.TestCase):
    def test_antiferromagnetic_order_by_selection(self):
        # B2-type antiferromagnet: corner atoms up, body centres down.
        result = run("set-magnetic-moments 2.2 | select-expression index % 2 == 1 | set-magnetic-moments -2.2 target=selected", iron())
        moments = [a.magmom for a in result.structure.atoms]
        self.assertEqual(moments[:4], [2.2, -2.2, 2.2, -2.2])
        self.assertAlmostEqual(sum(moments), 0.0, places=12)
        self.assertIn("MAGMOM = 1*2.2 1*-2.2", result.stages[-1]["notes"][0])
        # Expressions see the moments.
        down = run("set-magnetic-moments 2.2 | select-expression index % 2 == 1 | set-magnetic-moments -2.2 target=selected"
                   " | select-expression magmom < 0", iron())
        self.assertEqual(down.selected_count, 27)
        # Moments travel with atoms through edits (replicate, sort, delete).
        edited = run("set-magnetic-moments 1.5 | replicate 2 1 1 | sort-atoms key=x | select-random 0.5 | delete-selected", iron())
        self.assertTrue(all(a.magmom == 1.5 for a in edited.structure.atoms))

    def test_charges_report_neutrality(self):
        from_rocksalt = af.Structure()
        a = 4.2
        for i, (x, y, z) in enumerate([(0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)]):
            from_rocksalt.add_atom("Mg", x * a, y * a, z * a)
            from_rocksalt.add_atom("O", (x + .5) % 1 * a, y * a, z * a)
        from_rocksalt.cell = [[a, 0, 0], [0, a, 0], [0, 0, a]]
        result = run("select-element Mg | set-charges 2 target=selected | select-element O | set-charges -2 target=selected", from_rocksalt)
        self.assertEqual(sorted(a.charge for a in result.structure.atoms), [-2.0] * 4 + [2.0] * 4)
        self.assertIn("total charge 0 e (neutral)", result.stages[-1]["notes"][0])
        self.assertIn("not neutral", run("set-charges 1", from_rocksalt).stages[0]["notes"][0])

    def test_maxwell_boltzmann_velocities(self):
        big = iron(6)  # 432 atoms
        result = run("set-velocities 500 seed=4", big)
        v = np.array([a.velocity for a in result.structure.atoms])
        mass = 55.845
        # No net momentum; kinetic temperature exactly 500 K with 3N - 3 degrees of freedom.
        np.testing.assert_allclose(v.sum(axis=0), 0, atol=1e-12)
        kinetic = 0.5 * mass * (v ** 2).sum() * 103.6427
        temperature = 2 * kinetic / ((3 * len(v) - 3) * 8.617333262e-5)
        self.assertAlmostEqual(temperature, 500, places=6)
        # Components are Gaussian with sigma sqrt(kT/m).
        sigma = math.sqrt(8.617333262e-5 * 500 / (mass * 103.6427))
        self.assertAlmostEqual(v.std(), sigma, delta=0.05 * sigma)
        again = run("set-velocities 500 seed=4", big).structure
        np.testing.assert_allclose([a.velocity for a in again.atoms], v)
        cold = run("set-velocities 500 | set-velocities 0", big).structure
        self.assertTrue(all(a.velocity == (0.0, 0.0, 0.0) for a in cold.atoms))

    def test_python_xyz_round_trip_keeps_data_and_periodicity(self):
        s = run("set-magnetic-moments 2.2 | set-charges 0.5 | set-velocities 300 | set-pbc c=false", iron(2)).structure
        with tempfile.TemporaryDirectory() as folder:
            path = str(Path(folder) / "fe.xyz")
            af.save(s, path)
            back = af.load(path)
            self.assertEqual(back.pbc, [True, True, False])
            for a, b in zip(s.atoms, back.atoms):
                self.assertAlmostEqual(a.magmom, b.magmom, places=9)
                self.assertAlmostEqual(a.charge, b.charge, places=9)
                np.testing.assert_allclose(a.velocity, b.velocity, atol=1e-9)
            text = Path(path).read_text().splitlines()[1]
            self.assertIn("initial_magmoms:R:1:initial_charges:R:1:velo:R:3", text)
        if not missing("ase"):
            from ase.io import read
            with tempfile.TemporaryDirectory() as folder:
                path = str(Path(folder) / "fe.extxyz")
                af.save(s, path.replace(".extxyz", ".xyz"))
                atoms = read(path.replace(".extxyz", ".xyz"))
                np.testing.assert_allclose(atoms.get_initial_magnetic_moments(), [a.magmom for a in s.atoms])
                np.testing.assert_allclose(atoms.get_initial_charges(), [a.charge for a in s.atoms])

    @unittest.skipUnless(os.environ.get("ATOMFORGE_PATH"), "needs the AtomForge executable")
    def test_native_files_keep_the_data(self):
        exe = os.environ["ATOMFORGE_PATH"]
        s = run("set-magnetic-moments 2.2 | set-charges -0.25 | set-velocities 300 seed=2 | rotate axis=1,1,1 angle=17", iron(2)).structure
        with tempfile.TemporaryDirectory() as folder:
            source = str(Path(folder) / "in.xyz")
            af.save(s, source)
            # extended XYZ -> LAMMPS data (charge style, Velocities) -> extended XYZ, through the C++ readers and writers.
            data = str(Path(folder) / "fe.data")
            again = str(Path(folder) / "back.xyz")
            for args in ((source, data), (data, again)):
                done = subprocess.run([exe, "--convert", "--input", args[0], "--output", args[1]], capture_output=True, text=True, timeout=120)
                self.assertEqual(done.returncode, 0, done.stderr)
            text = Path(data).read_text()
            self.assertIn("Atoms  # charge", text)
            self.assertIn("Velocities", text)
            back = af.load(again)
            # LAMMPS rotates the box into its frame: compare frame-free quantities.
            speeds_in = sorted(np.linalg.norm(a.velocity) for a in s.atoms)
            speeds_out = sorted(np.linalg.norm(a.velocity) for a in back.atoms)
            np.testing.assert_allclose(speeds_out, speeds_in, rtol=1e-8)
            self.assertTrue(all(abs(a.charge + 0.25) < 1e-12 for a in back.atoms))
            # CASTEP keeps the moments (SPIN=), extended XYZ everything.
            cell = str(Path(folder) / "fe.cell")
            done = subprocess.run([exe, "--convert", "--input", source, "--output", cell], capture_output=True, text=True, timeout=120)
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertIn("SPIN=2.2", Path(cell).read_text())
            pwi = str(Path(folder) / "fe.pwi")
            subprocess.run([exe, "--convert", "--input", source, "--output", pwi], capture_output=True, text=True, timeout=120)
            text = Path(pwi).read_text()
            self.assertIn("nspin = 2", text)
            self.assertIn("starting_magnetization(1) = 0.44", text)


@unittest.skipIf(np is None, "numpy is required")
class Stacking(unittest.TestCase):
    def test_graphene_on_mos2_heterostructure(self):
        with tempfile.TemporaryDirectory() as folder:
            graphene = run("build-monolayer graphene na=4 nb=4").structure
            path = str(Path(folder) / "graphene.xyz")
            af.save(graphene, path)
            mos2 = run("build-monolayer mx2-2h na=3 nb=3").structure
            result = run(f"stack-layers {path} gap=3.4", mos2)
            s = result.structure
            self.assertEqual(len(s.atoms), 27 + 32)
            z = np.array([a.z for a in s.atoms])
            top_mos2 = z[:27].max()
            self.assertAlmostEqual(z[27:].min() - top_mos2, 3.4, places=9)
            self.assertLess(np.ptp(z[27:]), 1e-9)
            self.assertEqual(result.selected_indices, list(range(27, 59)))
            # Strained to the MoS2 cell: 3 x 3.18 against 4 x 2.4612.
            self.assertIn("strained by", " ".join(result.stages[0]["notes"]))
            strain = 3 * 3.18 / (4 * 2.4612) - 1
            self.assertIn(f"{100 * strain:.4g}", " ".join(result.stages[0]["notes"]))
            # c holds the stack and 15 A of vacuum.
            self.assertGreaterEqual(s.cell[2][2], z.max() - z.min() + 15 - 1e-9)
            # Unstrained: the graphene keeps its own in-plane coordinates.
            plain = run(f"stack-layers {path} gap=3.4 strain=false", mos2).structure
            np.testing.assert_allclose([a.x for a in plain.atoms[27:]], [a.x for a in graphene.atoms], atol=1e-9)


class Databases(unittest.TestCase):
    class Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    def test_cod_download(self):
        cif = """data_test
_cell_length_a 3.615
_cell_length_b 3.615
_cell_length_c 3.615
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Cu1 Cu 0 0 0
Cu2 Cu 0 0.5 0.5
Cu3 Cu 0.5 0 0.5
Cu4 Cu 0.5 0.5 0
"""
        seen = []

        def opener(request, timeout):
            seen.append(request.full_url)
            return self.Response(cif.encode())
        s = fetch_cod(9008468, opener=opener)
        self.assertEqual(seen, ["https://www.crystallography.net/cod/9008468.cif"])
        self.assertEqual(len(s.atoms), 4)
        self.assertAlmostEqual(s.cell[0][0], 3.615)
        with self.assertRaises(ValueError):
            fetch_cod("mp-13", opener=opener)
        with self.assertRaises(DatabaseError):
            fetch_cod(1, opener=lambda request, timeout: self.Response(b"data_empty\n"))

        def failing(request, timeout):
            raise OSError("offline")
        with self.assertRaises(DatabaseError):
            fetch_cod(1, opener=failing)

    def test_materials_project_with_moments(self):
        document = {"data": [{"material_id": "mp-13", "structure": {
            "lattice": {"matrix": [[2.87, 0, 0], [0, 2.87, 0], [0, 0, 2.87]]},
            "sites": [{"species": [{"element": "Fe", "occu": 1}], "xyz": [0, 0, 0], "properties": {"magmom": 2.3}},
                      {"species": [{"element": "Fe", "occu": 1}], "xyz": [1.435, 1.435, 1.435], "properties": {"magmom": -2.3}}]}}]}
        headers = {}

        def opener(request, timeout):
            headers.update(request.headers)
            self.assertIn("material_ids=mp-13", request.full_url)
            return self.Response(json.dumps(document).encode())
        s = fetch_materials_project("mp-13", api_key="secret-key", opener=opener)
        self.assertEqual(headers.get("X-api-key"), "secret-key")
        self.assertEqual([a.magmom for a in s.atoms], [2.3, -2.3])
        self.assertEqual(s.cell[1][1], 2.87)
        with mock.patch.dict(os.environ, {}, clear=True):
            with self.assertRaises(DatabaseError):
                fetch_materials_project("mp-13", opener=opener)
        with self.assertRaises(DatabaseError):
            fetch_materials_project("mp-0", api_key="k", opener=lambda request, timeout: self.Response(b'{"data": []}'))


if __name__ == "__main__":
    unittest.main()
