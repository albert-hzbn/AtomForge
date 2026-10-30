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
from atomforge.databases import DatabaseError, fetch_cod, fetch_materials_project, search_cod, search_materials_project
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
            self.assertIn("Strains: stacked layer", " ".join(result.stages[0]["notes"]))
            strain = 3 * 3.18 / (4 * 2.4612) - 1
            self.assertIn(f"{100 * strain:.4g}", " ".join(result.stages[0]["notes"]))
            # c holds the stack and 15 A of vacuum.
            self.assertGreaterEqual(s.cell[2][2], z.max() - z.min() + 15 - 1e-9)
            # Unstrained: the graphene keeps its own in-plane coordinates.
            plain = run(f"stack-layers {path} gap=3.4 strain=false", mos2).structure
            np.testing.assert_allclose([a.x for a in plain.atoms[27:]], [a.x for a in graphene.atoms], atol=1e-9)

    def test_stacking_options(self):
        with tempfile.TemporaryDirectory() as folder:
            graphene = run("build-monolayer graphene na=2 nb=2").structure
            hbn = run("build-monolayer hbn na=2 nb=2").structure
            path = str(Path(folder) / "hbn.xyz")
            af.save(hbn, path)
            a_g, a_h = graphene.cell[0][0], hbn.cell[0][0]
            for match, expected in (("layer", a_g), ("substrate", a_h), ("average", (a_g + a_h) / 2)):
                s = run(f"stack-layers {path} match={match}", graphene).structure
                self.assertAlmostEqual(s.cell[0][0], expected, places=9)
                # Both layers keep their fractional in-plane coordinates in the common cell.
                inverse = np.linalg.inv(np.array(s.cell))
                frac = np.array([[a.x, a.y, a.z] for a in s.atoms]) @ inverse
                ref = np.array([[a.x, a.y, a.z] for a in graphene.atoms]) @ np.linalg.inv(np.array(graphene.cell))
                np.testing.assert_allclose(frac[:8, :2] % 1, ref[:, :2] % 1, atol=1e-9)
            stack = run(f"stack-layers {path} copies=3 gap=3.3 shift=0.5,0,0", graphene).structure
            self.assertEqual(len(stack.atoms), 8 + 3 * 8)
            z = sorted({round(a.z, 6) for a in stack.atoms})
            np.testing.assert_allclose(np.diff(z), 3.3, atol=1e-6)
            # Each copy is shifted by half a lattice vector along a relative to the previous one.
            first = np.array([[a.x, a.y] for a in stack.atoms[8:16]])
            second = np.array([[a.x, a.y] for a in stack.atoms[16:24]])
            shift = (second - first) @ np.linalg.inv(np.array(stack.cell)[:2, :2])
            np.testing.assert_allclose((shift[:, 0] - 0.5 + 0.5) % 1 - 0.5, 0, atol=1e-9)
            twisted = run(f"stack-layers {path} twist=21.79 strain=false", graphene)
            self.assertIn("not periodic", " ".join(twisted.stages[0]["notes"]))


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


    def test_searches(self):
        seen = []

        def cod(request, timeout):
            seen.append(request.full_url)
            return self.Response(json.dumps([{"file": "9015662", "formula": "- O2 Ti -", "sg": "P 42/m n m", "a": "4.594", "b": "4.594",
                                              "c": "2.959", "alpha": "90", "beta": "90", "gamma": "90", "vol": "62.4",
                                              "mineral": "Rutile", "title": "Rutile", "year": "2001"}]).encode())
        hits = search_cod("TiO2", elements=["Ti", "O"], strict=True, opener=cod)
        self.assertIn("formula=O2+Ti", seen[0])
        self.assertIn("el1=Ti&el2=O&strictmin=2&strictmax=2", seen[0])
        self.assertEqual(hits[0]["id"], 9015662)
        self.assertEqual(hits[0]["formula"], "O2 Ti")
        self.assertAlmostEqual(hits[0]["a"], 4.594)
        with self.assertRaises(ValueError):
            search_cod(opener=cod)
        with self.assertRaises(ValueError):
            search_cod("Ti$O2", opener=cod)

        def mp(request, timeout):
            seen.append(request.full_url)
            return self.Response(json.dumps({"data": [
                {"material_id": "mp-390", "formula_pretty": "TiO2", "symmetry": {"symbol": "I4_1/amd", "crystal_system": "Tetragonal"},
                 "nsites": 6, "energy_above_hull": 0.006, "band_gap": 2.1, "is_stable": False, "is_magnetic": False},
                {"material_id": "mp-2657", "formula_pretty": "TiO2", "symmetry": {"symbol": "P4_2/mnm", "crystal_system": "Tetragonal"},
                 "nsites": 6, "energy_above_hull": 0.0, "band_gap": 1.8, "is_stable": True, "is_magnetic": False}]}).encode())
        hits = search_materials_project(formula="TiO2", max_energy_above_hull=0.05, api_key="k", opener=mp)
        self.assertIn("formula=TiO2", seen[-1])
        self.assertIn("energy_above_hull_max=0.05", seen[-1])
        self.assertEqual([h["id"] for h in hits], ["mp-2657", "mp-390"])  # most stable first
        self.assertEqual(hits[0]["spacegroup"], "P4_2/mnm")


class AdvancedPerAtomData(unittest.TestCase):
    def test_magnetic_orders(self):
        fe = iron(2)
        layered = run("set-magnetic-moments mode=layered value=2.2 axis=z", fe).structure
        for a in layered.atoms:
            layer = round(a.z / 1.435)
            self.assertEqual(a.magmom, 2.2 if layer % 2 == 0 else -2.2)
        doubled = run("set-magnetic-moments mode=layered value=1 axis=z period=2", fe).structure
        self.assertEqual(sorted({round(a.z / 1.435): a.magmom for a in doubled.atoms}.items())[:4], [(0, 1), (1, 1), (2, -1), (3, -1)])
        g = run("set-magnetic-moments mode=checkerboard value=2", fe).structure
        self.assertEqual(sum(a.magmom for a in g.atoms), 0)
        corner = [a.magmom for a in g.atoms if round(a.x / 1.435) % 2 == 0]
        self.assertEqual(len(set(corner)), 1)  # all corner atoms share a sign in bcc (G-type)
        spins = run("set-magnetic-moments mode=random value=1 seed=3", fe).structure
        self.assertEqual(sum(a.magmom for a in spins.atoms), 0)
        per_element = run("set-magnetic-moments mode=elements moments=Fe,3", fe).structure
        self.assertTrue(all(a.magmom == 3 for a in per_element.atoms))

    def test_charges_per_element_and_neutralising(self):
        salt = af.Structure()
        for k, symbol in enumerate(["Na", "Cl"] * 4):
            salt.add_atom(symbol, 2.8 * k, 0, 0)
        scaled = run("set-charges mode=elements charges=Na,0.9,Cl,-0.8 neutralize=scale", salt).structure
        self.assertAlmostEqual(sum(a.charge for a in scaled.atoms), 0, places=12)
        self.assertEqual({a.symbol: round(a.charge, 9) for a in scaled.atoms}, {"Na": 0.8, "Cl": -0.8})
        shifted = run("set-charges mode=elements charges=Na,0.9,Cl,-0.8 neutralize=shift", salt).structure
        self.assertEqual({a.symbol: round(a.charge, 9) for a in shifted.atoms}, {"Na": 0.85, "Cl": -0.85})

    def test_velocities_without_rotation(self):
        cluster = run("build-cluster element=Cu shells=4 | set-velocities temperature=400 zero_rotation=true seed=2").structure
        p = np.array([[a.x, a.y, a.z] for a in cluster.atoms])
        v = np.array([a.velocity for a in cluster.atoms])
        np.testing.assert_allclose(v.sum(axis=0), 0, atol=1e-12)
        np.testing.assert_allclose(np.cross(p - p.mean(axis=0), v).sum(axis=0), 0, atol=1e-10)
        kinetic = 0.5 * 63.546 * (v ** 2).sum() * 103.6427
        self.assertAlmostEqual(2 * kinetic / ((3 * len(v) - 6) * 8.617333262e-5), 400, places=6)


if __name__ == "__main__":
    unittest.main()
