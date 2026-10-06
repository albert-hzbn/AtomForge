"""atomforge.native: the compiled scientific tools called from Python."""

import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import atomforge as af
from atomforge import native


def copper(repeat=2, a=3.615):
    s = af.Structure()
    for i in range(repeat):
        for j in range(repeat):
            for k in range(repeat):
                for f in ((0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)):
                    s.add_atom("Cu", (i + f[0]) * a, (j + f[1]) * a, (k + f[2]) * a)
    s.cell = [[a * repeat, 0, 0], [0, a * repeat, 0], [0, 0, a * repeat]]
    return s


def rocksalt(a=5.64):
    s = af.Structure()
    for f in ((0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0)):
        s.add_atom("Na", f[0] * a, f[1] * a, f[2] * a)
        s.add_atom("Cl", (f[0] + .5) * a, f[1] * a, f[2] * a)
    s.cell = [[a, 0, 0], [0, a, 0], [0, 0, a]]
    return s


class NativeTests(unittest.TestCase):
    def test_catalog_lists_every_tool_with_parameters(self):
        tools = {tool["id"]: tool for tool in native.catalog()}
        for tool in ("structure-type", "dislocation-lines", "powder-xrd", "relax", "lammps-export"):
            self.assertIn(tool, tools)
        relax = {p["name"]: p for p in tools["relax"]["parameters"]}
        self.assertTrue(relax["structure"]["required"])
        self.assertEqual(relax["fmax"]["default"], 0.01)

    def test_structure_objects_in_and_out(self):
        result = native.structure_type(structure=copper())
        self.assertEqual(result["counts"]["fcc"], 32)
        self.assertEqual(result.tool, "structure-type")
        xrd = native.run("powder-xrd", structure=copper(1), two_theta_range_deg=[10, 100])
        self.assertEqual(len(xrd["peaks"]), 5)
        self.assertAlmostEqual(xrd["peaks"][0]["intensity"], 100)

    def test_symmetry_returns_structures(self):
        try:
            sym = native.symmetry(structure=rocksalt(), output_cell="primitive")
        except ValueError as error:
            self.assertIn("spglib", str(error))  # builds without spglib
            return
        self.assertEqual(sym["international_symbol"], "Fm-3m")
        self.assertEqual(sym["wyckoff_letters"][:2], ["a", "b"])
        primitive = sym.frames[-1]
        self.assertIsInstance(primitive, af.Structure)
        self.assertEqual(len(primitive), 2)
        self.assertIsNotNone(primitive.cell)

    def test_structure_lists_are_trajectory_frames(self):
        frames = []
        for step in range(4):
            s = copper(1)
            for atom in s.atoms:
                atom.x += 0.1 * step
            frames.append(s)
        msd = native.msd(positions=frames, timestep_fs=2.0)
        self.assertAlmostEqual(msd["msd_A2"][1], 0.01, places=12)
        self.assertEqual(msd["lag_fs"][1], 2.0)

    def test_dynamics_frames_and_velocities(self):
        md = native.nvt(structure=copper(), calculator={"potential": "EMT"}, steps=20, sample_interval=5, temperature_K=50.0)
        self.assertEqual(len(md.frames), len(md.velocities))
        self.assertEqual(len(md.velocities[0]), 32)
        self.assertEqual(md.times_fs[0], 0.0)

    def test_generated_files_and_file_references(self):
        with tempfile.TemporaryDirectory() as folder:
            exported = native.lammps_export(structure=copper(1), task="nvt")
            written = exported.write_files(folder)
            self.assertEqual(sorted(p.name for p in written), ["data.lammps", "in.lammps"])
            # Data files are referenced relative to base.
            Path(folder, "energies.csv").write_text("0.01,0.02\n", encoding="utf-8")
            thermo = native.run("harmonic-thermodynamics", energies_eV={"file": "energies.csv"}, temperatures_K=[300], base=folder)
            self.assertGreater(thermo["heat_capacity_eV_per_K"][0], 0)

    def test_errors_carry_the_engine_message(self):
        with self.assertRaisesRegex(ValueError, "relax"):
            native.run("relax", structure=copper(), steps=-3)
        with self.assertRaises(ValueError):
            native.run("no-such-tool")
        with self.assertRaises(AttributeError):
            native.not_a_tool


if __name__ == "__main__":
    unittest.main()
