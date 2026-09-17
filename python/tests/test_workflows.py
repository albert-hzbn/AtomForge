"""Cross-interface scientific reference cases and recipe regression tests."""

import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import atomforge as af
from atomforge.electronic import Grid, run_pipeline, batch_process, load_volume
from atomforge.electronic.__main__ import main


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

    def test_native_analyses(self):
        structure = af.Structure()
        structure.add_atom("O", 0, 0, 0)
        structure.add_atom("H", 1, 0, 0)
        structure.add_atom("C", 0, 1, 0)
        result = af.adf(structure, centre="O", neighbor_a="H", neighbor_b="C", cutoff=1.1)
        self.assertEqual(sum(row["raw_count"] for row in result), 1)
        self.assertEqual(result[90]["raw_count"], 1)
        result = af.rdf(structure, cutoff=1.1, bins=11, raw=True, no_pbc=True)
        self.assertEqual(sum(row["raw_count"] for row in result), 4)
        self.assertEqual(len(af.cna(structure, cutoff=1.1)), 3)
        with self.assertRaises(ValueError):
            af.rdf(structure, misspelled_parameter=True)

    def test_molecular_formats(self):
        source=af.Structure()
        source.add_atom("O",0,0,0); source.add_atom("H",.75,.58,0); source.add_atom("H",-.75,.58,0)
        for extension in ("sdf","mol","mol2"):
            with self.subTest(extension=extension):
                path=self.root/("water."+extension)
                source.save(path); restored=af.load(path)
                self.assertEqual(len(restored),3)
                for a,b in zip(source.atoms,restored.atoms):
                    self.assertEqual(a.symbol,b.symbol)
                    self.assertAlmostEqual(a.x,b.x,places=4)
                    self.assertAlmostEqual(a.y,b.y,places=4)

    def test_cli_rejects_unknown_options(self):
        source=self.root/"atom.xyz"
        source.write_text("1\nfixture\nCu 0 0 0\n")
        result=subprocess.run([os.environ["ATOMFORGE_PATH"],"--analyze","cna","--input",str(source),
            "--output",str(self.root/"result.csv"),"--cutof","3"],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn("Unsupported option",result.stderr)

    def test_metadata_round_trip_and_sculpt(self):
        source=af.Structure(); source.cell=[[4,0,0],[0,4,0],[0,0,4]]
        source.add_atom("Cu",1,1,1); source.add_atom("Ni",3,3,3)
        for index,atom in enumerate(source.atoms):
            atom.grain_color=(.1,.2,.3); atom.grain_region=10+index
        path=self.root/"source.vasp"; source.save(path)
        loaded=af.load(path)
        self.assertEqual([a.grain_region for a in loaded.atoms],[10,11])
        result=af.sculpt(source,[(1,0,0,-1.1,0),(0,1,0,-1.1,0),(0,0,1,-1.1,0)])
        self.assertEqual(len(result),1)
        self.assertEqual(result.atoms[0].grain_region,10)

    def test_stacking_fault_sequence(self):
        source=af.Structure(); source.cell=[[3.6,0,0],[0,3.6,0],[0,0,3.6]]
        for x,y,z in ((0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0)):
            source.add_atom("Cu",3.6*x,3.6*y,3.6*z)
        directory=self.root/"frames"
        result=af.build("stacking-fault",["--plane","1","--layers","6","--interval",".5",
            "--maximum","1","--sequence",str(directory)],source=source)
        paths=list(directory.glob("*.vasp"))
        self.assertGreaterEqual(len(paths),3)
        self.assertGreater(len(result.structure),0)
        self.assertTrue(all(len(af.load(path))==len(result.structure) for path in paths))

    def test_interface_same_lattice_is_identity(self):
        structure = af.Structure()
        structure.cell = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
        structure.add_atom("Cu", .6, .9, 0)
        source = self.root / "source.vasp"
        structure.save(source)
        result = af.build("interface", ["--layerA", str(source), "--layerB", str(source),
            "--nmax", "1", "--maxcells", "1", "--gap", "2", "--vacuum", "4"])
        self.assertEqual(len(result.structure), 2)
        atoms = result.structure.atoms
        self.assertAlmostEqual(atoms[0].x, atoms[1].x, places=5)
        self.assertAlmostEqual(atoms[0].y, atoms[1].y, places=5)
        self.assertAlmostEqual(atoms[1].z-atoms[0].z, 2, places=5)
        with self.assertRaises(subprocess.CalledProcessError):
            af.build("interface", ["--layerA", str(source), "--layerB", str(source), "--pick", "-1"])

    def test_recipe_conservation_and_batch(self):
        grid = Grid((4,4,4), ((2,0,0),(0,2,0),(0,0,2)), [2]*64,
                    periodic=True, unit="e/A^3")
        steps = [{"operation":"scale", "parameters":{"factor":3}, "name":"scaled"},
                 {"operation":"subtract", "reference":"input", "name":"difference"},
                 {"operation":"integrate", "name":"electrons"}]
        result = run_pipeline(grid, steps)
        self.assertAlmostEqual(result["electrons"], 32)
        self.assertAlmostEqual(grid.integrate(), 16)
        # Constant density has zero gradient, so the pipeline-dispatched
        # reduced_density_gradient operation is a well-defined zero field
        # (unlike dori(), which divides by the also-zero gradient norm).
        bonding = run_pipeline(grid, [{"operation":"reduced_density_gradient", "parameters":{"floor":1e-12}, "name":"rdg"}])
        self.assertEqual(bonding["rdg"].values, [0]*64)
        with self.assertRaises(ValueError):
            run_pipeline(grid, [{"operation":"__getattribute__"}])
        with self.assertRaises(ValueError):
            run_pipeline(grid, [{"operation":"scale", "parameters":{"factor":1}, "name":"input"}])
        paths = [self.root / "first.xsf", self.root / "second.xsf"]
        for path in paths:
            grid.save(path)
        report = batch_process(paths, steps, self.root / "batch")
        self.assertTrue(all(item["success"] for item in report), report)
        for path in paths:
            folder = self.root / "batch" / path.stem
            self.assertAlmostEqual(float((folder / "result.csv").read_text().strip()), 32)
            self.assertEqual(json.loads((folder / "recipe.json").read_text())["steps"], steps)
        repeat = batch_process(paths, steps, self.root / "batch")
        self.assertTrue(all(not item["success"] for item in repeat))

    def test_electronic_cli_reference(self):
        grid = Grid((3,3,3), ((1,0,0),(0,1,0),(0,0,1)), [2]*27, unit="e/A^3")
        path = self.root / "input.cube"
        grid.save(path)
        output = self.root / "difference.cube"
        self.assertEqual(main([str(path), "--operation", "subtract", "--reference", str(path),
                               "--output", str(output), "--quantity", "density"]), 0)
        self.assertAlmostEqual(load_volume(output).fields[0].integrate(), 0)

    def test_electronic_cli_bonding_topology_and_bader(self):
        # rho is an exact quadratic form; the underlying native operations
        # are already proven analytically correct in test_electronic.py, so
        # this exercises the CLI wiring itself -- argument parsing, operation
        # dispatch and output persistence -- by checking the CLI's output
        # matches the same operation called directly through the Python API
        # on the same (file-roundtripped) density.
        n = 6
        values = [5.0 + 0.01 * x * x + 0.02 * y * y + 0.03 * z * z
                  for z in range(n) for y in range(n) for x in range(n)]
        grid = Grid((n, n, n), ((4, 0, 0), (0, 4, 0), (0, 0, 4)), values, unit="e/A^3")
        path = self.root / "density.xsf"
        grid.save(path)
        reloaded = load_volume(path, "density").fields[0]

        for operation, expected in (
            ("dori", reloaded.dori()),
            ("reduced_density_gradient", reloaded.reduced_density_gradient()),
            ("signed_density", reloaded.signed_density()),
        ):
            output = self.root / (operation + ".xsf")
            self.assertEqual(main([str(path), "--operation", operation, "--output", str(output),
                                   "--quantity", "density"]), 0)
            restored = load_volume(output).fields[0]
            self.assertEqual(len(restored.values), len(expected.values))
            for actual, want in zip(restored.values, expected.values):
                self.assertAlmostEqual(actual, want, places=6)

        betti_output = self.root / "betti.csv"
        self.assertEqual(main([str(path), "--operation", "betti_curve",
                               "--parameters", json.dumps({"thresholds": [5.0, 5.5, 6.0]}),
                               "--output", str(betti_output), "--quantity", "density"]), 0)
        with open(betti_output, newline="", encoding="utf-8") as stream:
            rows = [tuple(map(float, row)) for row in csv.reader(stream)]
        expected_rows = reloaded.betti_curve([5.0, 5.5, 6.0])
        self.assertEqual(len(rows), len(expected_rows))
        for actual, want in zip(rows, expected_rows):
            for a, w in zip(actual, want):
                self.assertAlmostEqual(a, w, places=6)
        self.assertGreater(rows[0][1], 0)  # sanity: some threshold is non-trivially solid

        bader_output = self.root / "bader.csv"
        self.assertEqual(main([str(path), "--operation", "bader_partition", "--output", str(bader_output),
                               "--quantity", "density"]), 0)
        with open(bader_output, newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        partition = reloaded.bader_partition()
        self.assertEqual(len(rows), partition.num_basins)
        self.assertGreater(partition.num_basins, 0)
        for i, row in enumerate(rows):
            basin = partition.basin(i)
            self.assertEqual(int(row["basin"]), i)
            self.assertAlmostEqual(float(row["charge_e"]), basin["charge"], places=6)
            self.assertAlmostEqual(float(row["volume_A3"]), basin["volume"], places=6)
            for axis, key in enumerate(("max_x_A", "max_y_A", "max_z_A")):
                self.assertAlmostEqual(float(row[key]), basin["maximum"][axis], places=6)

        # Missing --output is a user error, not a crash: main() reports failure via its return code.
        self.assertEqual(main([str(path), "--operation", "bader_partition", "--quantity", "density"]), 1)


if __name__ == "__main__":
    unittest.main()
