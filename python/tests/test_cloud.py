"""atomforge.cloud: atom clouds of up to billions of atoms (the desktop's large-data viewer)."""

import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import atomforge as af
from atomforge import cloud


@unittest.skipUnless(os.environ.get("ATOMFORGE_PATH"), "needs the AtomForge executable")
class CloudTests(unittest.TestCase):
    def test_generate_info_and_read(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "fcc.afcloud"
            description = cloud.generate(path, "fcc", 3.615, (12, 12, 12), "Cu", second="Ni", chunk=2048)
            self.assertEqual(description.atoms, 4 * 12 ** 3)
            self.assertEqual([s["name"] for s in description.species], ["Cu", "Ni"])
            self.assertGreater(len(description.chunks), 1)
            self.assertAlmostEqual(description.cell[0][0], 12 * 3.615)
            positions, species = cloud.read(path, 0)
            self.assertEqual(len(positions), description.chunks[0]["count"])
            # Positions sit on lattice sites (half the lattice constant apart), within the quantization.
            residual = max(abs(v / 1.8075 - round(v / 1.8075)) * 1.8075 for p in positions for v in p)
            self.assertLess(residual, 2e-3)
            prefix, _ = cloud.read(path, 0, count=100)
            self.assertEqual(len(prefix), 100)
            self.assertTrue(all(abs(a - b) < 1e-9 for p, q in zip(prefix, positions[:100]) for a, b in zip(p, q)))
            self.assertTrue(set(int(s) for s in species) <= {0, 1})

    def test_build_from_dump_and_structure(self):
        with tempfile.TemporaryDirectory() as folder:
            dump = Path(folder) / "run.dump"
            dump.write_text("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n"
                            "ITEM: ATOMS id type x y z\n1 1 1 2 3\n2 2 4 5 6\n")
            description = cloud.build(dump, Path(folder) / "run.afcloud", types=["Fe", "Cr"])
            self.assertEqual(description.atoms, 2)
            self.assertEqual(sorted(s["name"] for s in description.species), ["Cr", "Fe"])
            s = af.Structure()
            s.add_atom("O", 0, 0, 0)
            s.add_atom("H", 0.96, 0, 0)
            s.add_atom("H", -0.24, 0.93, 0)
            self.assertEqual(cloud.from_structure(s, Path(folder) / "water.afcloud").atoms, 3)
            with self.assertRaises(ValueError):
                cloud.info(dump)
            with self.assertRaises(RuntimeError):
                cloud.generate(Path(folder) / "bad.afcloud", "quasicrystal", 3.0, (1, 1, 1), "Cu")


if __name__ == "__main__":
    unittest.main()
