"""Regression cases for file integrity and viewer failure cleanup."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import atomforge as af
from atomforge._io import _cart_to_frac


class FileIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def structure(self):
        s = af.Structure().set_cell(3, 4, 5, 80, 75, 70)
        s.add_atom("Zn", .5, 1, 1.5)
        s.add_atom("Fe", 1, 2, 3)
        return s

    def test_extensionless_vasp_roundtrip(self):
        for name in ("POSCAR", "CONTCAR"):
            path = self.root / name
            self.structure().save(path)
            loaded = af.load(path)
            self.assertEqual(len(loaded), 2)
            self.assertAlmostEqual(loaded.atoms[0].x, .5)

    def test_nonperiodic_cif_roundtrip(self):
        s = self.structure()
        s.cell = None
        path = self.root / "cluster.cif"
        s.save(path)
        loaded = af.load(path)
        self.assertIsNone(loaded.cell)
        self.assertEqual(len(loaded), 2)
        self.assertEqual(loaded.atoms[1].z, 3)

    def test_cif_blank_lines_between_atoms(self):
        path = self.root / "atoms.cif"
        self.structure().save(path)
        path.write_text(path.read_text().replace("Fe  Fe1", "\n# separator\nFe  Fe1"))
        self.assertEqual(len(af.load(path)), 2)

    def test_cif_wrapped_rows(self):
        path = self.root / "wrapped.cif"
        self.structure().save(path)
        path.write_text(path.read_text().replace("Zn  Zn1  ", "Zn\nZn1\n"))
        self.assertEqual(len(af.load(path)), 2)

    def test_invalid_cif_cell_rejected(self):
        path = self.root / "invalid.cif"
        self.structure().save(path)
        path.write_text(path.read_text().replace("70.0000", "0.0000"))
        with self.assertRaises(ValueError):
            af.load(path)

    def test_xyz_truncation_rejected(self):
        path = self.root / "bad.xyz"
        for text in ("2\ncomment\nFe 0 0 0\n", "-1\ncomment\n", "1\ncomment\nFe 0\n"):
            path.write_text(text)
            with self.assertRaises(ValueError):
                af.load(path)

    def test_vasp_count_mismatch_rejected(self):
        path = self.root / "bad.vasp"
        path.write_text("test\n1\n1 0 0\n0 1 0\n0 0 1\nFe Cu\n1\nDirect\n0 0 0\n")
        with self.assertRaises(ValueError):
            af.load(path)

    def test_singular_cell_rejected(self):
        with self.assertRaises(ValueError):
            _cart_to_frac(1, 2, 3, [[1, 0, 0], [2, 0, 0], [0, 0, 1]])

    def test_lammps_styles_and_image_flags(self):
        path = self.root / "styles.lmp"
        cases = [("atomic", "1 1 11 22 33 0 0 0 # Fe"),
                 ("full", "1 7 1 0.0 11 22 33 0 0 0"),
                 ("charge", "1 1 -0.5 11 22 33"),
                 ("molecular", "1 7 1 11 22 33")]
        for style, row in cases:
            with self.subTest(style=style):
                path.write_text(f"test\n\n1 atoms\n1 atom types\n10 20 xlo xhi\n20 30 ylo yhi\n30 40 zlo zhi\n\nMasses\n\n# comment\n1 55.845 # Fe\n\nAtoms # {style}\n\n{row}\n")
                s = af.load(path)
                self.assertEqual(len(s), 1)
                self.assertEqual(s.atoms[0].symbol, "Fe")
                self.assertEqual((s.atoms[0].x, s.atoms[0].y, s.atoms[0].z), (1, 2, 3))

    def test_triclinic_lammps_geometry_and_masses(self):
        s = self.structure()
        # Rotate cell and atoms away from the conventional basis.
        s.cell = [[-y, x, z] for x, y, z in s.cell]
        for a in s.atoms:
            a.x, a.y = -a.y, a.x
        path = self.root / "skew.lmp"
        s.save(path)
        loaded = af.load(path)
        self.assertIn("65.380", path.read_text())
        for before, after in zip(s.atoms, loaded.atoms):
            f = _cart_to_frac(before.x, before.y, before.z, s.cell)
            g = _cart_to_frac(after.x, after.y, after.z, loaded.cell)
            for a, b in zip(f, g):
                self.assertAlmostEqual(a, b, places=5)
        for a in range(3):
            for b in range(3):
                self.assertAlmostEqual(sum(x*y for x, y in zip(s.cell[a], s.cell[b])),
                                       sum(x*y for x, y in zip(loaded.cell[a], loaded.cell[b])), places=5)

    def test_cluster_lammps_box_encloses_atoms(self):
        s = self.structure()
        s.cell = None
        s.translate(-10, 20, 30)
        path = self.root / "cluster.lmp"
        s.save(path)
        loaded = af.load(path)
        for atom in loaded.atoms:
            for value, length in zip((atom.x, atom.y, atom.z), (loaded.cell[i][i] for i in range(3))):
                self.assertTrue(0 < value < length)

    def test_viewer_serialization_failure_removes_tempfile(self):
        from atomforge._viewer import view
        paths = []
        def fail(s, path):
            paths.append(path)
            raise OSError("write failed")
        with patch("atomforge._viewer._find_atomforge", return_value="viewer"), patch("atomforge._io._save_xyz", side_effect=fail):
            with self.assertRaises(OSError):
                view(self.structure())
        self.assertEqual(len(paths), 1)
        self.assertFalse(Path(paths[0]).exists())


if __name__ == "__main__":
    unittest.main()
