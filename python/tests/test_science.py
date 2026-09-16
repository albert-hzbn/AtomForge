"""Scientific integration tests with analytic fixtures and established engines."""

from pathlib import Path
import os
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.emt import EMT
import atomforge as af
from atomforge.science import (relax, molecular_dynamics, phonons, powder_diffraction,
    single_crystal, plot_diffraction, read_dos, plot_dos, read_bands, plot_bands,
    read_trajectory, write_trajectory, export_animation, hirshfeld, mode_frames)
from atomforge.electronic import Grid


class HarmonicBond(Calculator):
    implemented_properties = ["energy","forces"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms,properties,system_changes)
        delta=atoms.positions[1]-atoms.positions[0]
        length=np.linalg.norm(delta)
        force=2*(length-1)*delta/length
        self.results={"energy":(length-1)**2,"forces":np.array([force,-force])}


class ScienceTests(unittest.TestCase):
    def test_bandgrid_order_and_surface(self):
        from atomforge.science import read_bandgrid, plot_surface
        path=self.root/"bands.bxsf"
        values=[x/2+2*y/2+4*z/2 for x in range(3) for y in range(3) for z in range(3)]
        path.write_text("BEGIN_INFO\nFermi Energy: 3.5\nEND_INFO\nBEGIN_BLOCK_BANDGRID_3D\nfixture\n"
            "BEGIN_BANDGRID_3D\n1\n3 3 3\n0 0 0\n1 0 0\n0 1 0\n0 0 1\nBAND: 7\n"+
            " ".join(map(str,values))+"\nEND_BANDGRID_3D\nEND_BLOCK_BANDGRID_3D\n")
        grid=read_bandgrid(path)
        vertices=np.asarray(grid.surface(7).vertices)
        self.assertGreater(len(vertices),0)
        np.testing.assert_allclose(vertices[:,0]+2*vertices[:,1]+4*vertices[:,2],3.5,atol=1e-10)
        plot_surface(grid,7,self.root/"surface.png")
        with self.assertRaises(ValueError): grid.surface(7,100)
        path.write_text(path.read_text().replace("END_BANDGRID_3D",""))
        with self.assertRaises(ValueError): read_bandgrid(path)

    def test_symmetry_and_calculation_formats(self):
        path=self.root/"inversion.cif"
        path.write_text("data_test\n_cell_length_a 4\n_cell_length_b 4\n_cell_length_c 4\n"
            "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
            "_space_group_name_H-M_alt 'P -1'\nloop_\n_atom_site_label\n_atom_site_type_symbol\n"
            "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nCu1 Cu .1 .2 .3\n")
        structure=af.load(path)
        self.assertEqual(len(structure),2)
        np.testing.assert_allclose(sorted(a.x for a in structure.atoms),[.4,3.6])
        for extension in ("pwi","gjf","com"):
            with self.subTest(extension=extension):
                target=self.root/("structure."+extension)
                structure.save(target); restored=af.load(target)
                self.assertEqual(len(restored),2)
                np.testing.assert_allclose([[a.x,a.y,a.z] for a in restored.atoms],
                                           [[a.x,a.y,a.z] for a in structure.atoms],atol=1e-6)

    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)

    def dimer(self):
        result=af.Structure(); result.add_atom("H",0,0,0); result.add_atom("H",1.4,0,0)
        return result

    def copper(self):
        result=af.Structure(); result.cell=[[3.6,0,0],[0,3.6,0],[0,0,3.6]]
        for point in [(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0)]:
            result.add_atom("Cu",*(3.6*v for v in point))
        return result

    def test_relaxation_and_energy_conservation(self):
        source=self.dimer()
        result=relax(source,HarmonicBond(),fmax=1e-7)
        self.assertTrue(result.converged)
        self.assertLess(result.energy,1e-12)
        self.assertEqual(source.atoms[1].x,1.4)
        dynamics=molecular_dynamics(source,HarmonicBond(),steps=100,timestep_fs=.1,
                                     velocities=[[0,0,0],[0,0,0]])
        energy=np.asarray(dynamics["total_energies_eV"])
        self.assertLess(np.max(np.abs(energy-energy[0])),2e-5)
        self.assertEqual(len(dynamics["frames"]),101)

    def test_diffraction_extinctions_and_bragg_angle(self):
        source=self.copper()
        pattern=powder_diffraction(source,wavelength=1.5406)
        expected=2*np.degrees(np.arcsin(1.5406*np.sqrt(3)/(2*3.6)))
        self.assertAlmostEqual(pattern.x[0],expected,places=5)
        reflections=single_crystal(source,[(1,0,0),(1,1,1)],[[1]*4,[1]*4])
        self.assertLess(reflections["intensities"][0],1e-20)
        self.assertAlmostEqual(reflections["intensities"][1],16)
        plot_diffraction(pattern,self.root/"diffraction.svg")
        self.assertIn("<svg",(self.root/"diffraction.svg").read_text())

    def test_phonon_acoustic_modes_and_animation(self):
        source=self.copper()
        result=phonons(source,EMT(),[[0,0,0],[.5,0,0]],self.root/"phonons",supercell=(2,2,2))
        self.assertEqual(result["energies_eV"].shape,(2,12))
        self.assertLess(np.max(np.abs(result["energies_eV"][0,:3])),1e-5)
        self.assertTrue(np.isfinite(result["modes"]).all())
        frames=mode_frames(source,result["modes"][0,-1],frames=4)
        export_animation(frames,self.root/"mode.gif",fps=4)
        self.assertEqual((self.root/"mode.gif").read_bytes()[:3],b"GIF")
        write_trajectory(frames,self.root/"frames.extxyz")
        restored=read_trajectory(self.root/"frames.extxyz")
        self.assertEqual(len(restored),4)
        np.testing.assert_allclose(restored[0].cell,source.cell)

    def test_stockholder_conservation(self):
        geometry=((1,0,0),(0,1,0),(0,0,1))
        density=Grid((2,2,2),geometry,[4]*8,periodic=True,unit="e/A^3")
        first=Grid((2,2,2),geometry,[1]*8,periodic=True,unit="e/A^3")
        second=Grid((2,2,2),geometry,[3]*8,periodic=True,unit="e/A^3")
        result=hirshfeld(density,[first,second],reference_electrons=[2,2])
        np.testing.assert_allclose(result["electron_populations"],[1,3])
        np.testing.assert_allclose(result["net_atomic_charges"],[1,-1])
        self.assertAlmostEqual(sum(result["electron_populations"]),density.integrate())
        zero=Grid((2,2,2),geometry,[0]*8,periodic=True,unit="e/A^3")
        with self.assertRaises(ValueError):
            hirshfeld(density,[zero])

    def test_doscar_and_projected_data(self):
        path=self.root/"DOSCAR"
        path.write_text("1 1 1 0\nheader\nheader\nheader\nfixture\n1 -1 3 0 1\n"
            "-1 1 2 0 0\n0 2 3 1 2\n1 1 2 3 5\n1 -1 3 0 1\n"
            "-1 .1 .2\n0 .3 .4\n1 .5 .6\n")
        result=read_dos(path,projected_labels=["s up","s down"])
        np.testing.assert_allclose(result["dos"]["up"],[1,2,1])
        self.assertEqual(result["projected"][0].shape,(3,2))
        plot_dos(result,self.root/"dos.png",sites=[0])
        self.assertEqual((self.root/"dos.png").read_bytes()[:8],b"\x89PNG\r\n\x1a\n")
        path.write_text(path.read_text().rsplit("\n",2)[0])
        with self.assertRaises(ValueError):
            read_dos(path)

    def test_eigenval(self):
        path=self.root/"EIGENVAL"
        path.write_text("1 1 1 1\nheader\nheader\nheader\nfixture\n2 2 2\n\n"
                        "0 0 0 .5\n1 -1 2\n2 1 0\n\n.5 0 0 .5\n1 -.5 2\n2 1.5 0\n")
        result=read_bands(path)
        energies=next(iter(result["energies_eV"].values()))
        np.testing.assert_allclose(energies,[[-1,1],[-.5,1.5]])
        plot_bands(result,self.root/"bands.svg")


if __name__=="__main__":
    unittest.main()
