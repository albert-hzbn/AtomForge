"""Optional tests against official pymatgen reference data.

Set ATOMFORGE_SCIENCE_FIXTURES to a pymatgen-test-files checkout, revision
85005108e99e40cbd59e78f067459b0ca664efc6 (MIT licensed). No test downloads data.
These check real file integration; they do not claim to execute Chargemol.
"""

import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import atomforge as af
from atomforge.science import bader, ddec, orbital_density, read_projections, plot_bands
from atomforge.electronic import Grid, load_volume


@unittest.skipUnless(os.environ.get("ATOMFORGE_SCIENCE_FIXTURES"),"External fixtures not configured")
class ExternalScienceTests(unittest.TestCase):
    def setUp(self):
        self.fixtures=Path(os.environ["ATOMFORGE_SCIENCE_FIXTURES"])
        self.directory=tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root=Path(self.directory.name)

    @unittest.skipUnless(os.environ.get("ATOMFORGE_BADER_EXECUTABLE"),"Bader executable not configured")
    def test_actual_bader_solver_conservation(self):
        values=[]
        for z in range(32):
            for y in range(32):
                for x in range(32):
                    point=np.array([x,y,z])/4
                    rho=0.0
                    for center in ([2,4,4],[6,4,4]):
                        delta=point-center; delta-=8*np.round(delta/8)
                        rho+=np.exp(-np.dot(delta,delta))
                    values.append(rho)
        grid=Grid((32,32,32),((8,0,0),(0,8,0),(0,0,8)),values,periodic=True,unit="e/A^3")
        grid=(grid*(2/grid.integrate())).with_sites([(1,2,4,4),(1,6,4,4)])
        path=self.root/"CHGCAR"; grid.save(path,format="vasp")
        result=bader(path,executable=os.environ["ATOMFORGE_BADER_EXECUTABLE"])
        populations=[atom["charge"] for atom in result["atoms"]]
        self.assertEqual(len(populations),2)
        np.testing.assert_allclose(populations,[1,1],atol=.002)
        self.assertAlmostEqual(sum(populations)+result["vacuum_charge"],2,places=4)

    def test_wavecar_orbital_and_lattice_validation(self):
        from pymatgen.io.vasp.outputs import Wavecar
        path=self.fixtures/"io/vasp/outputs/WAVECAR.N2"
        wave=Wavecar(str(path))
        structure=af.Structure(); structure.cell=wave.a.tolist()
        structure.add_atom("N",5,5,4.5); structure.add_atom("N",5,5,5.5)
        poscar=self.root/"POSCAR"; structure.save(poscar)
        output=self.root/"PARCHG"
        result=orbital_density(path,poscar,output,band=0,spin=0)
        self.assertTrue(np.isfinite(result.data["total"]).all())
        self.assertGreaterEqual(float(result.data["total"].min()),0)
        grid=load_volume(output).fields[0]
        # Parseval: mean reconstructed |psi|^2 equals coefficient norm squared.
        expected=float(np.sum(np.abs(wave.coeffs[0][0])**2))
        self.assertAlmostEqual(grid.integrate(),expected,places=7)
        structure.cell[0][0]*=2; structure.save(poscar)
        with self.assertRaises(ValueError): orbital_density(path,poscar,output)

    def test_procar_fatbands(self):
        procar=read_projections(self.fixtures/"io/vasp/outputs/PROCAR.simple")
        energies={str(spin):values for spin,values in procar.eigenvalues.items()}
        weights={str(spin):values.sum(axis=(2,3)) for spin,values in procar.data.items()}
        self.assertTrue(energies)
        for key in energies:
            self.assertEqual(energies[key].shape,weights[key].shape)
            self.assertTrue(np.isfinite(weights[key]).all())
        plot_bands({"energies_eV":energies},self.root/"fatbands.png",projections=weights)
        self.assertGreater((self.root/"fatbands.png").stat().st_size,1000)

    def test_ddec_reference_output(self):
        folder=self.fixtures/"command_line/chargemol/spin_unpolarized"
        result=ddec(folder,run=False)
        lines=(folder/"DDEC6_even_tempered_net_atomic_charges.xyz").read_text().splitlines()
        count=int(lines[0]); expected=[float(line.split()[4]) for line in lines[2:2+count]]
        np.testing.assert_allclose(result.ddec_charges,expected,atol=1e-12)
        self.assertAlmostEqual(sum(result.ddec_charges),0,places=4)
        self.assertEqual(len(result.cm5_charges),count)


if __name__=="__main__": unittest.main()
