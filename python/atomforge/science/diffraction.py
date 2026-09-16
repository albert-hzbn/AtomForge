"""Powder X-ray/neutron patterns and explicit single-crystal structure factors."""

import math


def _structure(structure):
    from pymatgen.core import Lattice, Structure
    if not structure.cell:
        raise ValueError("Diffraction requires a periodic cell")
    return Structure(Lattice(structure.cell), [a.symbol for a in structure.atoms],
                     [(a.x,a.y,a.z) for a in structure.atoms], coords_are_cartesian=True)


def powder_diffraction(structure, *, radiation="xray", wavelength=1.5406, two_theta=(0,90)):
    """Return pymatgen diffraction peaks; wavelength A, scattering angle degrees."""
    if not math.isfinite(wavelength) or wavelength <= 0:
        raise ValueError("Wavelength must be positive and finite")
    if len(two_theta)!=2 or not 0 <= two_theta[0] < two_theta[1] <= 180:
        raise ValueError("Require 0 <= minimum < maximum <= 180 degrees")
    if radiation == "xray":
        from pymatgen.analysis.diffraction.xrd import XRDCalculator
        calculator = XRDCalculator(wavelength=wavelength)
    elif radiation == "neutron":
        from pymatgen.analysis.diffraction.neutron import NDCalculator
        calculator = NDCalculator(wavelength=wavelength)
    else:
        raise ValueError("Radiation must be xray or neutron")
    return calculator.get_pattern(_structure(structure), two_theta_range=tuple(two_theta))


def single_crystal(structure, indices, scattering):
    """Return hkl, reciprocal vectors (2*pi/A) and |F|^2 for explicit site factors.

    scattering is [reflection][atom], permitting complex anomalous factors and
    neutron isotope factors. Intensities omit instrument/polarization corrections.
    """
    import numpy as np
    crystal = _structure(structure)
    indices = np.asarray(indices, dtype=float)
    scattering = np.asarray(scattering, dtype=complex)
    if indices.ndim!=2 or indices.shape[1]!=3 or not np.isfinite(indices).all() or not np.equal(indices,np.round(indices)).all():
        raise ValueError("Miller indices must have shape (n,3) and be integers")
    if scattering.shape != (len(indices),len(structure)) or not np.isfinite(scattering).all():
        raise ValueError("One finite scattering factor per reflection and site is required")
    factors = (scattering * np.exp(2j*np.pi*(indices @ crystal.frac_coords.T))).sum(axis=1)
    return {"hkl":indices.astype(int), "q_1_A":indices @ crystal.lattice.reciprocal_lattice.matrix,
            "factors":factors, "intensities":np.abs(factors)**2}


def plot_diffraction(pattern, output, *, fwhm=0.1, points=4000):
    """Plot Gaussian-broadened powder peaks; FWHM is in degrees of 2 theta."""
    import numpy as np
    from .spectra import _figure
    if not math.isfinite(fwhm) or fwhm <= 0 or not isinstance(points,int) or points<2 or points>1000000:
        raise ValueError("Invalid peak width or point count")
    if len(pattern.x)==0:
        raise ValueError("No diffraction peaks in the selected range")
    x = np.linspace(max(0,min(pattern.x)-5*fwhm),min(180,max(pattern.x)+5*fwhm),points)
    y = np.zeros(points)
    for position, intensity in zip(pattern.x,pattern.y):
        y += intensity*np.exp(-4*np.log(2)*((x-position)/fwhm)**2)
    figure, axes = _figure(); axes.plot(x,y)
    axes.set_xlabel("2 theta (degrees)"); axes.set_ylabel("Relative intensity")
    figure.savefig(str(output),dpi=180)
    return figure
