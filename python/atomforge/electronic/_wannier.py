"""Wannier-interpolated bands, Berry curvature and Chern numbers.

Thin ctypes wrapper over the native engine's WannierHamiltonian (see
src/electronic/Wannier.h). Reads a Wannier90 seedname_hr.dat real-space
Hamiltonian (Pizzi et al., J. Phys.: Condens. Matter 2020; originally
Mostofi et al., Comput. Phys. Commun. 2014) and Fourier-interpolates it on
the Wigner-Seitz point set, matching Wannier90/WannierBerri's own convention:

    H(k) = sum_R exp(2*pi*i*k.R) * H(R) / ndegen(R)

with k in fractional reciprocal coordinates and ndegen(R) the Wigner-Seitz
degeneracy Wannier90 writes alongside each R (dividing by it is required;
omitting it double-counts boundary lattice points shared between equivalent
Wigner-Seitz images).

berry_curvature() is the Hamiltonian-gauge (Kubo) term only (Wang, Yates,
Souza and Vanderbilt, Phys. Rev. B 74, 195118, 2006): it omits the
position-operator corrections that require seedname_r.dat, so it is exact
only when the Wannier gauge and the Hamiltonian (band) gauge coincide, and
it diverges at exact band degeneracies. Curvature is with respect to the
fractional reciprocal coordinates above (not an actual-momentum/Cartesian
inverse-length convention): integrating it over the fractional Brillouin
zone [0,1)^2 and dividing by 2*pi reproduces the integer chern_number().
chern_number() uses the gauge-invariant Fukui-Hatsugai-Suzuki lattice method
(J. Phys. Soc. Jpn. 74, 1674, 2005) instead, which does not rely on that
approximation and stays well-defined arbitrarily close to (but not exactly
at) a degeneracy.
"""

import ctypes as ct

from ._native import doubles, library


class WannierHamiltonian:
    """A Wannier90 seedname_hr.dat Hamiltonian, interpolated by the native engine.

    Construct with read_hr(), not directly.
    """

    def __init__(self, pointer):
        self._lib = library()
        if not pointer:
            raise ValueError(self._lib.af_wannier_error().decode("utf-8", errors="replace"))
        self._pointer = pointer

    def __del__(self):
        pointer = getattr(self, "_pointer", None)
        if pointer:
            self._lib.af_wannier_free(pointer)
            self._pointer = None

    def _check(self, ok):
        if not ok:
            raise ValueError(self._lib.af_wannier_error().decode("utf-8", errors="replace"))

    @property
    def num_wann(self):
        return self._lib.af_wannier_num_wann(self._pointer)

    def hamiltonian(self, kpoint):
        """Return the interpolated Bloch Hamiltonian H(k); kpoint is fractional reciprocal."""
        n = self.num_wann
        buffer = (ct.c_double * (2 * n * n))()
        self._check(self._lib.af_wannier_hamiltonian(self._pointer, doubles(kpoint), buffer))
        return [complex(buffer[2 * i], buffer[2 * i + 1]) for i in range(n * n)]

    def bands(self, kpoints):
        """Return ascending-sorted eigenvalues per k point, as one list per k point."""
        kpoints = [tuple(k) for k in kpoints]
        n = self.num_wann
        buffer = (ct.c_double * (len(kpoints) * n))()
        flat = doubles([v for k in kpoints for v in k])
        self._check(self._lib.af_wannier_bands(self._pointer, flat, len(kpoints), buffer))
        return [list(buffer[i * n:(i + 1) * n]) for i in range(len(kpoints))]

    def berry_curvature(self, kpoint, plane=(0, 1)):
        """Hamiltonian-gauge (Kubo) Berry curvature of every band at kpoint, ascending order.

        plane selects two distinct fractional-reciprocal directions (0, 1, 2).
        See the module docstring for what this omits, and for exact
        degeneracies use chern_number() instead.
        """
        a, b = plane
        buffer = (ct.c_double * self.num_wann)()
        self._check(self._lib.af_wannier_berry_curvature(self._pointer, doubles(kpoint), a, b, buffer))
        return list(buffer)

    def chern_number(self, band, plane=(0, 1), grid=30, fixed=0.0):
        """Fukui-Hatsugai-Suzuki lattice Chern number of one non-degenerate band.

        Integrates the gauge-invariant plaquette flux over the full
        Brillouin zone spanned by the two fractional-reciprocal directions
        in plane, with the third fractional coordinate held at fixed. The
        result is an integer to numerical precision whenever band stays
        non-degenerate and gapped from its neighbors everywhere on the mesh;
        raise the mesh density or nudge fixed if a plaquette link collapses
        (a degeneracy on the mesh).
        """
        a, b = plane
        buffer = (ct.c_double * 1)()
        self._check(self._lib.af_wannier_chern_number(self._pointer, band, a, b, grid, fixed, buffer))
        return buffer[0]


def read_hr(path):
    """Parse a Wannier90 seedname_hr.dat file into a WannierHamiltonian."""
    return WannierHamiltonian(library().af_wannier_load(str(path).encode("utf-8")))
