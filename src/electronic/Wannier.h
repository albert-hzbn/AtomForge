#pragma once

#include <complex>
#include <string>
#include <vector>
#include <glm/glm.hpp>

// Wannier-interpolated bands, Berry curvature and Chern numbers from a
// Wannier90 seedname_hr.dat real-space Hamiltonian (Pizzi et al., J. Phys.:
// Condens. Matter 2020; originally Mostofi et al., Comput. Phys. Commun.
// 2014). H(k) = sum_R exp(2*pi*i*k.R) H(R) / ndegen(R), k fractional
// reciprocal, matching Wannier90/WannierBerri's own convention; dividing by
// the Wigner-Seitz degeneracy ndegen(R) is required, not optional (it avoids
// double-counting boundary lattice points shared between equivalent
// Wigner-Seitz images).
//
// berryCurvature() is the Hamiltonian-gauge (Kubo) term only (Wang, Yates,
// Souza and Vanderbilt, Phys. Rev. B 74, 195118, 2006): it omits the
// position-operator corrections that require seedname_r.dat, so it is exact
// only when the Wannier gauge and the Hamiltonian (band) gauge coincide, and
// it diverges at exact band degeneracies. Curvature is with respect to the
// fractional reciprocal coordinates above (not an actual-momentum/Cartesian
// inverse-length convention): integrating it over the fractional Brillouin
// zone [0,1)^2 and dividing by 2*pi reproduces the integer chernNumber().
// chernNumber() uses the gauge-invariant Fukui-Hatsugai-Suzuki lattice
// method (J. Phys. Soc. Jpn. 74, 1674, 2005) instead, which does not rely on
// that approximation and stays well-defined arbitrarily close to (but not
// exactly at) a degeneracy.
namespace atomforge::electronic
{
struct WannierHamiltonian
{
    int numWann = 0;
    std::vector<glm::ivec3> vectors;                              // R, file order
    std::vector<std::vector<std::complex<double>>> hamiltonians;  // H(R), row-major numWann x numWann, aligned with vectors
    std::vector<int> degeneracies;                                // Wigner-Seitz ndegen(R), aligned with vectors
};

WannierHamiltonian readWannierHr(const std::string& path);
// Row-major numWann x numWann H(k); kpoint is fractional reciprocal.
std::vector<std::complex<double>> blochHamiltonian(const WannierHamiltonian& hr, glm::dvec3 kpoint);
// Ascending eigenvalues at each kpoint, flattened (kpoints.size() x numWann).
std::vector<double> wannierBands(const WannierHamiltonian& hr, const std::vector<glm::dvec3>& kpoints);
// Hamiltonian-gauge (Kubo) Berry curvature of every band at kpoint, ascending order.
// plane selects two distinct fractional-reciprocal directions (0, 1 or 2).
std::vector<double> berryCurvature(const WannierHamiltonian& hr, glm::dvec3 kpoint, int planeA, int planeB);
// Fukui-Hatsugai-Suzuki lattice Chern number of one non-degenerate band over
// the grid x grid mesh spanning [0,1)^2 in (planeA, planeB), third fractional
// coordinate held at `fixed`.
double chernNumber(const WannierHamiltonian& hr, int band, int planeA, int planeB, int grid, double fixed);
}
