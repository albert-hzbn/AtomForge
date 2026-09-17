#pragma once

#include <vector>
#include <glm/glm.hpp>

// Native QTAIM/Bader charge partitioning: replaces AtomForge's dependence on
// the external Henkelman executable (see python/atomforge/science/
// population.py's bader()) for the on-grid variant of the algorithm
// (Henkelman, Sanville and Jonsson, Comput. Mater. Sci. 36, 354, 2006), the
// zero-flux atomic-basin definition underlying the Quantum Theory of Atoms
// in Molecules (Bader, "Atoms in Molecules: A Quantum Theory", 1990). This
// implements the on-grid step only, not the near-grid boundary refinement,
// so individual dividing-surface placement is accurate to about one grid
// spacing; integrated basin charges converge with grid resolution like the
// reference algorithm's own on-grid mode.
namespace atomforge::electronic
{
struct Grid;

struct BaderBasin { double charge = 0.0; double volume = 0.0; glm::dvec3 maximum{0.0}; };
struct BaderPartition { std::vector<int> basin; std::vector<BaderBasin> basins; };  // basin: one entry per voxel, x-fastest order

// Assigns every voxel to the atomic basin reached by repeatedly stepping to
// whichever of its 26 neighbors maximizes the directional density increase
// per unit Cartesian distance, until a local density maximum is reached.
// Requires nonnegative electron density in e/A^3.
BaderPartition baderOnGrid(const Grid& density);
// Sums each basin's charge into whichever of `sites` its maximum is nearest
// to (periodic minimum image when density.periodic); one population per site.
std::vector<double> baderPopulations(const Grid& density, const std::vector<glm::dvec3>& sites);
}
