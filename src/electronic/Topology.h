#pragma once

#include <vector>
#include <glm/glm.hpp>

// Betti curves: topological descriptors of a scalar field's superlevel sets
// (Betti curves for the electron density: Bartel et al., ACS Mater. Lett.
// 2025, 7, 2158; the underlying method is standard cubical-complex
// persistent-homology bookkeeping). At each threshold t, the "solid" set is
// {voxels with value >= t}; betti0(t)/betti1(t)/betti2(t) count its
// connected components, independent loops and enclosed cavities.
//
// Finite (non-periodic) grids only: the enclosed-cavity count (betti2) is
// computed as the number of complementary background components that do not
// touch the grid boundary, which assumes an R^3 (bounded-box) ambient space,
// not the 3-torus topology of a periodic cell. betti1 is then recovered from
// the Euler-Poincare identity chi = betti0 - betti1 + betti2, with chi
// counted exactly from the cubical complex formed by the union of solid
// voxels' closed unit cubes (no higher (betti3) term: a finite, non-periodic
// subcomplex of R^3 cannot enclose a genuine 3-cycle).
namespace atomforge::electronic
{
struct Grid;

struct BettiNumbers { int betti0 = 0; int betti1 = 0; int betti2 = 0; };

BettiNumbers bettiNumbers(const Grid& field, double threshold);
// One (threshold, betti0, betti1, betti2) row per requested threshold.
std::vector<glm::dvec4> bettiCurve(const Grid& field, const std::vector<double>& thresholds);
}
