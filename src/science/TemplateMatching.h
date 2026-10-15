#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// Polyhedral template matching (Larsen, Schmidt and Schiotz, 2016): each
// atom's nearest neighbours are fitted, by the best rotation and scale, to the
// ideal fcc, hcp, bcc, icosahedral, simple-cubic and diamond shells. The
// structure with the smallest RMSD below a cutoff wins; the fit gives a local
// lattice orientation (quaternion) and lattice constant.
Json polyhedralTemplateMatching(const Parameters& p);

// Grains of a polycrystal: crystalline neighbours with a small misorientation
// are joined (in increasing order of misorientation, while the mean
// orientations of the two grains also agree); non-crystalline atoms join the
// nearest grain.
Json grainSegmentation(const Parameters& p);

// Atomic elastic strain: the deformation gradient that maps the ideal template
// (scaled by the given lattice constant) onto the neighbours.
Json elasticStrain(const Parameters& p);

// Planar faults in fcc: hcp-coordinated atoms grouped by how many of their
// out-of-plane neighbours are hcp too (twin boundary, intrinsic stacking
// fault, multilayer fault).
Json fccPlanarFaults(const Parameters& p);
}
