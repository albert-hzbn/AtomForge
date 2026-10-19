#pragma once

#include "science/ScienceData.h"

namespace atomforge::science
{
// Chill+ (Nguyen and Molinero, 2015): water oxygens with four neighbours are
// hexagonal or cubic ice, interfacial ice, clathrate hydrate or interfacial
// hydrate according to the q3 correlations (staggered, eclipsed) of their bonds.
Json chillPlus(const Parameters& p);

// Local pair entropy fingerprint (Piaggi and Parrinello, 2017) from each
// atom's Gaussian-smoothed g(r), optionally averaged over its neighbours.
Json localEntropy(const Parameters& p);

// Shortest-path rings (Franzblau, Phys. Rev. B 44, 4925 (1991)) of the bond
// network (bonds shorter than the cutoff, periodic images included): closed
// paths in which every two atoms are joined around the ring by a shortest
// path of the network. Reports the ring-size distribution and rings per atom.
Json findRings(const Parameters& p);

// Displacement of each atom from a reference configuration (minimum image in
// periodic cells; optionally after mapping a deformed cell onto the reference).
Json displacementVectors(const Parameters& p);

// Distributions of bond lengths (per pair of species) and bond angles (per
// triplet, centre in the middle) for bonds shorter than the cutoff.
Json bondDistributions(const Parameters& p);
}
