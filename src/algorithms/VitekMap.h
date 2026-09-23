#pragma once

#include "model/Structure.h"

#include <glm/glm.hpp>

#include <string>
#include <vector>

// Vitek differential-displacement map: for each pair of near-neighbor atoms
// (in the undeformed reference lattice), the difference in their
// displacement (deformed minus reference) split into a component along the
// dislocation line/Burgers direction ("screw" component, classically drawn
// as arrow length in a differential-displacement plot of a screw
// dislocation core) and a component in the plane perpendicular to it
// ("edge" component). Follows the same algorithm as BABEL's vitek.f90
// (Vitek, Phil. Mag. 18 (1968) 773), generalized here to auto-detect the
// first-neighbor cutoff from the reference structure instead of BABEL's
// hardcoded bcc/fcc/hcp crystal-type table, so it isn't restricted to those
// specific cases.
namespace atomforge
{

struct VitekPair
{
    int indexI = -1;
    int indexJ = -1;
    glm::dvec3 positionI;   // reference (undeformed) position of atom i
    glm::dvec3 positionJ;   // reference (undeformed) position of atom j (possibly a periodic image)
    double screwComponent = 0.0;   // (u_j - u_i) . lineDir, wrapped to [-b/2, b/2)
    glm::dvec2 edgeComponent{0.0, 0.0}; // (u_j - u_i) projected onto the plane perpendicular to lineDir
};

struct VitekResult
{
    bool success = false;
    std::string message;
    std::vector<VitekPair> pairs;
};

// reference and deformed must have the same atom count and ordering.
// lineDirection is the dislocation line (and, for a pure screw
// dislocation, Burgers) direction; burgersMagnitude is used to wrap the
// screw component onto the elastic (not rigid-translation) part of the
// displacement, matching BABEL's convention. neighborCutoff should be
// roughly the midpoint between the first- and second-neighbor shell
// distances (auto-detected from `reference` if left at 0). usePbcXY also
// considers periodic images of atom j in the plane perpendicular to the
// line (needed when the two atoms are only neighbors through a periodic
// boundary).
VitekResult computeVitekMap(const Structure& reference,
                            const Structure& deformed,
                            const glm::dvec3& lineDirection,
                            double burgersMagnitude,
                            float neighborCutoff,
                            bool usePbcXY);

}
