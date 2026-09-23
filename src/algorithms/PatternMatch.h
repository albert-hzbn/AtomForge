#pragma once

#include "model/Structure.h"

#include <glm/glm.hpp>

#include <string>
#include <vector>

// Crystallographic pattern matching: build a reference "pattern" (the set
// of neighbor directions around an atom in a perfect crystal), then check
// which atoms in a possibly-defective structure still have a locally
// matching neighbor environment. Atoms that don't match are flagged as
// being near a defect (dislocation core, stacking fault, free surface, ...).
// Follows the same two-stage design as BABEL's patternInit/patternDetect
// programs, generalized here to auto-extract the pattern directly from any
// given perfect reference structure instead of BABEL's hardcoded bcc/fcc/hcp
// direction tables, so it isn't restricted to those specific crystals.
// Unlike BABEL's own algorithm, this does not search over trial rotations
// to accommodate small elastic lattice rotation near a defect (only a
// direct angular comparison against the pattern's fixed orientation); a
// real elastic rotation larger than angleThresholdDeg will itself register
// as "not matched", which is a coarser (but still meaningful) result than
// BABEL's rotation-tolerant matching.
namespace atomforge
{

struct CrystalPattern
{
    bool success = false;
    std::string message;

    // Unit vectors from the representative atom to each of its neighbors in
    // the reference structure.
    std::vector<glm::dvec3> directions;
};

// cutoff should enclose exactly the first coordination shell (a small
// tolerance is fine); the representative atom is the one nearest the
// reference structure's geometric center, to avoid a boundary atom giving
// an incomplete neighbor list.
CrystalPattern buildPattern(const Structure& reference, float cutoff, bool usePbc);

struct PatternMatchRow
{
    int index = -1;
    bool matched = false;
    int neighborCount = 0;          // actual neighbor count within cutoff
    double maxAngleDeviationDeg = 0.0; // for matched atoms: worst-case angle to its assigned pattern direction
};

struct PatternDetectResult
{
    bool success = false;
    std::string message;
    std::vector<PatternMatchRow> rows; // one per atom, `structure`-indexed
};

PatternDetectResult detectPattern(const Structure& structure,
                                  const CrystalPattern& pattern,
                                  float cutoff,
                                  double angleThresholdDeg,
                                  bool usePbc);

}
