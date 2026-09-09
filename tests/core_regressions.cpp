#include "graphics/Picking.h"
#include "math/StructureMath.h"
#include "algorithms/SubstitutionalSolidSolutionBuilder.h"
#include "algorithms/AmorphousBuilder.h"
#include <cstdlib>
#include <iostream>
#include <limits>

static void check(bool condition, const char* message)
{
    if (!condition) { std::cerr << message << '\n'; std::exit(1); }
}

int main()
{
    const glm::vec3 origin(0.0f), dir(0.0f, 0.0f, 1.0f);
    const std::vector<glm::vec3> positions{{0, 0, 5}, {0, 0, 6}};
    check(pickAtom(origin, dir, positions, {0.5f, 2.0f}, 1) == 1,
          "Picking must select the nearest surface, not center");
    check(pickAtom(origin, 3.0f * dir, positions, {0.5f, 2.0f}, 1) == 1,
          "Picking must accept a non-unit direction");
    check(pickAtom(origin, origin, positions, {}, 1) == -1, "Zero ray must miss");
    check(pickAtom(origin, dir, {{0, 0, -0.5f}}, {1}, 1) == 0,
          "Ray starting inside a sphere must hit");
    check(pickAtom(origin, dir, {{0, 0, -5}}, {1}, 1) == -1,
          "Sphere behind ray must miss");
    Structure base;
    base.atoms.resize(3);
    base.grainColors = {{{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}};
    base.grainRegionIds = {10, 20, 30};
    base.eraseAtom(1);
    check(base.atoms.size() == 2 && base.grainRegionIds == std::vector<int>({10, 30})
          && base.grainColors[1][2] == 1, "Deleting atoms must preserve metadata alignment");
    base.hasUnitCell = true;
    base.cellVectors = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    glm::mat3 cell, inverse;
    check(tryMakeCellMatrices(base, cell, inverse), "Valid cell rejected");
    base.cellVectors[0][0] = std::numeric_limits<double>::quiet_NaN();
    check(!tryMakeCellMatrices(base, cell, inverse), "Non-finite cell accepted");
    base.atoms.resize(101);
    SSSParams params;
    params.composition = {{26, 0.3f}, {28, 0.7f}};
    check(buildSubstitutionalSolidSolution(base, params).success, "Valid SSS failed");
    params.composition[0].fraction = std::numeric_limits<float>::quiet_NaN();
    check(!buildSubstitutionalSolidSolution(base, params).success, "NaN SSS accepted");
    params.composition[0].fraction = std::numeric_limits<float>::infinity();
    check(!buildSubstitutionalSolidSolution(base, params).success, "Infinite SSS accepted");
    AmorphousParams amorphous;
    amorphous.elements = {{26, std::numeric_limits<int>::max()}, {28, 1}};
    check(!buildAmorphousStructure(amorphous, {}, {}).success, "Overflowing counts accepted");
    std::cout << "Core regression tests passed\n";
}
