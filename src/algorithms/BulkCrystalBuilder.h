#pragma once

#include "model/Structure.h"

#include <glm/glm.hpp>
#include <string>
#include <vector>

// -- Crystal system and parameter types --------------------------------------

enum class CrystalSystem
{
    Triclinic = 0,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Trigonal,
    Hexagonal,
    Cubic,
};

struct BulkBuildResult
{
    bool success = false;
    std::string message;
    int generatedAtoms = 0;
    int spaceGroupNumber = 0;
    std::string spaceGroupSymbol;
};

struct SpaceGroupRange
{
    CrystalSystem system;
    const char* label;
    int first;
    int last;
};

struct LatticeParameters
{
    double a = 4.0;
    double b = 4.0;
    double c = 4.0;
    double alpha = 90.0;
    double beta = 90.0;
    double gamma = 90.0;
};

extern const SpaceGroupRange kSpaceGroupRanges[];
constexpr int kNumSpaceGroupRanges = 7;

// -- Helper functions --------------------------------------------------------

bool hasElementColor(int atomicNumber, const std::vector<glm::vec3>& elementColors);

void applyElementToAtom(AtomSite& atom,
                        int atomicNumber,
                        const std::vector<glm::vec3>& elementColors);

const char* crystalSystemLabel(CrystalSystem system);

const SpaceGroupRange& currentRange(CrystalSystem system);

glm::mat3 buildLatticeFromParameters(const LatticeParameters& params);

bool validateParameters(CrystalSystem system, const LatticeParameters& params,
                        std::string& error);

void applySystemConstraints(CrystalSystem system, LatticeParameters& params);

std::vector<int>& hallBySpaceGroup();

std::string spaceGroupLabel(int spaceGroupNumber);

glm::vec3 wrapFractional(const glm::vec3& frac);

bool sameFractionalPosition(const glm::vec3& lhs, const glm::vec3& rhs, float tol);

// -- Builder -----------------------------------------------------------------

[[nodiscard]] BulkBuildResult buildBulkCrystal(Structure& structure,
                                 CrystalSystem system,
                                 int spaceGroupNumber,
                                 const LatticeParameters& params,
                                 const std::vector<AtomSite>& asymmetricAtoms,
                                 const std::vector<glm::vec3>& elementColors);

void addDefaultAsymmetricAtom(std::vector<AtomSite>& atoms,
                              const std::vector<glm::vec3>& elementColors);
