#include "algorithms/BulkCrystalBuilder.h"

#include "ElementData.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#ifdef ATOMS_ENABLE_SPGLIB
#include <spglib.h>
#endif

// -- Constants ---------------------------------------------------------------

const SpaceGroupRange kSpaceGroupRanges[] = {
    { CrystalSystem::Triclinic,    "Triclinic",    1,   2 },
    { CrystalSystem::Monoclinic,   "Monoclinic",   3,  15 },
    { CrystalSystem::Orthorhombic, "Orthorhombic", 16, 74 },
    { CrystalSystem::Tetragonal,   "Tetragonal",   75, 142 },
    { CrystalSystem::Trigonal,     "Trigonal",     143, 167 },
    { CrystalSystem::Hexagonal,    "Hexagonal",    168, 194 },
    { CrystalSystem::Cubic,        "Cubic",        195, 230 },
};

// -- Helper functions --------------------------------------------------------

bool hasElementColor(int atomicNumber, const std::vector<glm::vec3>& elementColors)
{
    return atomicNumber >= 0 && atomicNumber < (int)elementColors.size();
}

void applyElementToAtom(AtomSite& atom,
                        int atomicNumber,
                        const std::vector<glm::vec3>& elementColors)
{
    atom.atomicNumber = atomicNumber;
    atom.symbol = elementSymbol(atomicNumber);

    if (hasElementColor(atomicNumber, elementColors))
    {
        atom.r = elementColors[atomicNumber].r;
        atom.g = elementColors[atomicNumber].g;
        atom.b = elementColors[atomicNumber].b;
    }
    else
    {
        getDefaultElementColor(atomicNumber, atom.r, atom.g, atom.b);
    }
}

const char* crystalSystemLabel(CrystalSystem system)
{
    for (int i = 0; i < kNumSpaceGroupRanges; ++i)
    {
        if (kSpaceGroupRanges[i].system == system)
            return kSpaceGroupRanges[i].label;
    }
    return "Unknown";
}

const SpaceGroupRange& currentRange(CrystalSystem system)
{
    for (int i = 0; i < kNumSpaceGroupRanges; ++i)
    {
        if (kSpaceGroupRanges[i].system == system)
            return kSpaceGroupRanges[i];
    }
    return kSpaceGroupRanges[0];
}

glm::dmat3 buildLatticeFromParameters(const LatticeParameters& params)
{
    const double degToRad = 3.14159265358979323846 / 180.0;
    const double alpha = params.alpha * degToRad;
    const double beta = params.beta * degToRad;
    const double gamma = params.gamma * degToRad;

    const double cosAlpha = std::cos(alpha);
    const double cosBeta = std::cos(beta);
    const double cosGamma = std::cos(gamma);
    const double sinGamma = std::sin(gamma);

    glm::dvec3 a(params.a, 0.0, 0.0);
    glm::dvec3 b((params.b * cosGamma),
                (params.b * sinGamma),
                0.0);

    double cx = params.c * cosBeta;
    double cy = 0.0;
    if (std::abs(sinGamma) > 1e-10)
        cy = params.c * (cosAlpha - cosBeta * cosGamma) / sinGamma;
    double cz2 = params.c * params.c - cx * cx - cy * cy;
    if (cz2 < 0.0 && cz2 > -1e-8)
        cz2 = 0.0;
    glm::dvec3 c(cx, cy, ((cz2 > 0.0) ? std::sqrt(cz2) : 0.0));

    return glm::dmat3(a, b, c);
}

bool validateParameters(CrystalSystem system, const LatticeParameters& params, std::string& error)
{
    auto near = [](double lhs, double rhs) {
        return std::abs(lhs - rhs) < 1e-6;
    };

    if (params.a <= 0.0 || params.b <= 0.0 || params.c <= 0.0)
    {
        error = "Lattice lengths must be positive.";
        return false;
    }

    if (params.alpha <= 0.0 || params.alpha >= 180.0 ||
        params.beta <= 0.0 || params.beta >= 180.0 ||
        params.gamma <= 0.0 || params.gamma >= 180.0)
    {
        error = "Angles must be in the range (0, 180).";
        return false;
    }

    switch (system)
    {
        case CrystalSystem::Monoclinic:
            if (!near(params.alpha, 90.0) || !near(params.gamma, 90.0))
            {
                error = "Monoclinic builder uses alpha=gamma=90 deg.";
                return false;
            }
            break;
        case CrystalSystem::Orthorhombic:
            if (!near(params.alpha, 90.0) || !near(params.beta, 90.0) || !near(params.gamma, 90.0))
            {
                error = "Orthorhombic builder uses alpha=beta=gamma=90 deg.";
                return false;
            }
            break;
        case CrystalSystem::Tetragonal:
            if (!near(params.a, params.b) || !near(params.alpha, 90.0) || !near(params.beta, 90.0) || !near(params.gamma, 90.0))
            {
                error = "Tetragonal builder uses a=b and alpha=beta=gamma=90 deg.";
                return false;
            }
            break;
        case CrystalSystem::Trigonal:
            if (!near(params.a, params.b) || !near(params.alpha, 90.0) || !near(params.beta, 90.0) || !near(params.gamma, 120.0))
            {
                error = "Trigonal builder uses the hexagonal setting: a=b, alpha=beta=90 deg, gamma=120 deg.";
                return false;
            }
            break;
        case CrystalSystem::Hexagonal:
            if (!near(params.a, params.b) || !near(params.alpha, 90.0) || !near(params.beta, 90.0) || !near(params.gamma, 120.0))
            {
                error = "Hexagonal builder uses a=b, alpha=beta=90 deg, gamma=120 deg.";
                return false;
            }
            break;
        case CrystalSystem::Cubic:
            if (!near(params.a, params.b) || !near(params.b, params.c) ||
                !near(params.alpha, 90.0) || !near(params.beta, 90.0) || !near(params.gamma, 90.0))
            {
                error = "Cubic builder uses a=b=c and alpha=beta=gamma=90 deg.";
                return false;
            }
            break;
        case CrystalSystem::Triclinic:
        default:
            break;
    }

    glm::dmat3 lattice = buildLatticeFromParameters(params);
    if (std::abs(glm::determinant(lattice)) < 1e-10)
    {
        error = "Lattice parameters produce a zero-volume cell.";
        return false;
    }

    return true;
}

void applySystemConstraints(CrystalSystem system, LatticeParameters& params)
{
    switch (system)
    {
        case CrystalSystem::Monoclinic:
            params.alpha = 90.0;
            params.gamma = 90.0;
            break;
        case CrystalSystem::Orthorhombic:
            params.alpha = 90.0;
            params.beta = 90.0;
            params.gamma = 90.0;
            break;
        case CrystalSystem::Tetragonal:
            params.b = params.a;
            params.alpha = 90.0;
            params.beta = 90.0;
            params.gamma = 90.0;
            break;
        case CrystalSystem::Trigonal:
            params.b = params.a;
            params.alpha = 90.0;
            params.beta = 90.0;
            params.gamma = 120.0;
            break;
        case CrystalSystem::Hexagonal:
            params.b = params.a;
            params.alpha = 90.0;
            params.beta = 90.0;
            params.gamma = 120.0;
            break;
        case CrystalSystem::Cubic:
            params.b = params.a;
            params.c = params.a;
            params.alpha = 90.0;
            params.beta = 90.0;
            params.gamma = 90.0;
            break;
        case CrystalSystem::Triclinic:
        default:
            break;
    }
}

std::vector<int>& hallBySpaceGroup()
{
    static std::vector<int> hallByNumber(231, 0);
    static bool initialized = false;
    if (!initialized)
    {
#ifdef ATOMS_ENABLE_SPGLIB
        for (int hall = 1; hall <= 530; ++hall)
        {
            SpglibSpacegroupType type = spg_get_spacegroup_type(hall);
            if (type.number >= 1 && type.number <= 230 && hallByNumber[type.number] == 0)
                hallByNumber[type.number] = hall;
        }
#endif
        initialized = true;
    }
    return hallByNumber;
}

std::string spaceGroupLabel(int spaceGroupNumber)
{
#ifdef ATOMS_ENABLE_SPGLIB
    const std::vector<int>& hallMap = hallBySpaceGroup();
    if (spaceGroupNumber >= 1 && spaceGroupNumber <= 230 && hallMap[spaceGroupNumber] > 0)
    {
        SpglibSpacegroupType type = spg_get_spacegroup_type(hallMap[spaceGroupNumber]);
        std::ostringstream out;
        out << type.number << " " << type.international_short;
        return out.str();
    }
#endif
    std::ostringstream out;
    out << spaceGroupNumber;
    return out.str();
}

// The Hall number of a space group's origin choice 2, or 0 when it has one origin.
static int secondOriginHall(int spaceGroupNumber)
{
#ifdef ATOMS_ENABLE_SPGLIB
    const std::vector<int>& hallMap = hallBySpaceGroup();
    if (spaceGroupNumber < 1 || spaceGroupNumber > 230 || hallMap[spaceGroupNumber] == 0) return 0;
    for (int hall = hallMap[spaceGroupNumber]; hall <= 530; ++hall)
    {
        const SpglibSpacegroupType type = spg_get_spacegroup_type(hall);
        if (type.number != spaceGroupNumber) break;
        if (std::string(type.choice) == "2") return hall;
    }
#else
    (void)spaceGroupNumber;
#endif
    return 0;
}

bool hasTwoOrigins(int spaceGroupNumber)
{
    return secondOriginHall(spaceGroupNumber) != 0;
}

glm::dvec3 wrapFractional(const glm::dvec3& frac)
{
    glm::dvec3 wrapped = frac;
    wrapped.x -= std::floor(wrapped.x);
    wrapped.y -= std::floor(wrapped.y);
    wrapped.z -= std::floor(wrapped.z);
    // Snap values that differ from a simple fraction only by rounding
    // (1/3 as 0.3333333) and ones a hair below 1.
    for (int k = 0; k < 3; ++k)
    {
        if (wrapped[k] >= 1.0 - 1e-9) wrapped[k] = 0.0;
        for (int d : {2, 3, 4, 6, 8, 12})
        {
            const double n = std::round(wrapped[k] * d);
            if (std::abs(wrapped[k] * d - n) < 1e-6 * d) { wrapped[k] = n / d; break; }
        }
        if (wrapped[k] >= 1.0) wrapped[k] = 0.0;
    }
    return wrapped;
}

bool sameFractionalPosition(const glm::dvec3& lhs, const glm::dvec3& rhs, double tol)
{
    glm::dvec3 delta = lhs - rhs;
    delta.x -= std::round(delta.x);
    delta.y -= std::round(delta.y);
    delta.z -= std::round(delta.z);
    return std::abs(delta.x) <= tol && std::abs(delta.y) <= tol && std::abs(delta.z) <= tol;
}

// -- Builder -----------------------------------------------------------------

BulkBuildResult buildBulkCrystal(Structure& structure,
                                 CrystalSystem system,
                                 int spaceGroupNumber,
                                 const LatticeParameters& params,
                                 const std::vector<AtomSite>& asymmetricAtoms,
                                 const std::vector<glm::vec3>& elementColors,
                                 int originChoice)
{
    BulkBuildResult result;
    result.spaceGroupNumber = spaceGroupNumber;
    result.spaceGroupSymbol = spaceGroupLabel(spaceGroupNumber);
    if (originChoice != 1 && originChoice != 2)
    {
        result.message = "The origin choice must be 1 or 2.";
        return result;
    }
    if (hasTwoOrigins(spaceGroupNumber))
        result.spaceGroupSymbol += ", origin choice " + std::to_string(originChoice);
    else if (originChoice == 2)
    {
        result.message = "Space group " + std::to_string(spaceGroupNumber) + " has only one origin choice.";
        return result;
    }

#ifndef ATOMS_ENABLE_SPGLIB
    (void)system;
    (void)params;
    (void)asymmetricAtoms;
    (void)elementColors;
    result.message = "Bulk crystal builder requires spglib at build time.";
    return result;
#else
    std::string validationError;
    if (!validateParameters(system, params, validationError))
    {
        result.message = validationError;
        return result;
    }
    if (asymmetricAtoms.empty())
    {
        result.message = "Add at least one asymmetric-unit atom before building.";
        return result;
    }

    const std::vector<int>& hallMap = hallBySpaceGroup();
    if (spaceGroupNumber < 1 || spaceGroupNumber > 230 || hallMap[spaceGroupNumber] == 0)
    {
        result.message = "Selected space group is not available in the symmetry database.";
        return result;
    }

    const int hallNumber = originChoice == 2 ? secondOriginHall(spaceGroupNumber) : hallMap[spaceGroupNumber];
    int rotations[192][3][3];
    double translations[192][3];
    int operationCount = spg_get_symmetry_from_database(rotations, translations, hallNumber);
    if (operationCount <= 0)
    {
        result.message = "Could not load symmetry operations for the selected space group.";
        return result;
    }

    const glm::dmat3 lattice = buildLatticeFromParameters(params);
    const glm::dvec3 a = lattice[0];
    const glm::dvec3 b = lattice[1];
    const glm::dvec3 c = lattice[2];

    std::vector<AtomSite> generatedAtoms;
    generatedAtoms.reserve(asymmetricAtoms.size() * (size_t)operationCount);
    std::vector<glm::dvec3> generatedFrac;
    generatedFrac.reserve(asymmetricAtoms.size() * (size_t)operationCount);

    for (size_t atomIndex = 0; atomIndex < asymmetricAtoms.size(); ++atomIndex)
    {
        const AtomSite& basisAtom = asymmetricAtoms[atomIndex];
        glm::dvec3 frac(basisAtom.x, basisAtom.y, basisAtom.z);
        frac = wrapFractional(frac);

        for (int op = 0; op < operationCount; ++op)
        {
            glm::dvec3 transformed(0.0);
            for (int row = 0; row < 3; ++row)
            {
                transformed[row] = translations[op][row];
                for (int col = 0; col < 3; ++col)
                    transformed[row] += rotations[op][row][col] * frac[col];
            }
            transformed = wrapFractional(transformed);

            bool duplicate = false;
            for (size_t existingIndex = 0; existingIndex < generatedAtoms.size(); ++existingIndex)
            {
                if (generatedAtoms[existingIndex].atomicNumber != basisAtom.atomicNumber)
                    continue;
                if (sameFractionalPosition(transformed, generatedFrac[existingIndex], 1e-5))
                {
                    duplicate = true;
                    break;
                }
            }
            if (duplicate)
                continue;

            const glm::dvec3 cart = a * transformed.x + b * transformed.y + c * transformed.z;

            AtomSite generated;
            generated.atomicNumber = basisAtom.atomicNumber;
            generated.symbol = elementSymbol(generated.atomicNumber);
            generated.x = transformed.x;
            generated.y = transformed.y;
            generated.z = transformed.z;
            if (generated.atomicNumber >= 0 && generated.atomicNumber < (int)elementColors.size())
            {
                generated.r = elementColors[generated.atomicNumber].r;
                generated.g = elementColors[generated.atomicNumber].g;
                generated.b = elementColors[generated.atomicNumber].b;
            }
            else
            {
                getDefaultElementColor(generated.atomicNumber, generated.r, generated.g, generated.b);
            }

            generated.x = cart.x;
            generated.y = cart.y;
            generated.z = cart.z;
            generatedFrac.push_back(transformed);
            generatedAtoms.push_back(generated);
        }
    }

    if (generatedAtoms.empty())
    {
        result.message = "Symmetry expansion produced no atoms.";
        return result;
    }

    structure.atoms.swap(generatedAtoms);
    structure.hasUnitCell = true;
    structure.cellOffset = {{0.0, 0.0, 0.0}};
    structure.cellVectors = {{
        {{ a.x, a.y, a.z }},
        {{ b.x, b.y, b.z }},
        {{ c.x, c.y, c.z }}
    }};

    result.generatedAtoms = (int)structure.atoms.size();
    result.success = true;
    result.message = "Bulk crystal generated from asymmetric-unit atoms and selected space-group symmetry.";
    return result;
#endif
}

void addDefaultAsymmetricAtom(std::vector<AtomSite>& atoms,
                             const std::vector<glm::vec3>& elementColors)
{
    AtomSite atom;
    applyElementToAtom(atom, 1, elementColors);
    atom.x = 0.0;
    atom.y = 0.0;
    atom.z = 0.0;
    atoms.push_back(atom);
}
