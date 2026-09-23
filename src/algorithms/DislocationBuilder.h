#pragma once

#include "model/Structure.h"

#include <glm/glm.hpp>

#include <string>
#include <vector>

enum class DislocationLatticeFamily
{
    Unknown = 0,
    Fcc,
    Hcp,
    Bcc,
};

enum class DislocationCharacter
{
    Edge = 0,
    Screw,
    Mixed,
};

enum class DislocationShape
{
    HalfPlane = 0,
    Cylinder,
    Sphere,
    Ellipsoid,
    Freeform2D,
};

struct DislocationDetectionResult
{
    bool success = false;
    std::string message;
    DislocationLatticeFamily family = DislocationLatticeFamily::Unknown;
    int fccCount = 0;
    int hcpCount = 0;
    int bccCount = 0;
    int recognizedCount = 0;
};

struct DislocationParams
{
    DislocationCharacter character = DislocationCharacter::Edge;
    DislocationShape shape = DislocationShape::HalfPlane;

    bool autoDirections = true;
    glm::ivec3 planeHkl = glm::ivec3(1, 1, 1);
    glm::ivec3 burgersUvw = glm::ivec3(1, -1, 0);
    glm::ivec3 lineUvw = glm::ivec3(1, 1, -2);

    bool useFractionalLinePoint = true;
    glm::vec3 linePointFractional = glm::vec3(0.5f, 0.5f, 0.5f);
    glm::vec3 linePointCartesianOffset = glm::vec3(0.0f);

    float burgersScale = 1.0f;
    float burgersOverrideMagnitude = -1.0f; // if > 0, use directly (Å) instead of auto-detecting
    float mixedCharacterAngleDeg = 45.0f;
    float poissonRatio = 0.33f;
    float coreRadius = 1.2f;
    float cutoffRadius = 20.0f;
    float lineHalfLength = 1000000.0f;

    // Anisotropic elasticity (Stroh sextic formalism; see
    // src/algorithms/AnisotropicDislocation.h), used in place of the
    // isotropic (poissonRatio-based) displacement field when enabled. The
    // crystal-frame elastic constants below are rotated into the
    // dislocation's local frame automatically; they are expected in the
    // SAME crystallographic axes as the input structure's unit cell
    // (x=a, y=b, z=c for an orthogonal cell), matching BABEL's
    // cubic_elasticity/hexagonal_elasticity/anisotropic_elasticity
    // conventions.
    bool anisotropicElasticity = false;
    enum class ElasticSymmetry { Cubic, Hexagonal, General } elasticSymmetry = ElasticSymmetry::Cubic;
    double elasticC11 = 0.0, elasticC12 = 0.0, elasticC13 = 0.0, elasticC33 = 0.0, elasticC44 = 0.0; // GPa
    double elasticVoigt[6][6] = {}; // GPa, used when elasticSymmetry == General
    // Relative random perturbation applied to the elastic constants before
    // solving, needed only for elastically isotropic or high-symmetry
    // orientations where the Stroh sextic is otherwise degenerate (matches
    // BABEL's documented CVoigt_noise workaround).
    double elasticNoiseAmplitude = 1e-4;
    unsigned elasticNoiseSeed = 1;

    // Dislocation dipole (BABEL's DDipole namelist): superposes a second
    // dislocation with the OPPOSITE Burgers vector, offset from the first
    // by dipoleOffset (in the local (e1, e2) slip-plane axes, Angstrom).
    // A dipole's net Burgers vector is zero, so -- unlike a single
    // dislocation -- it can be genuinely periodic: its long-range
    // displacement field is compatible with the original supercell's
    // periodicity, matching real periodic-array simulation cells.
    bool dipole = false;
    glm::vec2 dipoleOffset = glm::vec2(15.0f, 0.0f);

    float cylinderRadius = 12.0f;
    float sphereRadius = 12.0f;
    glm::vec3 ellipsoidRadii = glm::vec3(14.0f, 10.0f, 10.0f);

    std::vector<glm::vec2> freeformPoints = {
        glm::vec2(-8.0f, 0.0f),
        glm::vec2(8.0f, 0.0f),
        glm::vec2(8.0f, 8.0f),
        glm::vec2(-8.0f, 8.0f),
    };
};

struct DislocationValidationReport
{
    bool passed = false;
    bool finiteCoordinates = true;
    bool latticeFamilyPreserved = false;
    bool minDistanceSafe = true;

    DislocationLatticeFamily familyBefore = DislocationLatticeFamily::Unknown;
    DislocationLatticeFamily familyAfter = DislocationLatticeFamily::Unknown;

    float minInteratomicDistance = 0.0f;
    float rmsDisplacement = 0.0f;
    float maxDisplacement = 0.0f;

    std::string message;
};

struct DislocationResult
{
    bool success = false;
    std::string message;

    Structure output;
    DislocationLatticeFamily family = DislocationLatticeFamily::Unknown;

    int shiftedAtomCount = 0;
    float burgersMagnitude = 0.0f;
    glm::vec3 linePoint = glm::vec3(0.0f);
    glm::vec3 lineDirection = glm::vec3(0.0f, 0.0f, 1.0f);
    glm::vec3 burgersDirection = glm::vec3(1.0f, 0.0f, 0.0f);
    glm::vec3 e1 = glm::vec3(1.0f, 0.0f, 0.0f);  // Local frame e1 (edge displacement axis)
    glm::vec3 e2 = glm::vec3(0.0f, 1.0f, 0.0f);  // Local frame e2 (perpendicular to e1 and lineDir)
    std::vector<glm::vec3> loopPoints;  // Dislocation loop boundary (polygon vertices)

    DislocationValidationReport validation;
};

const char* dislocationLatticeFamilyName(DislocationLatticeFamily family);

DislocationDetectionResult detectDislocationLattice(const Structure& structure,
                                                    bool usePbcForDetection = true);

DislocationResult buildDislocation(const Structure& base,
                                   const DislocationParams& params);
