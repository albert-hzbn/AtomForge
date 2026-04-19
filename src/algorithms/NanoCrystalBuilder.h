#pragma once

#include "io/StructureLoader.h"

#include <glm/glm.hpp>
#include <string>
#include <vector>

// -- Shape definitions -------------------------------------------------------

enum class NanoShape
{
    Sphere = 0,
    Ellipsoid,
    Box,
    Cylinder,
    Octahedron,
    TruncatedOctahedron,
    Cuboctahedron,
    MeshModel,
};

constexpr int kNumShapes = 8;

enum class NanoGenerationMode
{
    Shape = 0,
    WulffConstruction,
};

struct WulffPlaneInput
{
    int   h = 1;
    int   k = 0;
    int   l = 0;
    float surfaceEnergy = 1.0f;
};

struct WulffPlane
{
    int familyIndex = -1;
    int h = 0;
    int k = 0;
    int l = 0;
    float distance = 0.0f;
    glm::vec3 normal = glm::vec3(0.0f, 0.0f, 1.0f);
};

struct WulffFace
{
    int familyIndex = -1;
    int h = 0;
    int k = 0;
    int l = 0;
    float distance = 0.0f;
    glm::vec3 normal = glm::vec3(0.0f, 0.0f, 1.0f);
    std::vector<glm::vec3> vertices;
};

struct WulffPreview
{
    bool success = false;
    std::string message;
    glm::vec3 center = glm::vec3(0.0f);
    glm::vec3 minPoint = glm::vec3(0.0f);
    glm::vec3 maxPoint = glm::vec3(0.0f);
    float maxPlaneDistance = 0.0f;
    float boundingRadius = 0.0f;
    std::vector<WulffPlane> planes;
    std::vector<WulffFace> faces;
};

const char* shapeLabel(NanoShape s);

struct NanoParams
{
    NanoGenerationMode generationMode = NanoGenerationMode::Shape;
    NanoShape shape = NanoShape::Sphere;

    float sphereRadius = 15.0f;

    float ellipRx = 15.0f;
    float ellipRy = 12.0f;
    float ellipRz = 10.0f;

    float boxHx = 15.0f;
    float boxHy = 15.0f;
    float boxHz = 15.0f;

    float cylRadius = 12.0f;
    float cylHeight = 30.0f;
    int   cylAxis   = 2;

    float octRadius = 15.0f;

    float truncOctRadius = 18.0f;
    float truncOctTrunc  = 12.0f;

    float cuboRadius = 15.0f;

    // Bounds of a loaded 3D model after centering at origin; maintained by UI.
    float modelHx = 15.0f;
    float modelHy = 15.0f;
    float modelHz = 15.0f;
    float modelScale = 1.0f;

    bool  autoCenterFromAtoms = true;
    float cx = 0.0f;
    float cy = 0.0f;
    float cz = 0.0f;

    // Optional crystal orientation applied before shape clipping.
    bool  applyCrystalOrientation = false;
    bool  useMillerOrientation = false;
    float orientXDeg = 0.0f;
    float orientYDeg = 0.0f;
    float orientZDeg = 0.0f;
    int   millerH = 1;
    int   millerK = 0;
    int   millerL = 0;

    bool autoReplicate = true;
    int  repA = 5;
    int  repB = 5;
    int  repC = 5;

    bool  setOutputCell = true;
    float vacuumPadding = 5.0f;

    std::vector<WulffPlaneInput> wulffPlanes;
    float wulffMaxRadius = 20.0f;
};

struct NanoBuildResult
{
    bool        success           = false;
    std::string message;
    NanoGenerationMode mode       = NanoGenerationMode::Shape;
    int         inputAtoms        = 0;
    int         outputAtoms       = 0;
    NanoShape   shape             = NanoShape::Sphere;
    float       estimatedDiameter = 0.0f;
    int         repA = 0, repB = 0, repC = 0;
    bool        tilingUsed  = false;
    bool        repClamped  = false;
    int         wulffFamilyCount = 0;
    int         wulffFaceCount = 0;
};

struct HalfExtents { float hx, hy, hz; };

// -- Geometry helpers --------------------------------------------------------

float computeBoundingRadius(const NanoParams& p);
HalfExtents computeShapeHalfExtents(const NanoParams& p);
bool isInsideShape(const glm::vec3& p, const NanoParams& params,
                   const std::vector<glm::vec3>& modelVertices = {},
                   const std::vector<unsigned int>& modelIndices = {});
glm::vec3 computeAtomCentroid(const std::vector<AtomSite>& atoms);
[[nodiscard]] WulffPreview computeWulffPreview(const Structure& reference,
                                               const NanoParams& params);

// -- Builder -----------------------------------------------------------------

[[nodiscard]] NanoBuildResult buildNanocrystal(Structure& structure,
                                 const Structure& reference,
                                 const NanoParams& params,
                                 const std::vector<glm::vec3>& elementColors,
                                 const std::vector<glm::vec3>& modelVertices,
                                 const std::vector<unsigned int>& modelIndices);
