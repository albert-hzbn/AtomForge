#include "StructureInstanceBuilder.h"

#include <glm/gtc/matrix_inverse.hpp>

#include <algorithm>
#include <cmath>

namespace
{
constexpr float kEpsilon = 1e-5f;
constexpr float kPbcBoundaryTol = 1e-4f;
constexpr float kAtomVisualScale = 0.60f;

float wrapAndSnapFractional(float value)
{
    value -= std::floor(value);

    if (std::abs(value) <= kPbcBoundaryTol)
        return 0.0f;
    if (std::abs(1.0f - value) <= kPbcBoundaryTol)
        return 0.0f;

    return value;
}

void appendBoxEdges(const glm::vec3 (&corners)[8], std::vector<glm::vec3>& boxLines)
{
    static const int edges[12][2] = {
        {0,1},{1,2},{2,3},{3,0},
        {4,5},{5,6},{6,7},{7,4},
        {0,4},{1,5},{2,6},{3,7}
    };

    for (const auto& edge : edges)
    {
        boxLines.push_back(corners[edge[0]]);
        boxLines.push_back(corners[edge[1]]);
    }
}

void appendPbcBoundaryImages(const Structure& structure,
                             std::vector<glm::vec3>& positions,
                             std::vector<glm::vec3>& colors,
                             std::vector<float>& scales,
                             std::vector<float>& bondRadii,
                             std::vector<float>& shininess,
                             std::vector<int>& atomicNumbers,
                             std::vector<int>& atomIndices)
{
    if (!structure.hasUnitCell || positions.empty())
        return;

    glm::vec3 origin((float)structure.cellOffset[0],
                     (float)structure.cellOffset[1],
                     (float)structure.cellOffset[2]);

    glm::vec3 a((float)structure.cellVectors[0][0], (float)structure.cellVectors[0][1], (float)structure.cellVectors[0][2]);
    glm::vec3 b((float)structure.cellVectors[1][0], (float)structure.cellVectors[1][1], (float)structure.cellVectors[1][2]);
    glm::vec3 c((float)structure.cellVectors[2][0], (float)structure.cellVectors[2][1], (float)structure.cellVectors[2][2]);

    glm::mat3 cellMat(a, b, c);
    float det = glm::determinant(cellMat);
    if (std::abs(det) <= kEpsilon)
        return;

    glm::mat3 invCellMat = glm::inverse(cellMat);

    const size_t baseCount = positions.size();
    for (size_t i = 0; i < baseCount; ++i)
    {
        const glm::vec3 pos = positions[i];
        glm::vec3 frac = invCellMat * (pos - origin);
        frac.x = wrapAndSnapFractional(frac.x);
        frac.y = wrapAndSnapFractional(frac.y);
        frac.z = wrapAndSnapFractional(frac.z);

        const glm::vec3 canonicalPos = origin + frac.x * a + frac.y * b + frac.z * c;
        positions[i] = canonicalPos;

        std::vector<int> shiftsX = {0};
        std::vector<int> shiftsY = {0};
        std::vector<int> shiftsZ = {0};

        // Mirror atoms on the low boundary to the high boundary.
        // After load-time wrapping, this handles face/edge/vertex continuity
        // without creating periodic images outside the displayed cell.
        if (std::abs(frac.x) <= kPbcBoundaryTol) shiftsX.push_back(1);
        if (std::abs(frac.y) <= kPbcBoundaryTol) shiftsY.push_back(1);
        if (std::abs(frac.z) <= kPbcBoundaryTol) shiftsZ.push_back(1);

        for (int sx : shiftsX)
        {
            for (int sy : shiftsY)
            {
                for (int sz : shiftsZ)
                {
                    if (sx == 0 && sy == 0 && sz == 0)
                        continue;

                    glm::vec3 imagePos = canonicalPos + (float)sx * a + (float)sy * b + (float)sz * c;
                    positions.push_back(imagePos);
                    colors.push_back(colors[i]);
                    scales.push_back(scales[i]);
                    bondRadii.push_back(bondRadii[i]);
                    shininess.push_back(shininess[i]);
                    atomicNumbers.push_back(atomicNumbers[i]);
                    atomIndices.push_back(atomIndices[i]);
                }
            }
        }
    }
}
}

StructureInstanceData buildStructureInstanceData(
    const Structure& structure,
    bool useTransformMatrix,
    const int (&transformMatrix)[3][3],
    const std::vector<float>& elementRadii,
    const std::vector<float>& elementShininess)
{
    StructureInstanceData data;

    data.positions.reserve(structure.atoms.size());
    data.colors.reserve(structure.atoms.size());
    data.scales.reserve(structure.atoms.size());
    data.bondRadii.reserve(structure.atoms.size());
    data.shininess.reserve(structure.atoms.size());
    data.atomicNumbers.reserve(structure.atoms.size());

    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        const auto& atom = structure.atoms[i];
        data.positions.emplace_back((float)atom.x, (float)atom.y, (float)atom.z);
        data.colors.emplace_back(atom.r, atom.g, atom.b);
        float radius = 1.0f;
        if (atom.atomicNumber >= 0 && atom.atomicNumber < (int)elementRadii.size() && elementRadii[atom.atomicNumber] > 0.0f)
            radius = elementRadii[atom.atomicNumber];
        data.scales.push_back(radius * kAtomVisualScale);
        data.bondRadii.push_back(radius);
        float atomShine = 32.0f;
        if (atom.atomicNumber >= 0 && atom.atomicNumber < (int)elementShininess.size())
            atomShine = elementShininess[atom.atomicNumber];
        data.shininess.push_back(atomShine);
        data.atomicNumbers.push_back(atom.atomicNumber);
        data.atomIndices.push_back(i);
    }

    if (useTransformMatrix && structure.hasUnitCell)
    {
        std::vector<glm::vec3> basePositions = data.positions;
        std::vector<glm::vec3> baseColors    = data.colors;
        std::vector<float>     baseScales    = data.scales;
        std::vector<float>     baseBondRadii = data.bondRadii;
        std::vector<float>     baseShininess = data.shininess;
        std::vector<int>       baseAtomicNumbers = data.atomicNumbers;
        std::vector<int>       baseIndices   = data.atomIndices;

        data.positions.clear();
        data.colors.clear();
        data.scales.clear();
        data.bondRadii.clear();
        data.shininess.clear();
        data.atomicNumbers.clear();
        data.atomIndices.clear();

        glm::vec3 origin((float)structure.cellOffset[0],
                         (float)structure.cellOffset[1],
                         (float)structure.cellOffset[2]);

        glm::vec3 a((float)structure.cellVectors[0][0], (float)structure.cellVectors[0][1], (float)structure.cellVectors[0][2]);
        glm::vec3 b((float)structure.cellVectors[1][0], (float)structure.cellVectors[1][1], (float)structure.cellVectors[1][2]);
        glm::vec3 c((float)structure.cellVectors[2][0], (float)structure.cellVectors[2][1], (float)structure.cellVectors[2][2]);

        glm::mat3 cellMat(a, b, c);
        glm::mat3 invCellMat = glm::inverse(cellMat);

        glm::mat3 transform(1.0f);
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                transform[col][row] = (float)transformMatrix[row][col];

        float det = glm::determinant(transform);
        if (std::abs(det) <= kEpsilon)
        {
            data.positions   = std::move(basePositions);
            data.colors      = std::move(baseColors);
            data.scales      = std::move(baseScales);
            data.bondRadii   = std::move(baseBondRadii);
            data.shininess   = std::move(baseShininess);
            data.atomicNumbers = std::move(baseAtomicNumbers);
            data.atomIndices = std::move(baseIndices);
        }
        else
        {
            glm::mat3 invTransform = glm::inverse(transform);

            int expectedCopies = std::max(1, (int)std::round(std::abs(det)));
            data.positions.reserve(basePositions.size() * (size_t)expectedCopies);
            data.colors.reserve(baseColors.size() * (size_t)expectedCopies);
            data.scales.reserve(baseScales.size() * (size_t)expectedCopies);
            data.bondRadii.reserve(baseBondRadii.size() * (size_t)expectedCopies);
            data.shininess.reserve(baseShininess.size() * (size_t)expectedCopies);
            data.atomicNumbers.reserve(baseAtomicNumbers.size() * (size_t)expectedCopies);

            int nMin[3] = {0, 0, 0};
            int nMax[3] = {0, 0, 0};
            for (int row = 0; row < 3; ++row)
            {
                int rowMin = 0;
                int rowMax = 0;
                for (int col = 0; col < 3; ++col)
                {
                    if (transformMatrix[row][col] < 0)
                        rowMin += transformMatrix[row][col];
                    if (transformMatrix[row][col] > 0)
                        rowMax += transformMatrix[row][col];
                }

                nMin[row] = rowMin - 1;
                nMax[row] = rowMax + 1;
            }

            for (size_t atomIndex = 0; atomIndex < basePositions.size(); ++atomIndex)
            {
                glm::vec3 frac = invCellMat * (basePositions[atomIndex] - origin);

                frac.x -= std::floor(frac.x);
                frac.y -= std::floor(frac.y);
                frac.z -= std::floor(frac.z);

                for (int ix = nMin[0]; ix <= nMax[0]; ++ix)
                {
                    for (int iy = nMin[1]; iy <= nMax[1]; ++iy)
                    {
                        for (int iz = nMin[2]; iz <= nMax[2]; ++iz)
                        {
                            glm::vec3 shiftedFrac = frac + glm::vec3((float)ix, (float)iy, (float)iz);
                            glm::vec3 mappedFrac = invTransform * shiftedFrac;

                            if (mappedFrac.x < -kEpsilon || mappedFrac.y < -kEpsilon || mappedFrac.z < -kEpsilon)
                                continue;
                            if (mappedFrac.x >= 1.0f - kEpsilon || mappedFrac.y >= 1.0f - kEpsilon || mappedFrac.z >= 1.0f - kEpsilon)
                                continue;

                            glm::vec3 worldPos = origin + shiftedFrac.x * a + shiftedFrac.y * b + shiftedFrac.z * c;
                            data.positions.push_back(worldPos);
                            data.colors.push_back(baseColors[atomIndex]);
                            data.scales.push_back(baseScales[atomIndex]);
                            data.bondRadii.push_back(baseBondRadii[atomIndex]);
                            data.shininess.push_back(baseShininess[atomIndex]);
                            data.atomicNumbers.push_back(baseAtomicNumbers[atomIndex]);
                            data.atomIndices.push_back(baseIndices[atomIndex]);
                        }
                    }
                }
            }
        }
    }

    // For periodic display, duplicate boundary atoms across neighboring cells
    // so atoms are visible on opposite faces/edges/vertices of the unit cell.
    if (!useTransformMatrix && structure.hasUnitCell)
        appendPbcBoundaryImages(structure, data.positions, data.colors, data.scales, data.bondRadii, data.shininess, data.atomicNumbers, data.atomIndices);

    if (data.positions.empty())
        return data;

    glm::vec3 minPos(1e9f), maxPos(-1e9f);
    for (const auto& p : data.positions)
    {
        minPos = glm::min(minPos, p);
        maxPos = glm::max(maxPos, p);
    }

    glm::vec3 corners[8];

    if (structure.hasUnitCell)
    {
        glm::vec3 origin((float)structure.cellOffset[0],
                         (float)structure.cellOffset[1],
                         (float)structure.cellOffset[2]);

        glm::vec3 a((float)structure.cellVectors[0][0], (float)structure.cellVectors[0][1], (float)structure.cellVectors[0][2]);
        glm::vec3 b((float)structure.cellVectors[1][0], (float)structure.cellVectors[1][1], (float)structure.cellVectors[1][2]);
        glm::vec3 c((float)structure.cellVectors[2][0], (float)structure.cellVectors[2][1], (float)structure.cellVectors[2][2]);

        if (useTransformMatrix)
        {
            glm::vec3 aT = (float)transformMatrix[0][0] * a + (float)transformMatrix[0][1] * b + (float)transformMatrix[0][2] * c;
            glm::vec3 bT = (float)transformMatrix[1][0] * a + (float)transformMatrix[1][1] * b + (float)transformMatrix[1][2] * c;
            glm::vec3 cT = (float)transformMatrix[2][0] * a + (float)transformMatrix[2][1] * b + (float)transformMatrix[2][2] * c;

            corners[0] = origin;
            corners[1] = origin + aT;
            corners[2] = origin + aT + bT;
            corners[3] = origin + bT;
            corners[4] = origin + cT;
            corners[5] = origin + aT + cT;
            corners[6] = origin + aT + bT + cT;
            corners[7] = origin + bT + cT;
        }
        else
        {
            corners[0] = origin;
            corners[1] = origin + a;
            corners[2] = origin + a + b;
            corners[3] = origin + b;
            corners[4] = origin + c;
            corners[5] = origin + a + c;
            corners[6] = origin + a + b + c;
            corners[7] = origin + b + c;
        }
    }
    else
    {
        corners[0] = {minPos.x, minPos.y, minPos.z};
        corners[1] = {maxPos.x, minPos.y, minPos.z};
        corners[2] = {maxPos.x, maxPos.y, minPos.z};
        corners[3] = {minPos.x, maxPos.y, minPos.z};
        corners[4] = {minPos.x, minPos.y, maxPos.z};
        corners[5] = {maxPos.x, minPos.y, maxPos.z};
        corners[6] = {maxPos.x, maxPos.y, maxPos.z};
        corners[7] = {minPos.x, maxPos.y, maxPos.z};
    }

    data.orbitCenter = 0.5f * (corners[0] + corners[6]);
    appendBoxEdges(corners, data.boxLines);

    return data;
}

    Structure buildSupercell(const Structure& structure, const int (&transformMatrix)[3][3])
    {
        if (!structure.hasUnitCell)
            return structure;

        glm::vec3 origin((float)structure.cellOffset[0],
                         (float)structure.cellOffset[1],
                         (float)structure.cellOffset[2]);

        glm::vec3 a((float)structure.cellVectors[0][0], (float)structure.cellVectors[0][1], (float)structure.cellVectors[0][2]);
        glm::vec3 b((float)structure.cellVectors[1][0], (float)structure.cellVectors[1][1], (float)structure.cellVectors[1][2]);
        glm::vec3 c((float)structure.cellVectors[2][0], (float)structure.cellVectors[2][1], (float)structure.cellVectors[2][2]);

        glm::mat3 cellMat(a, b, c);
        glm::mat3 invCellMat = glm::inverse(cellMat);

        glm::mat3 transform(1.0f);
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                transform[col][row] = (float)transformMatrix[row][col];

        float det = glm::determinant(transform);
        if (std::abs(det) <= kEpsilon)
            return structure;

        glm::mat3 invTransform = glm::inverse(transform);

        int nMin[3] = {0, 0, 0};
        int nMax[3] = {0, 0, 0};
        for (int row = 0; row < 3; ++row)
        {
            int rowMin = 0, rowMax = 0;
            for (int col = 0; col < 3; ++col)
            {
                if (transformMatrix[row][col] < 0) rowMin += transformMatrix[row][col];
                if (transformMatrix[row][col] > 0) rowMax += transformMatrix[row][col];
            }
            nMin[row] = rowMin - 1;
            nMax[row] = rowMax + 1;
        }

        Structure result;
        result.hasUnitCell = true;
        result.cellOffset  = structure.cellOffset;

        glm::vec3 aT = (float)transformMatrix[0][0]*a + (float)transformMatrix[0][1]*b + (float)transformMatrix[0][2]*c;
        glm::vec3 bT = (float)transformMatrix[1][0]*a + (float)transformMatrix[1][1]*b + (float)transformMatrix[1][2]*c;
        glm::vec3 cT = (float)transformMatrix[2][0]*a + (float)transformMatrix[2][1]*b + (float)transformMatrix[2][2]*c;
        result.cellVectors[0] = { (double)aT.x, (double)aT.y, (double)aT.z };
        result.cellVectors[1] = { (double)bT.x, (double)bT.y, (double)bT.z };
        result.cellVectors[2] = { (double)cT.x, (double)cT.y, (double)cT.z };

        for (const auto& atom : structure.atoms)
        {
            glm::vec3 worldPos((float)atom.x, (float)atom.y, (float)atom.z);
            glm::vec3 frac = invCellMat * (worldPos - origin);
            frac.x -= std::floor(frac.x);
            frac.y -= std::floor(frac.y);
            frac.z -= std::floor(frac.z);

            for (int ix = nMin[0]; ix <= nMax[0]; ++ix)
            for (int iy = nMin[1]; iy <= nMax[1]; ++iy)
            for (int iz = nMin[2]; iz <= nMax[2]; ++iz)
            {
                glm::vec3 shiftedFrac = frac + glm::vec3((float)ix, (float)iy, (float)iz);
                glm::vec3 mappedFrac  = invTransform * shiftedFrac;

                if (mappedFrac.x < -kEpsilon || mappedFrac.y < -kEpsilon || mappedFrac.z < -kEpsilon)
                    continue;
                if (mappedFrac.x >= 1.0f - kEpsilon || mappedFrac.y >= 1.0f - kEpsilon || mappedFrac.z >= 1.0f - kEpsilon)
                    continue;

                glm::vec3 newWorld = origin + shiftedFrac.x*a + shiftedFrac.y*b + shiftedFrac.z*c;
                AtomSite site = atom;
                site.x = (double)newWorld.x;
                site.y = (double)newWorld.y;
                site.z = (double)newWorld.z;
                result.atoms.push_back(site);
            }
        }

        return result;
    }
