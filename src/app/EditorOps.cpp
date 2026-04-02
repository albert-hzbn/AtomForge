#include "app/EditorOps.h"

#include "Camera.h"
#include "ElementData.h"
#include "graphics/StructureInstanceBuilder.h"
#include "math/StructureMath.h"

#include <algorithm>
#include <cmath>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace
{
struct Vec3iHash
{
    size_t operator()(const glm::ivec3& value) const
    {
        return ((size_t)value.x * 73856093u)
             ^ ((size_t)value.y * 19349663u)
             ^ ((size_t)value.z * 83492791u);
    }
};

struct Signature
{
    int common = 0;
    int bonds = 0;
    int chain = 0;

    bool operator<(const Signature& other) const
    {
        if (common != other.common) return common < other.common;
        if (bonds != other.bonds) return bonds < other.bonds;
        return chain < other.chain;
    }

    bool operator==(const Signature& other) const
    {
        return common == other.common && bonds == other.bonds && chain == other.chain;
    }
};

glm::ivec3 getGridCell(const glm::vec3& position)
{
    constexpr float kSpatialHashCellSize = 4.0f;
    return glm::ivec3((int)std::floor(position.x / kSpatialHashCellSize),
                      (int)std::floor(position.y / kSpatialHashCellSize),
                      (int)std::floor(position.z / kSpatialHashCellSize));
}

uint64_t edgeKey(int a, int b)
{
    if (a > b) std::swap(a, b);
    return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
}

int longestChainLength(const std::vector<int>& commonNodes,
                       const std::unordered_set<uint64_t>& edgeSet)
{
    if (commonNodes.empty())
        return 0;
    if (commonNodes.size() == 1)
        return 1;

    std::vector<std::vector<int>> adjacency(commonNodes.size());
    for (int i = 0; i < (int)commonNodes.size(); ++i)
    {
        for (int j = i + 1; j < (int)commonNodes.size(); ++j)
        {
            if (edgeSet.find(edgeKey(commonNodes[i], commonNodes[j])) != edgeSet.end())
            {
                adjacency[i].push_back(j);
                adjacency[j].push_back(i);
            }
        }
    }

    int best = 1;
    for (int src = 0; src < (int)commonNodes.size(); ++src)
    {
        std::vector<int> dist(commonNodes.size(), -1);
        std::queue<int> q;
        dist[src] = 0;
        q.push(src);

        while (!q.empty())
        {
            const int u = q.front();
            q.pop();
            for (int v : adjacency[u])
            {
                if (dist[v] >= 0)
                    continue;
                dist[v] = dist[u] + 1;
                best = std::max(best, dist[v] + 1);
                q.push(v);
            }
        }
    }

    return best;
}

bool isCrystalLikeEnvironment(const Signature& s)
{
    Signature fcc;  fcc.common = 4; fcc.bonds = 2; fcc.chain = 1;
    Signature hcp;  hcp.common = 4; hcp.bonds = 2; hcp.chain = 2;
    Signature bccA; bccA.common = 4; bccA.bonds = 4; bccA.chain = 1;
    Signature bccB; bccB.common = 6; bccB.bonds = 6; bccB.chain = 1;
    Signature ico;  ico.common = 5; ico.bonds = 5; ico.chain = 1;
    return s == fcc || s == hcp || s == bccA || s == bccB || s == ico;
}

void applyGrainBoundaryColors(Structure& structure)
{
    if (structure.atoms.empty())
        return;

    const glm::vec3 crystalColor(0.20f, 0.72f, 0.98f);
    const glm::vec3 gbColor(1.00f, 0.45f, 0.12f);

    if (structure.atoms.size() < 4)
    {
        for (auto& atom : structure.atoms)
        {
            atom.r = crystalColor.r;
            atom.g = crystalColor.g;
            atom.b = crystalColor.b;
        }
        return;
    }

    std::vector<float> radii = makeLiteratureCovalentRadii();
    glm::mat3 cell(1.0f), invCell(1.0f);
    const bool usePbc = structure.hasUnitCell && tryMakeCellMatrices(structure, cell, invCell);

    std::vector<glm::vec3> positions(structure.atoms.size());
    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        positions[i] = glm::vec3((float)structure.atoms[i].x,
                                 (float)structure.atoms[i].y,
                                 (float)structure.atoms[i].z);
    }

    std::vector<std::vector<int>> neighbors(structure.atoms.size());
    std::unordered_set<uint64_t> edgeSet;

    auto tryAddBond = [&](int i, int j)
    {
        if (j <= i)
            return;

        const int zi = structure.atoms[i].atomicNumber;
        const int zj = structure.atoms[j].atomicNumber;
        const float ri = (zi >= 0 && zi < (int)radii.size()) ? radii[zi] : 1.0f;
        const float rj = (zj >= 0 && zj < (int)radii.size()) ? radii[zj] : 1.0f;

        const glm::vec3 delta = minimumImageDelta(positions[j] - positions[i], usePbc, cell, invCell);
        const float d = glm::length(delta);
        if (d <= kMinBondDistance)
            return;

        const float cutoff = (ri + rj) * kBondToleranceFactor;
        if (d > cutoff)
            return;

        neighbors[i].push_back(j);
        neighbors[j].push_back(i);
        edgeSet.insert(edgeKey(i, j));
    };

    if (!usePbc)
    {
        std::unordered_map<glm::ivec3, std::vector<int>, Vec3iHash> grid;
        grid.reserve(positions.size());
        for (int i = 0; i < (int)positions.size(); ++i)
            grid[getGridCell(positions[i])].push_back(i);

        for (int i = 0; i < (int)positions.size(); ++i)
        {
            const glm::ivec3 cellCoord = getGridCell(positions[i]);
            for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz)
            {
                const glm::ivec3 neighborCell(cellCoord.x + dx, cellCoord.y + dy, cellCoord.z + dz);
                auto it = grid.find(neighborCell);
                if (it == grid.end())
                    continue;

                for (int idx : it->second)
                    tryAddBond(i, idx);
            }
        }
    }
    else
    {
        for (int i = 0; i < (int)structure.atoms.size(); ++i)
            for (int j = i + 1; j < (int)structure.atoms.size(); ++j)
                tryAddBond(i, j);
    }

    for (auto& nbrs : neighbors)
        std::sort(nbrs.begin(), nbrs.end());

    std::vector<bool> isBoundary(structure.atoms.size(), false);
    int boundarySeedCount = 0;

    if (structure.grainRegionIds.size() == structure.atoms.size())
    {
        for (int i = 0; i < (int)structure.atoms.size(); ++i)
        {
            const int regionId = structure.grainRegionIds[i];
            for (int ni : neighbors[i])
            {
                if (structure.grainRegionIds[ni] != regionId)
                {
                    // Skip bonds that cross a periodic boundary – these connect
                    // atoms on opposite faces of the box and are artefacts of
                    // the periodic Voronoi tessellation, not real interior GBs.
                    if (usePbc)
                    {
                        const glm::vec3 rawDelta = positions[ni] - positions[i];
                        const glm::vec3 minDelta = minimumImageDelta(rawDelta, true, cell, invCell);
                        if (glm::length(rawDelta - minDelta) > 0.1f)
                            continue;
                    }

                    if (!isBoundary[i])
                    {
                        isBoundary[i] = true;
                        ++boundarySeedCount;
                    }
                    if (!isBoundary[ni])
                    {
                        isBoundary[ni] = true;
                        ++boundarySeedCount;
                    }
                }
            }
        }

        if (boundarySeedCount > 0)
        {
            for (int i = 0; i < (int)structure.atoms.size(); ++i)
            {
                const glm::vec3 color = isBoundary[i] ? gbColor : crystalColor;
                structure.atoms[i].r = color.r;
                structure.atoms[i].g = color.g;
                structure.atoms[i].b = color.b;
            }
        }
        else
        {
            for (auto& atom : structure.atoms)
            {
                atom.r = crystalColor.r;
                atom.g = crystalColor.g;
                atom.b = crystalColor.b;
            }
        }
        return;
    }

    std::vector<int> coordination(structure.atoms.size(), 0);
    std::map<int, std::map<int, int>> coordinationHistogramByElement;
    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        coordination[i] = (int)neighbors[i].size();
        ++coordinationHistogramByElement[structure.atoms[i].atomicNumber][coordination[i]];
    }

    std::map<int, int> modalCoordinationByElement;
    for (std::map<int, std::map<int, int>>::const_iterator elementIt = coordinationHistogramByElement.begin();
         elementIt != coordinationHistogramByElement.end(); ++elementIt)
    {
        int modalCoordination = 0;
        int modalCount = -1;
        for (std::map<int, int>::const_iterator coordIt = elementIt->second.begin();
             coordIt != elementIt->second.end(); ++coordIt)
        {
            if (coordIt->second > modalCount)
            {
                modalCoordination = coordIt->first;
                modalCount = coordIt->second;
            }
        }
        modalCoordinationByElement[elementIt->first] = modalCoordination;
    }

    // Identify atoms near the bounding box edges.  Surface atoms naturally
    // have lower coordination and should not be marked as grain-boundary atoms.
    // For periodic structures the unit cell defines the box; for non-periodic
    // structures we fall back to the bounding box of all atom positions.
    std::vector<bool> nearEdge(structure.atoms.size(), false);
    {
        glm::vec3 bbMin, bbMax;
        if (usePbc)
        {
            bbMin = glm::vec3((float)structure.cellOffset[0],
                              (float)structure.cellOffset[1],
                              (float)structure.cellOffset[2]);
            bbMax = bbMin + glm::vec3(cell[0]) + glm::vec3(cell[1]) + glm::vec3(cell[2]);
        }
        else
        {
            bbMin = glm::vec3(std::numeric_limits<float>::max());
            bbMax = glm::vec3(std::numeric_limits<float>::lowest());
            for (const auto& pos : positions)
            {
                bbMin = glm::min(bbMin, pos);
                bbMax = glm::max(bbMax, pos);
            }
        }

        // Maximum possible bonding cutoff (two of the largest covalent radii).
        float maxRadius = 0.0f;
        for (int i = 0; i < (int)structure.atoms.size(); ++i)
        {
            const int z = structure.atoms[i].atomicNumber;
            const float r = (z >= 0 && z < (int)radii.size()) ? radii[z] : 1.0f;
            if (r > maxRadius) maxRadius = r;
        }
        const float edgeCutoff = maxRadius * 2.0f * kBondToleranceFactor;

        for (int i = 0; i < (int)structure.atoms.size(); ++i)
        {
            const glm::vec3& p = positions[i];
            if (p.x - bbMin.x < edgeCutoff || bbMax.x - p.x < edgeCutoff ||
                p.y - bbMin.y < edgeCutoff || bbMax.y - p.y < edgeCutoff ||
                p.z - bbMin.z < edgeCutoff || bbMax.z - p.z < edgeCutoff)
            {
                nearEdge[i] = true;
            }
        }
    }

    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        const int atomicNumber = structure.atoms[i].atomicNumber;
        std::map<int, int>::const_iterator modalIt = modalCoordinationByElement.find(atomicNumber);
        const int modalCoordination = (modalIt != modalCoordinationByElement.end()) ? modalIt->second : coordination[i];

        if (std::abs(coordination[i] - modalCoordination) >= 1 && !nearEdge[i])
        {
            isBoundary[i] = true;
            ++boundarySeedCount;
        }
    }

    if (boundarySeedCount > 0)
    {
        for (int i = 0; i < (int)structure.atoms.size(); ++i)
        {
            const glm::vec3 color = isBoundary[i] ? gbColor : crystalColor;
            structure.atoms[i].r = color.r;
            structure.atoms[i].g = color.g;
            structure.atoms[i].b = color.b;
        }
        return;
    }

    std::vector<std::map<Signature, int>> atomSignatures(structure.atoms.size());
    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        for (int t = 0; t < (int)neighbors[i].size(); ++t)
        {
            const int j = neighbors[i][t];
            if (j <= i)
                continue;

            std::vector<int> common;
            common.reserve(std::min(neighbors[i].size(), neighbors[j].size()));
            std::set_intersection(neighbors[i].begin(), neighbors[i].end(),
                                  neighbors[j].begin(), neighbors[j].end(),
                                  std::back_inserter(common));

            int bondCount = 0;
            for (int a = 0; a < (int)common.size(); ++a)
            {
                for (int b = a + 1; b < (int)common.size(); ++b)
                {
                    if (edgeSet.find(edgeKey(common[a], common[b])) != edgeSet.end())
                        ++bondCount;
                }
            }

            Signature sig;
            sig.common = (int)common.size();
            sig.bonds = bondCount;
            sig.chain = longestChainLength(common, edgeSet);
            ++atomSignatures[i][sig];
            ++atomSignatures[j][sig];
        }
    }

    std::vector<bool> fallbackBoundary(structure.atoms.size(), true);
    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        Signature dominant;
        int dominantCount = 0;
        for (std::map<Signature, int>::const_iterator it = atomSignatures[i].begin(); it != atomSignatures[i].end(); ++it)
        {
            if (it->second > dominantCount)
            {
                dominant = it->first;
                dominantCount = it->second;
            }
        }

        const bool crystalLike = dominantCount > 0 && isCrystalLikeEnvironment(dominant);
        fallbackBoundary[i] = !crystalLike;
    }

    for (int i = 0; i < (int)structure.atoms.size(); ++i)
    {
        const glm::vec3 color = fallbackBoundary[i] ? gbColor : crystalColor;
        structure.atoms[i].r = color.r;
        structure.atoms[i].g = color.g;
        structure.atoms[i].b = color.b;
    }
}

float estimateSceneRadius(const SceneBuffers& sceneBuffers)
{
    float maxRadius = 0.0f;

    if (!sceneBuffers.atomPositions.empty())
    {
        for (size_t i = 0; i < sceneBuffers.atomPositions.size(); ++i)
        {
            const float radius = (i < sceneBuffers.atomRadii.size()) ? sceneBuffers.atomRadii[i] : 0.0f;
            const float distance = glm::length(sceneBuffers.atomPositions[i] - sceneBuffers.orbitCenter) + radius;
            maxRadius = std::max(maxRadius, distance);
        }
    }
    else if (!sceneBuffers.boxLines.empty())
    {
        for (const glm::vec3& point : sceneBuffers.boxLines)
        {
            const float distance = glm::length(point - sceneBuffers.orbitCenter);
            maxRadius = std::max(maxRadius, distance);
        }
    }

    return std::max(maxRadius, 1.0f);
}
}

EditorSnapshot captureSnapshot(const EditorState& state)
{
    EditorSnapshot snapshot;
    snapshot.structure = state.structure;
    snapshot.elementRadii = state.editMenuDialogs.elementRadii;
    snapshot.elementColors = state.editMenuDialogs.elementColors;
    snapshot.elementShininess = state.editMenuDialogs.elementShininess;
    return snapshot;
}

void updateBuffers(EditorState& state)
{
    if (state.fileBrowser.isTransformMatrixEnabled() && state.structure.hasUnitCell)
    {
        const size_t inputAtomCount = state.structure.atoms.size();
        const int (&matrix)[3][3] = state.fileBrowser.getTransformMatrix();
        state.structure = buildSupercell(state.structure, state.fileBrowser.getTransformMatrix());
        std::cout << "[Operation] Applied supercell transform: "
                  << "matrix=[[" << matrix[0][0] << " " << matrix[0][1] << " " << matrix[0][2] << "]"
                  << ",[" << matrix[1][0] << " " << matrix[1][1] << " " << matrix[1][2] << "]"
                  << ",[" << matrix[2][0] << " " << matrix[2][1] << " " << matrix[2][2] << "]]"
                  << ", atoms=" << inputAtomCount << " -> " << state.structure.atoms.size()
                  << std::endl;
        state.fileBrowser.clearTransformMatrix();
    }

    for (auto& atom : state.structure.atoms)
    {
        int atomicNumber = atom.atomicNumber;
        if (atomicNumber >= 0 && atomicNumber < (int)state.editMenuDialogs.elementColors.size())
        {
            atom.r = state.editMenuDialogs.elementColors[atomicNumber].r;
            atom.g = state.editMenuDialogs.elementColors[atomicNumber].g;
            atom.b = state.editMenuDialogs.elementColors[atomicNumber].b;
        }
    }

    // Override with alternative color modes when selected.
    if (state.fileBrowser.getAtomColorMode() == AtomColorMode::CrystalOrientation)
    {
        if (state.structure.grainColors.size() == state.structure.atoms.size())
        {
            for (size_t i = 0; i < state.structure.atoms.size(); ++i)
            {
                state.structure.atoms[i].r = state.structure.grainColors[i][0];
                state.structure.atoms[i].g = state.structure.grainColors[i][1];
                state.structure.atoms[i].b = state.structure.grainColors[i][2];
            }
        }
        else
        {
            // No per-atom grain colors: assume identity orientation (IPF-Z [001] = red)
            for (size_t i = 0; i < state.structure.atoms.size(); ++i)
            {
                state.structure.atoms[i].r = 1.0f;
                state.structure.atoms[i].g = 0.0f;
                state.structure.atoms[i].b = 0.0f;
            }
        }
    }
    else if (state.fileBrowser.getAtomColorMode() == AtomColorMode::GrainBoundary)
    {
        applyGrainBoundaryColors(state.structure);
    }

    StructureInstanceData data = buildStructureInstanceData(
        state.structure,
        state.fileBrowser.isTransformMatrixEnabled(),
        state.fileBrowser.getTransformMatrix(),
        state.editMenuDialogs.elementRadii,
        state.editMenuDialogs.elementShininess);

    state.sceneBuffers.upload(
        data,
        state.fileBrowser.isBondElementFilterEnabled(),
        state.fileBrowser.getBondElementFilterMask());
    state.selectedInstanceIndices.clear();
    state.voronoiDirty = true;

    if (!state.suppressHistoryCommit)
        state.undoRedo.commit(captureSnapshot(state));
}

void updateBuffers(EditorState& state, Structure& structure)
{
    state.structure = structure;
    updateBuffers(state);
    structure = state.structure;
}

void applySnapshot(EditorState& state, const EditorSnapshot& snapshot)
{
    state.structure = snapshot.structure;
    state.editMenuDialogs.elementRadii = snapshot.elementRadii;
    state.editMenuDialogs.elementColors = snapshot.elementColors;
    state.editMenuDialogs.elementShininess = snapshot.elementShininess;
    state.suppressHistoryCommit = true;
    updateBuffers(state);
    state.suppressHistoryCommit = false;
    state.measurementState.clearVisuals();
}

void applyDefaultView(
    Camera& camera,
    const SceneBuffers& sceneBuffers,
    int viewportWidth,
    int viewportHeight,
    bool fitToStructure)
{
    constexpr float kIsoYawDeg = 45.0f;
    constexpr float kIsoPitchDeg = 35.2643897f;

    camera.yaw = kIsoYawDeg;
    camera.pitch = kIsoPitchDeg;
    camera.roll = 0.0f;
    camera.panOffset = glm::vec3(0.0f);

    if (!fitToStructure || sceneBuffers.atomCount == 0)
    {
        camera.distance = 10.0f;
        return;
    }

    float maxRadius = estimateSceneRadius(sceneBuffers);

    float aspect = (viewportHeight > 0) ? (float)viewportWidth / (float)viewportHeight : 1.0f;
    float verticalFov = glm::radians(45.0f);
    float horizontalFov = 2.0f * std::atan(std::tan(verticalFov * 0.5f) * aspect);
    float halfFov = 0.5f * std::min(verticalFov, horizontalFov);
    halfFov = std::max(halfFov, glm::radians(10.0f));

    float framedDistance = (maxRadius / std::sin(halfFov)) * 1.15f;
    camera.distance = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, framedDistance));
}

void clearSelection(EditorState& state)
{
    for (int index : state.selectedInstanceIndices)
        state.sceneBuffers.restoreAtomColor(index);

    state.selectedInstanceIndices.clear();
    state.measurementState.clearVisuals();
}

void deleteSelectedAtoms(EditorState& state)
{
    if (state.selectedInstanceIndices.empty())
        return;

    std::vector<int> baseIndicesToDelete;
    for (int selectedIndex : state.selectedInstanceIndices)
    {
        if (selectedIndex >= 0 && selectedIndex < (int)state.sceneBuffers.atomIndices.size())
        {
            int baseIndex = state.sceneBuffers.atomIndices[selectedIndex];
            if (baseIndex >= 0 && baseIndex < (int)state.structure.atoms.size())
                baseIndicesToDelete.push_back(baseIndex);
        }
    }

    std::sort(baseIndicesToDelete.begin(), baseIndicesToDelete.end(), std::greater<int>());
    baseIndicesToDelete.erase(
        std::unique(baseIndicesToDelete.begin(), baseIndicesToDelete.end()),
        baseIndicesToDelete.end());

    for (int baseIndex : baseIndicesToDelete)
    {
        if (baseIndex >= 0 && baseIndex < (int)state.structure.atoms.size())
            state.structure.atoms.erase(state.structure.atoms.begin() + baseIndex);
    }

    std::cout << "[Operation] Deleted atoms (selection): " << baseIndicesToDelete.size() << std::endl;

    updateBuffers(state);
}

void refreshSelectionHighlights(EditorState& state)
{
    for (int& index : state.selectedInstanceIndices)
    {
        if (index >= (int)state.sceneBuffers.atomCount)
        {
            state.sceneBuffers.restoreAtomColor(index);
            index = -1;
        }
        else if (index >= 0)
        {
            state.sceneBuffers.highlightAtom(index, glm::vec3(1.0f, 1.0f, 0.0f));
        }
    }

    state.selectedInstanceIndices.erase(
        std::remove(state.selectedInstanceIndices.begin(), state.selectedInstanceIndices.end(), -1),
        state.selectedInstanceIndices.end());
}

void rotateCrystalAroundAxis(Camera& camera, int axis, double angleDeg)
{
    // Rotating the crystal by +angleDeg visually equals orbiting the camera by -angleDeg.
    const float rad = glm::radians((float)-angleDeg);
    const glm::vec3 axisVec = (axis == 0) ? glm::vec3(1, 0, 0)
                            : (axis == 1) ? glm::vec3(0, 1, 0)
                                          : glm::vec3(0, 0, 1);
    const glm::mat3 rot = glm::mat3(glm::rotate(glm::mat4(1.0f), rad, axisVec));

    // Reconstruct full camera frame including current roll.
    const float yawR   = glm::radians(camera.yaw);
    const float pitchR = glm::radians(camera.pitch);
    const glm::vec3 forward(
        std::cos(pitchR) * std::sin(yawR),
        std::sin(pitchR),
        std::cos(pitchR) * std::cos(yawR));
    const glm::vec3 canonUp(
        -std::sin(pitchR) * std::sin(yawR),
         std::cos(pitchR),
        -std::sin(pitchR) * std::cos(yawR));
    // right = cross(orbit-forward, canonUp) — matches buildFrameView convention
    const glm::vec3 right = glm::normalize(glm::cross(forward, canonUp));
    const float rollRad   = glm::radians(camera.roll);
    const glm::vec3 upVec = std::cos(rollRad) * canonUp + std::sin(rollRad) * right;

    // Rotate both vectors.
    const glm::vec3 newForward = glm::normalize(rot * forward);
    const glm::vec3 newUp      = glm::normalize(rot * upVec);

    // Extract new yaw/pitch from rotated orbit-forward direction.
    const float newPitch = glm::degrees(std::asin(glm::clamp(newForward.y, -1.0f, 1.0f)));
    const float newYaw   = glm::degrees(std::atan2(newForward.x, newForward.z));

    // Compute roll as the signed angle from the canonical up of the new
    // orientation to the actual up, measured around the new forward axis.
    const float ny = glm::radians(newYaw);
    const float np = glm::radians(newPitch);
    const glm::vec3 newCanonUp(
        -std::sin(np) * std::sin(ny),
         std::cos(np),
        -std::sin(np) * std::cos(ny));
    const float sinRoll = glm::dot(glm::cross(newCanonUp, newUp), newForward);
    const float cosRoll = glm::dot(newCanonUp, newUp);
    const float newRoll = glm::degrees(std::atan2(sinRoll, cosRoll));

    camera.yaw   = newYaw;
    camera.pitch = newPitch;
    camera.roll  = newRoll;

    const char* axisNames[] = {"X", "Y", "Z"};
    std::cout << "[Operation] Camera orbited " << angleDeg << " deg around "
              << axisNames[axis] << " (visual crystal rotation)" << std::endl;
}