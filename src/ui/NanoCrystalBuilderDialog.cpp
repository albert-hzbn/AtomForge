#include "ui/NanoCrystalBuilderDialog.h"

#include "ElementData.h"
#include "graphics/CylinderMesh.h"
#include "graphics/Renderer.h"
#include "graphics/SphereMesh.h"
#include "graphics/StructureInstanceBuilder.h"
#include "graphics/ShadowMap.h"
#include "io/StructureLoader.h"
#include "app/SceneView.h"
#include "camera/Camera.h"
#include "imgui.h"

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef _WIN32
#include <dirent.h>
#include <sys/stat.h>
#else
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace
{

// ---------------------------------------------------------------------------
// Path / home-directory utilities
// ---------------------------------------------------------------------------

std::string nanoDetectHome()
{
#ifdef _WIN32
    if (const char* up = std::getenv("USERPROFILE"))
        return up;
    const char* hd = std::getenv("HOMEDRIVE");
    const char* hp = std::getenv("HOMEPATH");
    if (hd && hp) { return std::string(hd) + hp; }
    return "C:\\";
#else
    if (const char* home = std::getenv("HOME"))
        return home;
    return "/";
#endif
}

std::string nanoJoin(const std::string& base, const std::string& name)
{
    if (base.empty() || base == ".")
        return name;
    char last = base.back();
    if (last == '/' || last == '\\')
        return base + name;
    return base + "/" + name;
}

bool nanoIsSupportedFile(const std::string& name)
{
    static const char* exts[] = {
        ".cif",".mol",".pdb",".xyz",".sdf",".vasp",".mol2",".pwi",".gjf",nullptr
    };
    std::string low = name;
    for (char& c : low) c = (char)std::tolower((unsigned char)c);
    for (int i = 0; exts[i]; ++i) {
        const std::string ext(exts[i]);
        if (low.size() >= ext.size() &&
            low.compare(low.size() - ext.size(), ext.size(), ext) == 0)
            return true;
    }
    return false;
}

void nanoLoadDirEntries(const std::string& dir,
                        std::vector<std::pair<std::string, bool>>& entries)
{
    entries.clear();
#ifdef _WIN32
    std::string pat = dir;
    if (!pat.empty() && pat.back() != '/' && pat.back() != '\\')
        pat.push_back('\\');
    pat.push_back('*');
    WIN32_FIND_DATAA fd;
    HANDLE h = FindFirstFileA(pat.c_str(), &fd);
    if (h == INVALID_HANDLE_VALUE) return;
    do {
        std::string nm(fd.cFileName);
        if (nm == "." || nm == "..") continue;
        bool isDir = (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;
        if (!isDir && !nanoIsSupportedFile(nm)) continue;
        entries.emplace_back(nm, isDir);
    } while (FindNextFileA(h, &fd));
    FindClose(h);
#else
    DIR* d = opendir(dir.empty() ? "." : dir.c_str());
    if (!d) return;
    struct dirent* de;
    while ((de = readdir(d)) != nullptr) {
        std::string nm(de->d_name);
        if (nm == "." || nm == "..") continue;
        std::string full = nanoJoin(dir, nm);
        struct stat st;
        bool isDir = (stat(full.c_str(), &st) == 0 && S_ISDIR(st.st_mode));
        if (!isDir && !nanoIsSupportedFile(nm)) continue;
        entries.emplace_back(nm, isDir);
    }
    closedir(d);
#endif
    std::sort(entries.begin(), entries.end(),
              [](const std::pair<std::string,bool>& a,
                 const std::pair<std::string,bool>& b) {
        if (a.second != b.second) return a.second > b.second;
        return a.first < b.first;
    });
}

// ---------------------------------------------------------------------------
// Shape definitions
// ---------------------------------------------------------------------------

enum class NanoShape
{
    Sphere = 0,
    Ellipsoid,
    Box,
    Cylinder,
    Octahedron,
    TruncatedOctahedron,
    Cuboctahedron,
};

constexpr int kNumShapes = 7;

const char* shapeLabel(NanoShape s)
{
    switch (s) {
        case NanoShape::Sphere:              return "Sphere";
        case NanoShape::Ellipsoid:           return "Ellipsoid";
        case NanoShape::Box:                 return "Box";
        case NanoShape::Cylinder:            return "Cylinder";
        case NanoShape::Octahedron:          return "Octahedron";
        case NanoShape::TruncatedOctahedron: return "Truncated Octahedron";
        case NanoShape::Cuboctahedron:       return "Cuboctahedron";
    }
    return "Unknown";
}

struct NanoParams
{
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

    bool  autoCenterFromAtoms = true;
    float cx = 0.0f;
    float cy = 0.0f;
    float cz = 0.0f;

    bool autoReplicate = true;
    int  repA = 5;
    int  repB = 5;
    int  repC = 5;

    bool  setOutputCell = true;
    float vacuumPadding = 5.0f;
};

struct NanoBuildResult
{
    bool        success           = false;
    std::string message;
    int         inputAtoms        = 0;
    int         outputAtoms       = 0;
    NanoShape   shape             = NanoShape::Sphere;
    float       estimatedDiameter = 0.0f;
    int         repA = 0, repB = 0, repC = 0;
    bool        tilingUsed  = false;
    bool        repClamped  = false;
};

struct HalfExtents { float hx, hy, hz; };

float computeBoundingRadius(const NanoParams& p)
{
    switch (p.shape) {
        case NanoShape::Sphere:
            return p.sphereRadius;
        case NanoShape::Ellipsoid:
            return std::max({p.ellipRx, p.ellipRy, p.ellipRz});
        case NanoShape::Box:
            return std::sqrt(p.boxHx*p.boxHx + p.boxHy*p.boxHy + p.boxHz*p.boxHz);
        case NanoShape::Cylinder:
            return std::sqrt(p.cylRadius*p.cylRadius
                             + (p.cylHeight*0.5f)*(p.cylHeight*0.5f));
        case NanoShape::Octahedron:
            return p.octRadius;
        case NanoShape::TruncatedOctahedron:
            return std::max(p.truncOctRadius, p.truncOctTrunc * std::sqrt(3.0f));
        case NanoShape::Cuboctahedron:
            return p.cuboRadius;
    }
    return 10.0f;
}

HalfExtents computeShapeHalfExtents(const NanoParams& p)
{
    switch (p.shape) {
        case NanoShape::Sphere:
            return {p.sphereRadius, p.sphereRadius, p.sphereRadius};
        case NanoShape::Ellipsoid:
            return {p.ellipRx, p.ellipRy, p.ellipRz};
        case NanoShape::Box:
            return {p.boxHx, p.boxHy, p.boxHz};
        case NanoShape::Cylinder:
            if (p.cylAxis == 0) return {p.cylHeight*0.5f, p.cylRadius, p.cylRadius};
            if (p.cylAxis == 1) return {p.cylRadius, p.cylHeight*0.5f, p.cylRadius};
            return {p.cylRadius, p.cylRadius, p.cylHeight*0.5f};
        case NanoShape::Octahedron:
            return {p.octRadius, p.octRadius, p.octRadius};
        case NanoShape::TruncatedOctahedron: {
            float r = std::min(p.truncOctTrunc, p.truncOctRadius);
            return {r, r, r};
        }
        case NanoShape::Cuboctahedron:
            return {p.cuboRadius, p.cuboRadius, p.cuboRadius};
    }
    return {10.0f, 10.0f, 10.0f};
}

bool isInsideShape(const glm::vec3& p, const NanoParams& params)
{
    switch (params.shape) {
        case NanoShape::Sphere:
            return glm::dot(p, p) <= params.sphereRadius * params.sphereRadius;
        case NanoShape::Ellipsoid: {
            float fx = p.x / params.ellipRx;
            float fy = p.y / params.ellipRy;
            float fz = p.z / params.ellipRz;
            return fx*fx + fy*fy + fz*fz <= 1.0f;
        }
        case NanoShape::Box:
            return std::abs(p.x) <= params.boxHx &&
                   std::abs(p.y) <= params.boxHy &&
                   std::abs(p.z) <= params.boxHz;
        case NanoShape::Cylinder: {
            float r2, ax;
            if      (params.cylAxis == 0) { r2 = p.y*p.y + p.z*p.z; ax = p.x; }
            else if (params.cylAxis == 1) { r2 = p.x*p.x + p.z*p.z; ax = p.y; }
            else                          { r2 = p.x*p.x + p.y*p.y; ax = p.z; }
            return r2 <= params.cylRadius * params.cylRadius &&
                   std::abs(ax) <= params.cylHeight * 0.5f;
        }
        case NanoShape::Octahedron:
            return std::abs(p.x) + std::abs(p.y) + std::abs(p.z) <= params.octRadius;
        case NanoShape::TruncatedOctahedron:
            return (std::abs(p.x) + std::abs(p.y) + std::abs(p.z) <= params.truncOctRadius) &&
                   std::abs(p.x) <= params.truncOctTrunc &&
                   std::abs(p.y) <= params.truncOctTrunc &&
                   std::abs(p.z) <= params.truncOctTrunc;
        case NanoShape::Cuboctahedron:
            return std::abs(p.x) + std::abs(p.y) <= params.cuboRadius &&
                   std::abs(p.y) + std::abs(p.z) <= params.cuboRadius &&
                   std::abs(p.x) + std::abs(p.z) <= params.cuboRadius;
    }
    return false;
}

glm::vec3 computeAtomCentroid(const std::vector<AtomSite>& atoms)
{
    if (atoms.empty()) return glm::vec3(0.0f);
    glm::vec3 sum(0.0f);
    for (const AtomSite& a : atoms)
        sum += glm::vec3((float)a.x, (float)a.y, (float)a.z);
    return sum / (float)atoms.size();
}

float safeLen3(const glm::vec3& v)
{
    return std::sqrt(glm::dot(v, v));
}

// ---------------------------------------------------------------------------
// Nanocrystal builder core
// ---------------------------------------------------------------------------

NanoBuildResult buildNanocrystal(Structure& structure,
                                  const Structure& reference,
                                  const NanoParams& params,
                                  const std::vector<glm::vec3>& elementColors)
{
    NanoBuildResult result;
    result.shape      = params.shape;
    result.inputAtoms = (int)reference.atoms.size();

    if (reference.atoms.empty()) {
        result.message = "Reference structure has no atoms.";
        return result;
    }

    glm::vec3 center;
    if (params.autoCenterFromAtoms)
        center = computeAtomCentroid(reference.atoms);
    else
        center = glm::vec3(params.cx, params.cy, params.cz);

    const float maxR = computeBoundingRadius(params);
    std::vector<AtomSite> generatedAtoms;

    if (reference.hasUnitCell) {
        const auto& cv = reference.cellVectors;
        const glm::vec3 a((float)cv[0][0], (float)cv[0][1], (float)cv[0][2]);
        const glm::vec3 b((float)cv[1][0], (float)cv[1][1], (float)cv[1][2]);
        const glm::vec3 c((float)cv[2][0], (float)cv[2][1], (float)cv[2][2]);
        const float la = safeLen3(a), lb = safeLen3(b), lc = safeLen3(c);

        if (la < 1e-8f || lb < 1e-8f || lc < 1e-8f) {
            result.message = "Reference structure has degenerate lattice vectors.";
            return result;
        }

        const int kMaxReps = 40;
        bool clamped = false;
        auto safeRep = [&](float L) -> int {
            int n = (int)std::ceil(maxR / L) + 2;
            if (n > kMaxReps) { clamped = true; n = kMaxReps; }
            return n;
        };

        int nA, nB, nC;
        if (params.autoReplicate) {
            nA = safeRep(la); nB = safeRep(lb); nC = safeRep(lc);
        } else {
            nA = std::max(1, params.repA);
            nB = std::max(1, params.repB);
            nC = std::max(1, params.repC);
            clamped = false;
        }

        result.repA = nA; result.repB = nB; result.repC = nC;
        result.tilingUsed  = true;
        result.repClamped  = clamped;

        const long long total =
            (long long)(2*nA+1) * (long long)(2*nB+1) * (long long)(2*nC+1)
            * (long long)reference.atoms.size();

        if (total > 8000000LL) {
            std::ostringstream msg;
            msg << "Tiling would test " << total
                << " atoms (limit 8M). Reduce shape size or use manual replication.";
            result.message = msg.str();
            return result;
        }

        generatedAtoms.reserve((size_t)std::min(total, (long long)2000000));

        for (int ia = -nA; ia <= nA; ++ia)
        for (int ib = -nB; ib <= nB; ++ib)
        for (int ic = -nC; ic <= nC; ++ic) {
            const glm::vec3 offset = (float)ia*a + (float)ib*b + (float)ic*c;
            for (const AtomSite& atom : reference.atoms) {
                const glm::vec3 pos(
                    (float)atom.x + offset.x,
                    (float)atom.y + offset.y,
                    (float)atom.z + offset.z);
                if (!isInsideShape(pos - center, params)) continue;

                AtomSite out = atom;
                out.x = (double)pos.x;
                out.y = (double)pos.y;
                out.z = (double)pos.z;
                int z = out.atomicNumber;
                if (z >= 0 && z < (int)elementColors.size()) {
                    out.r = elementColors[z].r;
                    out.g = elementColors[z].g;
                    out.b = elementColors[z].b;
                } else {
                    getDefaultElementColor(z, out.r, out.g, out.b);
                }
                generatedAtoms.push_back(out);
            }
        }
    } else {
        result.tilingUsed = false;
        generatedAtoms.reserve(reference.atoms.size());
        for (const AtomSite& atom : reference.atoms) {
            const glm::vec3 pos((float)atom.x, (float)atom.y, (float)atom.z);
            if (!isInsideShape(pos - center, params)) continue;
            AtomSite out = atom;
            int z = out.atomicNumber;
            if (z >= 0 && z < (int)elementColors.size()) {
                out.r = elementColors[z].r;
                out.g = elementColors[z].g;
                out.b = elementColors[z].b;
            } else {
                getDefaultElementColor(z, out.r, out.g, out.b);
            }
            generatedAtoms.push_back(out);
        }
    }

    if (generatedAtoms.empty()) {
        result.message =
            "No atoms within the specified shape. "
            "Try increasing the size parameter(s).";
        return result;
    }

    result.estimatedDiameter = 2.0f * maxR;
    result.outputAtoms = (int)generatedAtoms.size();
    structure.atoms.swap(generatedAtoms);

    if (params.setOutputCell) {
        const HalfExtents he = computeShapeHalfExtents(params);
        const float pad = params.vacuumPadding;
        structure.hasUnitCell = true;
        structure.cellOffset  = {{
            (double)(center.x - he.hx - pad),
            (double)(center.y - he.hy - pad),
            (double)(center.z - he.hz - pad) }};
        structure.cellVectors = {{
            {{ 2.0*(he.hx+pad), 0.0, 0.0 }},
            {{ 0.0, 2.0*(he.hy+pad), 0.0 }},
            {{ 0.0, 0.0, 2.0*(he.hz+pad) }} }};
    } else {
        structure.hasUnitCell = false;
    }

    result.success = true;
    result.message = "Nanocrystal built successfully.";
    return result;
}

// ---------------------------------------------------------------------------
// Shape parameter UI helpers
// ---------------------------------------------------------------------------

void drawShapeParameters(NanoParams& params)
{
    switch (params.shape) {
        case NanoShape::Sphere:
            ImGui::SetNextItemWidth(160.0f);
            ImGui::InputFloat("Radius (A)##sph", &params.sphereRadius, 0.0f, 0.0f, "%.2f");
            if (params.sphereRadius < 0.1f) params.sphereRadius = 0.1f;
            break;

        case NanoShape::Ellipsoid:
            ImGui::SetNextItemWidth(160.0f);
            ImGui::InputFloat("Semi-axis X (A)##ell", &params.ellipRx, 0.0f, 0.0f, "%.2f");
            ImGui::SetNextItemWidth(160.0f);
            ImGui::InputFloat("Semi-axis Y (A)##ell", &params.ellipRy, 0.0f, 0.0f, "%.2f");
            ImGui::SetNextItemWidth(160.0f);
            ImGui::InputFloat("Semi-axis Z (A)##ell", &params.ellipRz, 0.0f, 0.0f, "%.2f");
            if (params.ellipRx < 0.1f) params.ellipRx = 0.1f;
            if (params.ellipRy < 0.1f) params.ellipRy = 0.1f;
            if (params.ellipRz < 0.1f) params.ellipRz = 0.1f;
            break;

        case NanoShape::Box:
            ImGui::SetNextItemWidth(160.0f);
            ImGui::InputFloat("Half-width X (A)##box", &params.boxHx, 0.0f, 0.0f, "%.2f");
            ImGui::SetNextItemWidth(160.0f);
            ImGui::InputFloat("Half-width Y (A)##box", &params.boxHy, 0.0f, 0.0f, "%.2f");
            ImGui::SetNextItemWidth(160.0f);
            ImGui::InputFloat("Half-width Z (A)##box", &params.boxHz, 0.0f, 0.0f, "%.2f");
            ImGui::TextDisabled("Full: %.2fx%.2fx%.2f A",
                                2.f*params.boxHx, 2.f*params.boxHy, 2.f*params.boxHz);
            if (params.boxHx < 0.1f) params.boxHx = 0.1f;
            if (params.boxHy < 0.1f) params.boxHy = 0.1f;
            if (params.boxHz < 0.1f) params.boxHz = 0.1f;
            break;

        case NanoShape::Cylinder: {
            const char* axisNames[] = {"X", "Y", "Z"};
            ImGui::SetNextItemWidth(160.0f);
            ImGui::InputFloat("Radius (A)##cyl", &params.cylRadius, 0.0f, 0.0f, "%.2f");
            ImGui::SetNextItemWidth(160.0f);
            ImGui::InputFloat("Height (A)##cyl", &params.cylHeight, 0.0f, 0.0f, "%.2f");
            ImGui::SetNextItemWidth(120.0f);
            ImGui::Combo("Cylinder axis##cyl", &params.cylAxis, axisNames, 3);
            if (params.cylRadius < 0.1f) params.cylRadius = 0.1f;
            if (params.cylHeight < 0.1f) params.cylHeight = 0.1f;
            break;
        }

        case NanoShape::Octahedron:
            ImGui::SetNextItemWidth(160.0f);
            ImGui::InputFloat("Radius (A)##oct", &params.octRadius, 0.0f, 0.0f, "%.2f");
            ImGui::TextDisabled("|x|+|y|+|z| <= R");
            if (params.octRadius < 0.1f) params.octRadius = 0.1f;
            break;

        case NanoShape::TruncatedOctahedron:
            ImGui::SetNextItemWidth(160.0f);
            ImGui::InputFloat("Octahedron R (A)##trunc", &params.truncOctRadius, 0.0f, 0.0f, "%.2f");
            ImGui::SetNextItemWidth(160.0f);
            ImGui::InputFloat("Truncation R (A)##trunc", &params.truncOctTrunc,  0.0f, 0.0f, "%.2f");
            ImGui::TextDisabled("|x|+|y|+|z|<=R_oct  AND  max(|x|,|y|,|z|)<=R_trunc");
            if (params.truncOctRadius < 0.1f) params.truncOctRadius = 0.1f;
            if (params.truncOctTrunc  < 0.1f) params.truncOctTrunc  = 0.1f;
            break;

        case NanoShape::Cuboctahedron:
            ImGui::SetNextItemWidth(160.0f);
            ImGui::InputFloat("Radius (A)##cubo", &params.cuboRadius, 0.0f, 0.0f, "%.2f");
            ImGui::TextDisabled("|x|+|y|, |y|+|z|, |x|+|z|  <= R");
            if (params.cuboRadius < 0.1f) params.cuboRadius = 0.1f;
            break;
    }
}

void drawReplicationPreview(const NanoParams& params, const Structure& reference)
{
    if (!reference.hasUnitCell) return;
    const auto& cv = reference.cellVectors;
    const glm::vec3 a((float)cv[0][0],(float)cv[0][1],(float)cv[0][2]);
    const glm::vec3 b((float)cv[1][0],(float)cv[1][1],(float)cv[1][2]);
    const glm::vec3 c((float)cv[2][0],(float)cv[2][1],(float)cv[2][2]);
    const float la = safeLen3(a), lb = safeLen3(b), lc = safeLen3(c);
    if (la < 1e-8f || lb < 1e-8f || lc < 1e-8f) return;

    const float maxR = computeBoundingRadius(params);
    const int kMax = 40;
    bool clamped = false;
    auto rep = [&](float L) -> int {
        int n = (int)std::ceil(maxR / L) + 2;
        if (n > kMax) { clamped = true; n = kMax; }
        return n;
    };
    int nA, nB, nC;
    if (params.autoReplicate) {
        nA = rep(la); nB = rep(lb); nC = rep(lc);
    } else {
        nA = std::max(1, params.repA);
        nB = std::max(1, params.repB);
        nC = std::max(1, params.repC);
    }
    const long long tested =
        (long long)(2*nA+1)*(long long)(2*nB+1)*(long long)(2*nC+1)
        * (long long)reference.atoms.size();
    ImGui::TextDisabled("Tiling: %dx%dx%d  (%lld atoms tested)",
                        2*nA+1, 2*nB+1, 2*nC+1, (long long)tested);
    if (clamped)
        ImGui::TextColored(ImVec4(1,0.7f,0,1), "Warning: replication clamped to %d.", kMax);
    if (tested > 8000000LL)
        ImGui::TextColored(ImVec4(1,0.3f,0.3f,1),
                           "Exceeds 8M limit -- build will be refused.");
}

void drawNanoBuildResult(const NanoBuildResult& result)
{
    if (result.message.empty()) return;
    ImGui::TextWrapped("Status: %s", result.message.c_str());
    if (!result.success) return;
    ImGui::Separator();
    ImGui::Text("Shape:            %s", shapeLabel(result.shape));
    ImGui::Text("Reference atoms:  %d", result.inputAtoms);
    ImGui::Text("Output atoms:     %d", result.outputAtoms);
    ImGui::Text("Est. diameter:    %.2f A", result.estimatedDiameter);
    if (result.tilingUsed)
        ImGui::Text("Tiling: %dx%dx%d (%s)",
                    2*result.repA+1, 2*result.repB+1, 2*result.repC+1,
                    result.repClamped ? "clamped" : "auto");
    else
        ImGui::Text("Tiling: none (no unit cell in reference)");
}

} // anonymous namespace

// ===========================================================================
// NanoCrystalBuilderDialog – public / private implementation
// ===========================================================================

NanoCrystalBuilderDialog::NanoCrystalBuilderDialog()
    : m_browsEntryDirty(true)
    , m_browsShowPanel(false)
    , m_glReady(false)
    , m_previewBufDirty(true)
    , m_camYaw(45.0f)
    , m_camPitch(35.0f)
    , m_camDistance(10.0f)
{
    m_browsFilename[0] = '\0';
    m_browsStatusMsg[0] = '\0';
    m_browsDir = nanoDetectHome();
}

NanoCrystalBuilderDialog::~NanoCrystalBuilderDialog()
{
    if (m_previewFBO)      { glDeleteFramebuffers(1,  &m_previewFBO);    m_previewFBO = 0; }
    if (m_previewColorTex) { glDeleteTextures(1,      &m_previewColorTex); m_previewColorTex = 0; }
    if (m_previewDepthRbo) { glDeleteRenderbuffers(1, &m_previewDepthRbo); m_previewDepthRbo = 0; }

    if (m_previewShadow.depthFBO)
        glDeleteFramebuffers(1, &m_previewShadow.depthFBO);
    if (m_previewShadow.depthTexture)
        glDeleteTextures(1, &m_previewShadow.depthTexture);

    delete m_previewSphere;
    delete m_previewCylinder;
}

void NanoCrystalBuilderDialog::initRenderResources(Renderer& renderer)
{
    m_renderer = &renderer;

    // Independent mesh objects so VAO state does not conflict with the main
    // SceneBuffers.
    m_previewSphere   = new SphereMesh(24, 24);
    m_previewCylinder = new CylinderMesh(16);
    m_previewBuffers.init(m_previewSphere->vao, m_previewCylinder->vao);

    // A 1×1 shadow map satisfies the sampler binding without drawing shadows.
    m_previewShadow = createShadowMap(1, 1);

    m_glReady = true;
}

void NanoCrystalBuilderDialog::feedDroppedFile(const std::string& path)
{
    m_pendingDropPath = path;
}

void NanoCrystalBuilderDialog::refreshBrowserEntries()
{
    nanoLoadDirEntries(m_browsDir, m_browsEntries);
    m_browsEntryDirty = false;
}

bool NanoCrystalBuilderDialog::tryLoadFile(const std::string& path,
                                            const std::vector<float>& radii,
                                            const std::vector<float>& shininess)
{
    Structure loaded;
    std::string err;
    if (!loadStructureFromFile(path, loaded, err)) {
        std::snprintf(m_browsStatusMsg, sizeof(m_browsStatusMsg),
                      "Load failed: %s", err.c_str());
        return false;
    }
    m_reference = std::move(loaded);
    std::snprintf(m_browsStatusMsg, sizeof(m_browsStatusMsg),
                  "Loaded: %d atoms", (int)m_reference.atoms.size());
    m_previewBufDirty = true;
    if (m_glReady) {
        rebuildPreviewBuffers(radii, shininess);
        autoFitPreviewCamera();
    }
    std::cout << "[NanoCrystalBuilder] Reference loaded: " << path
              << " (" << m_reference.atoms.size() << " atoms)" << std::endl;
    return true;
}

void NanoCrystalBuilderDialog::ensurePreviewFBO(int w, int h)
{
    if (w == m_previewW && h == m_previewH && m_previewFBO != 0)
        return;

    if (m_previewFBO)      { glDeleteFramebuffers(1,  &m_previewFBO);    m_previewFBO = 0; }
    if (m_previewColorTex) { glDeleteTextures(1,      &m_previewColorTex); m_previewColorTex = 0; }
    if (m_previewDepthRbo) { glDeleteRenderbuffers(1, &m_previewDepthRbo); m_previewDepthRbo = 0; }

    glGenTextures(1, &m_previewColorTex);
    glBindTexture(GL_TEXTURE_2D, m_previewColorTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenRenderbuffers(1, &m_previewDepthRbo);
    glBindRenderbuffer(GL_RENDERBUFFER, m_previewDepthRbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    glGenFramebuffers(1, &m_previewFBO);
    glBindFramebuffer(GL_FRAMEBUFFER, m_previewFBO);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                           GL_TEXTURE_2D, m_previewColorTex, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER, m_previewDepthRbo);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    m_previewW = w;
    m_previewH = h;
}

void NanoCrystalBuilderDialog::rebuildPreviewBuffers(const std::vector<float>& radii,
                                                      const std::vector<float>& shininess)
{
    if (!m_glReady || m_reference.atoms.empty())
        return;

    static const int kIdent[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    StructureInstanceData data = buildStructureInstanceData(
        m_reference, false, kIdent, radii, shininess);

    std::array<bool,119> noFilter = {};
    m_previewBuffers.upload(data, false, noFilter);
    m_previewBufDirty = false;
}

void NanoCrystalBuilderDialog::autoFitPreviewCamera()
{
    m_camYaw   = 45.0f;
    m_camPitch = 35.0f;

    if (m_previewBuffers.atomCount == 0) {
        m_camDistance = 10.0f;
        return;
    }

    float maxR = 0.0f;
    for (size_t i = 0; i < m_previewBuffers.atomPositions.size(); ++i) {
        float r = (i < m_previewBuffers.atomRadii.size())
                  ? m_previewBuffers.atomRadii[i] : 0.0f;
        float d = glm::length(m_previewBuffers.atomPositions[i]
                              - m_previewBuffers.orbitCenter) + r;
        maxR = std::max(maxR, d);
    }
    maxR = std::max(maxR, 1.0f);

    const float halfFov = glm::radians(22.5f); // half of 45 deg fov
    float dist = maxR / std::sin(halfFov) * 1.15f;
    dist = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, dist));
    m_camDistance = dist;
}

void NanoCrystalBuilderDialog::renderPreviewToFBO(int w, int h)
{
    if (!m_glReady || !m_renderer || m_previewBuffers.atomCount == 0)
        return;

    ensurePreviewFBO(w, h);

    // Save current GL state
    GLint prevFbo = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
    GLint prevVP[4];
    glGetIntegerv(GL_VIEWPORT, prevVP);

    // Build frame matrices using a temporary Camera
    Camera cam;
    cam.yaw      = m_camYaw;
    cam.pitch    = m_camPitch;
    cam.distance = m_camDistance;

    FrameView frame;
    frame.framebufferWidth  = w;
    frame.framebufferHeight = h;
    buildFrameView(cam, m_previewBuffers, false, frame);

    // Draw into preview FBO
    glBindFramebuffer(GL_FRAMEBUFFER, m_previewFBO);
    glViewport(0, 0, w, h);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.09f, 0.11f, 0.15f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    m_renderer->drawBonds(frame.projection, frame.view,
                          frame.lightPosition, frame.cameraPosition,
                          *m_previewCylinder, m_previewBuffers.bondCount);

    m_renderer->drawAtoms(frame.projection, frame.view,
                          frame.lightMVP, frame.lightPosition, frame.cameraPosition,
                          m_previewShadow, *m_previewSphere,
                          m_previewBuffers.atomCount);

    m_renderer->drawBoxLines(frame.projection, frame.view,
                             m_previewBuffers.lineVAO,
                             m_previewBuffers.boxLines.size());

    // Restore GL state
    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
    glViewport(prevVP[0], prevVP[1], prevVP[2], prevVP[3]);
}

// ---------------------------------------------------------------------------
// drawMenuItem
// ---------------------------------------------------------------------------

void NanoCrystalBuilderDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Nanocrystal...", NULL, false, enabled))
        m_openRequested = true;
}

// ---------------------------------------------------------------------------
// drawDialog  – the main ImGui popup
// ---------------------------------------------------------------------------

void NanoCrystalBuilderDialog::drawDialog(
    Structure& structure,
    const std::vector<glm::vec3>& elementColors,
    const std::vector<float>& elementRadii,
    const std::vector<float>& elementShininess,
    const std::function<void(Structure&)>& updateBuffers)
{
    static NanoParams      params;
    static NanoBuildResult lastResult;

    // Consume any pending file drop
    if (!m_pendingDropPath.empty()) {
        tryLoadFile(m_pendingDropPath, elementRadii, elementShininess);
        m_pendingDropPath.clear();
    }

    if (m_openRequested) {
        ImGui::OpenPopup("Build Nanocrystal");
        lastResult      = {};
        m_openRequested = false;
    }

    m_isOpen = ImGui::IsPopupOpen("Build Nanocrystal");

    ImGui::SetNextWindowSize(ImVec2(960.0f, 820.0f), ImGuiCond_FirstUseEver);
    bool dialogOpen = true;
    if (!ImGui::BeginPopupModal("Build Nanocrystal", &dialogOpen,
                                ImGuiWindowFlags_NoResize)) {
        m_isOpen = false;
        return;
    }
    m_isOpen = true;

    ImGui::TextWrapped(
        "Load a reference crystal structure via drag and drop, then choose a shape to carve from it. "
        "The reference unit cell is tiled automatically to fill the requested shape.");
    ImGui::Separator();

    // =========================================================================
    // Layout: Left = Structure view | Right = Builder options
    // =========================================================================

    constexpr float kPrevW   = 420.0f;
    constexpr float kPrevH   = 350.0f;
    constexpr float kSideH   = 400.0f;

    // ---- LEFT: Structure view with drag-and-drop ----
    ImGui::BeginChild("##nanoStructureView", ImVec2(kPrevW, kSideH), true);

    ImGui::Text("Reference Structure");
    ImGui::Separator();
    ImGui::Text("Status: %s", m_reference.atoms.empty() ? "(none)" : "loaded");

    ImGui::Spacing();

    // Drag-and-drop target using InvisibleButton
    ImGui::InvisibleButton("##nanoRefDropZone", ImVec2(-1.0f, kPrevH));

    // Draw drop zone border
    ImDrawList* dl = ImGui::GetWindowDrawList();
    ImVec2 dropMin = ImGui::GetItemRectMin();
    ImVec2 dropMax = ImGui::GetItemRectMax();
    const bool dropZoneHovered = ImGui::IsItemHovered();
    const bool dropZoneActive = ImGui::IsItemActive();
    dl->AddRect(dropMin, dropMax, ImGui::GetColorU32(ImGuiCol_Border), 2.0f);

    if (m_reference.atoms.empty()) {
        if (ImGui::BeginDragDropTarget()) {
            if (const ImGuiPayload* payload =
                    ImGui::AcceptDragDropPayload("DROP_FILE", ImGuiDragDropFlags_AcceptPeekOnly)) {
                const char* files = (const char*)payload->Data;
                tryLoadFile(files, elementRadii, elementShininess);
            }
            ImGui::EndDragDropTarget();
        }

        // Show drop hint
        ImVec2 dropMid = ImVec2((dropMin.x + dropMax.x) / 2.0f, (dropMin.y + dropMax.y) / 2.0f);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize(),
                    ImVec2(dropMid.x - 90.0f, dropMid.y - 20.0f),
                    ImGui::GetColorU32(ImGuiCol_TextDisabled),
                    "Drop a structure file here");
    } else {
        // Show 3D preview
        if (m_glReady) {
            if (m_previewBufDirty)
                rebuildPreviewBuffers(elementRadii, elementShininess);

            const float previewPadding = 5.0f;
            const ImVec2 previewSize(dropMax.x - dropMin.x - 2.0f * previewPadding,
                                     dropMax.y - dropMin.y - 2.0f * previewPadding);
            const int previewW = std::max(1, (int)previewSize.x);
            const int previewH = std::max(1, (int)previewSize.y);

            renderPreviewToFBO(previewW, previewH);

            const ImVec2 previewMin(dropMin.x + previewPadding, dropMin.y + previewPadding);
            const ImVec2 previewMax(previewMin.x + previewSize.x, previewMin.y + previewSize.y);

            dl->AddImage((ImTextureID)(intptr_t)m_previewColorTex,
                         previewMin,
                         previewMax,
                         ImVec2(0, 1), ImVec2(1, 0));

            if (ImGui::BeginDragDropTarget()) {
                if (const ImGuiPayload* payload =
                        ImGui::AcceptDragDropPayload("DROP_FILE", ImGuiDragDropFlags_AcceptPeekOnly)) {
                    const char* files = (const char*)payload->Data;
                    tryLoadFile(files, elementRadii, elementShininess);
                }
                ImGui::EndDragDropTarget();
            }

            if (dropZoneActive && ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f)) {
                ImVec2 delta = ImGui::GetIO().MouseDelta;
                m_camYaw   -= delta.x * 0.5f;
                m_camPitch += delta.y * 0.5f;
                if (m_camPitch >  89.0f) m_camPitch =  89.0f;
                if (m_camPitch < -89.0f) m_camPitch = -89.0f;
            }

            if (dropZoneHovered) {
                float wheel = ImGui::GetIO().MouseWheel;
                if (wheel != 0.0f) {
                    m_camDistance -= wheel * m_camDistance * 0.1f;
                    m_camDistance  = std::max(Camera::kMinDistance,
                                              std::min(Camera::kMaxDistance, m_camDistance));
                }
            }

            ImGui::TextDisabled("Left-drag = orbit   Scroll = zoom");
        }
    }

    // Reference info bar
    if (!m_reference.atoms.empty()) {
        ImGui::Spacing();
        ImGui::Separator();
        if (m_reference.hasUnitCell) {
            const auto& cv = m_reference.cellVectors;
            float la = safeLen3(glm::vec3((float)cv[0][0],(float)cv[0][1],(float)cv[0][2]));
            float lb = safeLen3(glm::vec3((float)cv[1][0],(float)cv[1][1],(float)cv[1][2]));
            float lc = safeLen3(glm::vec3((float)cv[2][0],(float)cv[2][1],(float)cv[2][2]));
            ImGui::TextDisabled("%d atoms", (int)m_reference.atoms.size());
            ImGui::TextDisabled("a=%.3f  b=%.3f  c=%.3f A", la, lb, lc);
            ImGui::TextDisabled("Unit cell: tiling enabled");
        } else {
            ImGui::TextDisabled("%d atoms (shape carving only)",
                        (int)m_reference.atoms.size());
        }
    }

    ImGui::EndChild(); // ##nanoStructureView

    ImGui::SameLine();

    // ---- RIGHT: Builder options ----
    ImGui::BeginChild("##nanoBuilderOptions", ImVec2(0, kSideH), true);

    ImGui::Text("Builder Options");
    ImGui::Separator();

    ImGui::Text("Shape");
    const char* shapeLabels[kNumShapes] = {
        "Sphere","Ellipsoid","Box","Cylinder",
        "Octahedron","Truncated Octahedron","Cuboctahedron"
    };
    int shapeInt = (int)params.shape;
    ImGui::SetNextItemWidth(-1.0f);
    if (ImGui::Combo("##nanoShapeCombo", &shapeInt, shapeLabels, kNumShapes))
        params.shape = (NanoShape)shapeInt;

    ImGui::Spacing();
    ImGui::Separator();

    // Use a child for scrollable parameters
    ImGui::BeginChild("##nanoParamsScroll", ImVec2(-1, -70), true, ImGuiWindowFlags_HorizontalScrollbar);
    
    drawShapeParameters(params);

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Text("Center of carving");
    ImGui::Checkbox("Auto-center##nano",
                    &params.autoCenterFromAtoms);
    if (!params.autoCenterFromAtoms) {
        ImGui::SetNextItemWidth(100.0f);
        ImGui::InputFloat("CX##cnano", &params.cx, 0.f, 0.f, "%.3f");
        ImGui::SetNextItemWidth(100.0f);
        ImGui::InputFloat("CY##cnano", &params.cy, 0.f, 0.f, "%.3f");
        ImGui::SetNextItemWidth(100.0f);
        ImGui::InputFloat("CZ##cnano", &params.cz, 0.f, 0.f, "%.3f");
    }

    if (m_reference.hasUnitCell) {
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("Supercell replication");
        ImGui::Checkbox("Auto-replicate##nano", &params.autoReplicate);
        if (!params.autoReplicate) {
            ImGui::SetNextItemWidth(70.0f);
            ImGui::InputInt("Reps a##nano", &params.repA);
            ImGui::SetNextItemWidth(70.0f);
            ImGui::InputInt("Reps b##nano", &params.repB);
            ImGui::SetNextItemWidth(70.0f);
            ImGui::InputInt("Reps c##nano", &params.repC);
            if (params.repA < 1) params.repA = 1;
            if (params.repB < 1) params.repB = 1;
            if (params.repC < 1) params.repC = 1;
        }
        drawReplicationPreview(params, m_reference);
    }

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Text("Output cell");
    ImGui::Checkbox("Rectangular cell##nano", &params.setOutputCell);
    if (params.setOutputCell) {
        ImGui::SetNextItemWidth(100.0f);
        ImGui::InputFloat("Vacuum (A)##nano",
                          &params.vacuumPadding, 0.f, 0.f, "%.2f");
        if (params.vacuumPadding < 0.f) params.vacuumPadding = 0.f;
        const HalfExtents he = computeShapeHalfExtents(params);
        ImGui::TextDisabled("%.2fx%.2fx%.2f A",
                            2.f*(he.hx+params.vacuumPadding),
                            2.f*(he.hy+params.vacuumPadding),
                            2.f*(he.hz+params.vacuumPadding));
    }

    ImGui::EndChild(); // ##nanoParamsScroll

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Text("Result");
    drawNanoBuildResult(lastResult);

    ImGui::EndChild(); // ##nanoBuilderOptions

    // =========================================================================
    // Action bar
    // =========================================================================

    ImGui::Separator();

    const bool canBuild = !m_reference.atoms.empty();
    if (!canBuild) ImGui::BeginDisabled();
    if (ImGui::Button("Build##nano", ImVec2(100.0f, 0.0f))) {
        lastResult = buildNanocrystal(structure, m_reference, params, elementColors);
        if (lastResult.success)
            updateBuffers(structure);
    }
    if (!canBuild) ImGui::EndDisabled();

    ImGui::SameLine();
    if (ImGui::Button("Close##nano", ImVec2(80.0f, 0.0f))) {
        lastResult  = {};
        m_isOpen    = false;
        ImGui::CloseCurrentPopup();
    }

    ImGui::EndPopup();
}
