#include "ui/NanoCrystalBuilderDialog.h"
#include "algorithms/NanoCrystalBuilder.h"
#include "util/PathUtils.h"

#include "ElementData.h"
#include "graphics/CylinderMesh.h"
#include "graphics/Renderer.h"
#include "graphics/Shader.h"
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
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <glm/gtc/type_ptr.hpp>

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
// File filter and directory loading (nano-specific)
// ---------------------------------------------------------------------------

bool nanoIsSupportedFile(const std::string& name)
{
    static const char* exts[] = {
        ".cif",".mol",".pdb",".xyz",".sdf",".vasp",".mol2",".pwi",".gjf",
        nullptr
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
        std::string full = joinPath(dir, nm);
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
// Local helper
// ---------------------------------------------------------------------------

static float safeLen3(const glm::vec3& v)
{
    return std::sqrt(glm::dot(v, v));
}

std::size_t hashCombine(std::size_t seed, std::size_t value)
{
    return seed ^ (value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2));
}

std::size_t structurePreviewHash(const Structure& structure)
{
    std::size_t seed = 0;
    seed = hashCombine(seed, std::hash<std::size_t>{}(structure.atoms.size()));
    seed = hashCombine(seed, std::hash<bool>{}(structure.hasUnitCell));
    if (structure.hasUnitCell)
    {
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                seed = hashCombine(seed, std::hash<double>{}(structure.cellVectors[row][col]));
    }
    const std::size_t sampleCount = std::min<std::size_t>(structure.atoms.size(), 32u);
    for (std::size_t index = 0; index < sampleCount; ++index)
    {
        const AtomSite& atom = structure.atoms[index];
        seed = hashCombine(seed, std::hash<int>{}(atom.atomicNumber));
        seed = hashCombine(seed, std::hash<double>{}(atom.x));
        seed = hashCombine(seed, std::hash<double>{}(atom.y));
        seed = hashCombine(seed, std::hash<double>{}(atom.z));
    }
    return seed;
}

std::size_t wulffParamsHash(const NanoParams& params)
{
    std::size_t seed = 0;
    seed = hashCombine(seed, std::hash<int>{}((int)params.generationMode));
    seed = hashCombine(seed, std::hash<float>{}(params.wulffMaxRadius));
    seed = hashCombine(seed, std::hash<bool>{}(params.autoCenterFromAtoms));
    seed = hashCombine(seed, std::hash<float>{}(params.cx));
    seed = hashCombine(seed, std::hash<float>{}(params.cy));
    seed = hashCombine(seed, std::hash<float>{}(params.cz));
    seed = hashCombine(seed, std::hash<bool>{}(params.applyCrystalOrientation));
    seed = hashCombine(seed, std::hash<bool>{}(params.useMillerOrientation));
    seed = hashCombine(seed, std::hash<float>{}(params.orientXDeg));
    seed = hashCombine(seed, std::hash<float>{}(params.orientYDeg));
    seed = hashCombine(seed, std::hash<float>{}(params.orientZDeg));
    seed = hashCombine(seed, std::hash<int>{}(params.millerH));
    seed = hashCombine(seed, std::hash<int>{}(params.millerK));
    seed = hashCombine(seed, std::hash<int>{}(params.millerL));
    seed = hashCombine(seed, std::hash<std::size_t>{}(params.wulffPlanes.size()));
    for (const WulffPlaneInput& plane : params.wulffPlanes)
    {
        seed = hashCombine(seed, std::hash<int>{}(plane.h));
        seed = hashCombine(seed, std::hash<int>{}(plane.k));
        seed = hashCombine(seed, std::hash<int>{}(plane.l));
        seed = hashCombine(seed, std::hash<float>{}(plane.surfaceEnergy));
    }
    return seed;
}

glm::vec3 wulffFamilyColor(int familyIndex)
{
    static const std::array<glm::vec3, 10> palette = {{
        glm::vec3(0.92f, 0.33f, 0.29f),
        glm::vec3(0.18f, 0.63f, 0.58f),
        glm::vec3(0.24f, 0.45f, 0.86f),
        glm::vec3(0.96f, 0.73f, 0.20f),
        glm::vec3(0.72f, 0.38f, 0.82f),
        glm::vec3(0.13f, 0.64f, 0.28f),
        glm::vec3(0.89f, 0.52f, 0.16f),
        glm::vec3(0.35f, 0.73f, 0.87f),
        glm::vec3(0.82f, 0.23f, 0.47f),
        glm::vec3(0.62f, 0.65f, 0.19f)
    }};
    return palette[(std::size_t)std::max(0, familyIndex) % palette.size()];
}

const char* generationModeLabel(NanoGenerationMode mode)
{
    switch (mode)
    {
        case NanoGenerationMode::Shape:
            return "Shape Cut";
        case NanoGenerationMode::WulffConstruction:
            return "Wulff Construction";
    }
    return "Unknown";
}

static const char* kWulffPlaneVS = R"(
    #version 130

    in vec3 position;

    uniform mat4 projection;
    uniform mat4 view;

    out vec3 fragWorldPos;

    void main()
    {
        fragWorldPos = position;
        gl_Position = projection * view * vec4(position, 1.0);
    }
)";

static const char* kWulffPlaneFS = R"(
    #version 130

    uniform vec3 faceColor;
    uniform vec3 faceNormal;
    uniform vec3 lightPos;
    uniform vec3 viewPos;

    in vec3 fragWorldPos;

    out vec4 color;

    void main()
    {
        vec3 norm = normalize(faceNormal);
        if (!gl_FrontFacing)
            norm = -norm;

        vec3 viewDir = normalize(viewPos - fragWorldPos);
        vec3 lightDir = normalize(lightPos - fragWorldPos);
        vec3 halfVec = normalize(lightDir + viewDir);

        float diffuse = max(dot(norm, lightDir), 0.0);
        diffuse = max(diffuse, 0.34);

        float specular = pow(max(dot(norm, halfVec), 0.0), 20.0);
        float fresnel = pow(1.0 - max(dot(norm, viewDir), 0.0), 2.6);

        vec3 reflected = reflect(-viewDir, norm);
        float skyMix = clamp(reflected.y * 0.5 + 0.5, 0.0, 1.0);
        vec3 envColor = mix(vec3(0.10, 0.11, 0.13), vec3(0.72, 0.78, 0.86), pow(skyMix, 1.35));

        vec3 baseLit = faceColor * (0.30 + 0.70 * diffuse);
        vec3 reflectedColor = mix(baseLit, envColor, 0.14 * fresnel);
        vec3 finalColor = reflectedColor + vec3(0.12 * specular + 0.08 * fresnel);
        color = vec4(finalColor, 1.0);
    }
)";

void drawWulffParameters(NanoParams& params)
{
    ImGui::SetNextItemWidth(180.0f);
    ImGui::DragFloat("Max radius (A)##wulffRadius", &params.wulffMaxRadius, 0.25f, 0.1f, 1000.0f, "%.2f");
    if (params.wulffMaxRadius < 0.1f)
        params.wulffMaxRadius = 0.1f;

    if (params.wulffPlanes.empty())
        params.wulffPlanes.push_back(WulffPlaneInput{});

    ImGui::Spacing();
    ImGui::Text("Facet families");
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextDisabled("Distances follow energy ratios; the highest-energy family sits at the max radius.");
    ImGui::PopTextWrapPos();

    for (std::size_t index = 0; index < params.wulffPlanes.size(); ++index)
    {
        WulffPlaneInput& plane = params.wulffPlanes[index];
        ImGui::PushID((int)index);
        ImGui::Separator();
        const glm::vec3 familyColor = wulffFamilyColor((int)index);
        ImGui::ColorButton("##familyColor", ImVec4(familyColor.r, familyColor.g, familyColor.b, 1.0f), ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoDragDrop, ImVec2(18.0f, 18.0f));
        ImGui::SameLine();
        ImGui::Text("Family %d", (int)index + 1);
        ImGui::SameLine();
        if (ImGui::SmallButton("Remove") && params.wulffPlanes.size() > 1)
        {
            params.wulffPlanes.erase(params.wulffPlanes.begin() + (long)index);
            ImGui::PopID();
            break;
        }

        ImGui::SetNextItemWidth(70.0f);
        ImGui::DragInt("h", &plane.h, 0.2f, -12, 12);
        ImGui::SameLine();
        ImGui::SetNextItemWidth(70.0f);
        ImGui::DragInt("k", &plane.k, 0.2f, -12, 12);
        ImGui::SameLine();
        ImGui::SetNextItemWidth(70.0f);
        ImGui::DragInt("l", &plane.l, 0.2f, -12, 12);
        ImGui::SameLine();
        ImGui::SetNextItemWidth(120.0f);
        ImGui::DragFloat("Energy", &plane.surfaceEnergy, 0.01f, 0.001f, 1000.0f, "%.4f");
        if (plane.surfaceEnergy < 0.001f)
            plane.surfaceEnergy = 0.001f;
        ImGui::PopID();
    }

    if (ImGui::Button("Add Facet Family##wulffAdd", ImVec2(150.0f, 0.0f)))
        params.wulffPlanes.push_back(WulffPlaneInput{});
}

// ---------------------------------------------------------------------------
// Shape parameter UI helpers
void drawShapeParameters(NanoParams& params)
{
    switch (params.shape) {
        case NanoShape::Sphere:
            ImGui::SetNextItemWidth(180.0f);
            ImGui::DragFloat("Radius (A)##sph", &params.sphereRadius, 0.5f, 0.1f, 500.0f, "%.2f");
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("Drag to adjust, Ctrl+click to type");
            if (params.sphereRadius < 0.1f) params.sphereRadius = 0.1f;
            break;

        case NanoShape::Ellipsoid:
            ImGui::SetNextItemWidth(180.0f);
            ImGui::DragFloat("Semi-axis a (A)##ell", &params.ellipRx, 0.5f, 0.1f, 500.0f, "%.2f");
            ImGui::SetNextItemWidth(180.0f);
            ImGui::DragFloat("Semi-axis b (A)##ell", &params.ellipRy, 0.5f, 0.1f, 500.0f, "%.2f");
            ImGui::SetNextItemWidth(180.0f);
            ImGui::DragFloat("Semi-axis c (A)##ell", &params.ellipRz, 0.5f, 0.1f, 500.0f, "%.2f");
            if (params.ellipRx < 0.1f) params.ellipRx = 0.1f;
            if (params.ellipRy < 0.1f) params.ellipRy = 0.1f;
            if (params.ellipRz < 0.1f) params.ellipRz = 0.1f;
            break;

        case NanoShape::Box:
            ImGui::SetNextItemWidth(180.0f);
            ImGui::DragFloat("Half-size X (A)##box", &params.boxHx, 0.5f, 0.1f, 500.0f, "%.2f");
            ImGui::SetNextItemWidth(180.0f);
            ImGui::DragFloat("Half-size Y (A)##box", &params.boxHy, 0.5f, 0.1f, 500.0f, "%.2f");
            ImGui::SetNextItemWidth(180.0f);
            ImGui::DragFloat("Half-size Z (A)##box", &params.boxHz, 0.5f, 0.1f, 500.0f, "%.2f");
            ImGui::TextDisabled("Full box: %.2f x %.2f x %.2f A",
                                2.f*params.boxHx, 2.f*params.boxHy, 2.f*params.boxHz);
            if (params.boxHx < 0.1f) params.boxHx = 0.1f;
            if (params.boxHy < 0.1f) params.boxHy = 0.1f;
            if (params.boxHz < 0.1f) params.boxHz = 0.1f;
            break;

        case NanoShape::Cylinder: {
            const char* axisNames[] = {"X", "Y", "Z"};
            ImGui::SetNextItemWidth(180.0f);
            ImGui::DragFloat("Radius (A)##cyl", &params.cylRadius, 0.5f, 0.1f, 500.0f, "%.2f");
            ImGui::SetNextItemWidth(180.0f);
            ImGui::DragFloat("Height (A)##cyl", &params.cylHeight, 0.5f, 0.1f, 1000.0f, "%.2f");
            ImGui::SetNextItemWidth(120.0f);
            ImGui::Combo("Axis##cyl", &params.cylAxis, axisNames, 3);
            if (params.cylRadius < 0.1f) params.cylRadius = 0.1f;
            if (params.cylHeight < 0.1f) params.cylHeight = 0.1f;
            break;
        }

        case NanoShape::Octahedron:
            ImGui::SetNextItemWidth(180.0f);
            ImGui::DragFloat("Radius (A)##oct", &params.octRadius, 0.5f, 0.1f, 500.0f, "%.2f");
            ImGui::TextDisabled("Condition: |x|+|y|+|z| <= R");
            if (params.octRadius < 0.1f) params.octRadius = 0.1f;
            break;

        case NanoShape::TruncatedOctahedron:
            ImGui::SetNextItemWidth(180.0f);
            ImGui::DragFloat("Octahedron R (A)##trunc", &params.truncOctRadius, 0.5f, 0.1f, 500.0f, "%.2f");
            ImGui::SetNextItemWidth(180.0f);
            ImGui::DragFloat("Truncation R (A)##trunc", &params.truncOctTrunc,  0.5f, 0.1f, 500.0f, "%.2f");
            ImGui::TextDisabled("Cond: |x|+|y|+|z|<=R_oct  AND  max(|x|,|y|,|z|)<=R_trunc");
            if (params.truncOctRadius < 0.1f) params.truncOctRadius = 0.1f;
            if (params.truncOctTrunc  < 0.1f) params.truncOctTrunc  = 0.1f;
            break;

        case NanoShape::Cuboctahedron:
            ImGui::SetNextItemWidth(180.0f);
            ImGui::DragFloat("Radius (A)##cubo", &params.cuboRadius, 0.5f, 0.1f, 500.0f, "%.2f");
            ImGui::TextDisabled("Cond: |x|+|y|, |y|+|z|, |x|+|z| <= R");
            if (params.cuboRadius < 0.1f) params.cuboRadius = 0.1f;
            break;

        default:
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
    if (ImGui::BeginTable("##nanoResultSummary", 2, ImGuiTableFlags_SizingStretchProp))
    {
        ImGui::TableSetupColumn("Label", ImGuiTableColumnFlags_WidthFixed, 120.0f);
        ImGui::TableSetupColumn("Value", ImGuiTableColumnFlags_WidthStretch);

        auto row = [](const char* label, const char* fmt, ...) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextDisabled("%s", label);
            ImGui::TableNextColumn();
            va_list args;
            va_start(args, fmt);
            ImGui::TextWrappedV(fmt, args);
            va_end(args);
        };

        row("Mode", "%s", generationModeLabel(result.mode));
        if (result.mode == NanoGenerationMode::Shape)
            row("Shape", "%s", shapeLabel(result.shape));
        row("Input atoms", "%d", result.inputAtoms);
        row("Output atoms", "%d", result.outputAtoms);
        row("Est. diameter", "%.2f A", result.estimatedDiameter);
        if (result.mode == NanoGenerationMode::WulffConstruction)
            row("Wulff faces", "%d", result.wulffFaceCount);
        if (result.tilingUsed)
            row("Tiling", "%dx%dx%d (%s)",
                2 * result.repA + 1, 2 * result.repB + 1, 2 * result.repC + 1,
                result.repClamped ? "clamped" : "auto");
        else
            row("Tiling", "%s", "none (no unit cell in reference)");

        ImGui::EndTable();
    }
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
    m_browsDir = detectHomePath();
}

NanoCrystalBuilderDialog::~NanoCrystalBuilderDialog()
{
    if (m_previewFBO)      { glDeleteFramebuffers(1,  &m_previewFBO);    m_previewFBO = 0; }
    if (m_previewColorTex) { glDeleteTextures(1,      &m_previewColorTex); m_previewColorTex = 0; }
    if (m_previewDepthRbo) { glDeleteRenderbuffers(1, &m_previewDepthRbo); m_previewDepthRbo = 0; }
    if (m_wulffPreviewFBO)      { glDeleteFramebuffers(1,  &m_wulffPreviewFBO);    m_wulffPreviewFBO = 0; }
    if (m_wulffPreviewColorTex) { glDeleteTextures(1,      &m_wulffPreviewColorTex); m_wulffPreviewColorTex = 0; }
    if (m_wulffPreviewDepthRbo) { glDeleteRenderbuffers(1, &m_wulffPreviewDepthRbo); m_wulffPreviewDepthRbo = 0; }
    if (m_wulffLineVAO) { glDeleteVertexArrays(1, &m_wulffLineVAO); m_wulffLineVAO = 0; }
    if (m_wulffLineVBO) { glDeleteBuffers(1, &m_wulffLineVBO); m_wulffLineVBO = 0; }
    if (m_wulffPlaneVAO) { glDeleteVertexArrays(1, &m_wulffPlaneVAO); m_wulffPlaneVAO = 0; }
    if (m_wulffPlaneVBO) { glDeleteBuffers(1, &m_wulffPlaneVBO); m_wulffPlaneVBO = 0; }
    if (m_wulffPlaneProgram) { glDeleteProgram(m_wulffPlaneProgram); m_wulffPlaneProgram = 0; }

    if (m_previewShadow.depthFBO)
        glDeleteFramebuffers(1, &m_previewShadow.depthFBO);
    if (m_previewShadow.depthTexture)
        glDeleteTextures(1, &m_previewShadow.depthTexture);
    delete m_previewSphere;
    delete m_previewCylinder;
}

void NanoCrystalBuilderDialog::ensureWulffFamilyColors(std::size_t familyCount)
{
    if (m_wulffFamilyColors.size() >= familyCount)
    {
        m_wulffFamilyColors.resize(familyCount);
        return;
    }

    const std::size_t oldSize = m_wulffFamilyColors.size();
    m_wulffFamilyColors.resize(familyCount);
    for (std::size_t index = oldSize; index < familyCount; ++index)
        m_wulffFamilyColors[index] = wulffFamilyColor((int)index);
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

    glGenVertexArrays(1, &m_wulffLineVAO);
    glGenBuffers(1, &m_wulffLineVBO);
    glBindVertexArray(m_wulffLineVAO);
    glBindBuffer(GL_ARRAY_BUFFER, m_wulffLineVBO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glBindVertexArray(0);

    m_wulffPlaneProgram = createProgram(kWulffPlaneVS, kWulffPlaneFS);
    glGenVertexArrays(1, &m_wulffPlaneVAO);
    glGenBuffers(1, &m_wulffPlaneVBO);
    glBindVertexArray(m_wulffPlaneVAO);
    glBindBuffer(GL_ARRAY_BUFFER, m_wulffPlaneVBO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glBindVertexArray(0);

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
    m_wulffPreviewDirty = true;
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

void NanoCrystalBuilderDialog::ensureWulffPreviewFBO(int w, int h)
{
    if (w == m_wulffPreviewW && h == m_wulffPreviewH && m_wulffPreviewFBO != 0)
        return;

    if (m_wulffPreviewFBO)      { glDeleteFramebuffers(1,  &m_wulffPreviewFBO);    m_wulffPreviewFBO = 0; }
    if (m_wulffPreviewColorTex) { glDeleteTextures(1,      &m_wulffPreviewColorTex); m_wulffPreviewColorTex = 0; }
    if (m_wulffPreviewDepthRbo) { glDeleteRenderbuffers(1, &m_wulffPreviewDepthRbo); m_wulffPreviewDepthRbo = 0; }

    glGenTextures(1, &m_wulffPreviewColorTex);
    glBindTexture(GL_TEXTURE_2D, m_wulffPreviewColorTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenRenderbuffers(1, &m_wulffPreviewDepthRbo);
    glBindRenderbuffer(GL_RENDERBUFFER, m_wulffPreviewDepthRbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    glGenFramebuffers(1, &m_wulffPreviewFBO);
    glBindFramebuffer(GL_FRAMEBUFFER, m_wulffPreviewFBO);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                           GL_TEXTURE_2D, m_wulffPreviewColorTex, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER, m_wulffPreviewDepthRbo);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    m_wulffPreviewW = w;
    m_wulffPreviewH = h;
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

void NanoCrystalBuilderDialog::rebuildWulffPreviewGeometry(const WulffPreview& preview)
{
    m_wulffPreviewBatches.clear();
    if (!preview.success)
        return;

    std::size_t familyCount = 0;
    for (const WulffFace& face : preview.faces)
        if (face.familyIndex >= 0)
            familyCount = std::max(familyCount, (std::size_t)face.familyIndex + 1);
    ensureWulffFamilyColors(familyCount);

    m_wulffPreviewBatches.resize(preview.faces.size());
    for (std::size_t faceIndex = 0; faceIndex < preview.faces.size(); ++faceIndex)
    {
        const WulffFace& face = preview.faces[faceIndex];
        WulffPreviewBatch& batch = m_wulffPreviewBatches[faceIndex];
        if (face.familyIndex >= 0 && (std::size_t)face.familyIndex < m_wulffFamilyColors.size())
            batch.color = m_wulffFamilyColors[(std::size_t)face.familyIndex];
        else
            batch.color = wulffFamilyColor(face.familyIndex);
        batch.normal = face.normal;
        batch.triangleVertices.clear();
        if (face.vertices.size() >= 3)
        {
            batch.triangleVertices.reserve((face.vertices.size() - 2) * 3);
            for (std::size_t vertexIndex = 1; vertexIndex + 1 < face.vertices.size(); ++vertexIndex)
            {
                batch.triangleVertices.push_back(face.vertices[0]);
                batch.triangleVertices.push_back(face.vertices[vertexIndex]);
                batch.triangleVertices.push_back(face.vertices[vertexIndex + 1]);
            }
        }
        batch.vertices.clear();
        if (face.vertices.size() < 2)
            continue;

        batch.vertices.reserve(face.vertices.size() * 2);
        for (std::size_t vertexIndex = 0; vertexIndex < face.vertices.size(); ++vertexIndex)
        {
            const glm::vec3& a = face.vertices[vertexIndex];
            const glm::vec3& b = face.vertices[(vertexIndex + 1) % face.vertices.size()];
            batch.vertices.push_back(a);
            batch.vertices.push_back(b);
        }
    }
}

void NanoCrystalBuilderDialog::autoFitWulffPreviewCamera()
{
    m_wulffCamYaw = 35.0f;
    m_wulffCamPitch = 30.0f;
    const float radius = std::max(m_wulffPreviewData.boundingRadius, 1.0f);
    const float halfFov = glm::radians(22.5f);
    float distance = radius / std::sin(halfFov) * 1.2f;
    distance = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, distance));
    m_wulffCamDistance = distance;
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
    buildFrameView(cam, m_previewBuffers, true, frame);

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

void NanoCrystalBuilderDialog::renderWulffPreviewToFBO(int w, int h)
{
    if (!m_glReady || !m_renderer || !m_wulffPreviewData.success)
        return;

    ensureWulffPreviewFBO(w, h);

    GLint prevFbo = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
    GLint prevVP[4];
    glGetIntegerv(GL_VIEWPORT, prevVP);

    SceneBuffers linesOnly;
    linesOnly.orbitCenter = m_wulffPreviewData.center;
    for (const WulffPreviewBatch& batch : m_wulffPreviewBatches)
        linesOnly.boxLines.insert(linesOnly.boxLines.end(), batch.vertices.begin(), batch.vertices.end());

    Camera cam;
    cam.yaw = m_wulffCamYaw;
    cam.pitch = m_wulffCamPitch;
    cam.distance = m_wulffCamDistance;

    FrameView frame;
    frame.framebufferWidth = w;
    frame.framebufferHeight = h;
    buildFrameView(cam, linesOnly, true, frame);

    glBindFramebuffer(GL_FRAMEBUFFER, m_wulffPreviewFBO);
    glViewport(0, 0, w, h);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.08f, 0.10f, 0.13f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (m_wulffPlaneProgram != 0)
    {
        glUseProgram(m_wulffPlaneProgram);
        glUniformMatrix4fv(glGetUniformLocation(m_wulffPlaneProgram, "projection"),
                           1, GL_FALSE, glm::value_ptr(frame.projection));
        glUniformMatrix4fv(glGetUniformLocation(m_wulffPlaneProgram, "view"),
                           1, GL_FALSE, glm::value_ptr(frame.view));
        glUniform3f(glGetUniformLocation(m_wulffPlaneProgram, "lightPos"),
                    frame.lightPosition.x, frame.lightPosition.y, frame.lightPosition.z);
        glUniform3f(glGetUniformLocation(m_wulffPlaneProgram, "viewPos"),
                    frame.cameraPosition.x, frame.cameraPosition.y, frame.cameraPosition.z);
        glDisable(GL_CULL_FACE);

        glBindVertexArray(m_wulffPlaneVAO);
        glBindBuffer(GL_ARRAY_BUFFER, m_wulffPlaneVBO);
        for (const WulffPreviewBatch& batch : m_wulffPreviewBatches)
        {
            if (batch.triangleVertices.empty())
                continue;

            glBufferData(GL_ARRAY_BUFFER,
                         (GLsizeiptr)(batch.triangleVertices.size() * sizeof(glm::vec3)),
                         batch.triangleVertices.data(),
                         GL_DYNAMIC_DRAW);
            glUniform3f(glGetUniformLocation(m_wulffPlaneProgram, "faceColor"),
                        batch.color.r, batch.color.g, batch.color.b);
            glUniform3f(glGetUniformLocation(m_wulffPlaneProgram, "faceNormal"),
                        batch.normal.x, batch.normal.y, batch.normal.z);
            glDrawArrays(GL_TRIANGLES, 0, (GLsizei)batch.triangleVertices.size());
        }
        glBindVertexArray(0);
        glUseProgram(0);
    }

    glBindVertexArray(m_wulffLineVAO);
    glBindBuffer(GL_ARRAY_BUFFER, m_wulffLineVBO);
    const glm::vec3 edgeColor(0.0f, 0.0f, 0.0f);
    for (const WulffPreviewBatch& batch : m_wulffPreviewBatches)
    {
        if (batch.vertices.empty())
            continue;
        glBufferData(GL_ARRAY_BUFFER,
                     (GLsizeiptr)(batch.vertices.size() * sizeof(glm::vec3)),
                     batch.vertices.data(),
                     GL_DYNAMIC_DRAW);
        m_renderer->drawBoxLines(frame.projection, frame.view,
                                 m_wulffLineVAO,
                                 batch.vertices.size(),
                                 edgeColor);
    }
    glBindVertexArray(0);

    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
    glViewport(prevVP[0], prevVP[1], prevVP[2], prevVP[3]);
}

// ---------------------------------------------------------------------------
// drawMenuItem
// ---------------------------------------------------------------------------

void NanoCrystalBuilderDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Nanocrystal", NULL, false, enabled))
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
    static NanoParams      params = []() {
        NanoParams defaults;
        defaults.generationMode = NanoGenerationMode::WulffConstruction;
        return defaults;
    }();
    static NanoBuildResult lastResult;

    if (params.wulffPlanes.empty())
        params.wulffPlanes.push_back(WulffPlaneInput{});

    ensureWulffFamilyColors(params.wulffPlanes.size());

    // Consume any pending file drop
    if (!m_pendingDropPath.empty()) {
        tryLoadFile(m_pendingDropPath, elementRadii, elementShininess);
        m_pendingDropPath.clear();
    }

    if (m_openRequested) {
        ImGui::OpenPopup("Build Nanocrystal");
        lastResult      = {};
        if (m_reference.atoms.empty() && !structure.atoms.empty()) {
            m_reference = structure;
            m_previewBufDirty = true;
            m_wulffPreviewDirty = true;
        }
        m_openRequested = false;
    }

    if (!m_reference.atoms.empty())
    {
        const std::size_t signature = hashCombine(structurePreviewHash(m_reference), wulffParamsHash(params));
        if (m_wulffPreviewDirty || signature != m_wulffPreviewSignature)
        {
            m_wulffPreviewData = computeWulffPreview(m_reference, params);
            rebuildWulffPreviewGeometry(m_wulffPreviewData);
            autoFitWulffPreviewCamera();
            m_wulffPreviewSignature = signature;
            m_wulffPreviewDirty = false;
        }
    }

    m_isOpen = ImGui::IsPopupOpen("Build Nanocrystal");

    ImGui::SetNextWindowSize(ImVec2(1160.0f, 900.0f), ImGuiCond_FirstUseEver);
    bool dialogOpen = true;
    if (!ImGui::BeginPopupModal("Build Nanocrystal", &dialogOpen, 0)) {
        m_isOpen = false;
        return;
    }
    m_isOpen = true;

    ImGui::TextWrapped("Drop a reference crystal into the left panel or use the currently loaded structure, then build the nanocrystal.");
    ImGui::Separator();

    // =========================================================================
    // Layout: Left = stacked previews (reference top, Wulff bottom)
    //         Right = Builder options
    // =========================================================================

    constexpr float kLeftW          = 460.0f;
    constexpr float kColumnH        = 700.0f;
    constexpr float kTopPanelH      = 360.0f;
    const float     kBottomPanelH   = kColumnH - kTopPanelH - ImGui::GetStyle().ItemSpacing.y;
    constexpr float kStructureViewH = 284.0f;
    constexpr float kWulffViewH     = 250.0f;

    ImGui::BeginGroup();

    // ---- TOP-LEFT: Structure view with drag-and-drop ----
    ImGui::BeginChild("##nanoStructureView", ImVec2(kLeftW, kTopPanelH), true);

    ImGui::Text("Reference Structure");
    ImGui::Separator();

    if (m_reference.atoms.empty())
    {
        ImGui::TextWrapped("Drop a supported structure file into this panel.");
    }
    else
    {
        if (m_browsStatusMsg[0] != '\0')
            ImGui::TextWrapped("%s", m_browsStatusMsg);
        else
            ImGui::TextDisabled("Using the current structure as reference.");

        ImGui::SameLine();
        if (ImGui::Button("Clear##nanoClearRef", ImVec2(70.0f, 0.0f))) {
            m_reference       = {};
            m_previewBufDirty = true;
            m_wulffPreviewDirty = true;
            m_wulffPreviewData = {};
            m_wulffPreviewBatches.clear();
            m_wulffPreviewSignature = 0;
            m_browsStatusMsg[0] = '\0';
        }
    }
    ImGui::Spacing();

    // Drag-and-drop target using InvisibleButton
    const float dropH = kStructureViewH - ImGui::GetCursorPosY() + ImGui::GetStyle().WindowPadding.y + 6.0f;
    ImGui::InvisibleButton("##nanoRefDropZone", ImVec2(-1.0f, dropH > 60.0f ? dropH : 60.0f));

    ImDrawList* dl = ImGui::GetWindowDrawList();
    ImVec2 dropMin = ImGui::GetItemRectMin();
    ImVec2 dropMax = ImGui::GetItemRectMax();
    const bool dropZoneHovered = ImGui::IsItemHovered();
    const bool dropZoneActive = ImGui::IsItemActive();

    if (m_reference.atoms.empty()) {
        if (ImGui::BeginDragDropTarget()) {
            if (const ImGuiPayload* payload =
                    ImGui::AcceptDragDropPayload("DROP_FILE", ImGuiDragDropFlags_AcceptPeekOnly)) {
                const char* files = (const char*)payload->Data;
                tryLoadFile(files, elementRadii, elementShininess);
            }
            ImGui::EndDragDropTarget();
        }

        // Show centered drop hint (two lines)
        const char* hint1 = "Drop a structure file here";
        const char* hint2 = "(.cif .xyz .pdb .sdf .mol .vasp ...)"; 
        const float lh = ImGui::GetTextLineHeight();
        const float w1 = ImGui::CalcTextSize(hint1).x;
        const float w2 = ImGui::CalcTextSize(hint2).x;
        const ImVec2 mid((dropMin.x + dropMax.x) * 0.5f, (dropMin.y + dropMax.y) * 0.5f);
        const ImU32 hintCol = ImGui::GetColorU32(ImGuiCol_TextDisabled);
        dl->AddText(ImVec2(mid.x - w1 * 0.5f, mid.y - lh * 1.1f), hintCol, hint1);
        dl->AddText(ImVec2(mid.x - w2 * 0.5f, mid.y + lh * 0.1f), hintCol, hint2);
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
            }

            if (dropZoneHovered) {
                float wheel = ImGui::GetIO().MouseWheel;
                if (wheel != 0.0f) {
                    m_camDistance -= wheel * m_camDistance * 0.1f;
                    m_camDistance  = std::max(Camera::kMinDistance,
                                              std::min(Camera::kMaxDistance, m_camDistance));
                }
            }

            ImGui::TextDisabled("Drag to orbit, scroll to zoom");
        }
    }

    // Reference info bar
    if (!m_reference.atoms.empty()) {
        ImGui::Separator();
        if (m_reference.hasUnitCell) {
            const auto& cv = m_reference.cellVectors;
            float la = safeLen3(glm::vec3((float)cv[0][0],(float)cv[0][1],(float)cv[0][2]));
            float lb = safeLen3(glm::vec3((float)cv[1][0],(float)cv[1][1],(float)cv[1][2]));
            float lc = safeLen3(glm::vec3((float)cv[2][0],(float)cv[2][1],(float)cv[2][2]));
            ImGui::TextDisabled("%d atoms | a=%.3f  b=%.3f  c=%.3f A | periodic",
                                (int)m_reference.atoms.size(), la, lb, lc);
        } else {
            ImGui::TextDisabled("%d atoms | no unit cell | shape carving only",
                                (int)m_reference.atoms.size());
        }
    }

    ImGui::EndChild(); // ##nanoStructureView

    // ---- BOTTOM-LEFT: Wulff plane preview ----
    ImGui::BeginChild("##nanoWulffPreview", ImVec2(kLeftW, kBottomPanelH), true);

    ImGui::Text("Wulff Plane View");
    ImGui::Separator();
    if (params.generationMode != NanoGenerationMode::WulffConstruction)
    {
        ImGui::TextWrapped("Switch generation mode to Wulff Construction to preview symmetry-expanded planes.");
    }
    else if (m_reference.atoms.empty())
    {
        ImGui::TextWrapped("Load a periodic reference structure to generate the Wulff plane preview.");
    }
    else if (!m_wulffPreviewData.success)
    {
        ImGui::TextWrapped("%s", m_wulffPreviewData.message.c_str());
    }
    else
    {
        ImGui::TextDisabled("%d faces from %d symmetry-expanded planes",
                            (int)m_wulffPreviewData.faces.size(),
                            (int)m_wulffPreviewData.planes.size());
        ImGui::Spacing();

        ImGui::InvisibleButton("##nanoWulffDropZone", ImVec2(-1.0f, kWulffViewH));
        ImDrawList* wulffDl = ImGui::GetWindowDrawList();
        ImVec2 previewMin = ImGui::GetItemRectMin();
        ImVec2 previewMax = ImGui::GetItemRectMax();
        wulffDl->AddRect(previewMin, previewMax, ImGui::GetColorU32(ImGuiCol_Border), 2.0f);

        const ImVec2 previewSize(previewMax.x - previewMin.x - 10.0f,
                                 previewMax.y - previewMin.y - 10.0f);
        const int previewW = std::max(1, (int)previewSize.x);
        const int previewH = std::max(1, (int)previewSize.y);
        renderWulffPreviewToFBO(previewW, previewH);

        wulffDl->AddImage((ImTextureID)(intptr_t)m_wulffPreviewColorTex,
                  ImVec2(previewMin.x + 5.0f, previewMin.y + 5.0f),
                  ImVec2(previewMin.x + 5.0f + previewSize.x, previewMin.y + 5.0f + previewSize.y),
                  ImVec2(0, 1), ImVec2(1, 0));

        if (ImGui::IsItemActive() && ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f))
        {
            const ImVec2 delta = ImGui::GetIO().MouseDelta;
            m_wulffCamYaw -= delta.x * 0.5f;
            m_wulffCamPitch += delta.y * 0.5f;
        }
        if (ImGui::IsItemHovered())
        {
            const float wheel = ImGui::GetIO().MouseWheel;
            if (wheel != 0.0f)
            {
                m_wulffCamDistance -= wheel * m_wulffCamDistance * 0.1f;
                m_wulffCamDistance = std::max(Camera::kMinDistance,
                                              std::min(Camera::kMaxDistance, m_wulffCamDistance));
            }
        }

        ImGui::TextDisabled("Left-drag = orbit   Scroll = zoom");
    }

    ImGui::EndChild(); // ##nanoWulffPreview

    ImGui::EndGroup();

    ImGui::SameLine();

    // ---- RIGHT: Builder options ----
    ImGui::BeginChild("##nanoBuilderOptions", ImVec2(0.0f, kColumnH), true);

    ImGui::Text("Builder Options");
    ImGui::Separator();

    int generationMode = (int)params.generationMode;
    const char* generationLabels[] = {"Shape Cut", "Wulff Construction"};
    ImGui::SetNextItemWidth(-1.0f);
    if (ImGui::Combo("##nanoGenMode", &generationMode, generationLabels, 2))
    {
        params.generationMode = (NanoGenerationMode)generationMode;
        m_wulffPreviewDirty = true;
        lastResult = {};
    }

    ImGui::Spacing();
    ImGui::Separator();

    if (params.generationMode == NanoGenerationMode::Shape)
    {
        ImGui::Text("Shape");
        static const char* kShapeLabels[] = {
            "Sphere","Ellipsoid","Box","Cylinder",
            "Octahedron","Truncated Octahedron","Cuboctahedron"
        };
        constexpr int kNumShapeOptions = 7;
        int shapeInt = (int)params.shape;
        if (shapeInt >= kNumShapeOptions) { shapeInt = 0; params.shape = NanoShape::Sphere; }
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::Combo("##nanoShapeCombo", &shapeInt, kShapeLabels, kNumShapeOptions))
            params.shape = (NanoShape)shapeInt;

        ImGui::Spacing();
        ImGui::Separator();
    }

    // Scrollable parameters area
    ImGui::BeginChild("##nanoParamsScroll", ImVec2(-1, -50), true);

    if (params.generationMode == NanoGenerationMode::Shape)
        drawShapeParameters(params);
    else
        drawWulffParameters(params);

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Text("Carving center");
    ImGui::Checkbox("Auto-center from atoms##nano", &params.autoCenterFromAtoms);
    if (!params.autoCenterFromAtoms) {
        ImGui::SetNextItemWidth(130.0f);
        ImGui::DragFloat("X (A)##cnano", &params.cx, 0.1f, -1000.f, 1000.f, "%.3f");
        ImGui::SetNextItemWidth(130.0f);
        ImGui::DragFloat("Y (A)##cnano", &params.cy, 0.1f, -1000.f, 1000.f, "%.3f");
        ImGui::SetNextItemWidth(130.0f);
        ImGui::DragFloat("Z (A)##cnano", &params.cz, 0.1f, -1000.f, 1000.f, "%.3f");
    }

    if (m_reference.hasUnitCell) {
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("Supercell replication");
        ImGui::Checkbox("Auto-replicate##nano", &params.autoReplicate);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Automatically tile the unit cell to fully cover the generated shape");
        if (!params.autoReplicate) {
            ImGui::SetNextItemWidth(80.0f);
            ImGui::DragInt("Along a##nano", &params.repA, 0.2f, 1, 40);
            ImGui::SetNextItemWidth(80.0f);
            ImGui::DragInt("Along b##nano", &params.repB, 0.2f, 1, 40);
            ImGui::SetNextItemWidth(80.0f);
            ImGui::DragInt("Along c##nano", &params.repC, 0.2f, 1, 40);
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

        if (params.generationMode == NanoGenerationMode::Shape)
        {
            const HalfExtents he = computeShapeHalfExtents(params);
            ImGui::TextDisabled("%.2fx%.2fx%.2f A",
                                2.f*(he.hx+params.vacuumPadding),
                                2.f*(he.hy+params.vacuumPadding),
                                2.f*(he.hz+params.vacuumPadding));
        }
        else if (m_wulffPreviewData.success)
        {
            const glm::vec3 extents = m_wulffPreviewData.maxPoint - m_wulffPreviewData.minPoint;
            ImGui::TextDisabled("%.2fx%.2fx%.2f A",
                                extents.x + 2.0f*params.vacuumPadding,
                                extents.y + 2.0f*params.vacuumPadding,
                                extents.z + 2.0f*params.vacuumPadding);
        }
    }

    if (params.generationMode == NanoGenerationMode::WulffConstruction)
    {
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("Facet colors");
        for (std::size_t index = 0; index < params.wulffPlanes.size(); ++index)
        {
            const WulffPlaneInput& plane = params.wulffPlanes[index];
            glm::vec3& familyColor = m_wulffFamilyColors[index];
            float color[3] = { familyColor.r, familyColor.g, familyColor.b };
            ImGui::PushID((int)index + 1000);
            ImGui::SetNextItemWidth(120.0f);
            if (ImGui::ColorEdit3("##wulffLegendColor", color,
                                  ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel))
            {
                familyColor = glm::vec3(color[0], color[1], color[2]);
                rebuildWulffPreviewGeometry(m_wulffPreviewData);
            }
            ImGui::SameLine();
            ImGui::PushTextWrapPos(0.0f);
            ImGui::TextWrapped("(%d %d %d)  E = %.4f", plane.h, plane.k, plane.l, plane.surfaceEnergy);
            ImGui::PopTextWrapPos();
            ImGui::PopID();
        }
        if (!m_wulffPreviewData.message.empty())
        {
            ImGui::Spacing();
            ImGui::TextWrapped("Preview: %s", m_wulffPreviewData.message.c_str());
        }
    }

    ImGui::EndChild(); // ##nanoParamsScroll

    ImGui::Separator();
    ImGui::Text("Result");
    drawNanoBuildResult(lastResult);

    ImGui::EndChild(); // ##nanoBuilderOptions

    // =========================================================================
    // Action bar
    // =========================================================================

    ImGui::Spacing();
    ImGui::Separator();
    const bool canBuild = !m_reference.atoms.empty()
                       && (params.generationMode != NanoGenerationMode::WulffConstruction
                           || m_wulffPreviewData.success);
    if (!canBuild) ImGui::BeginDisabled();
    if (ImGui::Button("Build Nanocrystal##nano", ImVec2(160.0f, 0.0f))) {
        static const std::vector<glm::vec3>       kNoVerts;
        static const std::vector<unsigned int>    kNoIdx;
        lastResult = buildNanocrystal(structure,
                                      m_reference,
                                      params,
                                      elementColors,
                                      kNoVerts,
                                      kNoIdx);
        if (lastResult.success)
            updateBuffers(structure);
    }
    if (!canBuild) ImGui::EndDisabled();

    ImGui::SameLine(0.0f, 8.0f);
    if (ImGui::Button("Close##nano", ImVec2(80.0f, 0.0f))) {
        lastResult  = {};
        m_isOpen    = false;
        ImGui::CloseCurrentPopup();
    }

    ImGui::EndPopup();
}
