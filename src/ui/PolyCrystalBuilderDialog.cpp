#include "ui/PolyCrystalBuilderDialog.h"
#include "algorithms/PolyCrystalBuilder.h"
#include "util/PathUtils.h"
#include "io/StructureLoader.h"

#include "graphics/CylinderMesh.h"
#include "graphics/Renderer.h"
#include "graphics/SphereMesh.h"
#include "graphics/StructureInstanceBuilder.h"
#include "graphics/ShadowMap.h"
#include "app/SceneView.h"
#include "camera/Camera.h"
#include "ui/ThemeUtils.h"
#include "imgui.h"

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// ===========================================================================
// PolyCrystalBuilderDialog implementation
// ===========================================================================

PolyCrystalBuilderDialog::PolyCrystalBuilderDialog()
    : m_glReady(false)
    , m_previewBufDirty(true)
    , m_camYaw(45.0f)
    , m_camPitch(35.0f)
    , m_camDistance(10.0f)
{
    m_statusMsg[0] = '\0';
}

PolyCrystalBuilderDialog::~PolyCrystalBuilderDialog()
{
    if (m_previewFBO)      { glDeleteFramebuffers(1,  &m_previewFBO);      m_previewFBO = 0; }
    if (m_previewColorTex) { glDeleteTextures(1,      &m_previewColorTex); m_previewColorTex = 0; }
    if (m_previewDepthRbo) { glDeleteRenderbuffers(1, &m_previewDepthRbo); m_previewDepthRbo = 0; }

    if (m_previewShadow.depthFBO)
        glDeleteFramebuffers(1, &m_previewShadow.depthFBO);
    if (m_previewShadow.depthTexture)
        glDeleteTextures(1, &m_previewShadow.depthTexture);

    delete m_previewSphere;
    delete m_previewCylinder;
}

void PolyCrystalBuilderDialog::initRenderResources(Renderer& renderer)
{
    m_renderer = &renderer;
    m_previewSphere   = new SphereMesh(24, 24);
    m_previewCylinder = new CylinderMesh(16);
    m_previewBuffers.init(m_previewSphere->vbo, m_previewSphere->ebo, m_previewSphere->indexCount,
                          m_previewCylinder->vbo, m_previewCylinder->vertexCount);
    m_previewShadow = createShadowMap(1, 1);
    m_glReady = true;
}

// ---------------------------------------------------------------------------

void PolyCrystalBuilderDialog::feedDroppedFile(const std::string& path)
{
    m_pendingDropPath = path;
}

bool PolyCrystalBuilderDialog::tryLoadFile(const std::string& path,
                                            const std::vector<float>& radii,
                                            const std::vector<float>& shininess)
{
    Structure loaded;
    std::string err;
    if (!loadStructureFromFile(path, loaded, err)) {
        std::snprintf(m_statusMsg, sizeof(m_statusMsg),
                      "Load failed: %s", err.c_str());
        return false;
    }
    if (!loaded.hasUnitCell) {
        std::snprintf(m_statusMsg, sizeof(m_statusMsg),
                      "Error: structure has no unit cell (required for polycrystal).");
        return false;
    }
    m_reference = std::move(loaded);

    // Extract basename for display
    std::string::size_type sep = path.find_last_of("/\\");
    m_referenceFilename = (sep != std::string::npos) ? path.substr(sep + 1) : path;

    std::snprintf(m_statusMsg, sizeof(m_statusMsg),
                  "Loaded: %d atoms", (int)m_reference.atoms.size());

    m_previewBufDirty = true;
    if (m_glReady) {
        rebuildPreviewBuffers(radii, shininess);
        autoFitPreviewCamera();
    }

    std::cout << "[PolyCrystalBuilder] Reference loaded: " << path
              << " (" << m_reference.atoms.size() << " atoms)" << std::endl;
    return true;
}

// ---------------------------------------------------------------------------
// Preview rendering helpers
// ---------------------------------------------------------------------------

void PolyCrystalBuilderDialog::ensurePreviewFBO(int w, int h)
{
    if (w == m_previewW && h == m_previewH && m_previewFBO != 0)
        return;

    if (m_previewFBO)      { glDeleteFramebuffers(1,  &m_previewFBO);      m_previewFBO = 0; }
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

void PolyCrystalBuilderDialog::rebuildPreviewBuffers(const std::vector<float>& radii,
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

void PolyCrystalBuilderDialog::autoFitPreviewCamera()
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

    const float halfFov = glm::radians(22.5f);
    float dist = maxR / std::sin(halfFov) * 1.15f;
    dist = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, dist));
    m_camDistance = dist;
}

void PolyCrystalBuilderDialog::renderPreviewToFBO(int w, int h)
{
    if (!m_glReady || !m_renderer || m_previewBuffers.atomCount == 0)
        return;

    ensurePreviewFBO(w, h);

    GLint prevFbo = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
    GLint prevVP[4];
    glGetIntegerv(GL_VIEWPORT, prevVP);

    Camera cam;
    cam.yaw      = m_camYaw;
    cam.pitch    = m_camPitch;
    cam.distance = m_camDistance;

    FrameView frame;
    frame.framebufferWidth  = w;
    frame.framebufferHeight = h;
    buildFrameView(cam, m_previewBuffers, true, frame);

    glBindFramebuffer(GL_FRAMEBUFFER, m_previewFBO);
    glViewport(0, 0, w, h);
    glEnable(GL_DEPTH_TEST);
    {
        const ImVec4& bg = ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
        glClearColor(bg.x, bg.y, bg.z, bg.w);
    }
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    m_renderer->drawBonds(frame.projection, frame.view,
                          frame.lightPosition, frame.cameraPosition,
                          m_previewBuffers.tabCylinderVAO, m_previewBuffers.tabCylinderVertexCount,
                          m_previewBuffers.bondCount);

    m_renderer->drawAtoms(frame.projection, frame.view,
                          frame.lightMVP, frame.lightPosition, frame.cameraPosition,
                          m_previewShadow,
                          m_previewBuffers.tabSphereVAO, m_previewBuffers.tabSphereIndexCount,
                          m_previewBuffers.atomCount);

    m_renderer->drawBoxLines(frame.projection, frame.view,
                             m_previewBuffers.lineVAO,
                             m_previewBuffers.boxLines.size());

    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
    glViewport(prevVP[0], prevVP[1], prevVP[2], prevVP[3]);
}

// ---------------------------------------------------------------------------
// Menu item
// ---------------------------------------------------------------------------

void PolyCrystalBuilderDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Polycrystal", NULL, false, enabled))
        m_openRequested = true;
}

// ---------------------------------------------------------------------------
// Main dialog
// ---------------------------------------------------------------------------

void PolyCrystalBuilderDialog::drawDialog(
    Structure& structure,
    const std::vector<glm::vec3>& elementColors,
    const std::vector<float>& elementRadii,
    const std::vector<float>& elementShininess,
    const std::function<void(Structure&)>& updateBuffers)
{
    static PolyParams      params;
    static PolyBuildResult lastResult;

    // Consume pending drop
    if (!m_pendingDropPath.empty()) {
        tryLoadFile(m_pendingDropPath, elementRadii, elementShininess);
        m_pendingDropPath.clear();
    }

    if (m_openRequested) {
        ImGui::OpenPopup("Build Polycrystal");
        lastResult      = {};
        m_openRequested = false;
        m_statusMsg[0]  = '\0';
    }

    m_isOpen = ImGui::IsPopupOpen("Build Polycrystal");

    ImGui::SetNextWindowSize(ImVec2(820.0f, 700.0f), ImGuiCond_FirstUseEver);
    bool dialogOpen = true;
    if (!ImGui::BeginPopupModal("Build Polycrystal", &dialogOpen, 0)) {
        m_isOpen = false;
        return;
    }
    m_isOpen = true;

    // =========================================================================
    // Left panel: reference structure preview with drag-and-drop
    // =========================================================================

    constexpr float kLeftW  = 340.0f;
    constexpr float kPanelH = 560.0f;
    constexpr float kPrevH  = 320.0f;

    ImGui::BeginChild("##polyRefPanel", ImVec2(kLeftW, kPanelH), true);

    ImGui::Text("Reference Single Crystal");
    ImGui::Separator();

    // Show filename if loaded
    if (!m_referenceFilename.empty()) {
        ImGui::Text("File: %s", m_referenceFilename.c_str());
        ImGui::Text("%d atoms", (int)m_reference.atoms.size());
        if (m_reference.hasUnitCell) {
            const auto& cv = m_reference.cellVectors;
            float la = glm::length(glm::vec3((float)cv[0][0],(float)cv[0][1],(float)cv[0][2]));
            float lb = glm::length(glm::vec3((float)cv[1][0],(float)cv[1][1],(float)cv[1][2]));
            float lc = glm::length(glm::vec3((float)cv[2][0],(float)cv[2][1],(float)cv[2][2]));
            ImGui::TextDisabled("a=%.3f  b=%.3f  c=%.3f A", la, lb, lc);
        }
    } else {
        ImGui::TextDisabled("(none loaded)");
    }

    ImGui::Spacing();

    // Drop zone / 3D preview area
    ImGui::InvisibleButton("##polyRefDropZone", ImVec2(-1.0f, kPrevH));
    ImDrawList* dl = ImGui::GetWindowDrawList();
    ImVec2 dropMin = ImGui::GetItemRectMin();
    ImVec2 dropMax = ImGui::GetItemRectMax();
    const bool dropZoneHovered = ImGui::IsItemHovered();
    const bool dropZoneActive  = ImGui::IsItemActive();
    dl->AddRect(dropMin, dropMax, ImGui::GetColorU32(ImGuiCol_Border), 2.0f);

    if (m_reference.atoms.empty()) {
        // Show drop hint
        ImVec2 dropMid((dropMin.x + dropMax.x) / 2.0f,
                       (dropMin.y + dropMax.y) / 2.0f);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize(),
                    ImVec2(dropMid.x - 90.0f, dropMid.y - 8.0f),
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
                         previewMin, previewMax,
                         ImVec2(0, 1), ImVec2(1, 0));

            // Orbit with left-drag
            if (dropZoneActive && ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f)) {
                ImVec2 delta = ImGui::GetIO().MouseDelta;
                m_camYaw   -= delta.x * 0.5f;
                m_camPitch += delta.y * 0.5f;
            }

            // Zoom with scroll
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

    ImGui::Spacing();
    if (m_statusMsg[0] != '\0')
        ImGui::TextWrapped("%s", m_statusMsg);

    ImGui::EndChild(); // ##polyRefPanel

    ImGui::SameLine();

    // =========================================================================
    // Right panel: builder options
    // =========================================================================

    ImGui::BeginChild("##polyBuilderOptions", ImVec2(0, kPanelH), true);

    ImGui::Text("Builder Options");
    ImGui::Separator();

    // --- Box size ---
    ImGui::Text("Simulation Box Size (Angstroms)");
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputFloat("X##polySize", &params.sizeX, 0.0f, 0.0f, "%.2f");
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputFloat("Y##polySize", &params.sizeY, 0.0f, 0.0f, "%.2f");
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputFloat("Z##polySize", &params.sizeZ, 0.0f, 0.0f, "%.2f");
    if (params.sizeX < 1.0f) params.sizeX = 1.0f;
    if (params.sizeY < 1.0f) params.sizeY = 1.0f;
    if (params.sizeZ < 1.0f) params.sizeZ = 1.0f;

    ImGui::Spacing();
    ImGui::Separator();

    // --- Number of grains ---
    ImGui::Text("Number of Grains");
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputInt("##polyNumGrains", &params.numGrains);
    if (params.numGrains < 1) params.numGrains = 1;
    if (params.numGrains > 500) params.numGrains = 500;

    ImGui::Spacing();
    ImGui::Separator();

    // --- Seed ---
    ImGui::Text("Random Seed");
    ImGui::SetNextItemWidth(120.0f);
    ImGui::InputInt("##polySeed", &params.seed);

    ImGui::Spacing();
    ImGui::Separator();

    // --- Orientation mode ---
    ImGui::Text("Grain Orientations");
    const char* modeLabels[] = { "All Random", "All Specified", "Partial (specify some, rest random)" };
    int modeInt = (int)params.orientationMode;
    ImGui::SetNextItemWidth(-1.0f);
    if (ImGui::Combo("##polyOrientMode", &modeInt, modeLabels, 3))
        params.orientationMode = (GrainOrientationMode)modeInt;

    ImGui::Spacing();

    // Orientation table (shown for Specified / Partial modes)
    if (params.orientationMode == GrainOrientationMode::AllSpecified)
    {
        if ((int)params.specifiedOrientations.size() != params.numGrains)
        {
            params.specifiedOrientations.resize(params.numGrains);
            for (int i = 0; i < params.numGrains; ++i)
                params.specifiedOrientations[i].grainIndex = i;
        }

        ImGui::TextDisabled("Bunge Euler angles (phi1, Phi, phi2) in degrees");
        ImGui::BeginChild("##polyOrientTable", ImVec2(-1, 180), true);

        for (int i = 0; i < params.numGrains; ++i)
        {
            auto& o = params.specifiedOrientations[i];
            ImGui::PushID(i);
            ImGui::Text("Grain %d:", i + 1);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(70.0f);
            ImGui::InputFloat("phi1##grain", &o.phi1, 0.0f, 0.0f, "%.1f");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(70.0f);
            ImGui::InputFloat("Phi##grain", &o.Phi, 0.0f, 0.0f, "%.1f");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(70.0f);
            ImGui::InputFloat("phi2##grain", &o.phi2, 0.0f, 0.0f, "%.1f");
            ImGui::PopID();
        }

        ImGui::EndChild();
    }
    else if (params.orientationMode == GrainOrientationMode::PartialSpecified)
    {
        ImGui::TextDisabled(
            "Specify orientations for selected grains. Unspecified grains get random orientations.");
        ImGui::TextDisabled("Bunge Euler angles (phi1, Phi, phi2) in degrees");

        static int newGrainIdx = 1;
        ImGui::SetNextItemWidth(70.0f);
        ImGui::InputInt("Grain #", &newGrainIdx);
        if (newGrainIdx < 1) newGrainIdx = 1;
        if (newGrainIdx > params.numGrains) newGrainIdx = params.numGrains;

        ImGui::SameLine();
        if (ImGui::Button("Add Grain Orientation##poly"))
        {
            bool found = false;
            for (const auto& o : params.specifiedOrientations) {
                if (o.grainIndex == newGrainIdx - 1) { found = true; break; }
            }
            if (!found) {
                GrainOrientation go;
                go.grainIndex = newGrainIdx - 1;
                params.specifiedOrientations.push_back(go);
            }
        }

        ImGui::BeginChild("##polyPartialOrientTable", ImVec2(-1, 180), true);

        int removeIdx = -1;
        for (int i = 0; i < (int)params.specifiedOrientations.size(); ++i)
        {
            auto& o = params.specifiedOrientations[i];
            ImGui::PushID(i);
            ImGui::Text("Grain %d:", o.grainIndex + 1);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(70.0f);
            ImGui::InputFloat("phi1##pgrain", &o.phi1, 0.0f, 0.0f, "%.1f");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(70.0f);
            ImGui::InputFloat("Phi##pgrain", &o.Phi, 0.0f, 0.0f, "%.1f");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(70.0f);
            ImGui::InputFloat("phi2##pgrain", &o.phi2, 0.0f, 0.0f, "%.1f");
            ImGui::SameLine();
            if (ImGui::Button("X##removePGrain"))
                removeIdx = i;
            ImGui::PopID();
        }
        if (removeIdx >= 0)
            params.specifiedOrientations.erase(
                params.specifiedOrientations.begin() + removeIdx);

        ImGui::EndChild();
    }
    else
    {
        ImGui::TextDisabled("All grains will receive uniformly random orientations.");
    }

    ImGui::Spacing();
    ImGui::Separator();

    // --- Result ---
    ImGui::Text("Result");
    if (lastResult.success)
    {
        ImGui::TextColored(ImVec4(0.2f, 0.9f, 0.3f, 1.0f), "%s", lastResult.message.c_str());
    }
    else if (!lastResult.message.empty())
    {
        ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.3f, 1.0f), "%s", lastResult.message.c_str());
    }
    else
    {
        ImGui::TextDisabled("(not built yet)");
    }

    ImGui::EndChild(); // ##polyBuilderOptions

    // =========================================================================
    // Action bar
    // =========================================================================

    const bool canBuild = !m_reference.atoms.empty() && m_reference.hasUnitCell;
    if (!canBuild) ImGui::BeginDisabled();
    if (ImGui::Button("Build##poly", ImVec2(100.0f, 0.0f)))
    {
        lastResult = buildPolycrystal(structure, m_reference, params, elementColors);
        if (lastResult.success)
            updateBuffers(structure);
    }
    if (!canBuild) ImGui::EndDisabled();

    ImGui::SameLine();
    if (ImGui::Button("Close##poly", ImVec2(80.0f, 0.0f)))
    {
        lastResult = {};
        m_isOpen   = false;
        ImGui::CloseCurrentPopup();
    }

    ImGui::EndPopup();
}
