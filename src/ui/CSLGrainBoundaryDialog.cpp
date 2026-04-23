#include "ui/CSLGrainBoundaryDialog.h"
#include "algorithms/CSLComputation.h"

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

#include <cstdio>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>

namespace
{

void drawBuildResultSummary(const CSLBuildResult& result)
{
    if (result.message.empty())
        return;

    if (result.success)
    {
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.4f, 1.0f, 0.5f, 1.0f));
        ImGui::TextWrapped("OK: %s", result.message.c_str());
        ImGui::PopStyleColor();
    }
    else
    {
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.45f, 0.35f, 1.0f));
        ImGui::TextWrapped("Error: %s", result.message.c_str());
        ImGui::PopStyleColor();
        return;
    }

    ImGui::Spacing();

    if (ImGui::BeginTable("##ResultSummary", 2, ImGuiTableFlags_SizingStretchProp))
    {
        ImGui::TableSetupColumn("Label", ImGuiTableColumnFlags_WidthFixed, 220.0f);
        ImGui::TableSetupColumn("Value");

        auto row = [](const char* label, const char* fmt, ...) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextDisabled("%s", label);
            ImGui::TableNextColumn();
            va_list args;
            va_start(args, fmt);
            ImGui::TextV(fmt, args);
            va_end(args);
        };

        row("Input atoms", "%d", result.inputAtoms);
        row("Output atoms", "%d", result.generatedAtoms);
        row("Removed overlaps", "%d", result.removedOverlap);
        row("Sigma", "%d", result.sigma);
        row("Misorientation angle", "%.4f deg", result.thetaDeg);
        row("Boundary type", "%s", result.boundaryType.empty() ? "unknown" : result.boundaryType.c_str());

        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::Separator();
        ImGui::TableNextColumn();
        ImGui::Separator();

        row("Axis [u v w]", "[%d %d %d]", result.axis[0], result.axis[1], result.axis[2]);
        row("GB plane (h k l)", "(%d %d %d)", result.gbPlane[0], result.gbPlane[1], result.gbPlane[2]);
        row("CSL matrix", "[%d %d %d]  [%d %d %d]  [%d %d %d]",
            result.csl[0][0], result.csl[0][1], result.csl[0][2],
            result.csl[1][0], result.csl[1][1], result.csl[1][2],
            result.csl[2][0], result.csl[2][1], result.csl[2][2]);
        row("Grain A / Grain B", "%d / %d unit cells", result.ucA, result.ucB);
        row("Vacuum padding", "%.3f A", result.vacuumPadding);
        row("Gap", "%.3f A", result.gap);
        row("Overlap distance", "%.3f A", result.overlapDist);

        ImGui::EndTable();
    }
}

} // namespace

// -- Dialog methods -------------------------------------------------
CSLGrainBoundaryDialog::CSLGrainBoundaryDialog() = default;

CSLGrainBoundaryDialog::~CSLGrainBoundaryDialog()
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

void CSLGrainBoundaryDialog::initRenderResources(Renderer& renderer)
{
    m_renderer        = &renderer;
    m_previewSphere   = new SphereMesh(24, 24);
    m_previewCylinder = new CylinderMesh(16);
    m_previewBuffers.init(m_previewSphere->vbo, m_previewSphere->ebo, m_previewSphere->indexCount,
                          m_previewCylinder->vbo, m_previewCylinder->vertexCount);
    m_previewShadow   = createShadowMap(1, 1);
    m_glReady         = true;
}

void CSLGrainBoundaryDialog::ensurePreviewFBO(int w, int h)
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

void CSLGrainBoundaryDialog::rebuildPreviewBuffers(const Structure& s,
                                                    const std::vector<float>& radii,
                                                    const std::vector<float>& shininess)
{
    if (!m_glReady || s.atoms.empty())
        return;

    static const int kIdent[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    StructureInstanceData data = buildStructureInstanceData(
        s, false, kIdent, radii, shininess);

    std::array<bool,119> noFilter = {};
    m_previewBuffers.upload(data, false, noFilter);
    m_previewBufDirty = false;
}

void CSLGrainBoundaryDialog::renderPreviewToFBO(int w, int h)
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
    glClearColor(0.09f, 0.11f, 0.15f, 1.0f);
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

void CSLGrainBoundaryDialog::autoFitPreviewCamera()
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

void CSLGrainBoundaryDialog::feedDroppedFile(const std::string& path)
{
    m_pendingDropPath = path;
}

void CSLGrainBoundaryDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("CSL Grain Boundary", NULL, false, enabled))
        m_openRequested = true;
}

void CSLGrainBoundaryDialog::drawDialog(Structure& structure,
                                        const std::vector<glm::vec3>& elementColors,
                                        const std::vector<float>& elementRadii,
                                        const std::vector<float>& elementShininess,
                                        const std::function<void(Structure&)>& updateBuffers)
{
    (void)elementColors;

// -- Static dialog state --------------------------------------------
    static Structure inputStructure;
    static char statusMsg[256] = "(no structure loaded)";
    static char loadedFileName[256] = "(none)";
    static int axis[3] = {0, 0, 1};
    static int sigmaMax = 200;
    static std::vector<SigmaCandidate> sigmaCandidates;
    static int sigmaSelection = 0;
    static int planeSelection = 0;
    static int lastAxisForSigma[3] = {0, 0, 0};
    static int lastSigmaMaxForSigma = 0;
    static int ucA = 1;
    static int ucB = 1;
    static float vacuumPadding = 0.0f;
    static float gapDist = 0.0f;
    static float overlapDist = 0.0f;
    static bool conventionalCell = false;
    static CSLBuildResult lastResult;

    // -- Handle pending drop ----------------------------------------
    if (!m_pendingDropPath.empty())
    {
        Structure loaded;
        std::string err;
        if (loadStructureFromFile(m_pendingDropPath, loaded, err))
        {
            if (!loaded.hasUnitCell)
            {
                std::snprintf(statusMsg, sizeof(statusMsg), "Error: structure has no unit cell");
                std::snprintf(loadedFileName, sizeof(loadedFileName), "(none)");
                std::cout << "[CSL] Dropped file has no unit cell: " << m_pendingDropPath << std::endl;
            }
            else
            {
                std::string fileName = m_pendingDropPath;
                const size_t slashPos = fileName.find_last_of("/\\");
                if (slashPos != std::string::npos)
                    fileName = fileName.substr(slashPos + 1);

                inputStructure = loaded;
                std::snprintf(statusMsg, sizeof(statusMsg), "Loaded: %d atoms", (int)inputStructure.atoms.size());
                std::snprintf(loadedFileName, sizeof(loadedFileName), "%s", fileName.c_str());
                std::cout << "[CSL] Loaded input structure: " << m_pendingDropPath
                          << " (" << inputStructure.atoms.size() << " atoms)" << std::endl;
                m_previewBufDirty = true;
                if (m_glReady) {
                    rebuildPreviewBuffers(inputStructure, elementRadii, elementShininess);
                    autoFitPreviewCamera();
                }
            }
        }
        else
        {
            std::snprintf(statusMsg, sizeof(statusMsg), "Error: %s", err.c_str());
            std::snprintf(loadedFileName, sizeof(loadedFileName), "(none)");
            std::cout << "[CSL] Drop-load failed: " << m_pendingDropPath
                      << " (" << err << ")" << std::endl;
        }
        m_pendingDropPath.clear();
    }

    // -- Open popup -------------------------------------------------
    if (m_openRequested)
    {
        ImGui::OpenPopup("CSL Grain Boundary Builder");
        m_openRequested = false;
        m_isOpen = true;
    }

    ImGui::SetNextWindowSize(ImVec2(1060.0f, 460.0f), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(1060.0f, 340.0f), ImVec2(FLT_MAX, FLT_MAX));
    bool dialogOpen = true;
    if (ImGui::BeginPopupModal("CSL Grain Boundary Builder", &dialogOpen, ImGuiWindowFlags_None))
    {
        constexpr float kLeftW  = 380.0f;
        const float kContentH = ImGui::GetContentRegionAvail().y;

        // =============== LEFT: Structure preview ===============
        ImGui::BeginChild("##cslLeft", ImVec2(kLeftW, kContentH), true);
        ImGui::Text("Input Structure");
        ImGui::Separator();
        ImGui::Text("Status: %s", statusMsg);
        ImGui::TextDisabled("File: %s", loadedFileName);
        ImGui::Spacing();

        {
            const float prevH = ImGui::GetContentRegionAvail().y - ImGui::GetFrameHeightWithSpacing() * 2.0f;
            ImGui::InvisibleButton("##cslDropZone", ImVec2(-1.0f, std::max(prevH, 80.0f)));
            ImDrawList* dl = ImGui::GetWindowDrawList();
            ImVec2 dropMin = ImGui::GetItemRectMin();
            ImVec2 dropMax = ImGui::GetItemRectMax();
            const bool dropZoneHovered = ImGui::IsItemHovered();
            const bool dropZoneActive  = ImGui::IsItemActive();
            dl->AddRect(dropMin, dropMax, ImGui::GetColorU32(ImGuiCol_Border), 2.0f);

            if (inputStructure.hasUnitCell && !inputStructure.atoms.empty())
            {
                if (m_glReady) {
                    if (m_previewBufDirty)
                        rebuildPreviewBuffers(inputStructure, elementRadii, elementShininess);

                    const float pad = 5.0f;
                    const ImVec2 prevSize(dropMax.x - dropMin.x - 2.0f * pad,
                                          dropMax.y - dropMin.y - 2.0f * pad);
                    const int pw = std::max(1, (int)prevSize.x);
                    const int ph = std::max(1, (int)prevSize.y);

                    renderPreviewToFBO(pw, ph);

                    const ImVec2 prevMin(dropMin.x + pad, dropMin.y + pad);
                    const ImVec2 prevMax(prevMin.x + prevSize.x, prevMin.y + prevSize.y);
                    dl->AddImage((ImTextureID)(intptr_t)m_previewColorTex,
                                 prevMin, prevMax,
                                 ImVec2(0, 1), ImVec2(1, 0));

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
                }
            }
            else
            {
                ImVec2 mid((dropMin.x + dropMax.x) * 0.5f, (dropMin.y + dropMax.y) * 0.5f);
                dl->AddText(ImGui::GetFont(), ImGui::GetFontSize(),
                            ImVec2(mid.x - 100.0f, mid.y - 10.0f),
                            ImGui::GetColorU32(ImGuiCol_TextDisabled),
                            "Drop a structure file here");
            }
        }

        if (inputStructure.hasUnitCell && !inputStructure.atoms.empty())
        {
            const auto& cv = inputStructure.cellVectors;
            double la = std::sqrt(cv[0][0]*cv[0][0] + cv[0][1]*cv[0][1] + cv[0][2]*cv[0][2]);
            double lb = std::sqrt(cv[1][0]*cv[1][0] + cv[1][1]*cv[1][1] + cv[1][2]*cv[1][2]);
            double lc = std::sqrt(cv[2][0]*cv[2][0] + cv[2][1]*cv[2][1] + cv[2][2]*cv[2][2]);
            ImGui::TextDisabled("%d atoms  a=%.2f b=%.2f c=%.2f",
                                (int)inputStructure.atoms.size(), la, lb, lc);
            ImGui::TextDisabled("Drag=orbit  Scroll=zoom");
        }

        ImGui::EndChild(); // ##cslLeft

        ImGui::SameLine();

        // -- Options -----------------------------------------------
        ImGui::BeginChild("##cslRight", ImVec2(0, kContentH), true);
        ImGui::SeparatorText("Misorientation");
        {
            ImGui::SetNextItemWidth(220.0f);
            ImGui::InputInt3("Axis [u v w]", axis);
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Crystallographic rotation axis in Miller indices.");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(110.0f);
            ImGui::InputInt("Max sigma", &sigmaMax);
            if (sigmaMax < 3) sigmaMax = 3;
            ImGui::SameLine();
            if (ImGui::Button("Find"))
            {
                sigmaCandidates = computeGBInfo(axis, sigmaMax);
                sigmaSelection = 0;
                planeSelection = 0;
            }

            // Clear sigma list when axis or max changes
            if (axis[0] != lastAxisForSigma[0] || axis[1] != lastAxisForSigma[1] ||
                axis[2] != lastAxisForSigma[2] || sigmaMax != lastSigmaMaxForSigma)
            {
                sigmaCandidates.clear();
                sigmaSelection = 0;
                planeSelection = 0;
                lastAxisForSigma[0] = axis[0];
                lastAxisForSigma[1] = axis[1];
                lastAxisForSigma[2] = axis[2];
                lastSigmaMaxForSigma = sigmaMax;
            }

            if (!sigmaCandidates.empty())
            {
                auto sigmaGetter = [](void* data, int idx) -> const char* {
                    static char label[96];
                    std::vector<SigmaCandidate>* vec = static_cast<std::vector<SigmaCandidate>*>(data);
                    if (idx < 0 || idx >= (int)vec->size()) return "";
                    const SigmaCandidate& c = (*vec)[idx];
                    std::snprintf(label, sizeof(label), "S%d  m=%d n=%d  %.2f deg",
                                  c.sigma, c.m, c.n, c.thetaDeg);
                    return label;
                };
                if (sigmaSelection >= (int)sigmaCandidates.size())
                    sigmaSelection = 0;
                ImGui::SetNextItemWidth(-1.0f);
                ImGui::Combo("##sigma", &sigmaSelection, sigmaGetter,
                             &sigmaCandidates, (int)sigmaCandidates.size());
            }
            else
            {
                ImGui::TextDisabled("Set axis and press Find.");
            }
        }

        // -- Boundary Geometry --------------------------------------
        ImGui::SeparatorText("Boundary Geometry");
        {
            if (!sigmaCandidates.empty() && sigmaSelection >= 0 && sigmaSelection < (int)sigmaCandidates.size())
            {
                const SigmaCandidate& sel = sigmaCandidates[sigmaSelection];

                // Plane dropdown: 3 planes (columns of CSL)
        // -- Build -------------------------------------------------
                char planeLabels[3][96];
                for (int i = 0; i < 3; i++)
                {
                    std::string type = classifyBoundaryType(axis, sel.plane[i].data());
                    std::snprintf(planeLabels[i], sizeof(planeLabels[i]), "(%d %d %d)  %s",
                                  sel.plane[i][0], sel.plane[i][1], sel.plane[i][2], type.c_str());
                }
                auto planeLabelGetter = [](void* data, int idx) -> const char* {
                    if (idx < 0 || idx >= 3) return "";
                    return static_cast<char(*)[96]>(data)[idx];
                };

                if (planeSelection >= 3) planeSelection = 0;
                ImGui::SetNextItemWidth(260.0f);
                ImGui::Combo("GB plane", &planeSelection, planeLabelGetter, planeLabels, 3);

                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("GB plane normals derived from CSL matrix columns.\n"
                                      "Twist = plane parallel to axis, Tilt = plane perpendicular to axis.");

                ImGui::SameLine();
                ImGui::AlignTextToFramePadding();
                std::string btype = classifyBoundaryType(axis, sel.plane[planeSelection].data());
                ImGui::TextDisabled("Type: %s", btype.c_str());

                // Show CSL matrix
                ImGui::TextDisabled("CSL: [%d %d %d] [%d %d %d] [%d %d %d]",
                                    sel.csl[0][0], sel.csl[0][1], sel.csl[0][2],
                                    sel.csl[1][0], sel.csl[1][1], sel.csl[1][2],
                                    sel.csl[2][0], sel.csl[2][1], sel.csl[2][2]);
            }
            else
            {
                ImGui::TextDisabled("Select a Sigma candidate first.");
            }

            ImGui::SetNextItemWidth(110.0f);
            ImGui::InputInt("Grain A (uc)", &ucA);
            if (ucA < 1) ucA = 1;
            ImGui::SameLine();
            ImGui::SetNextItemWidth(110.0f);
            ImGui::InputInt("Grain B (uc)", &ucB);
            if (ucB < 1) ucB = 1;
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Number of unit cells for each grain along the GB normal.");
        }

        // -- Options -----------------------------------------------
        ImGui::SeparatorText("Options");
        {
            ImGui::SetNextItemWidth(140.0f);
            ImGui::DragFloat("Vacuum (A)", &vacuumPadding, 0.1f, 0.0f, 50.0f, "%.2f");
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Vacuum space at the end of the bicrystal cell.");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(140.0f);
            ImGui::DragFloat("Gap (A)", &gapDist, 0.05f, 0.0f, 10.0f, "%.2f");
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Extra distance between the two grains at the interface.");

            ImGui::SetNextItemWidth(140.0f);
            ImGui::DragFloat("Overlap dist (A)", &overlapDist, 0.05f, 0.0f, 5.0f, "%.3f");
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Atoms closer than this distance will be removed (0 = off).");

            ImGui::Checkbox("Conventional cell", &conventionalCell);
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("If unchecked (default), reduce to the primitive cell.\nIf checked, keep the conventional (unreduced) cell.");
        }

        // -- Build -------------------------------------------------
        ImGui::Spacing();
        const float closeW = 110.0f;
        const float buildW = 110.0f;

        if (ImGui::Button("Build", ImVec2(buildW, 0.0f)))
        {
            if (!inputStructure.hasUnitCell || inputStructure.atoms.empty())
            {
                lastResult.success = false;
                lastResult.message = "Please load a structure with a unit cell first.";
            }
            else if (sigmaCandidates.empty())
            {
                lastResult.success = false;
                lastResult.message = "Generate and select a Sigma candidate first.";
            }
            else
            {
                const SigmaCandidate& sel = sigmaCandidates[sigmaSelection];
                const int* plane = sel.plane[planeSelection].data();
                int direction = planeSelection; // column index = stacking direction

        // -- Build -------------------------------------------------
                Grain inputGrain = structureToGrain(inputStructure);

                // Scaling matrix = CSL^T (transpose of CSL)
                int cslT[3][3];
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        cslT[i][j] = sel.csl[j][i];

                Grain grainA = makeSupercell(inputGrain, cslT);
                lastResult.inputAtoms = (int)inputStructure.atoms.size();

                // Force orthogonal if needed
                if (!isOrthogonal(grainA.cell))
                    grainA = setOrthogonalGrain(grainA, direction);

                // Scale grain A by ucB for creating grain B, then scale grain A by ucA
                int scaleB[3] = {1, 1, 1};
                scaleB[direction] = ucB;
                Grain tempA = makeSupercellDiag(grainA, scaleB[0], scaleB[1], scaleB[2]);

                Grain grainB = getBFromA(tempA);
                // Normalize grain B
                int ones[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
                grainB = makeSupercell(grainB, ones);

                int scaleA[3] = {1, 1, 1};
                scaleA[direction] = ucA;
                grainA = makeSupercellDiag(grainA, scaleA[0], scaleA[1], scaleA[2]);

                // Stack grains
                Grain gb = stackGrains(grainA, grainB, direction, vacuumPadding, gapDist);

                // Remove overlaps
                int removed = removeOverlaps(gb.atoms, overlapDist);

                // Convert to Structure
                structure = grainToStructure(gb);

                // Compute PBC boundary tolerance from minimum layer spacing
                // so that atoms slightly off frac 0.0 still get periodic images.
                {
                    double inv[3][3];
                    invertCell(gb.cell, inv);
                    std::vector<double> layers;
                    layers.reserve(structure.atoms.size());
                    for (const auto& a : structure.atoms)
                    {
                        double cart[3] = {a.x, a.y, a.z};
                        double frac[3];
                        cartToFrac(cart, inv, frac);
                        layers.push_back(wrapFrac(frac[direction]));
                    }
                    std::sort(layers.begin(), layers.end());
                    double minSpacing = 1.0;
                    for (size_t i = 1; i < layers.size(); i++)
                    {
                        double d = layers[i] - layers[i - 1];
                        if (d > 1e-8 && d < minSpacing) minSpacing = d;
                    }
                    structure.pbcBoundaryTol = (float)(minSpacing * 0.5 + 1e-6);
                }

                if (!conventionalCell)
                {
                    reduceToPrimitiveGB(structure, direction);
                }

                structure.grainRegionIds.assign(structure.atoms.size(), 0);
                {
                    double cell[3][3];
                    for (int i = 0; i < 3; ++i)
                        for (int j = 0; j < 3; ++j)
                            cell[i][j] = structure.cellVectors[i][j];

                    double inv[3][3];
                    invertCell(cell, inv);
                    for (size_t i = 0; i < structure.atoms.size(); ++i)
                    {
                        double cart[3] = {structure.atoms[i].x, structure.atoms[i].y, structure.atoms[i].z};
                        double frac[3];
                        cartToFrac(cart, inv, frac);
                        const double layer = wrapFrac(frac[direction]);
                        structure.grainRegionIds[i] = (layer >= 0.5) ? 1 : 0;
                    }
                }

                lastResult.success = true;
                lastResult.sigma = sel.sigma;
                lastResult.thetaDeg = sel.thetaDeg;
                lastResult.generatedAtoms = (int)structure.atoms.size();
                lastResult.removedOverlap = removed;
                lastResult.boundaryType = classifyBoundaryType(axis, plane);
                std::memcpy(lastResult.axis, axis, sizeof(lastResult.axis));
                lastResult.gbPlane[0] = plane[0];
                lastResult.gbPlane[1] = plane[1];
                lastResult.gbPlane[2] = plane[2];
                std::memcpy(lastResult.csl, sel.csl, sizeof(lastResult.csl));
                lastResult.ucA = ucA;
                lastResult.ucB = ucB;
                lastResult.vacuumPadding = vacuumPadding;
                lastResult.gap = gapDist;
                lastResult.overlapDist = overlapDist;

                lastResult.message = "CSL grain boundary generated.";

                updateBuffers(structure);

                std::cout << "[Operation] Built CSL grain boundary: "
                          << "Sigma=" << sel.sigma
                          << ", theta=" << sel.thetaDeg << " deg"
                          << ", plane=(" << plane[0] << " " << plane[1] << " " << plane[2] << ")"
                          << ", type=" << lastResult.boundaryType
                          << ", atoms=" << lastResult.generatedAtoms
                          << std::endl;
            }

            if (!lastResult.success)
            {
                std::cout << "[Operation] CSL grain boundary generation failed: "
                          << lastResult.message << std::endl;
            }
        }

        ImGui::SameLine();
        if (ImGui::Button("Close", ImVec2(closeW, 0.0f)))
        {
            dialogOpen = false;
            ImGui::CloseCurrentPopup();
        }

        // -- Result ------------------------------------------------
        ImGui::Spacing();
        if (!lastResult.message.empty())
        {
            ImGui::SeparatorText("Result");
            drawBuildResultSummary(lastResult);
        }

        ImGui::EndChild(); // ##cslRight

        ImGui::EndPopup();
    }

    if (!dialogOpen)
    {
        m_isOpen = false;
    }
}
