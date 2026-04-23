#include "ui/InterfaceBuilderDialog.h"

#include "algorithms/InterfaceBuilder.h"
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
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// Implemented from 
// Stradi, Daniele, et al. "Method for determining optimal supercell representation of interfaces."
// Journal of Physics: Condensed Matter 29.18 (2017): 185901.
// ============================================================================

// ============================================================================
// InterfaceBuilderDialog implementation
// ============================================================================

InterfaceBuilderDialog::InterfaceBuilderDialog()
{
    initDemoStiffness();
}

InterfaceBuilderDialog::~InterfaceBuilderDialog()
{
    delete m_previewSphereA;
    delete m_previewCylinderA;
    delete m_previewSphereB;
    delete m_previewCylinderB;
    delete m_previewSphereR;
    delete m_previewCylinderR;
    if (m_fboA)      glDeleteFramebuffers(1, &m_fboA);
    if (m_texA)      glDeleteTextures(1, &m_texA);
    if (m_rboA)      glDeleteRenderbuffers(1, &m_rboA);
    if (m_fboB)      glDeleteFramebuffers(1, &m_fboB);
    if (m_texB)      glDeleteTextures(1, &m_texB);
    if (m_rboB)      glDeleteRenderbuffers(1, &m_rboB);
    if (m_fboResult) glDeleteFramebuffers(1, &m_fboResult);
    if (m_texResult) glDeleteTextures(1, &m_texResult);
    if (m_rboResult) glDeleteRenderbuffers(1, &m_rboResult);
}

void InterfaceBuilderDialog::initDemoStiffness()
{
    // Default stiffness tensor in GPa
    std::memset(m_stiffness, 0, sizeof(m_stiffness));
    m_stiffness[0][0] = 198.0f; m_stiffness[0][1] =  81.0f; m_stiffness[0][2] =  53.0f;
    m_stiffness[1][0] =  81.0f; m_stiffness[1][1] = 198.0f; m_stiffness[1][2] =  53.0f;
    m_stiffness[2][0] =  53.0f; m_stiffness[2][1] =  53.0f; m_stiffness[2][2] = 246.0f;
    m_stiffness[3][3] =  54.0f; m_stiffness[4][4] =  54.0f; m_stiffness[5][5] =  59.0f;
}

void InterfaceBuilderDialog::initRenderResources(Renderer& renderer)
{
    m_renderer = &renderer;
    // Each SceneBuffers needs its own VAO pair because init() wires instance
    // VBOs into the VAO; sharing a VAO means only the last init'd buffer works.
    m_previewSphereA    = new SphereMesh(24, 24);
    m_previewCylinderA  = new CylinderMesh(16);
    m_previewSphereB    = new SphereMesh(24, 24);
    m_previewCylinderB  = new CylinderMesh(16);
    m_previewSphereR    = new SphereMesh(24, 24);
    m_previewCylinderR  = new CylinderMesh(16);
    m_previewBufA.init(m_previewSphereA->vbo, m_previewSphereA->ebo, m_previewSphereA->indexCount,
                       m_previewCylinderA->vbo, m_previewCylinderA->vertexCount);
    m_previewBufB.init(m_previewSphereB->vbo, m_previewSphereB->ebo, m_previewSphereB->indexCount,
                       m_previewCylinderB->vbo, m_previewCylinderB->vertexCount);
    m_previewBufResult.init(m_previewSphereR->vbo, m_previewSphereR->ebo, m_previewSphereR->indexCount,
                            m_previewCylinderR->vbo, m_previewCylinderR->vertexCount);
    m_previewShadow = createShadowMap(1, 1);
    m_glReady = true;
}

void InterfaceBuilderDialog::feedDroppedFile(const std::string& path, int slotHint)
{
    if (slotHint == 1)
        m_pendingDropPathA = path;
    else if (slotHint == 2)
        m_pendingDropPathB = path;
    else
    {
        // Auto: fill A first, then B
        if (m_structureA.atoms.empty())
            m_pendingDropPathA = path;
        else
            m_pendingDropPathB = path;
    }
}

void InterfaceBuilderDialog::ensureFBO(GLuint& fbo, GLuint& tex, GLuint& rbo,
                                        int& curW, int& curH, int w, int h)
{
    if (w == curW && h == curH && fbo != 0)
        return;

    if (fbo) { glDeleteFramebuffers(1, &fbo);    fbo = 0; }
    if (tex) { glDeleteTextures(1, &tex);         tex = 0; }
    if (rbo) { glDeleteRenderbuffers(1, &rbo);    rbo = 0; }

    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenRenderbuffers(1, &rbo);
    glBindRenderbuffer(GL_RENDERBUFFER, rbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rbo);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    curW = w;
    curH = h;
}

void InterfaceBuilderDialog::rebuildPreviewBuffers(
    SceneBuffers& buf, const Structure& src,
    const std::vector<float>& radii,
    const std::vector<float>& shininess)
{
    if (!m_glReady || src.atoms.empty())
        return;
    static const int kIdent[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    StructureInstanceData data = buildStructureInstanceData(src, false, kIdent, radii, shininess);
    std::array<bool, 119> noFilter = {};
    buf.upload(data, false, noFilter);
}

void InterfaceBuilderDialog::autoFitCamera(const SceneBuffers& buf, float& dist)
{
    if (buf.atomCount == 0) { dist = 10.0f; return; }
    float maxR = 0.0f;
    for (size_t i = 0; i < buf.atomPositions.size(); ++i)
    {
        float r = (i < buf.atomRadii.size()) ? buf.atomRadii[i] : 0.0f;
        float d = glm::length(buf.atomPositions[i] - buf.orbitCenter) + r;
        maxR = std::max(maxR, d);
    }
    maxR = std::max(maxR, 1.0f);
    float halfFov = glm::radians(22.5f);
    float newDist = maxR / std::sin(halfFov) * 1.15f;
    dist = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, newDist));
}

void InterfaceBuilderDialog::renderToFBO(GLuint fbo, int w, int h,
                                          SceneBuffers& buf,
                                          float yaw, float pitch, float dist)
{
    if (!m_glReady || !m_renderer || buf.atomCount == 0)
        return;

    GLint prevFbo = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
    GLint prevVP[4];
    glGetIntegerv(GL_VIEWPORT, prevVP);

    Camera cam;
    cam.yaw = yaw;
    cam.pitch = pitch;
    cam.distance = dist;

    FrameView frame;
    frame.framebufferWidth = w;
    frame.framebufferHeight = h;
    buildFrameView(cam, buf, true, frame);

    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glViewport(0, 0, w, h);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.09f, 0.11f, 0.15f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    m_renderer->drawBonds(frame.projection, frame.view,
                          frame.lightPosition, frame.cameraPosition,
                          buf.tabCylinderVAO, buf.tabCylinderVertexCount, buf.bondCount);
    m_renderer->drawAtoms(frame.projection, frame.view,
                          frame.lightMVP, frame.lightPosition, frame.cameraPosition,
                          m_previewShadow,
                          buf.tabSphereVAO, buf.tabSphereIndexCount, buf.atomCount);
    m_renderer->drawBoxLines(frame.projection, frame.view,
                             buf.lineVAO, buf.boxLines.size());

    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
    glViewport(prevVP[0], prevVP[1], prevVP[2], prevVP[3]);
}

bool InterfaceBuilderDialog::tryLoadFile(
    const std::string& path, Structure& dest, std::string& status,
    SceneBuffers& buf,
    const std::vector<float>& radii,
    const std::vector<float>& shininess,
    float& camDist)
{
    Structure loaded;
    std::string err;
    if (!loadStructureFromFile(path, loaded, err))
    {
        char msg[512];
        std::snprintf(msg, sizeof(msg), "Load failed: %s", err.c_str());
        status = msg;
        return false;
    }
    dest = std::move(loaded);
    // Extract filename from path for display
    std::string fname = path;
    auto slashPos = fname.find_last_of("/\\");
    if (slashPos != std::string::npos) fname = fname.substr(slashPos + 1);
    char msg[512];
    std::snprintf(msg, sizeof(msg), "%s  (%d atoms)", fname.c_str(), (int)dest.atoms.size());
    status = msg;

    if (m_glReady)
    {
        rebuildPreviewBuffers(buf, dest, radii, shininess);
        autoFitCamera(buf, camDist);
    }
    return true;
}

// ============================================================================
// runSearch - the core matching algorithm
// ============================================================================
void InterfaceBuilderDialog::runSearch()
{
    m_candidates.clear();
    m_interfaceStructures.clear();
    m_searchDone = false;
    m_selectedIdx = -1;
    m_lastRenderedResultIdx = -2;

    if (!m_structureA.hasUnitCell || !m_structureB.hasUnitCell)
    {
        std::cout << "[InterfaceBuilder] Both structures must have unit cells." << std::endl;
        return;
    }

    double basisA[2][2], basisB[2][2];
    get2DBasis(m_structureA, basisA);
    get2DBasis(m_structureB, basisB);

    auto superA = generateUniqueSupercells(basisA, m_nmax, m_maxCellsA);
    auto superB = generateUniqueSupercells(basisB, m_mmax, m_maxCellsB);

    std::vector<float> orAngles;
    if (m_useOrPlaneDir && m_structureA.hasUnitCell && m_structureB.hasUnitCell)
    {
        float phiA = orientationAngleFromPlaneDir(m_structureA, m_planeA, m_dirA);
        float phiB = orientationAngleFromPlaneDir(m_structureB, m_planeB, m_dirB);
        float orFromPlaneDir = phiA - phiB;
        orAngles.push_back(orFromPlaneDir);
    }
    else if (m_useOrAngle)
    {
        orAngles.push_back(m_orAngleDeg);
    }

    int idx = 0;
    for (const auto& sa : superA)
    {
        for (const auto& sb : superB)
        {
            double phiV1 = angleOf(sa.vecs[0][0], sa.vecs[0][1]);
            double phiU1 = angleOf(sb.vecs[0][0], sb.vecs[0][1]);
            double thetaDeg = wrapDegPm180(phiV1 * 180.0 / M_PI - phiU1 * 180.0 / M_PI);
            double theta = thetaDeg * M_PI / 180.0;

            if (std::abs(thetaDeg) > m_maxRotationDeg)
                continue;

            // OR check
            float orTargetDeg = 0.0f, orMisfitDeg = 0.0f;
            if (!orAngles.empty())
            {
                double bestTarget = orAngles[0];
                double bestMisfit = std::abs(wrapDegPm180(thetaDeg - bestTarget));
                for (size_t k = 1; k < orAngles.size(); ++k)
                {
                    double misfit = std::abs(wrapDegPm180(thetaDeg - orAngles[k]));
                    if (misfit < bestMisfit)
                    {
                        bestMisfit = misfit;
                        bestTarget = orAngles[k];
                    }
                }
                orTargetDeg = (float)bestTarget;
                orMisfitDeg = (float)bestMisfit;
                if (orMisfitDeg > m_orTolDeg)
                    continue;
            }

            // Rotate u by theta
            Mat2 rotTheta = rotation2D(theta);
            double uRot[2][2];
            for (int i = 0; i < 2; ++i)
            {
                uRot[i][0] = sb.vecs[i][0]*rotTheta.m[0][0] + sb.vecs[i][1]*rotTheta.m[1][0];
                uRot[i][1] = sb.vecs[i][0]*rotTheta.m[0][1] + sb.vecs[i][1]*rotTheta.m[1][1];
            }

            double exx, eyy, exy;
            if (!strainComponents(sa.vecs, uRot, exx, eyy, exy))
                continue;

            double ebar = meanAbsStrain(exx, eyy, exy);
            if (ebar > m_maxMeanStrain)
                continue;

            double elastic = cubicElasticDensity(exx, eyy, exy, m_stiffness);

            // Compute exact 2x2 deformation gradient F that maps B supercell
            // vectors (sb.vecs) to A supercell vectors (sa.vecs).
            // F @ u_col_i = v_col_i, so F = V_col @ inv(U_col)
            double detU = sb.vecs[0][0]*sb.vecs[1][1] - sb.vecs[1][0]*sb.vecs[0][1];
            double F2[2][2];
            F2[0][0] = ( sa.vecs[0][0]*sb.vecs[1][1] - sa.vecs[1][0]*sb.vecs[0][1]) / detU;
            F2[0][1] = (-sa.vecs[0][0]*sb.vecs[1][0] + sa.vecs[1][0]*sb.vecs[0][0]) / detU;
            F2[1][0] = ( sa.vecs[0][1]*sb.vecs[1][1] - sa.vecs[1][1]*sb.vecs[0][1]) / detU;
            F2[1][1] = (-sa.vecs[0][1]*sb.vecs[1][0] + sa.vecs[1][1]*sb.vecs[0][0]) / detU;

            // Build actual interface structure
            Structure aSuper2D = makeSupercell2D(m_structureA, sa.mat);
            Structure aSuper = repeatLayersZ(aSuper2D, m_layersA);
            Structure bSuper2D = makeSupercell2D(m_structureB, sb.mat);
            Structure bSuper = repeatLayersZ(bSuper2D, m_layersB);
            Structure bMatched = applyTransform2D(bSuper, F2);
            Structure iface = assembleInterface(aSuper, bMatched, m_zGap, m_vacuum);
            if (m_repeatX > 1 || m_repeatY > 1)
                iface = repeatInterfaceXY(iface, m_repeatX, m_repeatY);

            InterfaceCandidate cand;
            cand.idx = idx;
            std::memcpy(cand.nMatrix, sa.mat, sizeof(sa.mat));
            std::memcpy(cand.mMatrix, sb.mat, sizeof(sb.mat));
            cand.detN = sa.det;
            cand.detM = sb.det;
            cand.thetaDeg = (float)thetaDeg;
            cand.orTargetDeg = orTargetDeg;
            cand.orMisfitDeg = orMisfitDeg;
            cand.exx = (float)exx;
            cand.eyy = (float)eyy;
            cand.exy = (float)exy;
            cand.meanAbsStrain = (float)ebar;
            cand.elasticDensity = (float)elastic;
            cand.interfaceAtoms = (int)iface.atoms.size();

            m_candidates.push_back(cand);
            m_interfaceStructures.push_back(std::move(iface));
            ++idx;
        }
    }

    m_searchDone = true;
    std::cout << "[InterfaceBuilder] Found " << m_candidates.size()
              << " candidate interfaces." << std::endl;
}

void InterfaceBuilderDialog::buildInterfaceStructure(int candidateIdx)
{
    // Already built during search - just select it for preview
    m_selectedIdx = candidateIdx;
}

// ============================================================================
// Draw helpers
// ============================================================================

void InterfaceBuilderDialog::drawDropZone(
    const char* label, Structure& struc, const std::string& status,
    GLuint tex, SceneBuffers& buf,
    float& yaw, float& pitch, float& dist,
    const std::vector<float>& radii,
    const std::vector<float>& shininess,
    float width, float height, int slot)
{
    (void)tex;

    ImGui::BeginChild(label, ImVec2(width, height), true);
    ImGui::TextColored(ImVec4(0.7f, 0.85f, 1.0f, 1.0f), "%s", label);
    ImGui::Separator();
    if (!status.empty())
        ImGui::Text("%s", status.c_str());
    else
        ImGui::TextDisabled("No file loaded");
    ImGui::Spacing();

    float prevH = height - 80.0f;
    if (prevH < 50.0f) prevH = 50.0f;

    ImGui::InvisibleButton(
        (std::string("##dropZone") + label).c_str(),
        ImVec2(-1.0f, prevH));

    ImDrawList* dl = ImGui::GetWindowDrawList();
    ImVec2 dropMin = ImGui::GetItemRectMin();
    ImVec2 dropMax = ImGui::GetItemRectMax();
    bool zoneHovered = ImGui::IsItemHovered();
    bool zoneActive  = ImGui::IsItemActive();
    dl->AddRect(dropMin, dropMax, ImGui::GetColorU32(ImGuiCol_Border), 2.0f);

    if (struc.atoms.empty())
    {
        // Accept drag-and-drop
        if (ImGui::BeginDragDropTarget())
        {
            if (const ImGuiPayload* payload =
                    ImGui::AcceptDragDropPayload("DROP_FILE", ImGuiDragDropFlags_AcceptPeekOnly))
            {
                const char* files = (const char*)payload->Data;
                std::string statusOut;
                if (slot == 0)
                    tryLoadFile(files, m_structureA, m_statusA, m_previewBufA, radii, shininess, m_camADist);
                else
                    tryLoadFile(files, m_structureB, m_statusB, m_previewBufB, radii, shininess, m_camBDist);
            }
            ImGui::EndDragDropTarget();
        }

        ImVec2 mid((dropMin.x + dropMax.x) * 0.5f, (dropMin.y + dropMax.y) * 0.5f);
        dl->AddText(ImGui::GetFont(), ImGui::GetFontSize(),
                    ImVec2(mid.x - 90.0f, mid.y - 10.0f),
                    ImGui::GetColorU32(ImGuiCol_TextDisabled),
                    "Drop a structure file here");
    }
    else if (m_glReady)
    {
        float pad = 5.0f;
        ImVec2 prevSize(dropMax.x - dropMin.x - 2*pad, dropMax.y - dropMin.y - 2*pad);
        int pw = std::max(1, (int)prevSize.x);
        int ph = std::max(1, (int)prevSize.y);

        GLuint* fboPtr;
        GLuint* texPtr;
        GLuint* rboPtr;
        int *wPtr, *hPtr;
        if (slot == 0) { fboPtr = &m_fboA; texPtr = &m_texA; rboPtr = &m_rboA; wPtr = &m_fboAW; hPtr = &m_fboAH; }
        else           { fboPtr = &m_fboB; texPtr = &m_texB; rboPtr = &m_rboB; wPtr = &m_fboBW; hPtr = &m_fboBH; }

        ensureFBO(*fboPtr, *texPtr, *rboPtr, *wPtr, *hPtr, pw, ph);
        renderToFBO(*fboPtr, pw, ph, buf, yaw, pitch, dist);

        ImVec2 prevMin(dropMin.x + pad, dropMin.y + pad);
        ImVec2 prevMax(prevMin.x + prevSize.x, prevMin.y + prevSize.y);
        dl->AddImage((ImTextureID)(intptr_t)*texPtr, prevMin, prevMax,
                     ImVec2(0, 1), ImVec2(1, 0));

        // Accept new drag-drop even when loaded
        if (ImGui::BeginDragDropTarget())
        {
            if (const ImGuiPayload* payload =
                    ImGui::AcceptDragDropPayload("DROP_FILE", ImGuiDragDropFlags_AcceptPeekOnly))
            {
                const char* files = (const char*)payload->Data;
                if (slot == 0)
                    tryLoadFile(files, m_structureA, m_statusA, m_previewBufA, radii, shininess, m_camADist);
                else
                    tryLoadFile(files, m_structureB, m_statusB, m_previewBufB, radii, shininess, m_camBDist);
            }
            ImGui::EndDragDropTarget();
        }

        // Orbit interaction
        if (zoneActive && ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f))
        {
            ImVec2 delta = ImGui::GetIO().MouseDelta;
            yaw   -= delta.x * 0.5f;
            pitch += delta.y * 0.5f;
        }
        if (zoneHovered)
        {
            float wheel = ImGui::GetIO().MouseWheel;
            if (wheel != 0.0f)
            {
                dist -= wheel * dist * 0.1f;
                dist = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, dist));
            }
        }

        ImGui::TextDisabled("Left-drag = orbit   Scroll = zoom");
    }

    if (!struc.atoms.empty())
    {
        ImGui::Spacing();
        ImGui::Text("Atoms: %d  Cell: %s",
                    (int)struc.atoms.size(),
                    struc.hasUnitCell ? "yes" : "no");
    }

    ImGui::EndChild();
}

// ============================================================================
// 2D scatter plot (ImGui-drawn, no matplotlib)
// ============================================================================
void InterfaceBuilderDialog::draw2DPlot(float width, float height)
{
    if (m_candidates.empty())
    {
        ImGui::TextDisabled("No candidates. Run search first.");
        return;
    }

    ImGui::BeginChild("##ifacePlot", ImVec2(width, height), true);

    ImDrawList* dl = ImGui::GetWindowDrawList();

    // Margins for axis labels (left for Y, bottom for X)
    float leftMargin  = 48.0f;
    float bottomMargin = 28.0f;
    float topPad  = 4.0f;
    float rightPad = 10.0f;

    ImVec2 origin = ImGui::GetCursorScreenPos();
    ImVec2 canvasMin(origin.x + leftMargin, origin.y + topPad);
    ImVec2 canvasSize(width - 16.0f - leftMargin - rightPad,
                      height - 40.0f - bottomMargin);
    if (canvasSize.x < 80.0f) canvasSize.x = 80.0f;
    if (canvasSize.y < 60.0f) canvasSize.y = 60.0f;
    ImVec2 canvasMax(canvasMin.x + canvasSize.x, canvasMin.y + canvasSize.y);

    // Data ranges
    float xMin = 1e30f, xMax = -1e30f;
    float yMin = 1e30f, yMax = -1e30f;
    for (const auto& c : m_candidates)
    {
        float x = c.meanAbsStrain * 100.0f;
        float y = (float)c.interfaceAtoms;
        xMin = std::min(xMin, x); xMax = std::max(xMax, x);
        yMin = std::min(yMin, y); yMax = std::max(yMax, y);
    }
    float xPad = (xMax - xMin) * 0.05f + 0.01f;
    float yPad = (yMax - yMin) * 0.05f + 1.0f;
    xMin -= xPad; xMax += xPad;
    yMin -= yPad; yMax += yPad;
    if (xMax <= xMin) xMax = xMin + 1.0f;
    if (yMax <= yMin) yMax = yMin + 10.0f;

    auto mapX = [&](float v) -> float {
        return canvasMin.x + (v - xMin) / (xMax - xMin) * canvasSize.x;
    };
    auto mapY = [&](float v) -> float {
        return canvasMax.y - (v - yMin) / (yMax - yMin) * canvasSize.y;
    };

    // Background
    dl->AddRectFilled(canvasMin, canvasMax, IM_COL32(20, 25, 35, 255));

    // Grid lines and tick labels
    int nTicks = 5;
    for (int i = 0; i <= nTicks; ++i)
    {
        float t = (float)i / (float)nTicks;
        float gx = canvasMin.x + t * canvasSize.x;
        float gy = canvasMin.y + t * canvasSize.y;

        // Vertical grid line
        dl->AddLine(ImVec2(gx, canvasMin.y), ImVec2(gx, canvasMax.y),
                    IM_COL32(60, 60, 60, 80));
        // Horizontal grid line
        dl->AddLine(ImVec2(canvasMin.x, gy), ImVec2(canvasMax.x, gy),
                    IM_COL32(60, 60, 60, 80));

        char buf[64];

        // X-axis tick labels (bottom, centered)
        std::snprintf(buf, sizeof(buf), "%.2f", xMin + t * (xMax - xMin));
        ImVec2 txSz = ImGui::CalcTextSize(buf);
        dl->AddText(ImVec2(gx - txSz.x * 0.5f, canvasMax.y + 4.0f),
                    IM_COL32(180, 180, 180, 255), buf);

        // Y-axis tick labels (left, right-aligned)
        std::snprintf(buf, sizeof(buf), "%.0f", yMax - t * (yMax - yMin));
        txSz = ImGui::CalcTextSize(buf);
        dl->AddText(ImVec2(canvasMin.x - txSz.x - 6.0f, gy - txSz.y * 0.5f),
                    IM_COL32(180, 180, 180, 255), buf);
    }

    // Solid axis lines on left and bottom edges
    dl->AddLine(ImVec2(canvasMin.x, canvasMin.y),
                ImVec2(canvasMin.x, canvasMax.y),
                IM_COL32(180, 180, 180, 255), 2.0f);
    dl->AddLine(ImVec2(canvasMin.x, canvasMax.y),
                ImVec2(canvasMax.x, canvasMax.y),
                IM_COL32(180, 180, 180, 255), 2.0f);

    // Plot points
    ImVec2 mousePos = ImGui::GetIO().MousePos;
    m_hoveredCandidate = -1;

    for (size_t i = 0; i < m_candidates.size(); ++i)
    {
        const auto& c = m_candidates[i];
        float px = mapX(c.meanAbsStrain * 100.0f);
        float py = mapY((float)c.interfaceAtoms);
        float radius = 4.0f;

        // Color by elastic density (green = low, red = high)
        float eMax = 0.001f;
        for (const auto& cc : m_candidates)
            eMax = std::max(eMax, std::abs(cc.elasticDensity));
        float t = std::min(1.0f, std::abs(c.elasticDensity) / eMax);
        ImU32 col = IM_COL32(
            (int)(50 + 200 * t),
            (int)(200 - 150 * t),
            50,
            220);

        // Highlight selected
        if ((int)i == m_selectedIdx)
        {
            dl->AddCircleFilled(ImVec2(px, py), radius + 4.0f, IM_COL32(255, 50, 50, 180));
            col = IM_COL32(255, 255, 0, 255);
        }

        dl->AddCircleFilled(ImVec2(px, py), radius, col);

        // Hover detection
        float dx = mousePos.x - px;
        float dy = mousePos.y - py;
        if (dx*dx + dy*dy < (radius + 6.0f) * (radius + 6.0f))
            m_hoveredCandidate = (int)i;
    }

    // Tooltip on hover
    if (m_hoveredCandidate >= 0 && m_hoveredCandidate < (int)m_candidates.size())
    {
        const auto& c = m_candidates[m_hoveredCandidate];
        ImGui::BeginTooltip();
        ImGui::Text("Candidate #%d", c.idx);
        ImGui::Text("Atoms: %d", c.interfaceAtoms);
        ImGui::Text("Mean strain: %.3f%%", c.meanAbsStrain * 100.0f);
        ImGui::Text("exx=%.3f%%  eyy=%.3f%%  exy=%.3f%%",
                    c.exx*100.0f, c.eyy*100.0f, c.exy*100.0f);
        ImGui::Text("Theta: %.2f deg", c.thetaDeg);
        ImGui::Text("OR: (%.0f %.0f %.0f)A || (%.0f %.0f %.0f)B,  [%d %d 0]A || [%d %d 0]B",
                    m_planeA[0], m_planeA[1], m_planeA[2],
                    m_planeB[0], m_planeB[1], m_planeB[2],
                    c.nMatrix[0][0], c.nMatrix[0][1],
                    c.mMatrix[0][0], c.mMatrix[0][1]);
        ImGui::Text("Elastic density: %.5f", c.elasticDensity);
        ImGui::Text("Click to select");
        ImGui::EndTooltip();

        // Click to select
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left))
        {
            buildInterfaceStructure(m_hoveredCandidate);
        }
    }

    // Axis title
    ImGui::SetCursorPosY(ImGui::GetCursorPosY() + canvasSize.y + bottomMargin + topPad + 4.0f);
    ImGui::TextDisabled("X: Mean strain (%%)   Y: Atoms   Color: Elastic density");

    ImGui::EndChild();
}

// ============================================================================
// drawMenuItem
// ============================================================================
void InterfaceBuilderDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Interface Builder", NULL, false, enabled))
        m_openRequested = true;
}

// ============================================================================
// drawDialog - main ImGui popup
// ============================================================================
void InterfaceBuilderDialog::drawDialog(
    Structure& structure,
    const std::vector<glm::vec3>& elementColors,
    const std::vector<float>& elementRadii,
    const std::vector<float>& elementShininess,
    const std::function<void(Structure&)>& updateBuffers)
{
    (void)elementColors;

    // Consume pending file drops
    if (!m_pendingDropPathA.empty())
    {
        tryLoadFile(m_pendingDropPathA, m_structureA, m_statusA,
                    m_previewBufA, elementRadii, elementShininess, m_camADist);
        m_pendingDropPathA.clear();
        m_searchDone = false;
    }
    if (!m_pendingDropPathB.empty())
    {
        tryLoadFile(m_pendingDropPathB, m_structureB, m_statusB,
                    m_previewBufB, elementRadii, elementShininess, m_camBDist);
        m_pendingDropPathB.clear();
        m_searchDone = false;
    }

    if (m_openRequested)
    {
        ImGui::OpenPopup("Interface Builder");
        m_openRequested = false;
    }

    m_isOpen = ImGui::IsPopupOpen("Interface Builder");

    // Size the dialog to almost the full application window
    ImVec2 displaySize = ImGui::GetIO().DisplaySize;
    ImGui::SetNextWindowSize(ImVec2(displaySize.x * 0.95f, displaySize.y * 0.95f), ImGuiCond_Always);
    ImGui::SetNextWindowPos(ImVec2(displaySize.x * 0.025f, displaySize.y * 0.025f), ImGuiCond_Always);

    bool dialogOpen = true;
    if (!ImGui::BeginPopupModal("Interface Builder", &dialogOpen,
                                 ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize))
    {
        if (!dialogOpen) m_isOpen = false;
        return;
    }

    m_isOpen = true;

    float totalW = ImGui::GetContentRegionAvail().x;
    float gap = 4.0f;

    // =========================================================================
    // TOP ROW: Struct A | Struct B | OR | Params | Stiffness
    // =========================================================================
    float topH = 300.0f;
    float structW = totalW * 0.16f;
    float orW     = totalW * 0.17f;
    float paramW  = totalW * 0.26f;
    float stiffW  = totalW - 2.0f * structW - orW - paramW - 4.0f * gap;

    // --- Structure A ---
    drawDropZone("Structure A", m_structureA, m_statusA,
                 m_texA, m_previewBufA,
                 m_camAYaw, m_camAPitch, m_camADist,
                 elementRadii, elementShininess,
                 structW, topH, 0);

    ImGui::SameLine(0.0f, gap);

    // --- Structure B ---
    drawDropZone("Structure B", m_structureB, m_statusB,
                 m_texB, m_previewBufB,
                 m_camBYaw, m_camBPitch, m_camBDist,
                 elementRadii, elementShininess,
                 structW, topH, 1);

    ImGui::SameLine(0.0f, gap);

    // --- Orientation Relationship ---
    ImGui::BeginChild("##orPanel", ImVec2(orW, topH), true);
    {
        ImGui::TextColored(ImVec4(0.7f, 0.85f, 1.0f, 1.0f), "Orientation Relationship");
        ImGui::Separator();

        if (ImGui::RadioButton("None##or", !m_useOrAngle && !m_useOrPlaneDir))
        { m_useOrAngle = false; m_useOrPlaneDir = false; }
        if (ImGui::RadioButton("Manual angle##or", m_useOrAngle && !m_useOrPlaneDir))
        { m_useOrAngle = true; m_useOrPlaneDir = false; }
        if (ImGui::RadioButton("Plane + direction##or", m_useOrPlaneDir))
        { m_useOrPlaneDir = true; m_useOrAngle = false; }

        ImGui::Spacing();
        ImGui::AlignTextToFramePadding();
        ImGui::Text("Tolerance"); ImGui::SameLine();
        ImGui::SetNextItemWidth(50); ImGui::InputFloat("##ortol", &m_orTolDeg, 0, 0, "%.2f");
        if (m_orTolDeg < 0.0f) m_orTolDeg = 0.0f;
        ImGui::SameLine(); ImGui::TextDisabled("deg");

        if (m_useOrAngle)
        {
            ImGui::AlignTextToFramePadding();
            ImGui::Text("Angle"); ImGui::SameLine();
            ImGui::SetNextItemWidth(60);
            ImGui::InputFloat("##orangle", &m_orAngleDeg, 0, 0, "%.2f");
            ImGui::SameLine(); ImGui::TextDisabled("deg");
        }

        if (m_useOrPlaneDir)
        {
            ImGui::Spacing();
            float fw = 35.0f;
            if (ImGui::BeginTable("##orTbl", 5, ImGuiTableFlags_SizingFixedFit))
            {
                ImGui::TableSetupColumn("id",   0, 14.0f);
                ImGui::TableSetupColumn("type", 0, 72.0f);
                ImGui::TableSetupColumn("v1",   0, fw + 4);
                ImGui::TableSetupColumn("v2",   0, fw + 4);
                ImGui::TableSetupColumn("v3",   0, fw + 4);

                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("A");
                ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("Plane (hkl)");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fw); ImGui::InputFloat("##hA", &m_planeA[0], 0, 0, "%.0f");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fw); ImGui::InputFloat("##kA", &m_planeA[1], 0, 0, "%.0f");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fw); ImGui::InputFloat("##lA", &m_planeA[2], 0, 0, "%.0f");

                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("Dir [uvw]");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fw); ImGui::InputFloat("##uA", &m_dirA[0], 0, 0, "%.0f");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fw); ImGui::InputFloat("##vA", &m_dirA[1], 0, 0, "%.0f");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fw); ImGui::InputFloat("##wA", &m_dirA[2], 0, 0, "%.0f");

                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("B");
                ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("Plane (hkl)");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fw); ImGui::InputFloat("##hB", &m_planeB[0], 0, 0, "%.0f");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fw); ImGui::InputFloat("##kB", &m_planeB[1], 0, 0, "%.0f");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fw); ImGui::InputFloat("##lB", &m_planeB[2], 0, 0, "%.0f");

                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("Dir [uvw]");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fw); ImGui::InputFloat("##uB", &m_dirB[0], 0, 0, "%.0f");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fw); ImGui::InputFloat("##vB", &m_dirB[1], 0, 0, "%.0f");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fw); ImGui::InputFloat("##wB", &m_dirB[2], 0, 0, "%.0f");

                ImGui::EndTable();
            }

            if (m_structureA.hasUnitCell && m_structureB.hasUnitCell)
            {
                float phiA = orientationAngleFromPlaneDir(m_structureA, m_planeA, m_dirA);
                float phiB = orientationAngleFromPlaneDir(m_structureB, m_planeB, m_dirB);
                ImGui::TextDisabled("OR: %.2f deg", phiA - phiB);
            }
        }
    }
    ImGui::EndChild();

    ImGui::SameLine(0.0f, gap);

    // --- Matching Parameters ---
    ImGui::BeginChild("##paramsPanel", ImVec2(paramW, topH), true);
    {
        ImGui::TextColored(ImVec4(0.7f, 0.85f, 1.0f, 1.0f), "Matching Parameters");
        ImGui::Separator();

        float fieldW = 55.0f;
        float colLblW = 85.0f;

        ImGui::TextColored(ImVec4(0.6f, 0.75f, 0.9f, 1.0f), "Supercell");
        if (ImGui::BeginTable("##scTbl", 4, ImGuiTableFlags_SizingFixedFit))
        {
            ImGui::TableSetupColumn("L1", 0, colLblW);
            ImGui::TableSetupColumn("V1", 0, fieldW + 8);
            ImGui::TableSetupColumn("L2", 0, colLblW);
            ImGui::TableSetupColumn("V2", 0, fieldW + 8);

            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("N max");
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fieldW); ImGui::InputInt("##nmax", &m_nmax, 0, 0);
            if (m_nmax < 1) m_nmax = 1;
            if (m_nmax > 10) m_nmax = 10;
            ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("M max");
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fieldW); ImGui::InputInt("##mmax", &m_mmax, 0, 0);
            if (m_mmax < 1) m_mmax = 1;
            if (m_mmax > 10) m_mmax = 10;

            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("Max cells A");
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fieldW); ImGui::InputInt("##maxcellsA", &m_maxCellsA, 0, 0);
            if (m_maxCellsA < 1) m_maxCellsA = 1;
            ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("Max cells B");
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fieldW); ImGui::InputInt("##maxcellsB", &m_maxCellsB, 0, 0);
            if (m_maxCellsB < 1) m_maxCellsB = 1;

            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("Layers A");
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fieldW); ImGui::InputInt("##layersA", &m_layersA, 0, 0);
            if (m_layersA < 1) m_layersA = 1;
            ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("Layers B");
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(fieldW); ImGui::InputInt("##layersB", &m_layersB, 0, 0);
            if (m_layersB < 1) m_layersB = 1;

            ImGui::EndTable();
        }

        ImGui::Spacing();
        ImGui::TextColored(ImVec4(0.6f, 0.75f, 0.9f, 1.0f), "Thresholds");
        if (ImGui::BeginTable("##thTbl", 4, ImGuiTableFlags_SizingFixedFit))
        {
            ImGui::TableSetupColumn("L1", 0, 100.0f);
            ImGui::TableSetupColumn("V1", 0, 70.0f);
            ImGui::TableSetupColumn("L2", 0, 100.0f);
            ImGui::TableSetupColumn("V2", 0, 70.0f);

            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("Max strain");
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(60.0f); ImGui::InputFloat("##maxstrain", &m_maxMeanStrain, 0, 0, "%.4f");
            if (m_maxMeanStrain < 0.0001f) m_maxMeanStrain = 0.0001f;
            ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("Max rotation");
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(60.0f); ImGui::InputFloat("##maxrot", &m_maxRotationDeg, 0, 0, "%.1f");
            if (m_maxRotationDeg < 0.0f) m_maxRotationDeg = 0.0f;

            ImGui::EndTable();
        }

        ImGui::Spacing();
        ImGui::TextColored(ImVec4(0.6f, 0.75f, 0.9f, 1.0f), "Interface Geometry");
        if (ImGui::BeginTable("##geoTbl", 4, ImGuiTableFlags_SizingFixedFit))
        {
            ImGui::TableSetupColumn("L1", 0, colLblW);
            ImGui::TableSetupColumn("V1", 0, 65.0f);
            ImGui::TableSetupColumn("L2", 0, colLblW);
            ImGui::TableSetupColumn("V2", 0, 65.0f);

            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("Z gap");
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(55.0f); ImGui::InputFloat("##zgap", &m_zGap, 0, 0, "%.2f");
            if (m_zGap < 0.0f) m_zGap = 0.0f;
            ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("Vacuum");
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(55.0f); ImGui::InputFloat("##vacuum", &m_vacuum, 0, 0, "%.2f");
            if (m_vacuum < 0.0f) m_vacuum = 0.0f;

            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("Repeat X");
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(55.0f); ImGui::InputInt("##repeatX", &m_repeatX, 0, 0);
            if (m_repeatX < 1) m_repeatX = 1;
            ImGui::TableNextColumn(); ImGui::AlignTextToFramePadding(); ImGui::Text("Repeat Y");
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(55.0f); ImGui::InputInt("##repeatY", &m_repeatY, 0, 0);
            if (m_repeatY < 1) m_repeatY = 1;

            ImGui::EndTable();
        }
    }
    ImGui::EndChild();

    ImGui::SameLine(0.0f, gap);

    // --- Elastic Stiffness ---
    ImGui::BeginChild("##stiffPanel", ImVec2(stiffW, topH), true);
    {
        ImGui::TextColored(ImVec4(0.7f, 0.85f, 1.0f, 1.0f), "Stiffness (GPa)");
        ImGui::SameLine(0.0f, 10.0f);
        ImGui::RadioButton("Structure A##stiff", &m_stiffnessTarget, 0);
        ImGui::SameLine();
        ImGui::RadioButton("Structure B##stiff", &m_stiffnessTarget, 1);
        ImGui::Separator();

        float stiffInW = std::min(50.0f, std::max(35.0f, (stiffW - 30.0f) / 6.0f - 4.0f));
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(2.0f, 2.0f));
        for (int i = 0; i < 6; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                if (j > 0) ImGui::SameLine();
                char id[32];
                std::snprintf(id, sizeof(id), "##C%d%d", i + 1, j + 1);
                ImGui::SetNextItemWidth(stiffInW);
                ImGui::InputFloat(id, &m_stiffness[i][j], 0, 0, "%.1f");
            }
        }
        ImGui::PopStyleVar();
    }
    ImGui::EndChild();

    ImGui::Spacing();

    // Search button
    bool canSearch = !m_structureA.atoms.empty() && !m_structureB.atoms.empty()
                     && m_structureA.hasUnitCell && m_structureB.hasUnitCell;
    if (!canSearch) ImGui::BeginDisabled();
    if (ImGui::Button("Find Matches##iface", ImVec2(160.0f, 28.0f)))
        runSearch();
    if (!canSearch) ImGui::EndDisabled();

    ImGui::SameLine();
    if (m_searchDone)
        ImGui::Text("%d candidates found.", (int)m_candidates.size());
    else
        ImGui::TextDisabled("Load two structures with unit cells, then click Find Matches.");

    ImGui::Separator();

    // =========================================================================
    // BOTTOM ROW: Plot (left) | Preview (right)
    // =========================================================================
    float bottomH = ImGui::GetContentRegionAvail().y - 35.0f;
    float plotW = totalW * 0.45f;
    float previewW = totalW - plotW - gap;

    // --- Left: 2D scatter plot ---
    ImGui::BeginChild("##bottomLeftPlot", ImVec2(plotW, bottomH), false);
    {
        if (m_searchDone && !m_candidates.empty())
        {
            ImGui::Text("Candidate Selection (click a point)");
            draw2DPlot(plotW, bottomH - 20.0f);
        }
        else
        {
            ImGui::TextDisabled("Results will appear here after search.");
        }
    }
    ImGui::EndChild();

    ImGui::SameLine(0.0f, gap);

    // --- Right: Interface structure 3D preview ---
    ImGui::BeginChild("##bottomRightPreview", ImVec2(previewW, bottomH), false);
    {
        if (m_selectedIdx >= 0 && m_selectedIdx < (int)m_interfaceStructures.size())
        {
            const auto& selCand = m_candidates[m_selectedIdx];
            ImGui::Text("Interface Preview (#%d, %d atoms, strain %.3f%%)",
                        selCand.idx, selCand.interfaceAtoms,
                        selCand.meanAbsStrain * 100.0f);
            ImGui::Text("OR: (%.0f %.0f %.0f)A || (%.0f %.0f %.0f)B,  [%d %d 0]A || [%d %d 0]B   theta=%.2f deg",
                        m_planeA[0], m_planeA[1], m_planeA[2],
                        m_planeB[0], m_planeB[1], m_planeB[2],
                        selCand.nMatrix[0][0], selCand.nMatrix[0][1],
                        selCand.mMatrix[0][0], selCand.mMatrix[0][1],
                        selCand.thetaDeg);

            // Only rebuild buffers and auto-fit camera when selection changes
            if (m_selectedIdx != m_lastRenderedResultIdx)
            {
                const Structure& iface = m_interfaceStructures[m_selectedIdx];
                rebuildPreviewBuffers(m_previewBufResult, iface, elementRadii, elementShininess);
                autoFitCamera(m_previewBufResult, m_camRDist);
                m_lastRenderedResultIdx = m_selectedIdx;
            }

            float prevRegionW = previewW - 8.0f;
            float prevRegionH = bottomH - 80.0f;
            if (prevRegionH < 100.0f) prevRegionH = 100.0f;

            ImGui::InvisibleButton("##ifaceResultZone", ImVec2(prevRegionW, prevRegionH));
            ImDrawList* dl = ImGui::GetWindowDrawList();
            ImVec2 rMin = ImGui::GetItemRectMin();
            ImVec2 rMax = ImGui::GetItemRectMax();
            bool rHovered = ImGui::IsItemHovered();
            bool rActive  = ImGui::IsItemActive();
            dl->AddRect(rMin, rMax, ImGui::GetColorU32(ImGuiCol_Border), 2.0f);

            int rw = std::max(1, (int)(rMax.x - rMin.x - 10));
            int rh = std::max(1, (int)(rMax.y - rMin.y - 10));
            ensureFBO(m_fboResult, m_texResult, m_rboResult, m_fboRW, m_fboRH, rw, rh);
            renderToFBO(m_fboResult, rw, rh, m_previewBufResult,
                        m_camRYaw, m_camRPitch, m_camRDist);

            dl->AddImage((ImTextureID)(intptr_t)m_texResult,
                         ImVec2(rMin.x + 5, rMin.y + 5),
                         ImVec2(rMin.x + 5 + rw, rMin.y + 5 + rh),
                         ImVec2(0, 1), ImVec2(1, 0));

            if (rActive && ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f))
            {
                ImVec2 delta = ImGui::GetIO().MouseDelta;
                m_camRYaw   -= delta.x * 0.5f;
                m_camRPitch += delta.y * 0.5f;
            }
            if (rHovered)
            {
                float wheel = ImGui::GetIO().MouseWheel;
                if (wheel != 0.0f)
                {
                    m_camRDist -= wheel * m_camRDist * 0.1f;
                    m_camRDist = std::max(Camera::kMinDistance,
                                          std::min(Camera::kMaxDistance, m_camRDist));
                }
            }
            ImGui::TextDisabled("Left-drag = orbit   Scroll = zoom");
        }
        else
        {
            ImGui::TextDisabled("Select a candidate from the plot to preview.");
        }
    }
    ImGui::EndChild();

    // =========================================================================
    // Action bar
    // =========================================================================
    ImGui::Separator();

    bool canBuild = (m_selectedIdx >= 0 &&
                     m_selectedIdx < (int)m_interfaceStructures.size());
    if (!canBuild) ImGui::BeginDisabled();
    if (ImGui::Button("Build##iface", ImVec2(100.0f, 0.0f)))
    {
        structure = m_interfaceStructures[m_selectedIdx];
        updateBuffers(structure);
        ImGui::CloseCurrentPopup();
        m_isOpen = false;
    }
    if (!canBuild) ImGui::EndDisabled();

    ImGui::SameLine();
    if (ImGui::Button("Close##iface", ImVec2(80.0f, 0.0f)))
    {
        m_isOpen = false;
        ImGui::CloseCurrentPopup();
    }

    ImGui::EndPopup();
}
