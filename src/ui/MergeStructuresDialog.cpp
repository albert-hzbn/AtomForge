#include "ui/MergeStructuresDialog.h"

#include "app/SceneView.h"
#include "camera/Camera.h"
#include "graphics/CylinderMesh.h"
#include "graphics/Picking.h"
#include "graphics/Renderer.h"
#include "graphics/ShadowMap.h"
#include "graphics/SphereMesh.h"
#include "graphics/StructureInstanceBuilder.h"
#include "imgui.h"
#include "util/ElementData.h"

#include <glm/gtc/matrix_transform.hpp>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <utility>

namespace
{
glm::mat3 eulerRotationXYZ(const glm::vec3& rotationDeg)
{
    const float rx = glm::radians(rotationDeg.x);
    const float ry = glm::radians(rotationDeg.y);
    const float rz = glm::radians(rotationDeg.z);
    const glm::mat4 rot4 = glm::rotate(glm::mat4(1.0f), rz, glm::vec3(0.0f, 0.0f, 1.0f))
                         * glm::rotate(glm::mat4(1.0f), ry, glm::vec3(0.0f, 1.0f, 0.0f))
                         * glm::rotate(glm::mat4(1.0f), rx, glm::vec3(1.0f, 0.0f, 0.0f));
    return glm::mat3(rot4);
}

constexpr float kGizmoHitRadius = 8.0f;
constexpr int   kArcSegments    = 32;
constexpr float kArrowheadFrac  = 0.15f;

ImU32 axisColor(int axis, bool highlight)
{
    const ImU32 normal[3]  = { IM_COL32(220, 50, 50, 255),
                               IM_COL32(50, 200, 50, 255),
                               IM_COL32(50, 80, 220, 255) };
    const ImU32 bright[3]  = { IM_COL32(255, 130, 130, 255),
                               IM_COL32(130, 255, 130, 255),
                               IM_COL32(130, 160, 255, 255) };
    return highlight ? bright[axis] : normal[axis];
}

float pointToSegmentDist(const ImVec2& p, const ImVec2& a, const ImVec2& b)
{
    const ImVec2 ab(b.x - a.x, b.y - a.y);
    const ImVec2 ap(p.x - a.x, p.y - a.y);
    const float lenSq = ab.x * ab.x + ab.y * ab.y;
    if (lenSq < 1e-6f)
        return std::sqrt(ap.x * ap.x + ap.y * ap.y);
    float t = (ap.x * ab.x + ap.y * ab.y) / lenSq;
    t = std::max(0.0f, std::min(1.0f, t));
    const float dx = ap.x - t * ab.x;
    const float dy = ap.y - t * ab.y;
    return std::sqrt(dx * dx + dy * dy);
}
}

// ===========================================================================
// Lifecycle
// ===========================================================================
MergeStructuresDialog::MergeStructuresDialog() = default;

MergeStructuresDialog::~MergeStructuresDialog()
{
    if (m_previewFBO) glDeleteFramebuffers(1, &m_previewFBO);
    if (m_previewColorTex) glDeleteTextures(1, &m_previewColorTex);
    if (m_previewDepthRbo) glDeleteRenderbuffers(1, &m_previewDepthRbo);

    if (m_previewShadow.depthFBO)
        glDeleteFramebuffers(1, &m_previewShadow.depthFBO);
    if (m_previewShadow.depthTexture)
        glDeleteTextures(1, &m_previewShadow.depthTexture);

    delete m_previewSphere;
    delete m_previewCylinder;
}

void MergeStructuresDialog::initRenderResources(Renderer& renderer)
{
    m_renderer = &renderer;
    m_previewSphere = new SphereMesh(24, 24);
    m_previewCylinder = new CylinderMesh(16);
    m_previewBuffers.init(m_previewSphere->vbo, m_previewSphere->ebo, m_previewSphere->indexCount,
                          m_previewCylinder->vbo, m_previewCylinder->vertexCount);
    m_previewShadow = createShadowMap(1, 1);
    m_glReady = true;
}

// ===========================================================================
// Menu / drop
// ===========================================================================
void MergeStructuresDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Merge Structures", nullptr, false, enabled))
        m_openRequested = true;
}

void MergeStructuresDialog::feedDroppedFile(const std::string& path)
{
    m_pendingDropPaths.push_back(path);
}

// ===========================================================================
// Entry management
// ===========================================================================
bool MergeStructuresDialog::addStructureFromPath(const std::string& path)
{
    Structure loaded;
    std::string error;
    if (!loadStructureFromFile(path, loaded, error))
    {
        m_status = std::string("Failed to load: ") + (error.empty() ? path : error);
        return false;
    }

    MergeEntry entry;
    entry.pivot = computeCentroid(loaded);
    entry.structure = std::move(loaded);
    entry.name = baseName(path);
    m_entries.push_back(std::move(entry));
    if (m_selectedIndex < 0)
        m_selectedIndex = 0;

    std::ostringstream msg;
    msg << "Loaded " << m_entries.back().name
        << " (atoms=" << m_entries.back().structure.atoms.size() << ")";
    m_status = msg.str();
    m_previewDirty = true;
    m_autoFitOnRebuild = true;
    return true;
}

void MergeStructuresDialog::addStructureFromScene(const Structure& structure)
{
    if (structure.atoms.empty())
    {
        m_status = "Current scene is empty.";
        return;
    }

    MergeEntry entry;
    entry.structure = structure;
    entry.name = "scene";
    entry.pivot = computeCentroid(entry.structure);
    m_entries.push_back(std::move(entry));
    if (m_selectedIndex < 0)
        m_selectedIndex = 0;

    std::ostringstream msg;
    msg << "Added current scene (atoms=" << m_entries.back().structure.atoms.size() << ")";
    m_status = msg.str();
    m_previewDirty = true;
    m_autoFitOnRebuild = true;
}

glm::vec3 MergeStructuresDialog::computeCentroid(const Structure& structure)
{
    if (structure.atoms.empty())
        return glm::vec3(0.0f);

    glm::vec3 sum(0.0f);
    for (const AtomSite& atom : structure.atoms)
        sum += glm::vec3((float)atom.x, (float)atom.y, (float)atom.z);
    return sum / (float)structure.atoms.size();
}

std::string MergeStructuresDialog::baseName(const std::string& path)
{
    const std::string::size_type slash = path.find_last_of("/\\");
    return (slash == std::string::npos) ? path : path.substr(slash + 1);
}

// ===========================================================================
// Gizmo helpers
// ===========================================================================
glm::vec3 MergeStructuresDialog::getTransformedCentroid(int idx) const
{
    const MergeEntry& e = m_entries[idx];
    // Pivot rotated around itself stays at pivot; just apply translation.
    return e.pivot + e.translation;
}

glm::vec2 MergeStructuresDialog::projectToViewport(const glm::vec3& worldPos,
                                                    float vpX, float vpY,
                                                    float vpW, float vpH) const
{
    const glm::vec4 clip = m_lastFrame.projection * m_lastFrame.view * glm::vec4(worldPos, 1.0f);
    if (std::abs(clip.w) < 1e-6f)
        return glm::vec2(-1e5f, -1e5f);
    const glm::vec3 ndc = glm::vec3(clip) / clip.w;
    // NDC x,y in [-1,1]. FBO is Y-up; ImGui image is drawn flipped.
    const float sx = vpX + (ndc.x * 0.5f + 0.5f) * vpW;
    const float sy = vpY + (1.0f - (ndc.y * 0.5f + 0.5f)) * vpH;
    return glm::vec2(sx, sy);
}

void MergeStructuresDialog::drawGizmoOverlay(ImDrawList* dl,
                                             float vpX, float vpY,
                                             float vpW, float vpH,
                                             float gizmoScale,
                                             GizmoMode mode)
{
    if (m_selectedIndex < 0 || m_selectedIndex >= (int)m_entries.size())
        return;

    const glm::vec3 center = getTransformedCentroid(m_selectedIndex);
    const glm::vec2 c2d = projectToViewport(center, vpX, vpY, vpW, vpH);

    if (c2d.x < vpX - 50 || c2d.x > vpX + vpW + 50 ||
        c2d.y < vpY - 50 || c2d.y > vpY + vpH + 50)
        return;

    const ImVec2 cIm(c2d.x, c2d.y);

    // Draw a selection halo so the active structure/gizmo is easier to spot.
    const glm::vec2 haloPoint = projectToViewport(center + glm::vec3(gizmoScale * 1.1f, 0.0f, 0.0f),
                                                  vpX, vpY, vpW, vpH);
    const float haloRadius = std::max(14.0f, glm::length(haloPoint - c2d));
    dl->AddCircle(cIm, haloRadius, IM_COL32(255, 255, 255, 96), 0, 1.8f);
    dl->AddCircle(cIm, haloRadius * 0.72f, IM_COL32(255, 255, 255, 48), 0, 1.0f);

    const glm::vec3 axisDir[3] = {
        glm::vec3(1, 0, 0), glm::vec3(0, 1, 0), glm::vec3(0, 0, 1)
    };

    const bool dragHL[3] = {
        m_drag.active && m_drag.axis == GizmoAxis::X,
        m_drag.active && m_drag.axis == GizmoAxis::Y,
        m_drag.active && m_drag.axis == GizmoAxis::Z,
    };

    // Draw translation arrows in all modes. Active mode is emphasized.
    for (int a = 0; a < 3; ++a)
    {
        const glm::vec3 tip3d = center + axisDir[a] * gizmoScale;
        const glm::vec2 tip2d = projectToViewport(tip3d, vpX, vpY, vpW, vpH);
        const ImVec2 tipIm(tip2d.x, tip2d.y);
        const ImU32 col = axisColor(a, dragHL[a]);
        const float thickness = (mode == GizmoMode::Translate) ? 2.8f : 1.8f;
        dl->AddLine(cIm, tipIm, col, thickness);

        const ImVec2 dir(tipIm.x - cIm.x, tipIm.y - cIm.y);
        const float len = std::sqrt(dir.x * dir.x + dir.y * dir.y);
        if (len > 5.0f)
        {
            const float headLen = len * kArrowheadFrac;
            const ImVec2 nd(dir.x / len, dir.y / len);
            const ImVec2 perp(-nd.y, nd.x);
            const ImVec2 base(tipIm.x - nd.x * headLen, tipIm.y - nd.y * headLen);
            const ImVec2 left(base.x + perp.x * headLen * 0.5f,
                              base.y + perp.y * headLen * 0.5f);
            const ImVec2 right(base.x - perp.x * headLen * 0.5f,
                               base.y - perp.y * headLen * 0.5f);
            dl->AddTriangleFilled(tipIm, left, right, col);
        }

        const char* labels[3] = { "X", "Y", "Z" };
        dl->AddText(ImVec2(tipIm.x + 4, tipIm.y - 8), col, labels[a]);
    }

    // Draw rotation circles in all modes. Active mode is emphasized.
    glm::vec3 us[3] = { {0,1,0}, {1,0,0}, {1,0,0} };
    glm::vec3 vs[3] = { {0,0,1}, {0,0,1}, {0,1,0} };
    for (int a = 0; a < 3; ++a)
    {
        const ImU32 col = axisColor(a, dragHL[a]);
        const float radius = gizmoScale;
        const float thickness = (mode == GizmoMode::Rotate) ? 2.8f : 1.8f;
        ImVec2 prev;
        for (int s = 0; s <= kArcSegments; ++s)
        {
            const float angle = (float)s / (float)kArcSegments * 6.28318530718f;
            const glm::vec3 pt = center
                + (us[a] * std::cos(angle) + vs[a] * std::sin(angle)) * radius;
            const glm::vec2 p2 = projectToViewport(pt, vpX, vpY, vpW, vpH);
            const ImVec2 pIm(p2.x, p2.y);
            if (s > 0)
                dl->AddLine(prev, pIm, col, thickness);
            prev = pIm;
        }

        const glm::vec3 labelPt = center
            + (us[a] * 0.707f + vs[a] * 0.707f) * radius;
        const glm::vec2 lp = projectToViewport(labelPt, vpX, vpY, vpW, vpH);
        const char* labels[3] = { "Rx", "Ry", "Rz" };
        dl->AddText(ImVec2(lp.x + 4, lp.y - 8), col, labels[a]);
    }

    // Center dot
    dl->AddCircleFilled(cIm, 4.0f, IM_COL32(255, 255, 255, 200));
    dl->AddCircle(cIm, 4.0f, IM_COL32(0, 0, 0, 200), 0, 1.5f);
}

std::pair<MergeStructuresDialog::GizmoMode, MergeStructuresDialog::GizmoAxis>
MergeStructuresDialog::hitTestGizmo(float mx, float my,
                                     float vpX, float vpY,
                                     float vpW, float vpH,
                                     float gizmoScale) const
{
    if (m_selectedIndex < 0 || m_selectedIndex >= (int)m_entries.size())
        return { GizmoMode::Translate, GizmoAxis::None };

    const glm::vec3 center = getTransformedCentroid(m_selectedIndex);
    const glm::vec2 c2d = projectToViewport(center, vpX, vpY, vpW, vpH);
    const ImVec2 mousePos(mx, my);

    const glm::vec3 axisDir[3] = {
        glm::vec3(1, 0, 0), glm::vec3(0, 1, 0), glm::vec3(0, 0, 1)
    };

    float bestDist = kGizmoHitRadius;
    GizmoAxis bestAxis = GizmoAxis::None;
    GizmoMode bestMode = GizmoMode::Translate;

    // Test translation arrows
    for (int a = 0; a < 3; ++a)
    {
        const glm::vec3 tip3d = center + axisDir[a] * gizmoScale;
        const glm::vec2 tip2d = projectToViewport(tip3d, vpX, vpY, vpW, vpH);
        const float d = pointToSegmentDist(mousePos,
                                            ImVec2(c2d.x, c2d.y),
                                            ImVec2(tip2d.x, tip2d.y));
        if (d < bestDist)
        {
            bestDist = d;
            bestAxis = (GizmoAxis)a;
            bestMode = GizmoMode::Translate;
        }
    }

    // Test rotation arcs
    glm::vec3 us[3] = { {0,1,0}, {1,0,0}, {1,0,0} };
    glm::vec3 vs[3] = { {0,0,1}, {0,0,1}, {0,1,0} };
    for (int a = 0; a < 3; ++a)
    {
        const float radius = gizmoScale;
        for (int s = 0; s < kArcSegments; ++s)
        {
            const float a0 = (float)s / (float)kArcSegments * 6.28318530718f;
            const float a1 = (float)(s + 1) / (float)kArcSegments * 6.28318530718f;
            const glm::vec3 p0 = center + (us[a] * std::cos(a0) + vs[a] * std::sin(a0)) * radius;
            const glm::vec3 p1 = center + (us[a] * std::cos(a1) + vs[a] * std::sin(a1)) * radius;
            const glm::vec2 s0 = projectToViewport(p0, vpX, vpY, vpW, vpH);
            const glm::vec2 s1 = projectToViewport(p1, vpX, vpY, vpW, vpH);
            const float d = pointToSegmentDist(mousePos,
                                                ImVec2(s0.x, s0.y),
                                                ImVec2(s1.x, s1.y));
            if (d < bestDist)
            {
                bestDist = d;
                bestAxis = (GizmoAxis)a;
                bestMode = GizmoMode::Rotate;
            }
        }
    }

    return { bestMode, bestAxis };
}

int MergeStructuresDialog::pickEntryAtCursor(float lx, float ly,
                                              float vpW, float vpH) const
{
    if (m_entries.empty() || vpW < 1 || vpH < 1)
        return -1;

    const glm::vec3 origin = pickRayOrigin(lx, ly, (int)vpW, (int)vpH,
                                           m_lastFrame.projection, m_lastFrame.view);
    const glm::vec3 dir = pickRayDir(lx, ly, (int)vpW, (int)vpH,
                                     m_lastFrame.projection, m_lastFrame.view);

    static std::vector<float> radii = makeLiteratureCovalentRadii();

    int   bestEntry = -1;
    float bestT     = 1e30f;

    for (int ei = 0; ei < (int)m_entries.size(); ++ei)
    {
        const MergeEntry& entry = m_entries[ei];
        if (!entry.enabled)
            continue;

        const glm::mat3 rot = eulerRotationXYZ(entry.rotationDeg);
        for (const AtomSite& atom : entry.structure.atoms)
        {
            const glm::vec3 p((float)atom.x, (float)atom.y, (float)atom.z);
            const glm::vec3 tp = entry.pivot + rot * (p - entry.pivot) + entry.translation;

            float r = 1.0f;
            if (atom.atomicNumber >= 1 && atom.atomicNumber < (int)radii.size())
                r = radii[atom.atomicNumber];
            if (r <= 0.0f) r = 1.0f;

            const glm::vec3 oc = tp - origin;
            const float t = glm::dot(oc, dir);
            if (t < 0.0f) continue;

            const glm::vec3 closest = origin + dir * t;
            const float distSq = glm::dot(closest - tp, closest - tp);
            if (distSq <= r * r && t < bestT)
            {
                bestT = t;
                bestEntry = ei;
            }
        }
    }

    return bestEntry;
}

// ===========================================================================
// Combined preview structure
// ===========================================================================
Structure MergeStructuresDialog::buildCombinedPreviewStructure() const
{
    Structure combined;
    combined.hasUnitCell = false;

    for (size_t i = 0; i < m_entries.size(); ++i)
    {
        const MergeEntry& entry = m_entries[i];
        if (!entry.enabled)
            continue;

        const glm::mat3 rot = eulerRotationXYZ(entry.rotationDeg);
        for (const AtomSite& atom : entry.structure.atoms)
        {
            const glm::vec3 p((float)atom.x, (float)atom.y, (float)atom.z);
            const glm::vec3 transformed = entry.pivot + rot * (p - entry.pivot) + entry.translation;

            AtomSite out = atom;
            out.x = transformed.x;
            out.y = transformed.y;
            out.z = transformed.z;

            combined.atoms.push_back(out);
        }
    }

    return combined;
}

// ===========================================================================
// FBO management
// ===========================================================================
void MergeStructuresDialog::ensurePreviewFBO(int w, int h)
{
    if (w == m_previewW && h == m_previewH && m_previewFBO != 0)
        return;

    if (m_previewFBO)      { glDeleteFramebuffers(1, &m_previewFBO);      m_previewFBO = 0; }
    if (m_previewColorTex) { glDeleteTextures(1, &m_previewColorTex);     m_previewColorTex = 0; }
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
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_previewColorTex, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, m_previewDepthRbo);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    m_previewW = w;
    m_previewH = h;
}

void MergeStructuresDialog::rebuildPreviewBuffers()
{
    if (!m_glReady)
        return;

    Structure preview = buildCombinedPreviewStructure();
    if (preview.atoms.empty())
    {
        m_previewBuffers.atomCount = 0;
        m_previewDirty = false;
        return;
    }

    static const int kIdent[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    static std::vector<float> radii = makeLiteratureCovalentRadii();
    static std::vector<float> shininess(119, 32.0f);

    StructureInstanceData data = buildStructureInstanceData(preview, false, kIdent, radii, shininess);
    std::array<bool,119> noFilter = {};
    m_previewBuffers.upload(data, false, noFilter);
    if (m_autoFitOnRebuild)
    {
        autoFitPreviewCamera();
        m_autoFitOnRebuild = false;
    }
    m_previewDirty = false;
}

void MergeStructuresDialog::autoFitPreviewCamera()
{
    if (m_previewBuffers.atomCount == 0)
    {
        m_camDistance = 12.0f;
        return;
    }

    float maxR = 0.0f;
    for (size_t i = 0; i < m_previewBuffers.atomPositions.size(); ++i)
    {
        const float r = (i < m_previewBuffers.atomRadii.size()) ? m_previewBuffers.atomRadii[i] : 0.0f;
        const float d = glm::length(m_previewBuffers.atomPositions[i] - m_previewBuffers.orbitCenter) + r;
        maxR = std::max(maxR, d);
    }

    maxR = std::max(maxR, 1.0f);
    const float halfFov = glm::radians(22.5f);
    float dist = maxR / std::sin(halfFov) * 1.2f;
    dist = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, dist));
    m_camDistance = dist;
}

void MergeStructuresDialog::renderPreviewToFBO(int w, int h)
{
    if (!m_glReady || !m_renderer || m_previewBuffers.atomCount == 0)
        return;

    ensurePreviewFBO(w, h);

    GLint prevFbo = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
    GLint prevVP[4];
    glGetIntegerv(GL_VIEWPORT, prevVP);

    Camera cam;
    cam.yaw = m_camYaw;
    cam.pitch = m_camPitch;
    cam.distance = m_camDistance;

    FrameView frame;
    frame.framebufferWidth = w;
    frame.framebufferHeight = h;
    buildFrameView(cam, m_previewBuffers, true, frame);

    // Store for gizmo projection and picking
    m_lastFrame = frame;

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

    if (m_showPreviewBoundingBox)
    {
        m_renderer->drawBoxLines(frame.projection, frame.view,
                                 m_previewBuffers.lineVAO,
                                 m_previewBuffers.boxLines.size());
    }

    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
    glViewport(prevVP[0], prevVP[1], prevVP[2], prevVP[3]);
}

// ===========================================================================
// Main dialog
// ===========================================================================
void MergeStructuresDialog::drawDialog(Structure& structure,
                                       const std::function<void(Structure&)>& updateBuffers)
{
    static GizmoMode gizmoMode = GizmoMode::Translate;

    for (const std::string& path : m_pendingDropPaths)
        addStructureFromPath(path);
    m_pendingDropPaths.clear();

    if (m_openRequested)
    {
        ImGui::OpenPopup("Merge Structures");
        m_openRequested = false;
        m_status.clear();
    }

    m_isOpen = ImGui::IsPopupOpen("Merge Structures");

    ImGui::SetNextWindowSize(ImVec2(1600.0f, 960.0f), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(1100.0f, 680.0f), ImVec2(10000.0f, 10000.0f));
    bool keepOpen = true;
    if (!ImGui::BeginPopupModal("Merge Structures", &keepOpen, ImGuiWindowFlags_NoCollapse))
    {
        m_isOpen = false;
        return;
    }

    m_isOpen = true;

    // ---- Top bar ----
    ImGui::TextWrapped("Arrange multiple structures in 3D before merging. "
                       "Left-click atoms to select, drag gizmo handles to move/rotate, "
                       "and left-drag empty space to orbit the preview.");
    ImGui::TextDisabled("Shortcuts: T = Translate, R = Rotate, Del = Delete selected");

    ImGui::Spacing();

    if (ImGui::Button("Clear All"))
    {
        m_entries.clear();
        m_selectedIndex = -1;
        m_drag.active = false;
        m_status = "Cleared.";
        m_previewDirty = true;
        m_autoFitOnRebuild = true;
    }

    ImGui::SameLine();
    if (ImGui::Checkbox("Show Bounding Box", &m_showPreviewBoundingBox))
        m_previewDirty = true;

    // Keyboard shortcuts
    if (!ImGui::GetIO().WantTextInput)
    {
        if (ImGui::IsKeyPressed(ImGuiKey_T, false))
            gizmoMode = GizmoMode::Translate;
        if (ImGui::IsKeyPressed(ImGuiKey_R, false))
            gizmoMode = GizmoMode::Rotate;

        if ((ImGui::IsKeyPressed(ImGuiKey_Delete, false) ||
             ImGui::IsKeyPressed(ImGuiKey_Backspace, false))
            && m_selectedIndex >= 0 && m_selectedIndex < (int)m_entries.size())
        {
            m_entries.erase(m_entries.begin() + m_selectedIndex);
            if (m_entries.empty())
                m_selectedIndex = -1;
            else if (m_selectedIndex >= (int)m_entries.size())
                m_selectedIndex = (int)m_entries.size() - 1;
            m_drag.active = false;
            m_previewDirty = true;
            m_autoFitOnRebuild = true;
            m_status = "Deleted selected structure.";
        }
    }

    ImGui::Separator();

    // ---- Left panel: list + numeric controls ----
    const float panelH = std::max(520.0f, ImGui::GetContentRegionAvail().y - 80.0f);
    const float leftW = ImGui::GetContentRegionAvail().x * 0.24f;
    ImGui::BeginChild("##merge-left", ImVec2(leftW, panelH), true);
    {
        ImGui::Text("Loaded Structures (%d)", (int)m_entries.size());
        ImGui::Separator();

        int removeIndex = -1;
        for (size_t i = 0; i < m_entries.size(); ++i)
        {
            MergeEntry& entry = m_entries[i];
            ImGui::PushID((int)i);

            if (ImGui::Selectable((entry.name + "##sel").c_str(), m_selectedIndex == (int)i))
            {
                m_selectedIndex = (int)i;
            }
            ImGui::SameLine();
            ImGui::TextDisabled("atoms=%d", (int)entry.structure.atoms.size());

            bool changed = false;
            changed |= ImGui::DragFloat3("Pos (A)", &entry.translation.x,
                                          0.05f, -10000.0f, 10000.0f, "%.3f");
            changed |= ImGui::DragFloat3("Rot (deg)", &entry.rotationDeg.x,
                                          0.25f, -3600.0f, 3600.0f, "%.2f");

            if (ImGui::SmallButton("Reset"))
            {
                entry.translation = glm::vec3(0.0f);
                entry.rotationDeg = glm::vec3(0.0f);
                changed = true;
            }
            ImGui::SameLine();
            if (ImGui::SmallButton("Remove"))
                removeIndex = (int)i;

            if (changed)
                m_previewDirty = true;

            ImGui::Separator();
            ImGui::PopID();
        }

        if (removeIndex >= 0 && removeIndex < (int)m_entries.size())
        {
            m_entries.erase(m_entries.begin() + removeIndex);
            if (m_entries.empty())
                m_selectedIndex = -1;
            else if (m_selectedIndex >= (int)m_entries.size())
                m_selectedIndex = (int)m_entries.size() - 1;
            m_drag.active = false;
            m_previewDirty = true;
            m_autoFitOnRebuild = true;
        }

        if (m_entries.empty())
        {
            ImGui::TextDisabled("No structures loaded yet.");
            ImGui::TextDisabled("Drop CIF/XYZ/PDB/... files here.");
        }
    }
    ImGui::EndChild();

    // ---- Right panel: 3D preview + gizmo ----
    ImGui::SameLine();
    ImGui::BeginChild("##merge-right", ImVec2(0.0f, panelH), true);
    {
        const bool hasSelection = m_selectedIndex >= 0 && m_selectedIndex < (int)m_entries.size();
        if (m_previewDirty)
            rebuildPreviewBuffers();

        ImGui::Text("3D Arrangement View");
        ImGui::Separator();

        const ImVec2 avail = ImGui::GetContentRegionAvail();
        const float previewH = std::max(420.0f, avail.y - 70.0f);
        const float pw = avail.x;
        const float ph = previewH;

        if (m_previewBuffers.atomCount > 0)
        {
            const int iw = std::max(1, (int)pw);
            const int ih = std::max(1, (int)ph);
            renderPreviewToFBO(iw, ih);

            const ImVec2 cursor = ImGui::GetCursorScreenPos();
            ImGui::InvisibleButton("##merge-preview", ImVec2(pw, ph));
            const bool hovered = ImGui::IsItemHovered();
            const bool active  = ImGui::IsItemActive();

            ImDrawList* dl = ImGui::GetWindowDrawList();
            dl->AddImage((ImTextureID)(intptr_t)m_previewColorTex,
                         cursor, ImVec2(cursor.x + pw, cursor.y + ph),
                         ImVec2(0, 1), ImVec2(1, 0));

            const float gizmoWorldScale = m_camDistance * 0.15f;
            const ImVec2 mouseScreen = ImGui::GetIO().MousePos;
            const float mlx = mouseScreen.x - cursor.x;
            const float mly = mouseScreen.y - cursor.y;

            // -- Gizmo drag continuation --
            if (m_drag.active)
            {
                if (ImGui::IsMouseDown(ImGuiMouseButton_Left))
                {
                    const float dx = mlx - m_drag.startMouseX;
                    const float dy = mly - m_drag.startMouseY;
                    const float sensitivity = m_camDistance / std::max(pw, 1.0f) * 2.0f;

                    MergeEntry& sel = m_entries[m_selectedIndex];
                    if (m_drag.mode == GizmoMode::Translate)
                    {
                        const glm::vec3 center = getTransformedCentroid(m_selectedIndex);
                        const glm::vec3 axDirs[3] = { {1,0,0}, {0,1,0}, {0,0,1} };
                        const glm::vec2 c2 = projectToViewport(center, cursor.x, cursor.y, pw, ph);
                        const glm::vec2 a2 = projectToViewport(center + axDirs[(int)m_drag.axis],
                                                                cursor.x, cursor.y, pw, ph);
                        const glm::vec2 screenAxis = a2 - c2;
                        const float axLen = glm::length(screenAxis);
                        if (axLen > 1e-3f)
                        {
                            const glm::vec2 sa = screenAxis / axLen;
                            const float proj = sa.x * dx + sa.y * dy;
                            sel.translation = m_drag.startTranslation;
                            sel.translation[(int)m_drag.axis] += proj * sensitivity;
                        }
                    }
                    else
                    {
                        const float angleDelta = (dx - dy) * 0.5f;
                        sel.rotationDeg = m_drag.startRotationDeg;
                        sel.rotationDeg[(int)m_drag.axis] += angleDelta;
                    }
                    m_previewDirty = true;
                }
                else
                {
                    m_drag.active = false;
                }
            }
            else if (hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left))
            {
                // Check gizmo handles first
                auto hitResult = hitTestGizmo(mouseScreen.x, mouseScreen.y,
                                              cursor.x, cursor.y, pw, ph,
                                              gizmoWorldScale);
                if (hitResult.second != GizmoAxis::None && hasSelection)
                {
                    m_drag.active = true;
                    m_drag.mode   = hitResult.first;
                    m_drag.axis   = hitResult.second;
                    m_drag.startMouseX = mlx;
                    m_drag.startMouseY = mly;
                    m_drag.startTranslation = m_entries[m_selectedIndex].translation;
                    m_drag.startRotationDeg = m_entries[m_selectedIndex].rotationDeg;
                }
                else
                {
                    // Click-to-select via ray cast
                    const int picked = pickEntryAtCursor(mlx, mly, pw, ph);
                    if (picked >= 0)
                    {
                        m_selectedIndex = picked;
                    }
                }
            }

            // Camera orbit (left on empty, or middle/right anywhere in preview)
            if ((hovered || active) && !m_drag.active)
            {
                if (ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f) ||
                    ImGui::IsMouseDragging(ImGuiMouseButton_Middle, 0.0f) ||
                    ImGui::IsMouseDragging(ImGuiMouseButton_Right, 0.0f))
                {
                    const ImVec2 delta = ImGui::GetIO().MouseDelta;
                    m_camYaw   -= delta.x * 0.45f;
                    m_camPitch += delta.y * 0.45f;
                }
            }

            if (hovered)
            {
                const float wheel = ImGui::GetIO().MouseWheel;
                if (wheel != 0.0f)
                {
                    m_camDistance -= wheel * m_camDistance * 0.1f;
                    m_camDistance = std::max(Camera::kMinDistance,
                                             std::min(Camera::kMaxDistance, m_camDistance));
                }
            }

            // Draw gizmo overlay
            if (hasSelection && m_entries[m_selectedIndex].enabled)
            {
                drawGizmoOverlay(dl, cursor.x, cursor.y, pw, ph,
                                 gizmoWorldScale, gizmoMode);
            }
        }
        else
        {
            ImGui::Dummy(ImVec2(pw, ph));
            ImGui::TextDisabled("3D preview will appear after loading structures.");
        }

        // ---- Below-preview controls ----
        ImGui::SeparatorText("Selected Structure Controls");
        if (hasSelection)
        {
            MergeEntry& sel = m_entries[m_selectedIndex];
            ImGui::Text("%s  (%d atoms)", sel.name.c_str(), (int)sel.structure.atoms.size());
            ImGui::TextDisabled("Use gizmo handles in the 3D view to move/rotate this structure.");
            ImGui::Spacing();
            if (ImGui::Button("Delete Selected", ImVec2(150.0f, 0.0f)))
            {
                m_entries.erase(m_entries.begin() + m_selectedIndex);
                if (m_entries.empty())
                    m_selectedIndex = -1;
                else if (m_selectedIndex >= (int)m_entries.size())
                    m_selectedIndex = (int)m_entries.size() - 1;
                m_drag.active = false;
                m_previewDirty = true;
                m_autoFitOnRebuild = true;
            }
        }
        else
        {
            ImGui::TextDisabled("Click a structure in the 3D preview or list to select it.");
        }
    }
    ImGui::EndChild();

    // ---- Bottom bar ----
    const bool canMerge = !m_entries.empty();
    if (!canMerge)
        ImGui::BeginDisabled();

    if (ImGui::Button("Merge Structures", ImVec2(170.0f, 0.0f)))
    {
        Structure merged = buildCombinedPreviewStructure();
        if (merged.atoms.empty())
        {
            m_status = "No enabled structures to merge.";
        }
        else
        {
            structure.atoms.swap(merged.atoms);
            structure.hasUnitCell = false;
            structure.cellOffset = {{0.0, 0.0, 0.0}};
            structure.cellVectors = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
            structure.pbcBoundaryTol = 0.0f;
            structure.grainColors.clear();
            structure.grainRegionIds.clear();
            structure.ipfLoadStatus.clear();

            updateBuffers(structure);

            std::ostringstream msg;
            msg << "Merged " << structure.atoms.size() << " atoms.";
            m_status = msg.str();
        }
    }

    if (!canMerge)
        ImGui::EndDisabled();

    ImGui::SameLine();
    if (ImGui::Button("Close", ImVec2(100.0f, 0.0f)))
        ImGui::CloseCurrentPopup();

    if (!m_status.empty())
    {
        ImGui::SameLine();
        ImGui::TextWrapped("%s", m_status.c_str());
    }

    ImGui::EndPopup();
}
