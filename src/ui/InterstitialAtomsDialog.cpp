#include "ui/InterstitialAtomsDialog.h"

#include "algorithms/MeshLoader.h"
#include "algorithms/InterstitialVoidAnalysis.h"
#include "app/SceneView.h"
#include "camera/Camera.h"
#include "graphics/CustomStructureShaders.h"
#include "graphics/CylinderMesh.h"
#include "graphics/Renderer.h"
#include "graphics/Shader.h"
#include "graphics/SphereMesh.h"
#include "graphics/StructureInstanceBuilder.h"
#include "math/StructureMath.h"
#include "ui/PeriodicTableDialog.h"
#include "util/ElementData.h"

#include "imgui.h"

#include <GL/glew.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

namespace
{
bool rayIntersectsTriangle(const glm::vec3& orig,
                           const glm::vec3& dir,
                           const glm::vec3& v0,
                           const glm::vec3& v1,
                           const glm::vec3& v2)
{
    const glm::vec3 e1 = v1 - v0;
    const glm::vec3 e2 = v2 - v0;
    const glm::vec3 pvec = glm::cross(dir, e2);
    const float det = glm::dot(e1, pvec);

    if (std::abs(det) < 1e-8f)
        return false;

    const float invDet = 1.0f / det;
    const glm::vec3 tvec = orig - v0;
    const float u = glm::dot(tvec, pvec) * invDet;
    if (u < 0.0f || u > 1.0f)
        return false;

    const glm::vec3 qvec = glm::cross(tvec, e1);
    const float v = glm::dot(dir, qvec) * invDet;
    if (v < 0.0f || (u + v) > 1.0f)
        return false;

    const float t = glm::dot(e2, qvec) * invDet;
    return t > 1e-6f;
}

std::string toLower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return (char)std::tolower(c);
    });
    return s;
}

bool hasExt(const std::string& path, const char* ext)
{
    const std::string lower = toLower(path);
    const size_t n = std::strlen(ext);
    return lower.size() >= n && lower.compare(lower.size() - n, n, ext) == 0;
}

std::string baseName(const std::string& path)
{
    const size_t slash = path.find_last_of("\\/");
    return (slash == std::string::npos) ? path : path.substr(slash + 1);
}

glm::mat4 makeMeshTransform(float scale,
                            const glm::vec3& translation,
                            const glm::vec3& rotationDeg)
{
    glm::mat4 transform(1.0f);
    transform = glm::translate(transform, translation);
    transform = glm::rotate(transform, glm::radians(rotationDeg.x), glm::vec3(1.0f, 0.0f, 0.0f));
    transform = glm::rotate(transform, glm::radians(rotationDeg.y), glm::vec3(0.0f, 1.0f, 0.0f));
    transform = glm::rotate(transform, glm::radians(rotationDeg.z), glm::vec3(0.0f, 0.0f, 1.0f));
    transform = glm::scale(transform, glm::vec3(scale));
    return transform;
}

// ---- Gizmo helpers ----
constexpr float kMeshGizmoHitR  = 8.0f;
constexpr int   kMeshArcSegs    = 32;
constexpr float kMeshArrowFrac  = 0.15f;

ImU32 meshAxisCol(int axis, bool hi)
{
    const ImU32 n[3] = { IM_COL32(220, 50, 50, 255), IM_COL32(50, 200, 50, 255), IM_COL32(50, 80, 220, 255) };
    const ImU32 h[3] = { IM_COL32(255, 130, 130, 255), IM_COL32(130, 255, 130, 255), IM_COL32(130, 160, 255, 255) };
    return hi ? h[axis] : n[axis];
}

float meshSegDist(const ImVec2& p, const ImVec2& a, const ImVec2& b)
{
    const float dx = b.x - a.x, dy = b.y - a.y;
    const float len2 = dx * dx + dy * dy;
    if (len2 < 1e-6f) return std::hypot(p.x - a.x, p.y - a.y);
    const float t = std::max(0.0f, std::min(1.0f, ((p.x-a.x)*dx + (p.y-a.y)*dy) / len2));
    const float ex = a.x + t*dx - p.x, ey = a.y + t*dy - p.y;
    return std::sqrt(ex*ex + ey*ey);
}

} // namespace

InterstitialAtomsDialog::InterstitialAtomsDialog()
{
    rebuildVoidTypeOptions();
}

InterstitialAtomsDialog::~InterstitialAtomsDialog()
{
    delete m_previewSphere;
    delete m_previewCylinder;
    if (m_previewFBO)      glDeleteFramebuffers(1, &m_previewFBO);
    if (m_previewColorTex) glDeleteTextures(1, &m_previewColorTex);
    if (m_previewDepthRbo) glDeleteRenderbuffers(1, &m_previewDepthRbo);
    if (m_meshProgram) glDeleteProgram(m_meshProgram);
    if (m_meshVAO) glDeleteVertexArrays(1, &m_meshVAO);
    if (m_meshVBO) glDeleteBuffers(1, &m_meshVBO);
}

void InterstitialAtomsDialog::initRenderResources(Renderer& renderer)
{
    m_renderer        = &renderer;
    m_previewSphere   = new SphereMesh(24, 24);
    m_previewCylinder = new CylinderMesh(16);
    m_previewBuffers.init(m_previewSphere->vbo, m_previewSphere->ebo, m_previewSphere->indexCount,
                          m_previewCylinder->vbo, m_previewCylinder->vertexCount);
    m_previewShadow = createShadowMap(1, 1);
    m_meshProgram = createProgram(customStructureMeshVertexShader(),
                                  customStructureMeshFragmentShader());
    glBindAttribLocation(m_meshProgram, 0, "position");
    glBindAttribLocation(m_meshProgram, 1, "normal");
    glLinkProgram(m_meshProgram);
    m_glReady = true;
}

void InterstitialAtomsDialog::rebuildMeshSurface()
{
    if (m_meshVertices.empty() || m_meshIndices.empty())
    {
        m_meshTriCount = 0;
        return;
    }

    struct Vertex { glm::vec3 pos; glm::vec3 norm; };
    std::vector<Vertex> verts;
    verts.reserve(m_meshIndices.size());

    const glm::mat4 model = makeMeshTransform(m_meshScale, m_meshTranslation, m_meshRotationDeg);
    const glm::mat3 normalMat = glm::mat3(glm::transpose(glm::inverse(model)));

    for (size_t i = 0; i + 2 < m_meshIndices.size(); i += 3)
    {
        const unsigned int i0 = m_meshIndices[i];
        const unsigned int i1 = m_meshIndices[i + 1];
        const unsigned int i2 = m_meshIndices[i + 2];
        if (i0 >= m_meshVertices.size() || i1 >= m_meshVertices.size() || i2 >= m_meshVertices.size())
            continue;

        const glm::vec3 a = glm::vec3(model * glm::vec4(m_meshVertices[i0], 1.0f));
        const glm::vec3 b = glm::vec3(model * glm::vec4(m_meshVertices[i1], 1.0f));
        const glm::vec3 c = glm::vec3(model * glm::vec4(m_meshVertices[i2], 1.0f));
        glm::vec3 n = glm::cross(b - a, c - a);
        const float len = glm::length(n);
        if (len > 1e-8f)
            n /= len;
        else
            n = glm::normalize(normalMat * glm::vec3(0.0f, 0.0f, 1.0f));

        verts.push_back({a, n});
        verts.push_back({b, n});
        verts.push_back({c, n});
    }

    m_meshTriCount = (int)(verts.size() / 3);
    if (m_meshTriCount <= 0)
        return;

    if (!m_meshVAO)
    {
        glGenVertexArrays(1, &m_meshVAO);
        glGenBuffers(1, &m_meshVBO);
    }

    glBindVertexArray(m_meshVAO);
    glBindBuffer(GL_ARRAY_BUFFER, m_meshVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 (GLsizeiptr)(verts.size() * sizeof(Vertex)),
                 verts.data(), GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                          (void*)offsetof(Vertex, pos));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                          (void*)offsetof(Vertex, norm));
    glBindVertexArray(0);
}

void InterstitialAtomsDialog::renderMeshSurfaceToPreview()
{
    if (!m_meshLoaded || !m_meshProgram || !m_meshVAO || m_meshTriCount <= 0)
        return;

    glUseProgram(m_meshProgram);
    glUniformMatrix4fv(glGetUniformLocation(m_meshProgram, "projection"),
                       1, GL_FALSE, &m_lastPreviewFrame.projection[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(m_meshProgram, "view"),
                       1, GL_FALSE, &m_lastPreviewFrame.view[0][0]);
    const glm::vec3 lightDir = glm::normalize(m_lastPreviewFrame.cameraPosition + glm::vec3(0.3f, 0.8f, 0.2f));
    glUniform3fv(glGetUniformLocation(m_meshProgram, "uLightDir"), 1, &lightDir[0]);
    glUniform3fv(glGetUniformLocation(m_meshProgram, "uViewPos"), 1, &m_lastPreviewFrame.cameraPosition[0]);
    glUniform3f(glGetUniformLocation(m_meshProgram, "uColor"), 0.45f, 0.65f, 0.85f);

    glBindVertexArray(m_meshVAO);
    glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(m_meshTriCount * 3));
    glBindVertexArray(0);
    glUseProgram(0);
}

void InterstitialAtomsDialog::ensurePreviewFBO(int w, int h)
{
    if (w == m_previewFboW && h == m_previewFboH && m_previewFBO != 0)
        return;

    if (m_previewFBO)      { glDeleteFramebuffers(1,  &m_previewFBO);      m_previewFBO      = 0; }
    if (m_previewColorTex) { glDeleteTextures(1,       &m_previewColorTex); m_previewColorTex = 0; }
    if (m_previewDepthRbo) { glDeleteRenderbuffers(1,  &m_previewDepthRbo); m_previewDepthRbo = 0; }

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

    m_previewFboW = w;
    m_previewFboH = h;
}

void InterstitialAtomsDialog::rebuildPreviewBuffers(const std::vector<float>& radii,
                                                    const std::vector<float>& shininess)
{
    if (!m_glReady || m_source.atoms.empty())
        return;

    static const int kIdent[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    StructureInstanceData data = buildStructureInstanceData(m_source, false, kIdent, radii, shininess);
    std::array<bool, 119> noFilter = {};
    m_previewBuffers.upload(data, false, noFilter);
    m_previewBufDirty = false;
}

void InterstitialAtomsDialog::renderPreviewToFBO(int w, int h)
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

    m_lastPreviewFrame.framebufferWidth  = w;
    m_lastPreviewFrame.framebufferHeight = h;
    buildFrameView(cam, m_previewBuffers, true, m_lastPreviewFrame);

    glBindFramebuffer(GL_FRAMEBUFFER, m_previewFBO);
    glViewport(0, 0, w, h);
    glEnable(GL_DEPTH_TEST);

    const ImVec4& bg = ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
    glClearColor(bg.x, bg.y, bg.z, bg.w);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    m_renderer->drawBonds(m_lastPreviewFrame.projection, m_lastPreviewFrame.view,
                          m_lastPreviewFrame.lightPosition, m_lastPreviewFrame.cameraPosition,
                          m_previewBuffers.tabCylinderVAO, m_previewBuffers.tabCylinderVertexCount,
                          m_previewBuffers.bondCount);

    m_renderer->drawAtoms(m_lastPreviewFrame.projection, m_lastPreviewFrame.view,
                          m_lastPreviewFrame.lightMVP, m_lastPreviewFrame.lightPosition,
                          m_lastPreviewFrame.cameraPosition,
                          m_previewShadow,
                          m_previewBuffers.tabSphereVAO, m_previewBuffers.tabSphereIndexCount,
                          m_previewBuffers.atomCount);

    m_renderer->drawBoxLines(m_lastPreviewFrame.projection, m_lastPreviewFrame.view,
                             m_previewBuffers.lineVAO,
                             m_previewBuffers.boxLines.size());

    if (m_meshLoaded)
    {
        rebuildMeshSurface();
        renderMeshSurfaceToPreview();
    }

    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
    glViewport(prevVP[0], prevVP[1], prevVP[2], prevVP[3]);
}

void InterstitialAtomsDialog::autoFitPreviewCamera()
{
    m_camYaw   = 45.0f;
    m_camPitch = 35.0f;

    if (m_previewBuffers.atomCount == 0)
    {
        m_camDistance = 10.0f;
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
    float dist = maxR / std::sin(halfFov) * 1.15f;
    m_camDistance = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, dist));
}

// ---------------------------------------------------------------------------
// Mesh gizmo helper methods
// ---------------------------------------------------------------------------

glm::vec3 InterstitialAtomsDialog::getMeshCenter() const
{
    if (!m_meshLoaded || m_meshVertices.empty())
        return m_meshTranslation;

    glm::vec3 sum(0.0f);
    for (const glm::vec3& v : m_meshVertices)
        sum += v;
    const glm::vec3 localCenter = sum / (float)m_meshVertices.size();
    const glm::mat4 t = makeMeshTransform(m_meshScale, m_meshTranslation, m_meshRotationDeg);
    return glm::vec3(t * glm::vec4(localCenter, 1.0f));
}

glm::vec2 InterstitialAtomsDialog::projectToCanvas(const glm::vec3& wp,
                                                    const ImVec2& canvasMin,
                                                    const ImVec2& canvasSize) const
{
    const glm::vec4 clip = m_lastPreviewFrame.projection
                         * m_lastPreviewFrame.view
                         * glm::vec4(wp, 1.0f);
    if (std::abs(clip.w) < 1e-6f)
        return glm::vec2(-1e5f, -1e5f);
    const glm::vec3 ndc = glm::vec3(clip) / clip.w;
    return glm::vec2(canvasMin.x + (ndc.x * 0.5f + 0.5f) * canvasSize.x,
                     canvasMin.y + (1.0f - (ndc.y * 0.5f + 0.5f)) * canvasSize.y);
}

void InterstitialAtomsDialog::drawMeshGizmo(ImDrawList* dl,
                                             const ImVec2& canvasMin,
                                             const ImVec2& canvasSize,
                                             float gizmoScale)
{
    if (!dl || !m_meshLoaded || m_meshVertices.empty())
        return;

    const glm::vec3 center = getMeshCenter();
    const glm::vec2 c2     = projectToCanvas(center, canvasMin, canvasSize);
    const ImVec2 cIm(c2.x, c2.y);

    const bool tDrag = m_meshGizmoDrag.active && m_meshGizmoDrag.mode == MeshGizmoMode::Translate;
    const bool rDrag = m_meshGizmoDrag.active && m_meshGizmoDrag.mode == MeshGizmoMode::Rotate;
    const bool sDrag = m_meshGizmoDrag.active && m_meshGizmoDrag.mode == MeshGizmoMode::Scale;

    const glm::vec3 axisDir[3] = { {1,0,0}, {0,1,0}, {0,0,1} };

    // Translation arrows
    for (int a = 0; a < 3; ++a)
    {
        const bool hl = tDrag && (int)m_meshGizmoDrag.axis == a;
        const glm::vec2 tip2 = projectToCanvas(center + axisDir[a] * gizmoScale, canvasMin, canvasSize);
        const ImVec2 tipIm(tip2.x, tip2.y);
        const ImU32 col = meshAxisCol(a, hl);
        const float thick = (m_meshGizmoMode == MeshGizmoMode::Translate) ? 2.8f : 1.4f;
        dl->AddLine(cIm, tipIm, col, thick);

        const float dx = tipIm.x - cIm.x, dy = tipIm.y - cIm.y;
        const float len = std::sqrt(dx*dx + dy*dy);
        if (len > 5.0f)
        {
            const float hl2 = len * kMeshArrowFrac;
            const float nx = dx/len, ny = dy/len;
            const ImVec2 base(tipIm.x - nx*hl2, tipIm.y - ny*hl2);
            const ImVec2 L(base.x + (-ny)*hl2*0.5f, base.y + nx*hl2*0.5f);
            const ImVec2 R(base.x - (-ny)*hl2*0.5f, base.y - nx*hl2*0.5f);
            dl->AddTriangleFilled(tipIm, L, R, col);
        }
        const char* lbl[3] = {"X","Y","Z"};
        dl->AddText(ImVec2(tipIm.x+4, tipIm.y-8), col, lbl[a]);
    }

    // Rotation arcs
    const glm::vec3 us[3] = { {0,1,0}, {1,0,0}, {1,0,0} };
    const glm::vec3 vs[3] = { {0,0,1}, {0,0,1}, {0,1,0} };
    for (int a = 0; a < 3; ++a)
    {
        const bool hl = rDrag && (int)m_meshGizmoDrag.axis == a;
        const ImU32 col = meshAxisCol(a, hl);
        const float thick = (m_meshGizmoMode == MeshGizmoMode::Rotate) ? 2.8f : 1.4f;
        ImVec2 prev;
        for (int s = 0; s <= kMeshArcSegs; ++s)
        {
            const float ang = (float)s / kMeshArcSegs * 6.28318530718f;
            const glm::vec3 pt = center + (us[a]*std::cos(ang) + vs[a]*std::sin(ang)) * gizmoScale;
            const glm::vec2 p2 = projectToCanvas(pt, canvasMin, canvasSize);
            const ImVec2 pIm(p2.x, p2.y);
            if (s > 0) dl->AddLine(prev, pIm, col, thick);
            prev = pIm;
        }
        const glm::vec3 lPt = center + (us[a]*0.707f + vs[a]*0.707f) * gizmoScale;
        const glm::vec2 lp  = projectToCanvas(lPt, canvasMin, canvasSize);
        const char* lbl[3] = {"Rx","Ry","Rz"};
        dl->AddText(ImVec2(lp.x+4, lp.y-8), col, lbl[a]);
    }

    // Scale handle (uniform, along XYZ diagonal)
    {
        const glm::vec3 diag = glm::normalize(glm::vec3(1.0f));
        const glm::vec2 h2   = projectToCanvas(center + diag * gizmoScale * 1.15f, canvasMin, canvasSize);
        const ImVec2 hIm(h2.x, h2.y);
        const bool inScale = (m_meshGizmoMode == MeshGizmoMode::Scale);
        const ImU32 scaleCol = sDrag ? IM_COL32(255, 255, 100, 255)
                             : inScale ? IM_COL32(235, 220, 150, 240)
                             : IM_COL32(200, 190, 130, 170);
        const float thick = inScale ? 2.0f : 1.2f;
        dl->AddLine(cIm, hIm, scaleCol, thick);
        const float hs = inScale ? 7.0f : 5.0f;
        dl->AddRectFilled(ImVec2(hIm.x-hs, hIm.y-hs), ImVec2(hIm.x+hs, hIm.y+hs), scaleCol);
        dl->AddRect(ImVec2(hIm.x-hs, hIm.y-hs), ImVec2(hIm.x+hs, hIm.y+hs), IM_COL32(0,0,0,200), 0, 0, 1.5f);
        dl->AddText(ImVec2(hIm.x+hs+2, hIm.y-7), scaleCol, "S");
    }

    // Center dot
    dl->AddCircleFilled(cIm, 4.0f, IM_COL32(255,255,255,200));
    dl->AddCircle(cIm, 4.0f, IM_COL32(0,0,0,200), 0, 1.5f);
}

std::pair<InterstitialAtomsDialog::MeshGizmoMode, InterstitialAtomsDialog::MeshGizmoAxis>
InterstitialAtomsDialog::hitTestMeshGizmo(float mx, float my,
                                           const ImVec2& canvasMin,
                                           const ImVec2& canvasSize,
                                           float gizmoScale) const
{
    if (!m_meshLoaded || m_meshVertices.empty())
        return { MeshGizmoMode::Translate, MeshGizmoAxis::None };

    const glm::vec3 center = getMeshCenter();
    const glm::vec2 c2     = projectToCanvas(center, canvasMin, canvasSize);
    const ImVec2    mp(mx, my);

    float         bestDist = kMeshGizmoHitR;
    MeshGizmoAxis bestAxis = MeshGizmoAxis::None;
    MeshGizmoMode bestMode = MeshGizmoMode::Translate;

    const glm::vec3 axisDir[3] = { {1,0,0}, {0,1,0}, {0,0,1} };

    // Test translation arrows
    for (int a = 0; a < 3; ++a)
    {
        const glm::vec2 tip2 = projectToCanvas(center + axisDir[a] * gizmoScale, canvasMin, canvasSize);
        const float d = meshSegDist(mp, ImVec2(c2.x, c2.y), ImVec2(tip2.x, tip2.y));
        if (d < bestDist)
        {
            bestDist = d;
            bestAxis = (MeshGizmoAxis)a;
            bestMode = MeshGizmoMode::Translate;
        }
    }

    // Test rotation arcs
    const glm::vec3 us[3] = { {0,1,0}, {1,0,0}, {1,0,0} };
    const glm::vec3 vs[3] = { {0,0,1}, {0,0,1}, {0,1,0} };
    for (int a = 0; a < 3; ++a)
    {
        for (int s = 0; s < kMeshArcSegs; ++s)
        {
            const float a0 = (float)s / kMeshArcSegs * 6.28318530718f;
            const float a1 = (float)(s+1) / kMeshArcSegs * 6.28318530718f;
            const glm::vec2 s0 = projectToCanvas(center + (us[a]*std::cos(a0) + vs[a]*std::sin(a0))*gizmoScale, canvasMin, canvasSize);
            const glm::vec2 s1 = projectToCanvas(center + (us[a]*std::cos(a1) + vs[a]*std::sin(a1))*gizmoScale, canvasMin, canvasSize);
            const float d = meshSegDist(mp, ImVec2(s0.x, s0.y), ImVec2(s1.x, s1.y));
            if (d < bestDist)
            {
                bestDist = d;
                bestAxis = (MeshGizmoAxis)a;
                bestMode = MeshGizmoMode::Rotate;
            }
        }
    }

    // Test scale handle
    {
        const glm::vec3 diag = glm::normalize(glm::vec3(1.0f));
        const glm::vec2 h2   = projectToCanvas(center + diag * gizmoScale * 1.15f, canvasMin, canvasSize);
        const float dx = mx - h2.x, dy = my - h2.y;
        const float d  = std::sqrt(dx*dx + dy*dy);
        if (d < bestDist + 3.0f)
        {
            bestAxis = MeshGizmoAxis::None;
            bestMode = MeshGizmoMode::Scale;
        }
    }

    return { bestMode, bestAxis };
}

void InterstitialAtomsDialog::drawPreviewOverlays(ImDrawList* dl,
                                                  const ImVec2& canvasMin,
                                                  const ImVec2& canvasSize)
{
    if (!dl)
        return;

    const bool hasVoids = !m_voids.empty();
    const bool hasMesh  = m_meshLoaded && !m_meshVertices.empty();
    if (!hasVoids && !hasMesh)
        return;

    // ---------------------------------------------------------------
    // Gizmo interaction (before void picking so it can consume clicks)
    // ---------------------------------------------------------------
    const ImVec2 mp          = ImGui::GetIO().MousePos;
    const bool   mouseClick  = ImGui::IsItemHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Left);
    const bool   mouseDown   = ImGui::IsMouseDown(ImGuiMouseButton_Left);
    bool         gizmoConsumed = false;

    if (hasMesh)
    {
        const float gizmoScale = m_camDistance * 0.14f;
        const glm::vec3 meshCenter = getMeshCenter();

        if (m_meshGizmoDrag.active)
        {
            if (mouseDown)
            {
                const float dx = mp.x - m_meshGizmoDrag.startMouseX;
                const float dy = mp.y - m_meshGizmoDrag.startMouseY;

                if (m_meshGizmoDrag.mode == MeshGizmoMode::Translate)
                {
                    const glm::vec3 axDirs[3] = { {1,0,0}, {0,1,0}, {0,0,1} };
                    const int ai = (int)m_meshGizmoDrag.axis;
                    const glm::vec2 c2 = projectToCanvas(meshCenter, canvasMin, canvasSize);
                    const glm::vec2 a2 = projectToCanvas(meshCenter + axDirs[ai], canvasMin, canvasSize);
                    const glm::vec2 sa = a2 - c2;
                    const float saLen = glm::length(sa);
                    if (saLen > 1e-3f)
                    {
                        const glm::vec2 saN = sa / saLen;
                        const float proj = saN.x * dx + saN.y * dy;
                        const float wpx  = m_camDistance * 2.0f / std::max(canvasSize.x, 1.0f);
                        m_meshTranslation = m_meshGizmoDrag.startTranslation;
                        m_meshTranslation[ai] += proj * wpx;
                    }
                }
                else if (m_meshGizmoDrag.mode == MeshGizmoMode::Rotate)
                {
                    const float delta = (dx - dy) * 0.5f;
                    m_meshRotationDeg = m_meshGizmoDrag.startRotationDeg;
                    m_meshRotationDeg[(int)m_meshGizmoDrag.axis] += delta;
                }
                else // Scale (uniform)
                {
                    const glm::vec3 diag = glm::normalize(glm::vec3(1.0f));
                    const glm::vec2 c2   = projectToCanvas(meshCenter, canvasMin, canvasSize);
                    const glm::vec2 h2   = projectToCanvas(meshCenter + diag * gizmoScale * 1.15f, canvasMin, canvasSize);
                    const glm::vec2 sa   = h2 - c2;
                    const float saLen    = glm::length(sa);
                    if (saLen > 1e-3f)
                    {
                        const glm::vec2 saN = sa / saLen;
                        const float proj    = saN.x * dx + saN.y * dy;
                        const float relD    = proj / 100.0f;
                        m_meshScale = std::max(1e-4f, m_meshGizmoDrag.startScale * (1.0f + relD));
                    }
                }
            }
            else
            {
                m_meshGizmoDrag.active = false;
            }
            gizmoConsumed = true;
        }
        else if (mouseClick)
        {
            auto [mode, axis] = hitTestMeshGizmo(mp.x, mp.y, canvasMin, canvasSize, gizmoScale);
            if (axis != MeshGizmoAxis::None || mode == MeshGizmoMode::Scale)
            {
                m_meshGizmoDrag.active           = true;
                m_meshGizmoDrag.mode             = mode;
                m_meshGizmoDrag.axis             = axis;
                m_meshGizmoDrag.startMouseX      = mp.x;
                m_meshGizmoDrag.startMouseY      = mp.y;
                m_meshGizmoDrag.startTranslation = m_meshTranslation;
                m_meshGizmoDrag.startRotationDeg = m_meshRotationDeg;
                m_meshGizmoDrag.startScale       = m_meshScale;
                gizmoConsumed = true;
            }
        }
    }

    // ---------------------------------------------------------------
    // Shared projection helper
    // ---------------------------------------------------------------
    auto project = [&](const glm::vec3& wp) -> std::pair<ImVec2, bool>
    {
        const glm::vec4 clip = m_lastPreviewFrame.projection
                             * m_lastPreviewFrame.view
                             * glm::vec4(wp, 1.0f);
        if (clip.w <= 1e-5f)
            return { {}, false };
        const float nx = clip.x / clip.w;
        const float ny = clip.y / clip.w;
        return { { canvasMin.x + (nx * 0.5f + 0.5f) * canvasSize.x,
                   canvasMin.y + (1.0f - (ny * 0.5f + 0.5f)) * canvasSize.y }, true };
    };

    // ---------------------------------------------------------------
    // Void center rendering and click-selection
    // ---------------------------------------------------------------
    if (hasVoids)
    {
        if (m_selectedVoidIndex < 0 || m_selectedVoidIndex >= (int)m_voids.size())
            m_selectedVoidIndex = 0;
        if (m_selectedVoidIndices.empty())
            m_selectedVoidIndices.push_back(m_selectedVoidIndex);

        const std::array<ImU32, 3> fillCol = {
            ImColor(0.20f, 0.80f, 1.00f, 0.18f),
            ImColor(1.00f, 0.68f, 0.20f, 0.20f),
            ImColor(0.78f, 0.40f, 1.00f, 0.16f)
        };
        const std::array<ImU32, 3> edgeCol = {
            ImColor(0.20f, 0.80f, 1.00f, 0.95f),
            ImColor(1.00f, 0.68f, 0.20f, 0.95f),
            ImColor(0.78f, 0.40f, 1.00f, 0.92f)
        };

        const bool canPick = !gizmoConsumed
                           && ImGui::IsItemHovered()
                           && ImGui::IsMouseClicked(ImGuiMouseButton_Left);
        const bool ctrlHeld = ImGui::GetIO().KeyCtrl;
        float bestPickDist2 = 64.0f;
        int   bestPick      = -1;

        auto isVoidSelected = [&](int idx) -> bool
        {
            return std::find(m_selectedVoidIndices.begin(), m_selectedVoidIndices.end(), idx)
                != m_selectedVoidIndices.end();
        };

        for (int vi = 0; vi < (int)m_voids.size(); ++vi)
        {
            const VoidCandidate& vc = m_voids[(size_t)vi];
            const int ci = (vc.type == VoidType::Tetrahedral) ? 0
                        : (vc.type == VoidType::Octahedral)   ? 1 : 2;

            auto [sp, spOk] = project(vc.position);
            if (!spOk) continue;

            const bool  selected  = isVoidSelected(vi);
            const bool  primary   = (vi == m_selectedVoidIndex);
            const float radius    = primary ? 5.4f : (selected ? 4.6f : 3.6f);
            const float thickness = primary ? 2.3f : (selected ? 1.9f : 1.2f);
            dl->AddCircle(sp, radius, edgeCol[(size_t)ci], 0, thickness);
            dl->AddCircleFilled(sp, selected ? 2.0f : 1.4f, edgeCol[(size_t)ci]);

            if (canPick)
            {
                const float ddx = mp.x - sp.x, ddy = mp.y - sp.y;
                const float d2  = ddx*ddx + ddy*ddy;
                if (d2 < bestPickDist2) { bestPickDist2 = d2; bestPick = vi; }
            }
        }

        if (bestPick >= 0)
        {
            m_selectedVoidIndex = bestPick;
            if (ctrlHeld)
            {
                const auto it = std::find(m_selectedVoidIndices.begin(), m_selectedVoidIndices.end(), bestPick);
                if (it == m_selectedVoidIndices.end())
                    m_selectedVoidIndices.push_back(bestPick);
                else if (m_selectedVoidIndices.size() > 1)
                    m_selectedVoidIndices.erase(it);
            }
            else
            {
                m_selectedVoidIndices.clear();
                m_selectedVoidIndices.push_back(bestPick);
            }
        }

        // Draw selected void polyhedron
        const VoidCandidate& v     = m_voids[(size_t)m_selectedVoidIndex];
        const int typeIdx = (v.type == VoidType::Tetrahedral) ? 0
                          : (v.type == VoidType::Octahedral)  ? 1 : 2;
        for (const std::array<glm::vec3, 3>& tri3 : v.faceVertices)
        {
            auto [pa, paOk] = project(tri3[0]);
            auto [pb, pbOk] = project(tri3[1]);
            auto [pc, pcOk] = project(tri3[2]);
            if (!paOk || !pbOk || !pcOk) continue;
            const ImVec2 tri[3] = { pa, pb, pc };
            dl->AddConvexPolyFilled(tri, 3, fillCol[(size_t)typeIdx]);
            dl->AddPolyline(tri, 3, edgeCol[(size_t)typeIdx], ImDrawFlags_Closed, 1.0f);
        }
    }

    // ---------------------------------------------------------------
    // Gizmo overlay
    // ---------------------------------------------------------------
    if (hasMesh)
    {
        const float gizmoScale = m_camDistance * 0.14f;
        drawMeshGizmo(dl, canvasMin, canvasSize, gizmoScale);
    }
}

void InterstitialAtomsDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Add Interstitial Atoms", nullptr, false, enabled))
        m_openRequested = true;
}

void InterstitialAtomsDialog::feedDroppedFile(const std::string& path)
{
    m_pendingDropPaths.push_back(path);
}

bool InterstitialAtomsDialog::acceptsStructurePath(const std::string& path) const
{
    return hasExt(path, ".cif") || hasExt(path, ".xyz") || hasExt(path, ".pdb")
        || hasExt(path, ".mol") || hasExt(path, ".sdf") || hasExt(path, ".mol2")
        || hasExt(path, ".vasp") || hasExt(path, ".pwi") || hasExt(path, ".gjf");
}

bool InterstitialAtomsDialog::acceptsMeshPath(const std::string& path) const
{
    return hasExt(path, ".obj") || hasExt(path, ".stl");
}

bool InterstitialAtomsDialog::loadStructureFromPath(const std::string& path,
                                                    const std::vector<float>& radii,
                                                    const std::vector<float>& shininess)
{
    Structure loaded;
    std::string err;
    if (!loadStructureFromFile(path, loaded, err))
    {
        m_statusIsError = true;
        m_statusMsg = std::string("Failed to load structure: ") + (err.empty() ? path : err);
        return false;
    }

    if (loaded.atoms.empty())
    {
        m_statusIsError = true;
        m_statusMsg = "Loaded structure has no atoms.";
        return false;
    }

    m_source = std::move(loaded);
    m_sourceLoaded = true;
    m_sourceLabel = baseName(path);
    m_voids.clear();
    m_selectedVoidIndices.clear();
    m_selectedVoidIndex = -1;
    rebuildVoidTypeOptions();
    m_statusIsError = false;
    m_statusMsg = "Structure loaded: " + m_sourceLabel;

    rebuildPreviewBuffers(radii, shininess);
    autoFitPreviewCamera();
    return true;
}

bool InterstitialAtomsDialog::loadMeshFromPath(const std::string& path)
{
    std::vector<glm::vec3> vertices;
    std::vector<unsigned int> indices;
    std::string err;
    bool ok = false;

    if (acceptsMeshPath(path) && hasExt(path, ".obj"))
        ok = parseObjMesh(path, vertices, indices, err);
    else if (acceptsMeshPath(path))
        ok = parseStlMesh(path, vertices, indices, err);

    if (!ok)
    {
        m_statusIsError = true;
        m_statusMsg = std::string("Failed to load mesh: ") + (err.empty() ? path : err);
        return false;
    }

    m_meshVertices.swap(vertices);
    m_meshIndices.swap(indices);
    m_meshLoaded = true;
    m_meshLabel = baseName(path);
    rebuildMeshSurface();
    m_statusIsError = false;
    m_statusMsg = "Mesh loaded: " + m_meshLabel;
    return true;
}

bool InterstitialAtomsDialog::detectVoids()
{
    m_voids.clear();
    m_selectedVoidIndices.clear();
    m_selectedVoidIndex = -1;
    rebuildVoidTypeOptions();
    if (!m_sourceLoaded || m_source.atoms.empty())
        return false;

    InterstitialVoidDetectionParams params;
    params.gridResolution = m_gridResolution;
    params.maxVoids = m_maxVoids;
    params.minClearance = m_minClearance;
    params.minSeparation = m_minSeparation;

    InterstitialVoidDetectionResult detection = detectInterstitialVoidRegions(m_source, params);
    if (!detection.success)
        return false;

    m_voids.reserve(detection.regions.size());
    for (const InterstitialVoidRegion& region : detection.regions)
    {
        VoidCandidate v;
        v.position = region.position;
        v.clearance = region.clearance;
        v.coordination = region.coordination;
        v.volume = region.volume;
        v.volumeBin = region.volumeBin;
        if (region.kind == InterstitialVoidKind::Tetrahedral)
            v.type = VoidType::Tetrahedral;
        else if (region.kind == InterstitialVoidKind::Octahedral)
            v.type = VoidType::Octahedral;
        else
            v.type = VoidType::Irregular;

        v.neighborAtomIndices = region.neighborAtomIndices;
        v.faces = region.faces;
        v.faceVertices = region.faceVertices;
        m_voids.push_back(std::move(v));
    }
    if (!m_voids.empty())
    {
        m_selectedVoidIndex = 0;
        m_selectedVoidIndices.push_back(0);
    }
    rebuildVoidTypeOptions();
    return !m_voids.empty();
}

void InterstitialAtomsDialog::rebuildVoidTypeOptions()
{
    std::map<int, bool> previousSelections;
    for (size_t i = 0; i < m_voidTypeOptions.size() && i < m_voidTypeSelections.size(); ++i)
        previousSelections[m_voidTypeOptions[i].key] = (m_voidTypeSelections[i] != 0);

    std::map<int, std::pair<int, double>> volumeGroups;
    for (const VoidCandidate& v : m_voids)
    {
        auto& g = volumeGroups[v.volumeBin];
        g.first += 1;
        g.second += (double)v.volume;
    }

    m_voidTypeOptions.clear();
    m_voidTypeSelections.clear();

    int voidIndex = 1;
    for (const auto& kv : volumeGroups)
    {
        const int bin = kv.first;
        const int count = kv.second.first;
        const double avgVolume = (count > 0) ? (kv.second.second / (double)count) : 0.0;

        std::ostringstream label;
        label.setf(std::ios::fixed, std::ios::floatfield);
        label.precision(3);
        label << "Void " << voidIndex << " (" << count << ") [V " << avgVolume << " A^3]";
        m_voidTypeOptions.push_back({label.str(), 1000 + bin, count});
        const auto it = previousSelections.find(1000 + bin);
        m_voidTypeSelections.push_back((it == previousSelections.end() || it->second) ? 1 : 0);
        ++voidIndex;
    }
}

bool InterstitialAtomsDialog::matchesSelectedVoidFilter(const VoidCandidate& v) const
{
    if (m_voidTypeOptions.empty())
        return true;

    if (m_voidTypeSelections.size() != m_voidTypeOptions.size())
        return true;

    for (int optionIndex = 0; optionIndex < (int)m_voidTypeOptions.size(); ++optionIndex)
    {
        if (m_voidTypeSelections[(size_t)optionIndex] == 0)
            continue;

        const int key = m_voidTypeOptions[(size_t)optionIndex].key;
        if (key >= 1000 && v.volumeBin == (key - 1000))
            return true;
    }

    return false;
}

bool InterstitialAtomsDialog::pointInsideTransformedMesh(const glm::vec3& point) const
{
    if (!m_meshLoaded || m_meshVertices.empty() || m_meshIndices.empty() || m_meshScale <= 0.0f)
        return false;

    const glm::mat4 transform = makeMeshTransform(m_meshScale, m_meshTranslation, m_meshRotationDeg);

    const glm::mat4 inv = glm::inverse(transform);
    const glm::vec3 pLocal = glm::vec3(inv * glm::vec4(point, 1.0f));

    static const glm::vec3 rayDirs[3] = {
        glm::normalize(glm::vec3(1.0f, 0.137f, 0.319f)),
        glm::normalize(glm::vec3(0.173f, 1.0f, 0.271f)),
        glm::normalize(glm::vec3(0.259f, 0.347f, 1.0f)),
    };

    int insideVotes = 0;
    for (int r = 0; r < 3; ++r)
    {
        int hits = 0;
        for (size_t i = 0; i + 2 < m_meshIndices.size(); i += 3)
        {
            const unsigned int i0 = m_meshIndices[i];
            const unsigned int i1 = m_meshIndices[i + 1];
            const unsigned int i2 = m_meshIndices[i + 2];
            if (i0 >= m_meshVertices.size() || i1 >= m_meshVertices.size() || i2 >= m_meshVertices.size())
                continue;
            if (rayIntersectsTriangle(pLocal, rayDirs[r],
                                      m_meshVertices[i0], m_meshVertices[i1], m_meshVertices[i2]))
            {
                ++hits;
            }
        }
        if ((hits % 2) == 1)
            ++insideVotes;
    }

    return insideVotes >= 2;
}

int InterstitialAtomsDialog::applyPlacement(Structure& structure,
                                            const std::vector<glm::vec3>& elementColors)
{
    if (!m_sourceLoaded)
        return 0;

    std::vector<int> pool;
    pool.reserve(m_voids.size());

    if (m_placementMode == PlacementMode::ManualSelection)
    {
        std::vector<char> seen(m_voids.size(), 0);
        for (int idx : m_selectedVoidIndices)
        {
            if (idx < 0 || idx >= (int)m_voids.size())
                continue;
            if (seen[(size_t)idx])
                continue;
            seen[(size_t)idx] = 1;
            pool.push_back(idx);
        }
    }
    else
    {
        for (int i = 0; i < (int)m_voids.size(); ++i)
        {
            const VoidCandidate& v = m_voids[(size_t)i];
            if (!matchesSelectedVoidFilter(v))
                continue;

            if (m_placementMode == PlacementMode::MeshRegion && !pointInsideTransformedMesh(v.position))
                continue;

            pool.push_back(i);
        }
    }

    if (pool.empty())
        return 0;

    int desired = 0;
    if (m_placementMode == PlacementMode::ManualSelection)
    {
        desired = (int)pool.size();
    }
    else if (m_targetMode == TargetMode::Percent)
    {
        const float p = std::max(0.0f, std::min(100.0f, m_targetPercent));
        desired = (int)std::lround((double)pool.size() * (double)p / 100.0);
    }
    else
    {
        desired = std::max(0, m_targetCount);
    }
    desired = std::min(desired, (int)pool.size());
    if (desired <= 0)
        return 0;

    if (m_placementMode != PlacementMode::ManualSelection)
    {
        std::mt19937 rng((unsigned int)std::max(0, m_seed));
        std::shuffle(pool.begin(), pool.end(), rng);
    }

    structure = m_source;

    glm::mat3 cell(1.0f);
    glm::mat3 invCell(1.0f);
    const bool usePbc = tryMakeCellMatrices(structure, cell, invCell);

    std::vector<glm::vec3> existing;
    existing.reserve(structure.atoms.size() + (size_t)desired);
    for (const AtomSite& a : structure.atoms)
        existing.push_back(glm::vec3((float)a.x, (float)a.y, (float)a.z));

    const int z = std::max(1, std::min(118, m_selectedElement));
    const glm::vec3 c = (z >= 0 && z < (int)elementColors.size())
                      ? elementColors[(size_t)z]
                      : glm::vec3(0.8f, 0.8f, 0.2f);

    const float minInsert = std::max(0.0f, m_minInsertDistance);
    int added = 0;

    for (int n = 0; n < desired; ++n)
    {
        const VoidCandidate& v = m_voids[(size_t)pool[(size_t)n]];

        bool clash = false;
        for (const glm::vec3& ep : existing)
        {
            glm::vec3 d = v.position - ep;
            d = minimumImageDelta(d, usePbc, cell, invCell);
            if (glm::dot(d, d) < minInsert * minInsert)
            {
                clash = true;
                break;
            }
        }
        if (clash)
            continue;

        AtomSite site;
        site.atomicNumber = z;
        site.symbol = elementSymbol(z);
        site.x = (double)v.position.x;
        site.y = (double)v.position.y;
        site.z = (double)v.position.z;
        site.r = c.r;
        site.g = c.g;
        site.b = c.b;

        structure.atoms.push_back(site);
        existing.push_back(v.position);
        ++added;
    }

    return added;
}

void InterstitialAtomsDialog::drawDialog(Structure& structure,
                                         const std::vector<glm::vec3>& elementColors,
                                         const std::vector<float>& elementRadii,
                                         const std::vector<float>& elementShininess,
                                         const std::function<void(Structure&)>& updateBuffers)
{
    for (const std::string& p : m_pendingDropPaths)
    {
        if (acceptsMeshPath(p))
            loadMeshFromPath(p);
        else if (acceptsStructurePath(p))
            loadStructureFromPath(p, elementRadii, elementShininess);
    }
    m_pendingDropPaths.clear();

    if (m_openRequested)
    {
        ImGui::OpenPopup("Add Interstitial Atoms");
        m_openRequested = false;
    }

    m_isOpen = ImGui::IsPopupOpen("Add Interstitial Atoms");
    ImGui::SetNextWindowSize(ImVec2(1100.0f, 720.0f), ImGuiCond_FirstUseEver);

    bool keepOpen = true;
    if (!ImGui::BeginPopupModal("Add Interstitial Atoms", &keepOpen, ImGuiWindowFlags_NoCollapse))
    {
        m_isOpen = false;
        return;
    }
    m_isOpen = true;

    // Rebuild preview buffers if dirty
    if (m_previewBufDirty)
        rebuildPreviewBuffers(elementRadii, elementShininess);

    // ---- Layout: left preview panel | right controls, bottom bar ----
    const float kPreviewW   = 560.0f;
    const float kBottomBarH = 54.0f;
    const float contentH    = ImGui::GetContentRegionAvail().y - kBottomBarH;

    // ---- Left: GL FBO preview / drop zone ----
    ImGui::BeginChild("##interstitial-preview-panel", ImVec2(kPreviewW, contentH), false,
                      ImGuiWindowFlags_NoScrollbar);
    {
        const ImVec2 avail    = ImGui::GetContentRegionAvail();
        const int    previewW = std::max(4, (int)avail.x);
        const int    previewH = std::max(4, (int)avail.y);

        if (m_glReady && m_sourceLoaded && m_previewBuffers.atomCount > 0)
            renderPreviewToFBO(previewW, previewH);

        const ImVec2 canvasMin  = ImGui::GetCursorScreenPos();
        const ImVec2 canvasSize((float)previewW, (float)previewH);
        const ImVec2 canvasMax(canvasMin.x + canvasSize.x, canvasMin.y + canvasSize.y);

        ImGui::InvisibleButton("##interstitial-orbit-zone",
                               ImVec2((float)previewW, (float)previewH),
                               ImGuiButtonFlags_MouseButtonLeft);

        const ImGuiIO& io = ImGui::GetIO();

        // Orbit — blocked while a mesh gizmo drag is in progress
        if (ImGui::IsItemActive() && !m_meshGizmoDrag.active
            && ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f))
        {
            const ImVec2 delta = io.MouseDelta;
            m_camYaw   -= delta.x * 0.5f;
            m_camPitch += delta.y * 0.5f;
            m_camPitch  = std::max(-85.0f, std::min(85.0f, m_camPitch));
        }
        // Zoom
        if (ImGui::IsItemHovered() && io.MouseWheel != 0.0f)
        {
            m_camDistance -= io.MouseWheel * m_camDistance * 0.08f;
            m_camDistance  = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, m_camDistance));
        }

        ImDrawList* dl = ImGui::GetWindowDrawList();

        if (m_glReady && m_previewFBO != 0 && m_previewColorTex != 0
            && m_sourceLoaded && m_previewBuffers.atomCount > 0)
        {
            dl->AddImage((ImTextureID)(intptr_t)m_previewColorTex,
                         canvasMin, canvasMax,
                         ImVec2(0.0f, 1.0f), ImVec2(1.0f, 0.0f));
            drawPreviewOverlays(dl, canvasMin, canvasSize);
        }
        else
        {
            dl->AddRectFilled(canvasMin, canvasMax, ImColor(0.06f, 0.07f, 0.09f, 1.0f));
            dl->AddRect(canvasMin, canvasMax, ImColor(0.25f, 0.27f, 0.31f, 1.0f));

            const char* hint = m_sourceLoaded
                ? (m_glReady ? "No atoms to preview." : "GL resources not ready.")
                : "Drop a structure file here\nor press 'Use Current Scene'.";
            const ImVec2 textSz = ImGui::CalcTextSize(hint);
            dl->AddText(
                ImVec2(canvasMin.x + (canvasSize.x - textSz.x) * 0.5f,
                       canvasMin.y + (canvasSize.y - textSz.y) * 0.5f),
                ImColor(0.65f, 0.67f, 0.72f, 1.0f), hint);
        }
    }
    ImGui::EndChild();

    ImGui::SameLine();

    // ---- Right: controls ----
    ImGui::BeginChild("##interstitial-controls", ImVec2(0.0f, contentH), false);
    {
        // Source
        ImGui::Text("Source");
        ImGui::Separator();

        if (ImGui::Button("Use Current Scene"))
        {
            if (structure.atoms.empty())
            {
                m_statusIsError = true;
                m_statusMsg = "Current scene has no atoms.";
            }
            else
            {
                m_source        = structure;
                m_sourceLoaded  = true;
                m_sourceLabel   = "scene";
                m_voids.clear();
                m_selectedVoidIndices.clear();
                m_selectedVoidIndex = -1;
                rebuildVoidTypeOptions();
                rebuildPreviewBuffers(elementRadii, elementShininess);
                autoFitPreviewCamera();
                m_statusIsError = false;
                m_statusMsg     = "Source set from current scene.";
            }
        }
        ImGui::SameLine();
        ImGui::TextDisabled("or drop a file onto the preview.");

        if (m_sourceLoaded)
        {
            ImGui::Text("Loaded: %s", m_sourceLabel.c_str());
            ImGui::TextDisabled("Atoms: %d", (int)m_source.atoms.size());
        }
        else
        {
            ImGui::TextDisabled("No source loaded.");
        }

        ImGui::Spacing();

        // Void Detection
        ImGui::Text("Void Detection");
        ImGui::Separator();

        ImGui::SetNextItemWidth(180.0f);
        ImGui::SliderInt("Grid", &m_gridResolution, 8, 30);
        ImGui::SetNextItemWidth(180.0f);
        ImGui::InputInt("Max Voids", &m_maxVoids);
        m_maxVoids = std::max(0, m_maxVoids);
        ImGui::TextDisabled("0 = detect all voids");
        ImGui::SetNextItemWidth(180.0f);
        ImGui::DragFloat("Min Clearance (A)", &m_minClearance, 0.02f, 0.05f, 5.0f, "%.2f");
        ImGui::SetNextItemWidth(180.0f);
        ImGui::DragFloat("Min Separation (A)", &m_minSeparation, 0.02f, 0.05f, 5.0f, "%.2f");

        if (!m_sourceLoaded)
            ImGui::BeginDisabled();
        if (ImGui::Button("Detect Voids", ImVec2(140.0f, 0.0f)))
        {
            if (detectVoids())
            {
                m_statusIsError = false;
                if (m_maxVoids > 0)
                    m_statusMsg = "Detected " + std::to_string(m_voids.size()) + " void candidates (capped).";
                else
                    m_statusMsg = "Detected " + std::to_string(m_voids.size()) + " void candidates.";
            }
            else
            {
                m_statusIsError = true;
                m_statusMsg = "No voids detected. Increase grid or lower clearance threshold.";
            }
        }
        if (!m_sourceLoaded)
            ImGui::EndDisabled();

        for (size_t gi = 0; gi < m_voidTypeOptions.size(); ++gi)
            ImGui::BulletText("%s", m_voidTypeOptions[gi].label.c_str());

        ImGui::Spacing();

        // Placement
        ImGui::Text("Placement");
        ImGui::Separator();

        if (!m_voids.empty())
        {
            if (m_selectedVoidIndex < 0 || m_selectedVoidIndex >= (int)m_voids.size())
                m_selectedVoidIndex = 0;
            if (m_selectedVoidIndices.empty())
                m_selectedVoidIndices.push_back(m_selectedVoidIndex);

            const VoidCandidate& sv = m_voids[(size_t)m_selectedVoidIndex];
            const char* typeText = (sv.type == VoidType::Tetrahedral)
                                 ? "Tetrahedral"
                                 : ((sv.type == VoidType::Octahedral) ? "Octahedral" : "Irregular");
            ImGui::Text("Selected Void Info");
            ImGui::TextDisabled("Click to select one, Ctrl+click to add/remove multiple.");
            ImGui::BulletText("Index: %d / %d", m_selectedVoidIndex + 1, (int)m_voids.size());
            ImGui::BulletText("Manual selection count: %d", (int)m_selectedVoidIndices.size());
            ImGui::BulletText("Type: %s", typeText);
            ImGui::BulletText("Volume: %.3f A^3", sv.volume);
            ImGui::BulletText("Clearance: %.3f A", sv.clearance);
            ImGui::BulletText("Coordination: %d", sv.coordination);
            ImGui::BulletText("Center: (%.3f, %.3f, %.3f)", sv.position.x, sv.position.y, sv.position.z);
            ImGui::TextDisabled("Only selected void polygon is shown.");
        }

        if (m_voidTypeSelections.size() != m_voidTypeOptions.size())
            m_voidTypeSelections.assign(m_voidTypeOptions.size(), 1);

        ImGui::Text("Void Types For Distribution");
        if (!m_voidTypeOptions.empty())
        {
            if (ImGui::Button("Select All Types"))
                std::fill(m_voidTypeSelections.begin(), m_voidTypeSelections.end(), (char)1);
            ImGui::SameLine();
            if (ImGui::Button("Clear All Types"))
                std::fill(m_voidTypeSelections.begin(), m_voidTypeSelections.end(), (char)0);

            ImGui::TextDisabled("Use checkboxes to choose one or multiple void types.");
            for (size_t i = 0; i < m_voidTypeOptions.size(); ++i)
            {
                bool enabled = m_voidTypeSelections[i] != 0;
                if (ImGui::Checkbox(m_voidTypeOptions[i].label.c_str(), &enabled))
                    m_voidTypeSelections[i] = enabled ? 1 : 0;
            }
        }
        else
        {
            ImGui::TextDisabled("Detect voids to populate void types.");
        }

        bool anyVoidTypeEnabled = false;
        for (char v : m_voidTypeSelections)
            anyVoidTypeEnabled = anyVoidTypeEnabled || (v != 0);
        if (!anyVoidTypeEnabled && !m_voidTypeSelections.empty())
            ImGui::TextColored(ImVec4(0.95f, 0.42f, 0.42f, 1.0f), "No void type selected. Placement pool will be empty.");

        ImGui::Checkbox("Show Void Polyhedra in Scene", &m_showVoidOverlay);

        if (ImGui::Button("Select Element"))
            openPeriodicTable();
        std::vector<ElementSelection> pick;
        if (drawPeriodicTable(pick) && !pick.empty())
            m_selectedElement = pick.front().atomicNumber;
        ImGui::SameLine();
        ImGui::Text("%s", elementSymbol(std::max(1, std::min(118, m_selectedElement))));

        const char* placementItems[] = {
            "Random distribution",
            "Inside OBJ/STL mesh region",
            "Manual (selected voids only)"
        };
        int placement = (m_placementMode == PlacementMode::Random) ? 0
                      : (m_placementMode == PlacementMode::MeshRegion) ? 1 : 2;
        ImGui::SetNextItemWidth(280.0f);
        if (ImGui::Combo("Mode", &placement, placementItems, 3))
        {
            if (placement == 0)
                m_placementMode = PlacementMode::Random;
            else if (placement == 1)
                m_placementMode = PlacementMode::MeshRegion;
            else
                m_placementMode = PlacementMode::ManualSelection;
        }

        const bool manualMode = (m_placementMode == PlacementMode::ManualSelection);
        if (manualMode)
        {
            ImGui::TextDisabled("Manual mode inserts atoms into the currently selected void centers.");
            ImGui::TextDisabled("Selected centers: %d", (int)m_selectedVoidIndices.size());
        }
        else
        {
            const char* targetItems[] = {"By percentage", "By atom count"};
            int target = (m_targetMode == TargetMode::Percent) ? 0 : 1;
            ImGui::SetNextItemWidth(160.0f);
            if (ImGui::Combo("Target", &target, targetItems, 2))
                m_targetMode = (target == 0) ? TargetMode::Percent : TargetMode::Count;

            if (m_targetMode == TargetMode::Percent)
            {
                ImGui::SetNextItemWidth(160.0f);
                ImGui::SliderFloat("Percent", &m_targetPercent, 0.0f, 100.0f, "%.1f%%");
            }
            else
            {
                ImGui::SetNextItemWidth(160.0f);
                ImGui::InputInt("Atom Count", &m_targetCount);
                m_targetCount = std::max(0, m_targetCount);
            }
        }

        ImGui::SetNextItemWidth(160.0f);
        ImGui::InputInt("Seed", &m_seed);
        ImGui::SetNextItemWidth(180.0f);
        ImGui::DragFloat("Min Insert Distance (A)", &m_minInsertDistance, 0.02f, 0.0f, 5.0f, "%.2f");

        // Mesh Region
        if (m_placementMode == PlacementMode::MeshRegion)
        {
            ImGui::Spacing();
            ImGui::Text("3D Region (OBJ/STL)");
            ImGui::Separator();
            ImGui::TextDisabled("Drop .obj or .stl onto preview, then use the gizmo to position it.");

            if (m_meshLoaded)
                ImGui::Text("Loaded: %s", m_meshLabel.c_str());
            else
                ImGui::TextDisabled("Drop an OBJ/STL mesh onto the preview.");

            // Gizmo mode toggle
            ImGui::Text("Gizmo Mode");
            if (ImGui::RadioButton("Translate (T)", m_meshGizmoMode == MeshGizmoMode::Translate))
                m_meshGizmoMode = MeshGizmoMode::Translate;
            ImGui::SameLine();
            if (ImGui::RadioButton("Rotate (R)", m_meshGizmoMode == MeshGizmoMode::Rotate))
                m_meshGizmoMode = MeshGizmoMode::Rotate;
            ImGui::SameLine();
            if (ImGui::RadioButton("Scale (S)", m_meshGizmoMode == MeshGizmoMode::Scale))
                m_meshGizmoMode = MeshGizmoMode::Scale;

            // Keyboard shortcuts for gizmo mode
            if (!ImGui::GetIO().WantTextInput)
            {
                if (ImGui::IsKeyPressed(ImGuiKey_T, false)) m_meshGizmoMode = MeshGizmoMode::Translate;
                if (ImGui::IsKeyPressed(ImGuiKey_R, false)) m_meshGizmoMode = MeshGizmoMode::Rotate;
                if (ImGui::IsKeyPressed(ImGuiKey_S, false)) m_meshGizmoMode = MeshGizmoMode::Scale;
            }

            ImGui::Spacing();
            ImGui::TextDisabled("Fine-tune with numeric controls:");
            ImGui::SetNextItemWidth(160.0f);
            ImGui::DragFloat("Scale", &m_meshScale, 0.01f, 1e-4f, 1000.0f, "%.3f");
            ImGui::DragFloat3("Translate", &m_meshTranslation.x, 0.05f, -10000.0f, 10000.0f, "%.2f");
            ImGui::DragFloat3("Rotate (deg)", &m_meshRotationDeg.x, 0.5f, -360.0f, 360.0f, "%.1f");

            if (ImGui::Button("Reset Transform"))
            {
                m_meshTranslation = glm::vec3(0.0f);
                m_meshRotationDeg = glm::vec3(0.0f);
                m_meshScale       = 1.0f;
            }
        }
    }
    ImGui::EndChild();

    // ---- Bottom bar ----
    ImGui::Separator();

    const bool hasManualSelection = !m_selectedVoidIndices.empty();
    const bool canApply = m_sourceLoaded && !m_voids.empty() &&
                          (m_placementMode != PlacementMode::MeshRegion || m_meshLoaded) &&
                          (m_placementMode != PlacementMode::ManualSelection || hasManualSelection);
    if (!canApply)
        ImGui::BeginDisabled();

    if (ImGui::Button("Apply Interstitials", ImVec2(180.0f, 0.0f)))
    {
        Structure out;
        const int added = applyPlacement(out, elementColors);
        if (added > 0)
        {
            structure = std::move(out);
            updateBuffers(structure);

            m_source       = structure;
            m_sourceLoaded = true;
            m_previewBufDirty = true;
            rebuildPreviewBuffers(elementRadii, elementShininess);
            m_statusIsError = false;
            m_statusMsg = "Added " + std::to_string(added) + " interstitial atom(s).";
        }
        else
        {
            m_statusIsError = true;
            m_statusMsg = "No interstitials placed. Check target, type filter, or mesh transform.";
        }
    }
    if (!canApply)
        ImGui::EndDisabled();

    ImGui::SameLine();
    if (ImGui::Button("Close", ImVec2(120.0f, 0.0f)))
    {
        ImGui::CloseCurrentPopup();
        m_isOpen = false;
    }

    if (!m_statusMsg.empty())
    {
        ImGui::SameLine();
        if (m_statusIsError)
            ImGui::TextColored(ImVec4(0.95f, 0.42f, 0.42f, 1.0f), "%s", m_statusMsg.c_str());
        else
            ImGui::TextColored(ImVec4(0.50f, 0.92f, 0.58f, 1.0f), "%s", m_statusMsg.c_str());
    }

    ImGui::EndPopup();
}





void InterstitialAtomsDialog::drawVoidOverlay(ImDrawList* drawList,
                                              const glm::mat4& projection,
                                              const glm::mat4& view,
                                              int framebufferWidth,
                                              int framebufferHeight,
                                              const Structure& activeStructure) const
{
    (void)drawList;
    (void)projection;
    (void)view;
    (void)framebufferWidth;
    (void)framebufferHeight;
    (void)activeStructure;

    // Intentionally disabled: avoid drawing interstitial void centers/polygons in main scene.
    return;
}
