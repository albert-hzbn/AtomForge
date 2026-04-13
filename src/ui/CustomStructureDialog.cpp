#include "ui/CustomStructureDialog.h"

#include "algorithms/MeshLoader.h"
#include "algorithms/NanoCrystalBuilder.h"
#include "graphics/CylinderMesh.h"
#include "graphics/CustomStructureShaders.h"
#include "graphics/Renderer.h"
#include "graphics/SphereMesh.h"
#include "graphics/StructureInstanceBuilder.h"
#include "graphics/ShadowMap.h"
#include "graphics/Shader.h"
#include "app/SceneView.h"
#include "camera/Camera.h"
#include "util/PathUtils.h"

#include "imgui.h"

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>

namespace
{
void drawDropRegion(const char* id, const char* title,
                    const char* bodyLine1, const char* bodyLine2,
                    float height, bool loaded, bool highlight)
{
    ImGui::TextUnformatted(title);
    ImGui::InvisibleButton(id, ImVec2(-1.0f, height));
    ImDrawList* draw = ImGui::GetWindowDrawList();
    const ImVec2 pMin = ImGui::GetItemRectMin();
    const ImVec2 pMax = ImGui::GetItemRectMax();
    const bool hovered = ImGui::IsItemHovered();

    // Background fill
    if (loaded)
        draw->AddRectFilled(pMin, pMax, IM_COL32(30, 80, 30, 80), 4.0f);
    else if (hovered || highlight)
        draw->AddRectFilled(pMin, pMax, IM_COL32(60, 60, 100, 60), 4.0f);

    // Border
    const ImU32 borderCol = loaded
        ? IM_COL32(80, 200, 80, 200)
        : (hovered ? IM_COL32(120, 120, 200, 200) : ImGui::GetColorU32(ImGuiCol_Border));
    draw->AddRect(pMin, pMax, borderCol, 4.0f, 0, loaded ? 2.0f : 1.0f);

    // Centered text
    const float midY = (pMin.y + pMax.y) * 0.5f;
    const float midX = (pMin.x + pMax.x) * 0.5f;
    const ImVec2 ts1 = ImGui::CalcTextSize(bodyLine1);
    draw->AddText(ImVec2(midX - ts1.x * 0.5f, midY - ts1.y - 2.0f),
                  loaded ? IM_COL32(100, 230, 100, 255) : ImGui::GetColorU32(ImGuiCol_TextDisabled),
                  bodyLine1);
    if (bodyLine2 && bodyLine2[0])
    {
        const ImVec2 ts2 = ImGui::CalcTextSize(bodyLine2);
        draw->AddText(ImVec2(midX - ts2.x * 0.5f, midY + 4.0f),
                      ImGui::GetColorU32(ImGuiCol_TextDisabled), bodyLine2);
    }
}
}

// ---------------------------------------------------------------------------
// Construction / Destruction / GL init
// ---------------------------------------------------------------------------

CustomStructureDialog::CustomStructureDialog() = default;

CustomStructureDialog::~CustomStructureDialog()
{
    delete m_previewSphere;
    delete m_previewCylinder;

    if (m_refFBO)      glDeleteFramebuffers(1, &m_refFBO);
    if (m_refColorTex) glDeleteTextures(1, &m_refColorTex);
    if (m_refDepthRbo) glDeleteRenderbuffers(1, &m_refDepthRbo);

    if (m_modelFBO)      glDeleteFramebuffers(1, &m_modelFBO);
    if (m_modelColorTex) glDeleteTextures(1, &m_modelColorTex);
    if (m_modelDepthRbo) glDeleteRenderbuffers(1, &m_modelDepthRbo);

    if (m_meshProgram) glDeleteProgram(m_meshProgram);
    if (m_meshVAO) glDeleteVertexArrays(1, &m_meshVAO);
    if (m_meshVBO) glDeleteBuffers(1, &m_meshVBO);
    if (m_meshIBO) glDeleteBuffers(1, &m_meshIBO);
}

void CustomStructureDialog::initRenderResources(Renderer& renderer)
{
    m_renderer = &renderer;
    m_previewSphere   = new SphereMesh(24, 24);
    m_previewCylinder = new CylinderMesh(16);
    m_previewBuffers.init(m_previewSphere->vao, m_previewCylinder->vao);
    m_previewShadow = createShadowMap(1, 1);
    m_meshProgram = createProgram(customStructureMeshVertexShader(),
                                  customStructureMeshFragmentShader());
    glBindAttribLocation(m_meshProgram, 0, "position");
    glBindAttribLocation(m_meshProgram, 1, "normal");
    glLinkProgram(m_meshProgram);
    m_glReady = true;
}

// ---------------------------------------------------------------------------
// FBO helpers
// ---------------------------------------------------------------------------

void CustomStructureDialog::ensureRefFBO(int w, int h)
{
    if (w == m_refFboW && h == m_refFboH && m_refFBO != 0)
        return;

    if (m_refFBO)      { glDeleteFramebuffers(1, &m_refFBO);      m_refFBO = 0; }
    if (m_refColorTex) { glDeleteTextures(1, &m_refColorTex);     m_refColorTex = 0; }
    if (m_refDepthRbo) { glDeleteRenderbuffers(1, &m_refDepthRbo); m_refDepthRbo = 0; }

    glGenTextures(1, &m_refColorTex);
    glBindTexture(GL_TEXTURE_2D, m_refColorTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenRenderbuffers(1, &m_refDepthRbo);
    glBindRenderbuffer(GL_RENDERBUFFER, m_refDepthRbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    glGenFramebuffers(1, &m_refFBO);
    glBindFramebuffer(GL_FRAMEBUFFER, m_refFBO);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                           GL_TEXTURE_2D, m_refColorTex, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER, m_refDepthRbo);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    m_refFboW = w;
    m_refFboH = h;
}

void CustomStructureDialog::ensureModelFBO(int w, int h)
{
    if (w == m_modelFboW && h == m_modelFboH && m_modelFBO != 0)
        return;

    if (m_modelFBO)      { glDeleteFramebuffers(1, &m_modelFBO);      m_modelFBO = 0; }
    if (m_modelColorTex) { glDeleteTextures(1, &m_modelColorTex);     m_modelColorTex = 0; }
    if (m_modelDepthRbo) { glDeleteRenderbuffers(1, &m_modelDepthRbo); m_modelDepthRbo = 0; }

    glGenTextures(1, &m_modelColorTex);
    glBindTexture(GL_TEXTURE_2D, m_modelColorTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenRenderbuffers(1, &m_modelDepthRbo);
    glBindRenderbuffer(GL_RENDERBUFFER, m_modelDepthRbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    glGenFramebuffers(1, &m_modelFBO);
    glBindFramebuffer(GL_FRAMEBUFFER, m_modelFBO);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                           GL_TEXTURE_2D, m_modelColorTex, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER, m_modelDepthRbo);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    m_modelFboW = w;
    m_modelFboH = h;
}

// ---------------------------------------------------------------------------
// Structure preview
// ---------------------------------------------------------------------------

void CustomStructureDialog::rebuildPreviewBuffers(const std::vector<float>& radii,
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

void CustomStructureDialog::autoFitRefCamera()
{
    m_refCamYaw   = 45.0f;
    m_refCamPitch = 35.0f;

    if (m_previewBuffers.atomCount == 0) {
        m_refCamDistance = 10.0f;
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
    m_refCamDistance = dist;
}

void CustomStructureDialog::renderRefPreviewToFBO(int w, int h)
{
    if (!m_glReady || !m_renderer || m_previewBuffers.atomCount == 0)
        return;

    ensureRefFBO(w, h);

    GLint prevFbo = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
    GLint prevVP[4];
    glGetIntegerv(GL_VIEWPORT, prevVP);

    Camera cam;
    cam.yaw      = m_refCamYaw;
    cam.pitch    = m_refCamPitch;
    cam.distance = m_refCamDistance;

    FrameView frame;
    frame.framebufferWidth  = w;
    frame.framebufferHeight = h;
    buildFrameView(cam, m_previewBuffers, true, frame);

    glBindFramebuffer(GL_FRAMEBUFFER, m_refFBO);
    glViewport(0, 0, w, h);
    glEnable(GL_DEPTH_TEST);
    {
        const ImVec4& bg = ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
        glClearColor(bg.x, bg.y, bg.z, bg.w);
    }
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

    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
    glViewport(prevVP[0], prevVP[1], prevVP[2], prevVP[3]);
}

// ---------------------------------------------------------------------------
// Model wireframe preview
// ---------------------------------------------------------------------------

void CustomStructureDialog::rebuildMeshSurface()
{
    if (m_modelVertices.empty() || m_modelIndices.empty())
        return;

    // Build interleaved position+normal buffer with per-face normals
    struct Vertex { glm::vec3 pos; glm::vec3 norm; };
    std::vector<Vertex> verts;
    verts.reserve(m_modelIndices.size());

    for (size_t i = 0; i + 2 < m_modelIndices.size(); i += 3)
    {
        const glm::vec3& a = m_modelVertices[m_modelIndices[i]];
        const glm::vec3& b = m_modelVertices[m_modelIndices[i + 1]];
        const glm::vec3& c = m_modelVertices[m_modelIndices[i + 2]];
        glm::vec3 n = glm::normalize(glm::cross(b - a, c - a));
        verts.push_back({a, n});
        verts.push_back({b, n});
        verts.push_back({c, n});
    }
    m_meshTriCount = (int)(verts.size() / 3);

    if (!m_meshVAO)
    {
        glGenVertexArrays(1, &m_meshVAO);
        glGenBuffers(1, &m_meshVBO);
    }

    glBindVertexArray(m_meshVAO);
    glBindBuffer(GL_ARRAY_BUFFER, m_meshVBO);
    glBufferData(GL_ARRAY_BUFFER,
                 (GLsizeiptr)(verts.size() * sizeof(Vertex)),
                 verts.data(), GL_STATIC_DRAW);
    // position
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                          (void*)offsetof(Vertex, pos));
    // normal
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                          (void*)offsetof(Vertex, norm));
    glBindVertexArray(0);
}

void CustomStructureDialog::autoFitModelCamera()
{
    m_modelCamYaw   = 45.0f;
    m_modelCamPitch = 35.0f;

    float maxR = glm::length(m_modelHalfExtents);
    maxR = std::max(maxR, 1.0f);

    const float halfFov = glm::radians(22.5f);
    float dist = maxR / std::sin(halfFov) * 1.15f;
    dist = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, dist));
    m_modelCamDistance = dist;
}

void CustomStructureDialog::renderModelPreviewToFBO(int w, int h)
{
    if (!m_glReady || !m_renderer || m_meshTriCount == 0)
        return;

    ensureModelFBO(w, h);

    GLint prevFbo = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
    GLint prevVP[4];
    glGetIntegerv(GL_VIEWPORT, prevVP);

    // Build view/projection manually for the model
    float aspect = (h > 0) ? (float)w / (float)h : 1.0f;
    const float verticalFov = glm::radians(45.0f);
    const float sceneRadius = std::max(1.0f, glm::length(m_modelHalfExtents));
    const float depthPadding = std::max(10.0f, sceneRadius * 0.25f);
    const float halfHeight = std::max(0.1f, m_modelCamDistance * std::tan(verticalFov * 0.5f));
    const float halfWidth = halfHeight * aspect;
    const float depthRange = std::max(1000.0f, m_modelCamDistance + sceneRadius + depthPadding);
    glm::mat4 projection = glm::ortho(
        -halfWidth,
        halfWidth,
        -halfHeight,
        halfHeight,
        -depthRange,
        depthRange);

    float yawRad   = glm::radians(m_modelCamYaw);
    float pitchRad = glm::radians(m_modelCamPitch);
    glm::vec3 eye(
        m_modelCamDistance * std::cos(pitchRad) * std::cos(yawRad),
        m_modelCamDistance * std::sin(pitchRad),
        m_modelCamDistance * std::cos(pitchRad) * std::sin(yawRad)
    );
    glm::mat4 view = glm::lookAt(eye, glm::vec3(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));

    glBindFramebuffer(GL_FRAMEBUFFER, m_modelFBO);
    glViewport(0, 0, w, h);
    glEnable(GL_DEPTH_TEST);
    {
        const ImVec4& bg = ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
        glClearColor(bg.x, bg.y, bg.z, bg.w);
    }
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glUseProgram(m_meshProgram);
    glUniformMatrix4fv(glGetUniformLocation(m_meshProgram, "projection"),
                       1, GL_FALSE, &projection[0][0]);
    glUniformMatrix4fv(glGetUniformLocation(m_meshProgram, "view"),
                       1, GL_FALSE, &view[0][0]);
    glm::vec3 lightDir = glm::normalize(eye + glm::vec3(0.3f, 0.8f, 0.2f));
    glUniform3fv(glGetUniformLocation(m_meshProgram, "uLightDir"), 1, &lightDir[0]);
    glUniform3fv(glGetUniformLocation(m_meshProgram, "uViewPos"), 1, &eye[0]);
    glUniform3f(glGetUniformLocation(m_meshProgram, "uColor"), 0.45f, 0.65f, 0.85f);

    glBindVertexArray(m_meshVAO);
    glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(m_meshTriCount * 3));
    glBindVertexArray(0);
    glUseProgram(0);

    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
    glViewport(prevVP[0], prevVP[1], prevVP[2], prevVP[3]);
}

void CustomStructureDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Custom Structure", NULL, false, enabled))
        m_openRequested = true;
}

void CustomStructureDialog::feedDroppedFile(const std::string& path)
{
    m_pendingDropPaths.push_back(path);
    std::cout << "[CustomStructure] File queued: " << path << std::endl;
}

bool CustomStructureDialog::loadReferenceStructure(const std::string& path)
{
    Structure loaded;
    std::string error;
    if (!loadStructureFromFile(path, loaded, error))
    {
        m_status = std::string("Failed to load crystal structure: ") + error;
        std::cout << "[CustomStructure] " << m_status << std::endl;
        return false;
    }

    m_reference = std::move(loaded);
    const size_t slash = path.find_last_of("/\\");
    m_refFileName = (slash == std::string::npos) ? path : path.substr(slash + 1);
    m_status = std::string("Crystal reference loaded: ") + m_refFileName
             + " (atoms=" + std::to_string(m_reference.atoms.size())
             + ", unitCell=" + (m_reference.hasUnitCell ? "yes" : "no") + ")";
    std::cout << "[CustomStructure] " << m_status << std::endl;
    m_previewBufDirty = true;
    return true;
}

bool CustomStructureDialog::loadModelMesh(const std::string& path)
{
    std::vector<glm::vec3> verts;
    std::vector<unsigned int> idx;
    std::string error;

    bool ok = false;
    if (hasExtension(toLowerStr(path), ".obj"))
        ok = parseObjMesh(path, verts, idx, error);
    else if (hasExtension(toLowerStr(path), ".stl"))
        ok = parseStlMesh(path, verts, idx, error);

    if (!ok)
    {
        if (error.empty())
            error = "Unsupported model format. Use OBJ or STL.";
        m_status = std::string("Failed to load model: ") + error;
        std::cout << "[CustomStructure] " << m_status << std::endl;
        return false;
    }

    glm::vec3 minP(1e30f), maxP(-1e30f);
    for (const glm::vec3& v : verts)
    {
        minP = glm::min(minP, v);
        maxP = glm::max(maxP, v);
    }
    const glm::vec3 center = 0.5f * (minP + maxP);
    m_modelHalfExtents = 0.5f * (maxP - minP);

    for (size_t i = 0; i < verts.size(); ++i)
        verts[i] -= center;

    m_modelVertices.swap(verts);
    m_modelIndices.swap(idx);

    const size_t slash = path.find_last_of("/\\");
    m_modelName = (slash == std::string::npos) ? path : path.substr(slash + 1);
    m_status = std::string("Model loaded: ") + m_modelName
             + " (triangles=" + std::to_string(m_modelIndices.size() / 3)
             + ", bbox=" + std::to_string((int)(2.0f * m_modelHalfExtents.x))
             + "x" + std::to_string((int)(2.0f * m_modelHalfExtents.y))
             + "x" + std::to_string((int)(2.0f * m_modelHalfExtents.z)) + ")";
    std::cout << "[CustomStructure] " << m_status << std::endl;
    rebuildMeshSurface();
    autoFitModelCamera();
    return true;
}

void CustomStructureDialog::drawDialog(Structure& structure,
                                         const std::vector<glm::vec3>& elementColors,
                                         const std::vector<float>& elementRadii,
                                         const std::vector<float>& elementShininess,
                                         const std::function<void(Structure&)>& updateBuffers)
{
    static NanoParams params;
    static NanoBuildResult lastResult;

    params.shape = NanoShape::MeshModel;
    if (params.modelScale <= 0.0f)
        params.modelScale = 1.0f;

    // Process all queued file drops (supports multi-file drop).
    bool justDroppedRef = false;
    bool justDroppedModel = false;
    for (const auto& dropPath : m_pendingDropPaths)
    {
        if (hasExtension(toLowerStr(dropPath), ".obj") || hasExtension(toLowerStr(dropPath), ".stl"))
        {
            if (loadModelMesh(dropPath))
                justDroppedModel = true;
        }
        else
        {
            if (loadReferenceStructure(dropPath))
                justDroppedRef = true;
        }
    }
    m_pendingDropPaths.clear();

    if (m_openRequested)
    {
        ImGui::OpenPopup("Custom Structure");
        m_openRequested = false;
        lastResult = NanoBuildResult();
        m_status.clear();
    }

    m_isOpen = ImGui::IsPopupOpen("Custom Structure");

    ImGui::SetNextWindowSize(ImVec2(1080.0f, 760.0f), ImGuiCond_FirstUseEver);
    bool keepOpen = true;
    if (!ImGui::BeginPopupModal("Custom Structure", &keepOpen, 0))
    {
        m_isOpen = false;
        return;
    }

    m_isOpen = true;

    ImGui::TextWrapped("Drag-and-drop files from your file manager onto this window. "
                       "Crystal structure files (CIF, XYZ, VASP, PDB, ...) go to the left panel. "
                       "3D model files (OBJ, STL) go to the right panel automatically.");
    ImGui::Separator();

    const float panelH = 340.0f;
    const float panelW = (ImGui::GetContentRegionAvail().x - 10.0f) * 0.5f;
    const float previewH = 220.0f;

    // ========== LEFT PANEL: Crystal Structure ==========
    ImGui::BeginChild("##scf-left", ImVec2(panelW, panelH), true);
    {
        const bool refLoaded = !m_reference.atoms.empty();
        ImGui::TextUnformatted("Crystal Structure");

        if (refLoaded && m_glReady)
        {
            if (m_previewBufDirty)
            {
                rebuildPreviewBuffers(elementRadii, elementShininess);
                autoFitRefCamera();
            }

            const ImVec2 avail = ImGui::GetContentRegionAvail();
            const float pw = avail.x;
            const float ph = std::min(avail.y - 60.0f, previewH);
            const int iw = std::max(1, (int)pw);
            const int ih = std::max(1, (int)ph);

            renderRefPreviewToFBO(iw, ih);

            ImVec2 cursor = ImGui::GetCursorScreenPos();
            ImGui::InvisibleButton("##scf-refpreview", ImVec2(pw, ph));
            bool previewHovered = ImGui::IsItemHovered();
            bool previewActive  = ImGui::IsItemActive();

            ImDrawList* dl = ImGui::GetWindowDrawList();
            dl->AddImage((ImTextureID)(intptr_t)m_refColorTex,
                         cursor, ImVec2(cursor.x + pw, cursor.y + ph),
                         ImVec2(0, 1), ImVec2(1, 0));

            if (previewActive && ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f))
            {
                ImVec2 delta = ImGui::GetIO().MouseDelta;
                m_refCamYaw   -= delta.x * 0.5f;
                m_refCamPitch += delta.y * 0.5f;
            }
            if (previewHovered)
            {
                float wheel = ImGui::GetIO().MouseWheel;
                if (wheel != 0.0f)
                {
                    m_refCamDistance -= wheel * m_refCamDistance * 0.1f;
                    m_refCamDistance = std::max(Camera::kMinDistance,
                                               std::min(Camera::kMaxDistance, m_refCamDistance));
                }
            }

            ImGui::TextDisabled("%s  |  %d atoms  |  cell: %s",
                                m_refFileName.c_str(),
                                (int)m_reference.atoms.size(),
                                m_reference.hasUnitCell ? "yes" : "no");
            if (m_reference.hasUnitCell)
            {
                const auto& cv = m_reference.cellVectors;
                float la = std::sqrt((float)(cv[0][0]*cv[0][0]+cv[0][1]*cv[0][1]+cv[0][2]*cv[0][2]));
                float lb = std::sqrt((float)(cv[1][0]*cv[1][0]+cv[1][1]*cv[1][1]+cv[1][2]*cv[1][2]));
                float lc = std::sqrt((float)(cv[2][0]*cv[2][0]+cv[2][1]*cv[2][1]+cv[2][2]*cv[2][2]));
                ImGui::TextDisabled("Cell: a=%.2f  b=%.2f  c=%.2f A", la, lb, lc);
            }
        }
        else
        {
            drawDropRegion("##scf-drop-structure",
                           "",
                           refLoaded ? m_refFileName.c_str() : "Drop crystal structure file here",
                           refLoaded ? "" : "(CIF, XYZ, PDB, VASP, ...)",
                           previewH, refLoaded, justDroppedRef);
        }

        if (refLoaded)
        {
            if (ImGui::SmallButton("Clear##ref"))
            {
                m_reference = Structure();
                m_refFileName.clear();
                m_previewBufDirty = true;
                m_status = "Crystal reference cleared.";
            }
            ImGui::SameLine();
        }
        if (ImGui::SmallButton("Use Current Scene##ref"))
        {
            if (!structure.atoms.empty() && structure.hasUnitCell)
            {
                m_reference = structure;
                m_refFileName = "scene (unit cell)";
                m_previewBufDirty = true;
                m_status = "Current scene copied as crystal reference.";
            }
            else if (!structure.atoms.empty())
            {
                m_reference = structure;
                m_refFileName = "scene (no unit cell)";
                m_previewBufDirty = true;
                m_status = "Warning: scene has no unit cell - tiling disabled.";
            }
            else
            {
                m_status = "Scene is empty. Load a structure first.";
            }
        }
    }
    ImGui::EndChild();

    ImGui::SameLine();

    // ========== RIGHT PANEL: 3D Model ==========
    ImGui::BeginChild("##scf-right", ImVec2(0.0f, panelH), true);
    {
        const bool modelLoaded = !m_modelVertices.empty();
        ImGui::TextUnformatted("3D Object");

        if (modelLoaded && m_glReady && m_meshTriCount > 0)
        {
            const ImVec2 avail = ImGui::GetContentRegionAvail();
            const float pw = avail.x;
            const float ph = std::min(avail.y - 60.0f, previewH);
            const int iw = std::max(1, (int)pw);
            const int ih = std::max(1, (int)ph);

            renderModelPreviewToFBO(iw, ih);

            ImVec2 cursor = ImGui::GetCursorScreenPos();
            ImGui::InvisibleButton("##scf-modelpreview", ImVec2(pw, ph));
            bool previewHovered = ImGui::IsItemHovered();
            bool previewActive  = ImGui::IsItemActive();

            ImDrawList* dl = ImGui::GetWindowDrawList();
            dl->AddImage((ImTextureID)(intptr_t)m_modelColorTex,
                         cursor, ImVec2(cursor.x + pw, cursor.y + ph),
                         ImVec2(0, 1), ImVec2(1, 0));

            if (previewActive && ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f))
            {
                ImVec2 delta = ImGui::GetIO().MouseDelta;
                m_modelCamYaw   += delta.x * 0.5f;
                m_modelCamPitch += delta.y * 0.5f;
            }
            if (previewHovered)
            {
                float wheel = ImGui::GetIO().MouseWheel;
                if (wheel != 0.0f)
                {
                    m_modelCamDistance -= wheel * m_modelCamDistance * 0.1f;
                    m_modelCamDistance = std::max(Camera::kMinDistance,
                                                 std::min(Camera::kMaxDistance, m_modelCamDistance));
                }
            }

            ImGui::TextDisabled("%s  |  %d triangles",
                                m_modelName.c_str(),
                                (int)(m_modelIndices.size() / 3));
            ImGui::TextDisabled("Model size: %.1f x %.1f x %.1f (model units)",
                                2.0f * m_modelHalfExtents.x,
                                2.0f * m_modelHalfExtents.y,
                                2.0f * m_modelHalfExtents.z);

            const float sx = 2.0f * m_modelHalfExtents.x * params.modelScale;
            const float sy = 2.0f * m_modelHalfExtents.y * params.modelScale;
            const float sz = 2.0f * m_modelHalfExtents.z * params.modelScale;
            ImGui::TextColored(ImVec4(0.5f, 0.9f, 1.0f, 1.0f),
                               "= %.1f x %.1f x %.1f Angstrom", sx, sy, sz);
        }
        else
        {
            drawDropRegion("##scf-drop-model",
                           "",
                           modelLoaded ? m_modelName.c_str() : "Drop 3D model file here",
                           modelLoaded ? "" : "(OBJ or STL)",
                           previewH, modelLoaded, justDroppedModel);
        }

        if (modelLoaded)
        {
            if (ImGui::SmallButton("Clear##model"))
            {
                m_modelVertices.clear();
                m_modelIndices.clear();
                m_modelName.clear();
                m_modelHalfExtents = glm::vec3(15.0f);
                m_meshTriCount = 0;
                m_status = "Model cleared.";
            }
        }
    }
    ImGui::EndChild();

    ImGui::Spacing();
    ImGui::SeparatorText("Fill Options");

    // ========== Model Scale with Angstrom reference ==========
    ImGui::SetNextItemWidth(140.0f);
    ImGui::InputFloat("Scale (Angstrom per model unit)", &params.modelScale, 0.1f, 1.0f, "%.3f");
    if (params.modelScale <= 1e-4f)
        params.modelScale = 1e-4f;
    if (!m_modelVertices.empty())
    {
        const float sx = 2.0f * m_modelHalfExtents.x * params.modelScale;
        const float sy = 2.0f * m_modelHalfExtents.y * params.modelScale;
        const float sz = 2.0f * m_modelHalfExtents.z * params.modelScale;
        ImGui::SameLine();
        ImGui::TextDisabled("-> %.1f x %.1f x %.1f A", sx, sy, sz);
    }

    ImGui::Checkbox("Auto-center from reference atoms", &params.autoCenterFromAtoms);
    if (!params.autoCenterFromAtoms)
    {
        ImGui::SetNextItemWidth(110.0f);
        ImGui::InputFloat("Center X", &params.cx, 0.0f, 0.0f, "%.3f");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(110.0f);
        ImGui::InputFloat("Center Y", &params.cy, 0.0f, 0.0f, "%.3f");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(110.0f);
        ImGui::InputFloat("Center Z", &params.cz, 0.0f, 0.0f, "%.3f");
    }

    if (m_reference.hasUnitCell)
    {
        ImGui::Checkbox("Auto-replicate reference cell", &params.autoReplicate);
        if (!params.autoReplicate)
        {
            ImGui::SetNextItemWidth(80.0f);
            ImGui::InputInt("Reps a", &params.repA);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(80.0f);
            ImGui::InputInt("Reps b", &params.repB);
            ImGui::SameLine();
            ImGui::SetNextItemWidth(80.0f);
            ImGui::InputInt("Reps c", &params.repC);
            params.repA = std::max(1, params.repA);
            params.repB = std::max(1, params.repB);
            params.repC = std::max(1, params.repC);
        }
    }

    ImGui::Checkbox("Set rectangular output cell", &params.setOutputCell);
    if (params.setOutputCell)
    {
        ImGui::SameLine();
        ImGui::SetNextItemWidth(120.0f);
        ImGui::InputFloat("Vacuum (A)", &params.vacuumPadding, 0.0f, 0.0f, "%.2f");
        if (params.vacuumPadding < 0.0f)
            params.vacuumPadding = 0.0f;
    }

    params.modelHx = m_modelHalfExtents.x;
    params.modelHy = m_modelHalfExtents.y;
    params.modelHz = m_modelHalfExtents.z;

    ImGui::Spacing();

    // ========== Build / Close ==========
    const bool canBuild = !m_reference.atoms.empty() && !m_modelVertices.empty() && !m_modelIndices.empty();
    if (!canBuild)
    {
        ImGui::BeginDisabled();
        if (m_reference.atoms.empty() && m_modelVertices.empty())
            ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.2f, 1.0f),
                               "Drop a crystal structure and a 3D model to begin.");
        else if (m_reference.atoms.empty())
            ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.2f, 1.0f),
                               "Drop a crystal structure file (CIF, XYZ, ...) to continue.");
        else
            ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.2f, 1.0f),
                               "Drop a 3D model (OBJ/STL) to continue.");
    }
    if (ImGui::Button("Build Fill", ImVec2(140.0f, 0.0f)))
    {
        std::cout << "[CustomStructure] Building: refAtoms=" << m_reference.atoms.size()
                  << ", unitCell=" << (m_reference.hasUnitCell ? "yes" : "no")
                  << ", modelTris=" << (m_modelIndices.size() / 3)
                  << ", modelScale=" << params.modelScale
                  << ", fillVolume=" << (2.0f * m_modelHalfExtents.x * params.modelScale)
                  << "x" << (2.0f * m_modelHalfExtents.y * params.modelScale)
                  << "x" << (2.0f * m_modelHalfExtents.z * params.modelScale)
                  << " A" << std::endl;

        lastResult = buildNanocrystal(structure,
                                      m_reference,
                                      params,
                                      elementColors,
                                      m_modelVertices,
                                      m_modelIndices);
        if (lastResult.success)
        {
            updateBuffers(structure);
            m_status = lastResult.message;
        }
        else
        {
            m_status = lastResult.message;
        }
        std::cout << "[CustomStructure] Result: " << m_status << std::endl;
    }
    if (!canBuild)
        ImGui::EndDisabled();

    ImGui::SameLine();
    if (ImGui::Button("Close", ImVec2(100.0f, 0.0f)))
        ImGui::CloseCurrentPopup();

    if (!m_status.empty())
    {
        ImGui::Spacing();
        ImGui::TextWrapped("Status: %s", m_status.c_str());
    }

    ImGui::EndPopup();
}
