#include "ui/DislocationBuilderDialog.h"

#include "app/SceneView.h"
#include "camera/Camera.h"
#include "graphics/CylinderMesh.h"
#include "graphics/Renderer.h"
#include "graphics/SphereMesh.h"
#include "graphics/StructureInstanceBuilder.h"
#include "io/StructureLoader.h"

#include "imgui.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <future>

#include <glm/gtc/type_ptr.hpp>

namespace
{
struct CharacterOption
{
    const char* label;
    DislocationCharacter value;
};

struct ShapeOption
{
    const char* label;
    DislocationShape value;
};

constexpr CharacterOption kCharacterOptions[] = {
    {"Edge", DislocationCharacter::Edge},
    {"Screw", DislocationCharacter::Screw},
    {"Mixed", DislocationCharacter::Mixed},
};

constexpr ShapeOption kShapeOptions[] = {
    {"Half-plane", DislocationShape::HalfPlane},
    {"Cylinder", DislocationShape::Cylinder},
    {"Sphere", DislocationShape::Sphere},
    {"Ellipsoid", DislocationShape::Ellipsoid},
    {"Freeform 2D", DislocationShape::Freeform2D},
};

const char* characterLabel(DislocationCharacter character)
{
    for (const CharacterOption& option : kCharacterOptions)
    {
        if (option.value == character)
            return option.label;
    }
    return "Edge";
}

const char* shapeLabel(DislocationShape shape)
{
    for (const ShapeOption& option : kShapeOptions)
    {
        if (option.value == shape)
            return option.label;
    }
    return "Half-plane";
}

ImVec4 validationStateColor(const DislocationResult& result)
{
    if (!result.success)
        return ImVec4(0.70f, 0.70f, 0.70f, 1.0f);
    if (result.validation.passed)
        return ImVec4(0.50f, 0.90f, 0.50f, 1.0f);
    return ImVec4(1.00f, 0.75f, 0.35f, 1.0f);
}

const char* validationStateLabel(const DislocationResult& result)
{
    if (!result.success)
        return "Not Generated";
    return result.validation.passed ? "Validated" : "Validation Warning";
}

void drawPreviewHint(ImDrawList* drawList,
                     const ImVec2& minCorner,
                     const ImVec2& maxCorner,
                     const char* line1,
                     const char* line2)
{
    const float lineHeight = ImGui::GetTextLineHeight();
    const ImVec2 mid((minCorner.x + maxCorner.x) * 0.5f,
                     (minCorner.y + maxCorner.y) * 0.5f);
    const float w1 = ImGui::CalcTextSize(line1).x;
    const float w2 = ImGui::CalcTextSize(line2).x;
    const ImU32 color = ImGui::GetColorU32(ImGuiCol_TextDisabled);
    drawList->AddText(ImVec2(mid.x - w1 * 0.5f, mid.y - lineHeight * 1.1f), color, line1);
    drawList->AddText(ImVec2(mid.x - w2 * 0.5f, mid.y + lineHeight * 0.1f), color, line2);
}
}

DislocationBuilderDialog::DislocationBuilderDialog() = default;

DislocationBuilderDialog::~DislocationBuilderDialog()
{
    if (m_generationRunning && m_generationFuture.valid())
        m_generationFuture.wait();

    releasePreview(m_sourcePreview);
    releasePreview(m_outputPreview);
}

void DislocationBuilderDialog::initRenderResources(Renderer& renderer)
{
    m_renderer = &renderer;
    m_sourcePreview.sphere = new SphereMesh(24, 24);
    m_sourcePreview.cylinder = new CylinderMesh(16);
    m_outputPreview.sphere = new SphereMesh(24, 24);
    m_outputPreview.cylinder = new CylinderMesh(16);

    m_sourcePreview.buffers.init(m_sourcePreview.sphere->vbo,
                                 m_sourcePreview.sphere->ebo,
                                 m_sourcePreview.sphere->indexCount,
                                 m_sourcePreview.cylinder->vbo,
                                 m_sourcePreview.cylinder->vertexCount);
    m_outputPreview.buffers.init(m_outputPreview.sphere->vbo,
                                 m_outputPreview.sphere->ebo,
                                 m_outputPreview.sphere->indexCount,
                                 m_outputPreview.cylinder->vbo,
                                 m_outputPreview.cylinder->vertexCount);

    m_sourcePreview.shadow = createShadowMap(1, 1);
    m_outputPreview.shadow = createShadowMap(1, 1);
    m_glReady = true;
}

void DislocationBuilderDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Insert Dislocation", nullptr, false, enabled))
        m_openRequested = true;
}

void DislocationBuilderDialog::feedDroppedFile(const std::string& path)
{
    m_pendingDropPath = path;
}

bool DislocationBuilderDialog::tryLoadFile(const std::string& path,
                                           const std::vector<float>& radii,
                                           const std::vector<float>& shininess)
{
    Structure loaded;
    std::string error;
    if (!loadStructureFromFile(path, loaded, error))
    {
        m_statusMsg = std::string("Load failed: ") + (error.empty() ? path : error);
        m_statusIsError = true;
        return false;
    }

    if (!loaded.hasUnitCell)
    {
        m_statusMsg = "Dislocation tool requires a structure with a unit cell.";
        m_statusIsError = true;
        return false;
    }

    m_source = std::move(loaded);
    m_sourceLoaded = true;
    m_useCurrentSceneSource = false;

    const std::string::size_type slash = path.find_last_of("\\/");
    m_sourceLabel = (slash == std::string::npos) ? path : path.substr(slash + 1);

    analyzeSource(m_source);
    m_sourcePreview.dirty = true;
    rebuildPreviewBuffers(m_sourcePreview, m_source, radii, shininess);
    autoFitPreviewCamera(m_sourcePreview);

    m_result = {};
    m_outputPreview.dirty = true;
    return true;
}

void DislocationBuilderDialog::loadFromScene(const Structure& scene,
                                             const std::vector<float>& radii,
                                             const std::vector<float>& shininess)
{
    if (scene.atoms.empty())
    {
        m_statusMsg = "Current scene has no atoms.";
        m_statusIsError = true;
        return;
    }

    if (!scene.hasUnitCell)
    {
        m_statusMsg = "Current scene has no unit cell.";
        m_statusIsError = true;
        return;
    }

    m_sourceLoaded = true;
    m_useCurrentSceneSource = true;
    m_sourceLabel = "scene";

    analyzeSource(scene);
    rebuildPreviewBuffers(m_sourcePreview, scene, radii, shininess);
    autoFitPreviewCamera(m_sourcePreview);

    m_result = {};
    m_outputPreview.dirty = true;
}

void DislocationBuilderDialog::analyzeSource(const Structure& source)
{
    m_detection = detectDislocationLattice(source, true);
    m_statusMsg = m_detection.message;
    m_statusIsError = !m_detection.success;
}

void DislocationBuilderDialog::generateDislocation(const Structure& source)
{
    if (!m_sourceLoaded)
    {
        m_statusMsg = "Load a structure first.";
        m_statusIsError = true;
        return;
    }

    if (m_generationRunning)
        return;

    const Structure sourceSnapshot = source;
    const DislocationParams paramsSnapshot = m_params;

    m_generationFuture = std::async(std::launch::async, [sourceSnapshot, paramsSnapshot]() {
        return buildDislocation(sourceSnapshot, paramsSnapshot);
    });
    m_generationRunning = true;
    m_statusMsg = "Generating dislocation...";
    m_statusIsError = false;
}

void DislocationBuilderDialog::pollGenerationResult(const std::vector<float>& radii,
                                                    const std::vector<float>& shininess)
{
    if (!m_generationRunning)
        return;

    const auto status = m_generationFuture.wait_for(std::chrono::milliseconds(0));
    if (status != std::future_status::ready)
        return;

    m_result = m_generationFuture.get();
    m_generationRunning = false;
    m_statusMsg = m_result.message;
    m_statusIsError = !m_result.success;

    if (m_result.success)
    {
        rebuildPreviewBuffers(m_outputPreview, m_result.output, radii, shininess, &m_result);
        autoFitPreviewCamera(m_outputPreview);
    }
    else
    {
        m_outputPreview.dirty = true;
    }
}

void DislocationBuilderDialog::releasePreview(PreviewState& preview)
{
    if (preview.fbo) glDeleteFramebuffers(1, &preview.fbo);
    if (preview.colorTex) glDeleteTextures(1, &preview.colorTex);
    if (preview.depthRbo) glDeleteRenderbuffers(1, &preview.depthRbo);
    if (preview.shadow.depthFBO) glDeleteFramebuffers(1, &preview.shadow.depthFBO);
    if (preview.shadow.depthTexture) glDeleteTextures(1, &preview.shadow.depthTexture);
    if (preview.dislocationLineVAO) glDeleteVertexArrays(1, &preview.dislocationLineVAO);
    if (preview.dislocationLineVBO) glDeleteBuffers(1, &preview.dislocationLineVBO);
    delete preview.sphere;
    delete preview.cylinder;
    preview = {};
}

void DislocationBuilderDialog::ensurePreviewFBO(PreviewState& preview, int width, int height)
{
    if (preview.fbo != 0 && preview.width == width && preview.height == height)
        return;

    if (preview.fbo) glDeleteFramebuffers(1, &preview.fbo);
    if (preview.colorTex) glDeleteTextures(1, &preview.colorTex);
    if (preview.depthRbo) glDeleteRenderbuffers(1, &preview.depthRbo);

    glGenTextures(1, &preview.colorTex);
    glBindTexture(GL_TEXTURE_2D, preview.colorTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenRenderbuffers(1, &preview.depthRbo);
    glBindRenderbuffer(GL_RENDERBUFFER, preview.depthRbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width, height);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    glGenFramebuffers(1, &preview.fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, preview.fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                           GL_TEXTURE_2D, preview.colorTex, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER, preview.depthRbo);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    preview.width = width;
    preview.height = height;
}

void DislocationBuilderDialog::rebuildPreviewBuffers(PreviewState& preview,
                                                     const Structure& structure,
                                                     const std::vector<float>& radii,
                                                     const std::vector<float>& shininess,
                                                     const DislocationResult* dislocationOverlay)
{
    if (!m_glReady || structure.atoms.empty())
        return;

    static const int kIdentity[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    StructureInstanceData data = buildStructureInstanceData(structure,
                                                            false,
                                                            kIdentity,
                                                            radii,
                                                            shininess);

    preview.hasDislocationLine = false;
    if (dislocationOverlay != nullptr && dislocationOverlay->success && !dislocationOverlay->loopPoints.empty())
    {
        preview.hasDislocationLine = true;

        // Update dislocation loop VAO/VBO with the loop points
        if (preview.dislocationLineVAO == 0)
        {
            glGenVertexArrays(1, &preview.dislocationLineVAO);
            glGenBuffers(1, &preview.dislocationLineVBO);
        }

        glBindVertexArray(preview.dislocationLineVAO);
        glBindBuffer(GL_ARRAY_BUFFER, preview.dislocationLineVBO);
        glBufferData(GL_ARRAY_BUFFER,
                     dislocationOverlay->loopPoints.size() * sizeof(glm::vec3),
                     dislocationOverlay->loopPoints.data(),
                     GL_DYNAMIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

        preview.dislocationLoopPointCount = dislocationOverlay->loopPoints.size();
    }

    std::array<bool, 119> noFilter = {};
    preview.buffers.upload(data, false, noFilter);
    preview.dirty = false;
}

void DislocationBuilderDialog::autoFitPreviewCamera(PreviewState& preview)
{
    preview.yaw = 45.0f;
    preview.pitch = 35.0f;
    if (preview.buffers.atomCount == 0)
    {
        preview.distance = 10.0f;
        return;
    }

    float maxRadius = 0.0f;
    for (size_t i = 0; i < preview.buffers.atomPositions.size(); ++i)
    {
        const float radius = (i < preview.buffers.atomRadii.size())
            ? preview.buffers.atomRadii[i] : 0.0f;
        const float distance = glm::length(preview.buffers.atomPositions[i]
                                         - preview.buffers.orbitCenter) + radius;
        maxRadius = std::max(maxRadius, distance);
    }

    const float halfFov = glm::radians(22.5f);
    const float distance = maxRadius / std::sin(halfFov) * 1.15f;
    preview.distance = std::max(Camera::kMinDistance,
                                std::min(Camera::kMaxDistance, distance));
}

void DislocationBuilderDialog::renderPreviewToFBO(PreviewState& preview, int width, int height)
{
    if (!m_glReady || !m_renderer || !preview.sphere || !preview.cylinder || preview.buffers.atomCount == 0)
        return;

    ensurePreviewFBO(preview, width, height);

    GLint previousFbo = 0;
    GLint previousViewport[4];
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &previousFbo);
    glGetIntegerv(GL_VIEWPORT, previousViewport);

    Camera camera;
    camera.yaw = preview.yaw;
    camera.pitch = preview.pitch;
    camera.distance = preview.distance;

    FrameView frame;
    frame.framebufferWidth = width;
    frame.framebufferHeight = height;
    buildFrameView(camera, preview.buffers, true, frame);

    glBindFramebuffer(GL_FRAMEBUFFER, preview.fbo);
    glViewport(0, 0, width, height);
    glEnable(GL_DEPTH_TEST);
    {
        const ImVec4& bg = ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
        glClearColor(bg.x, bg.y, bg.z, bg.w);
    }
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    m_renderer->drawBonds(frame.projection, frame.view,
                          frame.lightPosition, frame.cameraPosition,
                          preview.buffers.tabCylinderVAO, preview.buffers.tabCylinderVertexCount,
                          preview.buffers.bondCount);
    m_renderer->drawAtoms(frame.projection, frame.view,
                          frame.lightMVP, frame.lightPosition, frame.cameraPosition,
                          preview.shadow,
                          preview.buffers.tabSphereVAO, preview.buffers.tabSphereIndexCount,
                          preview.buffers.atomCount);
    m_renderer->drawBoxLines(frame.projection, frame.view,
                             preview.buffers.lineVAO,
                             preview.buffers.boxLines.size());

    // Draw dislocation loop in bright red
    if (preview.hasDislocationLine && preview.dislocationLineVAO != 0 && preview.dislocationLoopPointCount > 0)
    {
        // Use the line program to render the dislocation loop
        glUseProgram(m_renderer->lineProgram);
        glUniformMatrix4fv(glGetUniformLocation(m_renderer->lineProgram, "projection"),
                           1, GL_FALSE, glm::value_ptr(frame.projection));
        glUniformMatrix4fv(glGetUniformLocation(m_renderer->lineProgram, "view"),
                           1, GL_FALSE, glm::value_ptr(frame.view));
        glUniform3f(glGetUniformLocation(m_renderer->lineProgram, "uColor"), 1.0f, 0.0f, 0.0f);  // Bright red

        glLineWidth(3.5f);
        glBindVertexArray(preview.dislocationLineVAO);
        glDrawArrays(GL_LINE_LOOP, 0, (GLsizei)preview.dislocationLoopPointCount);
        glBindVertexArray(0);
    }

    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)previousFbo);
    glViewport(previousViewport[0], previousViewport[1], previousViewport[2], previousViewport[3]);
}

void DislocationBuilderDialog::drawDialog(
    Structure& structure,
    const std::vector<glm::vec3>& /*elementColors*/,
    const std::vector<float>& elementRadii,
    const std::vector<float>& elementShininess,
    const std::function<void(Structure&)>& updateBuffers)
{
    pollGenerationResult(elementRadii, elementShininess);

    if (!m_pendingDropPath.empty())
    {
        tryLoadFile(m_pendingDropPath, elementRadii, elementShininess);
        m_pendingDropPath.clear();
    }

    if (m_openRequested)
    {
        ImGui::OpenPopup("Insert Dislocation");
        m_openRequested = false;
    }

    m_isOpen = ImGui::IsPopupOpen("Insert Dislocation");

    ImGui::SetNextWindowSize(ImVec2(1240.0f, 780.0f), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(960.0f, 600.0f), ImVec2(3200.0f, 3200.0f));

    bool keepOpen = true;
    if (!ImGui::BeginPopupModal("Insert Dislocation", &keepOpen, ImGuiWindowFlags_NoCollapse))
    {
        m_isOpen = false;
        return;
    }
    m_isOpen = true;

    ImGui::TextDisabled("Isotropic-elastic dislocation insertion with automatic FCC/HCP/BCC detection and shape-controlled application region.");
    ImGui::SameLine();
    ImGui::TextColored(validationStateColor(m_result), "%s", validationStateLabel(m_result));
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    const Structure* activeSource = nullptr;
    if (m_sourceLoaded)
        activeSource = m_useCurrentSceneSource ? &structure : &m_source;

    const float contentHeight = ImGui::GetContentRegionAvail().y - 58.0f;
    const float previewPanelWidth = ImGui::GetContentRegionAvail().x * 0.62f;

    ImGui::BeginChild("##disloc_previews", ImVec2(previewPanelWidth, contentHeight), false);
    {
        const float splitGap = ImGui::GetStyle().ItemSpacing.y;
        const float upperHeight = (ImGui::GetContentRegionAvail().y - splitGap) * 0.5f;
        const float lowerHeight = ImGui::GetContentRegionAvail().y - upperHeight - splitGap;

        ImGui::BeginChild("##disloc_source", ImVec2(0.0f, upperHeight), false);
        {
            ImGui::Text("Input Preview");
            ImGui::SameLine();
            ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(10.0f, 6.0f));
            if (ImGui::Button("Use Current Scene"))
            {
                loadFromScene(structure, elementRadii, elementShininess);
                activeSource = &structure;
            }
            ImGui::PopStyleVar();

            if (activeSource != nullptr)
            {
                ImGui::TextColored(ImVec4(0.40f, 0.90f, 0.40f, 1.0f), "%s", m_sourceLabel.c_str());
                ImGui::SameLine();
                ImGui::TextDisabled("(%d atoms)", (int)activeSource->atoms.size());
            }
            else
            {
                ImGui::TextDisabled("Drop a structure here or use the current scene.");
            }

            const float previewHeight = ImGui::GetContentRegionAvail().y - 2.0f;
            const float previewWidth = ImGui::GetContentRegionAvail().x;
            const ImVec2 minCorner = ImGui::GetCursorScreenPos();
            const ImVec2 maxCorner(minCorner.x + previewWidth, minCorner.y + previewHeight);

            ImGui::InvisibleButton("##disloc_source_preview", ImVec2(previewWidth, previewHeight));
            const bool hovered = ImGui::IsItemHovered();
            const bool active = ImGui::IsItemActive();
            ImDrawList* drawList = ImGui::GetWindowDrawList();
            const ImVec4& bg = ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
            drawList->AddRectFilled(minCorner, maxCorner,
                                    IM_COL32((int)(bg.x * 255), (int)(bg.y * 255), (int)(bg.z * 255), 255));
            drawList->AddRect(minCorner, maxCorner, ImGui::GetColorU32(ImGuiCol_Separator), 3.0f);

            if (activeSource == nullptr)
            {
                drawPreviewHint(drawList, minCorner, maxCorner,
                                "Drop a structure file here",
                                "or use the current scene");
            }
            else
            {
                if (m_sourcePreview.dirty)
                    rebuildPreviewBuffers(m_sourcePreview, *activeSource, elementRadii, elementShininess);

                renderPreviewToFBO(m_sourcePreview, std::max(1, (int)previewWidth), std::max(1, (int)previewHeight));
                drawList->AddImage((ImTextureID)(intptr_t)m_sourcePreview.colorTex,
                                   minCorner, maxCorner, ImVec2(0, 1), ImVec2(1, 0));

                if (active && ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f))
                {
                    const ImVec2 delta = ImGui::GetIO().MouseDelta;
                    m_sourcePreview.yaw -= delta.x * 0.5f;
                    m_sourcePreview.pitch += delta.y * 0.5f;
                }
                if (hovered)
                {
                    const float wheel = ImGui::GetIO().MouseWheel;
                    if (wheel != 0.0f)
                    {
                        m_sourcePreview.distance -= wheel * m_sourcePreview.distance * 0.1f;
                        m_sourcePreview.distance = std::max(Camera::kMinDistance,
                                                            std::min(Camera::kMaxDistance, m_sourcePreview.distance));
                    }
                }
            }
        }
        ImGui::EndChild();

        ImGui::BeginChild("##disloc_output", ImVec2(0.0f, lowerHeight), false);
        {
            ImGui::Text("Output Preview");
            if (m_result.success)
            {
                ImGui::SameLine();
                ImGui::TextColored(validationStateColor(m_result), "%s", validationStateLabel(m_result));
            }
            else if (m_generationRunning)
            {
                ImGui::SameLine();
                ImGui::TextDisabled("(building...)");
            }

            const float previewHeight = ImGui::GetContentRegionAvail().y - 2.0f;
            const float previewWidth = ImGui::GetContentRegionAvail().x;
            const ImVec2 minCorner = ImGui::GetCursorScreenPos();
            const ImVec2 maxCorner(minCorner.x + previewWidth, minCorner.y + previewHeight);

            ImGui::InvisibleButton("##disloc_output_preview", ImVec2(previewWidth, previewHeight));
            const bool hovered = ImGui::IsItemHovered();
            const bool active = ImGui::IsItemActive();
            ImDrawList* drawList = ImGui::GetWindowDrawList();
            const ImVec4& bg = ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
            drawList->AddRectFilled(minCorner, maxCorner,
                                    IM_COL32((int)(bg.x * 255), (int)(bg.y * 255), (int)(bg.z * 255), 255));
            drawList->AddRect(minCorner, maxCorner, ImGui::GetColorU32(ImGuiCol_Separator), 3.0f);

            if (!m_result.success)
            {
                drawPreviewHint(drawList, minCorner, maxCorner,
                                "Generate dislocation structure",
                                "to preview and validate");
            }
            else
            {
                if (m_outputPreview.dirty)
                    rebuildPreviewBuffers(m_outputPreview,
                                          m_result.output,
                                          elementRadii,
                                          elementShininess,
                                          &m_result);

                renderPreviewToFBO(m_outputPreview, std::max(1, (int)previewWidth), std::max(1, (int)previewHeight));
                drawList->AddImage((ImTextureID)(intptr_t)m_outputPreview.colorTex,
                                   minCorner, maxCorner, ImVec2(0, 1), ImVec2(1, 0));

                if (active && ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f))
                {
                    const ImVec2 delta = ImGui::GetIO().MouseDelta;
                    m_outputPreview.yaw -= delta.x * 0.5f;
                    m_outputPreview.pitch += delta.y * 0.5f;
                }
                if (hovered)
                {
                    const float wheel = ImGui::GetIO().MouseWheel;
                    if (wheel != 0.0f)
                    {
                        m_outputPreview.distance -= wheel * m_outputPreview.distance * 0.1f;
                        m_outputPreview.distance = std::max(Camera::kMinDistance,
                                                            std::min(Camera::kMaxDistance, m_outputPreview.distance));
                    }
                }
            }
        }
        ImGui::EndChild();
    }
    ImGui::EndChild();

    ImGui::SameLine();

    ImGui::BeginChild("##disloc_controls", ImVec2(0, contentHeight), false);
    {
        ImGui::Text("Dislocation Options");
        ImGui::Separator();

        const bool canGenerate = (activeSource != nullptr) && !m_generationRunning;
        if (m_generationRunning)
            ImGui::TextDisabled("Generation running in background...");

        if (ImGui::Button("Detect Lattice", ImVec2(-1.0f, 0.0f)))
        {
            if (!m_sourceLoaded)
            {
                m_statusMsg = "Load a structure first.";
                m_statusIsError = true;
            }
            else
            {
                analyzeSource(*activeSource);
            }
        }

        if (m_detection.success)
        {
            ImGui::Text("Detected lattice: %s", dislocationLatticeFamilyName(m_detection.family));
            ImGui::TextDisabled("Recognized atoms: %d", m_detection.recognizedCount);
        }
        else
        {
            ImGui::TextDisabled("Detected lattice: not identified yet");
        }

        ImGui::Separator();
        ImGui::TextDisabled("Character / Region");

        if (ImGui::BeginCombo("Character", characterLabel(m_params.character)))
        {
            for (const CharacterOption& option : kCharacterOptions)
            {
                const bool selected = option.value == m_params.character;
                if (ImGui::Selectable(option.label, selected))
                    m_params.character = option.value;
                if (selected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        if (ImGui::BeginCombo("Shape", shapeLabel(m_params.shape)))
        {
            for (const ShapeOption& option : kShapeOptions)
            {
                const bool selected = option.value == m_params.shape;
                if (ImGui::Selectable(option.label, selected))
                    m_params.shape = option.value;
                if (selected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        ImGui::Separator();
        ImGui::TextDisabled("Crystallographic Directions");
        ImGui::Checkbox("Auto lattice vectors", &m_params.autoDirections);
        if (!m_params.autoDirections)
        {
            int plane[3] = {m_params.planeHkl.x, m_params.planeHkl.y, m_params.planeHkl.z};
            int burgers[3] = {m_params.burgersUvw.x, m_params.burgersUvw.y, m_params.burgersUvw.z};
            int line[3] = {m_params.lineUvw.x, m_params.lineUvw.y, m_params.lineUvw.z};

            ImGui::DragInt3("Plane hkl", plane, 0.2f, -6, 6);
            ImGui::DragInt3("Burgers uvw", burgers, 0.2f, -6, 6);
            ImGui::DragInt3("Line uvw", line, 0.2f, -6, 6);

            m_params.planeHkl = glm::ivec3(plane[0], plane[1], plane[2]);
            m_params.burgersUvw = glm::ivec3(burgers[0], burgers[1], burgers[2]);
            m_params.lineUvw = glm::ivec3(line[0], line[1], line[2]);
        }

        ImGui::Checkbox("Line point in fractional coordinates", &m_params.useFractionalLinePoint);
        if (m_params.useFractionalLinePoint)
            ImGui::DragFloat3("Line point frac", &m_params.linePointFractional.x, 0.01f, -1.0f, 2.0f, "%.3f");

        ImGui::DragFloat3("Line point offset (A)", &m_params.linePointCartesianOffset.x, 0.05f, -100.0f, 100.0f, "%.3f");

        ImGui::Separator();
        ImGui::TextDisabled("Elastic Parameters");
        ImGui::DragFloat("Burgers scale", &m_params.burgersScale, 0.02f, 0.05f, 3.0f, "%.3f");
        ImGui::DragFloat("Poisson ratio", &m_params.poissonRatio, 0.005f, 0.05f, 0.49f, "%.3f");
        ImGui::DragFloat("Core radius (A)", &m_params.coreRadius, 0.05f, 0.05f, 20.0f, "%.3f");
        ImGui::DragFloat("Cutoff radius (A)", &m_params.cutoffRadius, 0.2f, 0.0f, 100.0f, "%.2f");
        ImGui::DragFloat("Line half length (A)", &m_params.lineHalfLength, 0.5f, 0.1f, 1000000.0f, "%.1f");

        if (m_params.character == DislocationCharacter::Mixed)
            ImGui::SliderFloat("Mixed angle (deg)", &m_params.mixedCharacterAngleDeg, 0.0f, 90.0f, "%.1f");

        ImGui::Separator();
        ImGui::TextDisabled("Shape Parameters");
        if (m_params.shape == DislocationShape::Cylinder)
        {
            ImGui::DragFloat("Cylinder radius (A)", &m_params.cylinderRadius, 0.2f, 0.1f, 200.0f, "%.2f");
        }
        else if (m_params.shape == DislocationShape::Sphere)
        {
            ImGui::DragFloat("Sphere radius (A)", &m_params.sphereRadius, 0.2f, 0.1f, 200.0f, "%.2f");
        }
        else if (m_params.shape == DislocationShape::Ellipsoid)
        {
            ImGui::DragFloat3("Ellipsoid radii (A)", &m_params.ellipsoidRadii.x, 0.2f, 0.1f, 200.0f, "%.2f");
        }
        else if (m_params.shape == DislocationShape::Freeform2D)
        {
            ImGui::TextDisabled("Freeform polygon in local (x, y) around the dislocation line.");
            if (ImGui::Button("Add Point"))
            {
                if (m_params.freeformPoints.empty())
                    m_params.freeformPoints.push_back(glm::vec2(0.0f, 0.0f));
                else
                    m_params.freeformPoints.push_back(m_params.freeformPoints.back() + glm::vec2(1.0f, 0.0f));
            }
            ImGui::SameLine();
            if (ImGui::Button("Remove Last") && !m_params.freeformPoints.empty())
                m_params.freeformPoints.pop_back();

            const float pointsHeight = 112.0f;
            if (ImGui::BeginChild("##freeform_points", ImVec2(0.0f, pointsHeight), true))
            {
                for (int i = 0; i < (int)m_params.freeformPoints.size(); ++i)
                {
                    ImGui::PushID(i + 24000);
                    ImGui::SetNextItemWidth(-1.0f);
                    ImGui::DragFloat2("##pt", &m_params.freeformPoints[i].x, 0.05f, -500.0f, 500.0f, "%.2f");
                    ImGui::PopID();
                }
                ImGui::EndChild();
            }
        }

        ImGui::Spacing();
        if (!canGenerate)
            ImGui::BeginDisabled();
        if (ImGui::Button("Generate Dislocation", ImVec2(-1.0f, 0.0f)))
            generateDislocation(*activeSource);
        if (!canGenerate)
            ImGui::EndDisabled();

        if (m_result.success)
        {
            ImGui::Separator();
            ImGui::Text("Result summary");
            ImGui::TextDisabled("Shifted atoms: %d", m_result.shiftedAtomCount);
            ImGui::TextDisabled("|b|: %.4f A", m_result.burgersMagnitude);
            ImGui::TextDisabled("RMS displacement: %.4f A", m_result.validation.rmsDisplacement);
            ImGui::TextDisabled("Max displacement: %.4f A", m_result.validation.maxDisplacement);
            ImGui::TextDisabled("Min distance: %.4f A", m_result.validation.minInteratomicDistance);
            ImGui::TextDisabled("Family before/after: %s -> %s",
                                dislocationLatticeFamilyName(m_result.validation.familyBefore),
                                dislocationLatticeFamilyName(m_result.validation.familyAfter));

            if (ImGui::Button("Replace Main Scene", ImVec2(-1.0f, 0.0f)))
            {
                structure = m_result.output;
                structure.dislocationLoopPoints = m_result.loopPoints;
                updateBuffers(structure);
                m_statusMsg = "Main scene replaced with generated dislocation structure.";
                m_statusIsError = false;
            }
        }
    }
    ImGui::EndChild();

    ImGui::Separator();
    if (!m_statusMsg.empty())
    {
        const ImVec4 color = m_statusIsError
            ? ImVec4(1.0f, 0.4f, 0.4f, 1.0f)
            : ImVec4(0.5f, 0.9f, 0.5f, 1.0f);
        ImGui::TextColored(color, "%s", m_statusMsg.c_str());
        ImGui::SameLine();
    }

    const float buttonWidth = 100.0f;
    ImGui::SetCursorPosX(ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x - buttonWidth);
    if (ImGui::Button("Close", ImVec2(buttonWidth, 0.0f)))
        ImGui::CloseCurrentPopup();

    ImGui::EndPopup();
}
