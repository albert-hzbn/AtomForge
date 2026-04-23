#include "ui/StackingFaultBuilderDialog.h"

#include "app/SceneView.h"
#include "camera/Camera.h"
#include "graphics/CylinderMesh.h"
#include "graphics/Renderer.h"
#include "graphics/SphereMesh.h"
#include "graphics/StructureInstanceBuilder.h"
#include "io/StructureLoader.h"
#include "util/PathUtils.h"

#include "imgui.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <sstream>

namespace
{
struct ExportFormatOption
{
    const char* label;
    const char* extension;
    const char* format;
};

constexpr ExportFormatOption kExportFormats[] = {
    {"XYZ (.xyz)", ".xyz", "xyz"},
    {"CIF (.cif)", ".cif", "cif"},
    {"VASP (.vasp)", ".vasp", "vasp"},
    {"Quantum ESPRESSO (.pwi)", ".pwi", "pwscf"},
};

struct PlaneOption
{
    const char* label;
    StackingFaultPlane value;
};

constexpr PlaneOption kAutoPlaneOptions[] = {
    {"Auto", StackingFaultPlane::Auto},
};

constexpr PlaneOption kFccPlaneOptions[] = {
    {"(112) <111>", StackingFaultPlane::Fcc111},
};

constexpr PlaneOption kHcpPlaneOptions[] = {
    {"Basal", StackingFaultPlane::HcpBasal},
    {"Prismatic", StackingFaultPlane::HcpPrismatic},
    {"Pyramidal", StackingFaultPlane::HcpPyramidal},
};

constexpr PlaneOption kBccPlaneOptions[] = {
    {"(110) <111>", StackingFaultPlane::Bcc110},
    {"(112) <111>", StackingFaultPlane::Bcc112},
};

const PlaneOption* getPlaneOptions(StackingFaultFamily family, int& count)
{
    switch (family)
    {
    case StackingFaultFamily::Fcc:
        count = (int)(sizeof(kFccPlaneOptions) / sizeof(kFccPlaneOptions[0]));
        return kFccPlaneOptions;
    case StackingFaultFamily::Hcp:
        count = (int)(sizeof(kHcpPlaneOptions) / sizeof(kHcpPlaneOptions[0]));
        return kHcpPlaneOptions;
    case StackingFaultFamily::Bcc:
        count = (int)(sizeof(kBccPlaneOptions) / sizeof(kBccPlaneOptions[0]));
        return kBccPlaneOptions;
    default:
        count = 1;
        return kAutoPlaneOptions;
    }
}

const char* cellModeLabel(StackingFaultCellMode mode)
{
    switch (mode)
    {
    case StackingFaultCellMode::OrthogonalCell:
        return "Orthogonal Cell";
    case StackingFaultCellMode::SmallestUnitCell:
    default:
        return "Smallest Unit Cell";
    }
}

std::string basenameFromPath(const std::string& path)
{
    const std::string::size_type slash = path.find_last_of("\\/");
    return (slash == std::string::npos) ? path : path.substr(slash + 1);
}

std::string stemFromPath(const std::string& path)
{
    const std::string base = basenameFromPath(path);
    const std::string::size_type dot = base.find_last_of('.');
    return (dot == std::string::npos) ? base : base.substr(0, dot);
}

std::string sanitizeName(const std::string& value)
{
    std::string out;
    out.reserve(value.size());
    for (char ch : value)
    {
        if ((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z')
            || (ch >= '0' && ch <= '9'))
            out.push_back(ch);
        else if (ch == '-' || ch == '_')
            out.push_back(ch);
        else if (ch == ' ')
            out.push_back('_');
    }
    if (out.empty())
        out = "stacking_fault";
    return out;
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

StackingFaultBuilderDialog::StackingFaultBuilderDialog()
{
    std::snprintf(m_exportDir, sizeof(m_exportDir), "%s", ".");
    std::snprintf(m_exportFolderName, sizeof(m_exportFolderName), "%s", "stacking_fault_sequence");
}

StackingFaultBuilderDialog::~StackingFaultBuilderDialog()
{
    releasePreview(m_sourcePreview);
    releasePreview(m_outputPreview);
}

void StackingFaultBuilderDialog::initRenderResources(Renderer& renderer)
{
    m_renderer = &renderer;
    m_sourcePreview.sphere   = new SphereMesh(24, 24);
    m_sourcePreview.cylinder = new CylinderMesh(16);
    m_outputPreview.sphere   = new SphereMesh(24, 24);
    m_outputPreview.cylinder = new CylinderMesh(16);
    m_sourcePreview.buffers.init(m_sourcePreview.sphere->vbo, m_sourcePreview.sphere->ebo, m_sourcePreview.sphere->indexCount,
                                  m_sourcePreview.cylinder->vbo, m_sourcePreview.cylinder->vertexCount);
    m_outputPreview.buffers.init(m_outputPreview.sphere->vbo, m_outputPreview.sphere->ebo, m_outputPreview.sphere->indexCount,
                                  m_outputPreview.cylinder->vbo, m_outputPreview.cylinder->vertexCount);
    m_sourcePreview.shadow = createShadowMap(1, 1);
    m_outputPreview.shadow = createShadowMap(1, 1);
    m_glReady = true;
}

void StackingFaultBuilderDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Stacking Faults", nullptr, false, enabled))
        m_openRequested = true;
}

void StackingFaultBuilderDialog::feedDroppedFile(const std::string& path)
{
    m_pendingDropPath = path;
}

bool StackingFaultBuilderDialog::tryLoadFile(const std::string& path,
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
        m_statusMsg = "Stacking faults builder requires a structure with a unit cell.";
        m_statusIsError = true;
        return false;
    }

    m_source = std::move(loaded);
    m_sourceLoaded = true;
    m_sourcePath = path;
    m_sourceLabel = basenameFromPath(path);
    std::snprintf(m_exportDir, sizeof(m_exportDir), "%s", normalizePathSeparators(parentPath(path)).c_str());
    std::snprintf(m_exportFolderName, sizeof(m_exportFolderName), "%s",
                  (sanitizeName(stemFromPath(path)) + "_stacking_fault").c_str());

    analyzeSource();
    m_sourcePreview.dirty = true;
    rebuildPreviewBuffers(m_sourcePreview, m_source, radii, shininess);
    autoFitPreviewCamera(m_sourcePreview);

    m_result = {};
    m_selectedStructureIndex = 0;
    m_outputPreview.dirty = true;
    return true;
}

void StackingFaultBuilderDialog::loadFromScene(const Structure& scene,
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

    m_source = scene;
    m_sourceLoaded = true;
    m_sourcePath.clear();
    m_sourceLabel = "scene";
    std::snprintf(m_exportDir, sizeof(m_exportDir), "%s", ".");
    std::snprintf(m_exportFolderName, sizeof(m_exportFolderName), "%s", "scene_stacking_fault");

    analyzeSource();
    rebuildPreviewBuffers(m_sourcePreview, m_source, radii, shininess);
    autoFitPreviewCamera(m_sourcePreview);

    m_result = {};
    m_selectedStructureIndex = 0;
    m_outputPreview.dirty = true;
}

void StackingFaultBuilderDialog::analyzeSource()
{
    m_detection = detectStackingFaultFamily(m_source, true);
    m_statusMsg = m_detection.message;
    m_statusIsError = !m_detection.success;

    if (!m_detection.success)
        return;

    int planeCount = 0;
    const PlaneOption* options = getPlaneOptions(m_detection.family, planeCount);
    if (planeCount > 0)
        m_params.plane = options[0].value;
}

void StackingFaultBuilderDialog::regenerateSequence(const std::vector<float>& radii,
                                                    const std::vector<float>& shininess)
{
    if (!m_sourceLoaded)
    {
        m_statusMsg = "Load a structure first.";
        m_statusIsError = true;
        return;
    }

    m_result = buildStackingFaultSequence(m_source, m_params);
    m_statusMsg = m_result.message;
    m_statusIsError = !m_result.success;
    m_selectedStructureIndex = 0;

    if (m_result.success && !m_result.sequence.empty())
    {
        rebuildPreviewBuffers(m_outputPreview,
                              m_result.sequence[m_selectedStructureIndex].structure,
                              radii,
                              shininess);
        autoFitPreviewCamera(m_outputPreview);
    }
    else
    {
        m_outputPreview.dirty = true;
    }
}

bool StackingFaultBuilderDialog::exportSequence()
{
    if (!m_result.success || m_result.sequence.empty())
    {
        m_statusMsg = "Generate a stacking-fault sequence before exporting.";
        m_statusIsError = true;
        return false;
    }

    const ExportFormatOption& format = kExportFormats[m_exportFormatIndex];
    std::filesystem::path baseDir(normalizePathSeparators(std::string(m_exportDir)));
    std::filesystem::path targetDir = baseDir / sanitizeName(m_exportFolderName);

    std::error_code ec;
    std::filesystem::create_directories(targetDir, ec);
    if (ec)
    {
        m_statusMsg = std::string("Could not create export folder: ") + ec.message();
        m_statusIsError = true;
        return false;
    }

    for (size_t i = 0; i < m_result.sequence.size(); ++i)
    {
        char indexBuf[16];
        std::snprintf(indexBuf, sizeof(indexBuf), "%03d", (int)i);
        const std::string filename = std::string("sf_") + indexBuf + "_"
                                   + sanitizeName(m_result.sequence[i].label) + format.extension;
        const std::filesystem::path filePath = targetDir / filename;
        if (!saveStructure(m_result.sequence[i].structure,
                           normalizePathSeparators(filePath.string()),
                           format.format))
        {
            m_statusMsg = std::string("Failed to export: ") + filename;
            m_statusIsError = true;
            return false;
        }
    }

    m_statusMsg = std::string("Exported ") + std::to_string(m_result.sequence.size())
                + " structures to " + normalizePathSeparators(targetDir.string());
    m_statusIsError = false;
    return true;
}

void StackingFaultBuilderDialog::releasePreview(PreviewState& preview)
{
    if (preview.fbo) glDeleteFramebuffers(1, &preview.fbo);
    if (preview.colorTex) glDeleteTextures(1, &preview.colorTex);
    if (preview.depthRbo) glDeleteRenderbuffers(1, &preview.depthRbo);
    if (preview.shadow.depthFBO) glDeleteFramebuffers(1, &preview.shadow.depthFBO);
    if (preview.shadow.depthTexture) glDeleteTextures(1, &preview.shadow.depthTexture);
    delete preview.sphere;
    delete preview.cylinder;
    preview = {};
}

void StackingFaultBuilderDialog::ensurePreviewFBO(PreviewState& preview, int width, int height)
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

void StackingFaultBuilderDialog::rebuildPreviewBuffers(PreviewState& preview,
                                                       const Structure& structure,
                                                       const std::vector<float>& radii,
                                                       const std::vector<float>& shininess)
{
    if (!m_glReady || structure.atoms.empty())
        return;

    static const int kIdentity[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    const StructureInstanceData data = buildStructureInstanceData(structure,
                                                                  false,
                                                                  kIdentity,
                                                                  radii,
                                                                  shininess);
    std::array<bool, 119> noFilter = {};
    preview.buffers.upload(data, false, noFilter);
    preview.dirty = false;
}

void StackingFaultBuilderDialog::autoFitPreviewCamera(PreviewState& preview)
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
    float distance = maxRadius / std::sin(halfFov) * 1.15f;
    preview.distance = std::max(Camera::kMinDistance,
                                std::min(Camera::kMaxDistance, distance));
}

void StackingFaultBuilderDialog::renderPreviewToFBO(PreviewState& preview, int width, int height)
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

    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)previousFbo);
    glViewport(previousViewport[0], previousViewport[1], previousViewport[2], previousViewport[3]);
}

void StackingFaultBuilderDialog::drawDialog(
    Structure& structure,
    const std::vector<glm::vec3>& /*elementColors*/,
    const std::vector<float>& elementRadii,
    const std::vector<float>& elementShininess,
    const std::function<void(Structure&)>& updateBuffers)
{
    if (!m_pendingDropPath.empty())
    {
        tryLoadFile(m_pendingDropPath, elementRadii, elementShininess);
        m_pendingDropPath.clear();
    }

    if (m_openRequested)
    {
        ImGui::OpenPopup("Stacking Faults");
        m_openRequested = false;
    }

    m_isOpen = ImGui::IsPopupOpen("Stacking Faults");

    ImGui::SetNextWindowSize(ImVec2(1220.0f, 760.0f), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(900.0f, 560.0f), ImVec2(3000.0f, 3000.0f));

    bool keepOpen = true;
    if (!ImGui::BeginPopupModal("Stacking Faults", &keepOpen, ImGuiWindowFlags_NoCollapse))
    {
        m_isOpen = false;
        return;
    }
    m_isOpen = true;

    ImGui::TextDisabled("Load a crystal, identify whether it is FCC-like, HCP-like, or BCC-like, and build a sliding sequence across the detected stacking-fault family.");
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    const float contentHeight = ImGui::GetContentRegionAvail().y - 58.0f;
    const float leftWidth = ImGui::GetContentRegionAvail().x * 0.34f;
    const float centerWidth = ImGui::GetContentRegionAvail().x * 0.31f;

    ImGui::BeginChild("##sf_source", ImVec2(leftWidth, contentHeight), false);
    {
        ImGui::Text("Input Structure");
        ImGui::SameLine();
        if (ImGui::SmallButton("Use Current Scene"))
            loadFromScene(structure, elementRadii, elementShininess);

        if (m_sourceLoaded)
        {
            ImGui::TextColored(ImVec4(0.40f, 0.90f, 0.40f, 1.0f), "%s", m_sourceLabel.c_str());
            ImGui::SameLine();
            ImGui::TextDisabled("(%d atoms)", (int)m_source.atoms.size());
        }
        else
        {
            ImGui::TextDisabled("Drop a structure here or use the current scene.");
        }

        const float previewHeight = ImGui::GetContentRegionAvail().y - 2.0f;
        const float previewWidth = ImGui::GetContentRegionAvail().x;
        const ImVec2 minCorner = ImGui::GetCursorScreenPos();
        const ImVec2 maxCorner(minCorner.x + previewWidth, minCorner.y + previewHeight);

        ImGui::InvisibleButton("##sf_source_preview", ImVec2(previewWidth, previewHeight));
        const bool hovered = ImGui::IsItemHovered();
        const bool active = ImGui::IsItemActive();
        ImDrawList* drawList = ImGui::GetWindowDrawList();
        const ImVec4& bg = ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
        drawList->AddRectFilled(minCorner, maxCorner,
                                IM_COL32((int)(bg.x * 255), (int)(bg.y * 255), (int)(bg.z * 255), 255));
        drawList->AddRect(minCorner, maxCorner, ImGui::GetColorU32(ImGuiCol_Separator), 3.0f);

        if (!m_sourceLoaded)
        {
            drawPreviewHint(drawList, minCorner, maxCorner,
                            "Drop a structure file here",
                            "or use the current scene");
        }
        else
        {
            if (m_sourcePreview.dirty)
                rebuildPreviewBuffers(m_sourcePreview, m_source, elementRadii, elementShininess);

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

    ImGui::SameLine();

    ImGui::BeginChild("##sf_controls", ImVec2(centerWidth, contentHeight), false);
    {
        ImGui::Text("Builder");
        ImGui::Separator();
        if (m_detection.success)
        {
            ImGui::Text("Detected phase: %s", stackingFaultFamilyName(m_detection.family));
            ImGui::TextDisabled("Recognized atoms: %d", m_detection.recognizedCount);
        }
        else
        {
            ImGui::TextDisabled("Detected phase: not identified yet");
        }

        if (ImGui::Button("Detect Structure", ImVec2(-1.0f, 0.0f)))
        {
            if (!m_sourceLoaded)
            {
                m_statusMsg = "Load a structure first.";
                m_statusIsError = true;
            }
            else
            {
                analyzeSource();
            }
        }

        ImGui::Spacing();
        int planeCount = 0;
        const PlaneOption* planeOptions = getPlaneOptions(m_detection.family, planeCount);
        const char* selectedPlane = planeOptions[0].label;
        for (int i = 0; i < planeCount; ++i)
        {
            if (planeOptions[i].value == m_params.plane)
            {
                selectedPlane = planeOptions[i].label;
                break;
            }
        }
        if (ImGui::BeginCombo("Plane", selectedPlane))
        {
            for (int i = 0; i < planeCount; ++i)
            {
                const bool selected = planeOptions[i].value == m_params.plane;
                if (ImGui::Selectable(planeOptions[i].label, selected))
                    m_params.plane = planeOptions[i].value;
                if (selected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        const char* currentMode = cellModeLabel(m_params.cellMode);
        if (ImGui::BeginCombo("Cell Mode", currentMode))
        {
            const bool smallest = m_params.cellMode == StackingFaultCellMode::SmallestUnitCell;
            if (ImGui::Selectable("Smallest Unit Cell", smallest))
                m_params.cellMode = StackingFaultCellMode::SmallestUnitCell;
            if (smallest)
                ImGui::SetItemDefaultFocus();

            const bool orthogonal = m_params.cellMode == StackingFaultCellMode::OrthogonalCell;
            if (ImGui::Selectable("Orthogonal Cell", orthogonal))
                m_params.cellMode = StackingFaultCellMode::OrthogonalCell;
            ImGui::EndCombo();
        }

        ImGui::SliderInt("Layers", &m_params.layerCount, 2, 96);
        ImGui::DragFloat("Interval", &m_params.interval, 0.01f, 0.01f, 1.0f, "%.2f");
        ImGui::DragFloat("Max Displacement", &m_params.maxDisplacementFactor, 0.05f, 0.10f, 4.0f, "%.2f");

        ImGui::Spacing();
        if (ImGui::Button("Generate Sequence", ImVec2(-1.0f, 0.0f)))
            regenerateSequence(elementRadii, elementShininess);

        if (m_result.success && !m_result.sequence.empty())
        {
            ImGui::Separator();
            ImGui::Text("Sequence");
            int previewIndex = m_selectedStructureIndex;
            if (ImGui::SliderInt("Structure", &previewIndex, 0, (int)m_result.sequence.size() - 1))
            {
                m_selectedStructureIndex = previewIndex;
                rebuildPreviewBuffers(m_outputPreview,
                                      m_result.sequence[m_selectedStructureIndex].structure,
                                      elementRadii,
                                      elementShininess);
            }
            ImGui::Text("Label: %s", m_result.sequence[m_selectedStructureIndex].label.c_str());
            ImGui::TextDisabled("Partial displacement: %.4f A", m_result.partialDisplacement);
            ImGui::TextDisabled("Shifted atoms: %d", m_result.shiftedAtomCount);

            if (ImGui::Button("Use Selected Structure", ImVec2(-1.0f, 0.0f)))
            {
                structure = m_result.sequence[m_selectedStructureIndex].structure;
                updateBuffers(structure);
                m_statusMsg = std::string("Applied ") + m_result.sequence[m_selectedStructureIndex].label + " to the main scene.";
                m_statusIsError = false;
                ImGui::CloseCurrentPopup();
            }

            ImGui::Separator();
            ImGui::Text("Export");
            ImGui::InputText("Folder", m_exportDir, sizeof(m_exportDir));
            ImGui::InputText("Subfolder", m_exportFolderName, sizeof(m_exportFolderName));
            if (m_sourceLoaded && !m_sourcePath.empty())
            {
                if (ImGui::SmallButton("Use Source Folder"))
                    std::snprintf(m_exportDir, sizeof(m_exportDir), "%s",
                                  normalizePathSeparators(parentPath(m_sourcePath)).c_str());
            }

            const ExportFormatOption& currentFormat = kExportFormats[m_exportFormatIndex];
            if (ImGui::BeginCombo("Format", currentFormat.label))
            {
                for (int i = 0; i < (int)(sizeof(kExportFormats) / sizeof(kExportFormats[0])); ++i)
                {
                    const bool selected = i == m_exportFormatIndex;
                    if (ImGui::Selectable(kExportFormats[i].label, selected))
                        m_exportFormatIndex = i;
                    if (selected)
                        ImGui::SetItemDefaultFocus();
                }
                ImGui::EndCombo();
            }

            if (ImGui::Button("Export All Structures", ImVec2(-1.0f, 0.0f)))
                exportSequence();
        }
    }
    ImGui::EndChild();

    ImGui::SameLine();

    ImGui::BeginChild("##sf_output", ImVec2(0, contentHeight), false);
    {
        ImGui::Text("Generated Structures");
        if (m_result.success && !m_result.sequence.empty())
        {
            ImGui::SameLine();
            ImGui::TextDisabled("%s", m_result.sequence[m_selectedStructureIndex].label.c_str());
        }

        const float previewHeight = ImGui::GetContentRegionAvail().y - 2.0f;
        const float previewWidth = ImGui::GetContentRegionAvail().x;
        const ImVec2 minCorner = ImGui::GetCursorScreenPos();
        const ImVec2 maxCorner(minCorner.x + previewWidth, minCorner.y + previewHeight);

        ImGui::InvisibleButton("##sf_output_preview", ImVec2(previewWidth, previewHeight));
        const bool hovered = ImGui::IsItemHovered();
        const bool active = ImGui::IsItemActive();
        ImDrawList* drawList = ImGui::GetWindowDrawList();
        const ImVec4& bg = ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
        drawList->AddRectFilled(minCorner, maxCorner,
                                IM_COL32((int)(bg.x * 255), (int)(bg.y * 255), (int)(bg.z * 255), 255));
        drawList->AddRect(minCorner, maxCorner, ImGui::GetColorU32(ImGuiCol_Separator), 3.0f);

        if (!m_result.success || m_result.sequence.empty())
        {
            drawPreviewHint(drawList, minCorner, maxCorner,
                            "Generate a stacking-fault sequence",
                            "to preview intermediate structures");
        }
        else
        {
            if (m_outputPreview.dirty)
            {
                rebuildPreviewBuffers(m_outputPreview,
                                      m_result.sequence[m_selectedStructureIndex].structure,
                                      elementRadii,
                                      elementShininess);
            }

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
