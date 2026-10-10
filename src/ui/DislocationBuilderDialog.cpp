#include "ui/ResponsiveLayout.h"
#include "ui/DislocationBuilderDialog.h"

#include "app/SceneView.h"
#include "camera/Camera.h"
#include "graphics/CylinderMesh.h"
#include "graphics/Renderer.h"
#include "graphics/SphereMesh.h"
#include "graphics/StructureInstanceBuilder.h"
#include "io/StructureLoader.h"
#include "pipeline/Options.h"

#include "imgui.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cctype>
#include <cmath>
#include <future>
#include <sstream>

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
    {
        m_step.finish();  // from the menu: build normally
        m_stepLoadInput = false;
        m_openRequested = true;
    }
}

bool DislocationBuilderDialog::editStep(StepEdit edit)
{
    m_step = std::move(edit);
    const auto options = m_step.reader();
    const auto lower = [](std::string text) {
        for (char& c : text) c = (char)std::tolower((unsigned char)c);
        return text;
    };
    const auto ivec3 = [&](const char* flag, glm::ivec3 fallback) {
        const auto v = options.numbers(flag);
        return v.size() == 3 ? glm::ivec3((int)std::lround(v[0]), (int)std::lround(v[1]), (int)std::lround(v[2])) : fallback;
    };
    const auto vec3 = [&](const char* flag, glm::vec3 fallback) {
        const auto v = options.numbers(flag);
        return v.size() == 3 ? glm::vec3((float)v[0], (float)v[1], (float)v[2]) : fallback;
    };

    // Missing flags take the defaults of AtomForge --build dislocation.
    DislocationParams p;
    const std::string character = lower(options.text("--character", "edge"));
    p.character = character == "screw" ? DislocationCharacter::Screw
                : character == "mixed" ? DislocationCharacter::Mixed : DislocationCharacter::Edge;
    const std::string shape = lower(options.text("--shape", "halfplane"));
    p.shape = shape == "cylinder" ? DislocationShape::Cylinder
            : shape == "sphere" ? DislocationShape::Sphere
            : shape == "ellipsoid" ? DislocationShape::Ellipsoid
            : shape == "freeform" ? DislocationShape::Freeform2D : DislocationShape::HalfPlane;
    p.autoDirections = !options.has("--manual-vectors");
    p.planeHkl = ivec3("--plane", p.planeHkl);
    p.burgersUvw = ivec3("--burgers", p.burgersUvw);
    p.lineUvw = ivec3("--line", p.lineUvw);
    p.linePointFractional = vec3("--line-frac", p.linePointFractional);
    p.useFractionalLinePoint = !options.has("--line-center");
    p.linePointCartesianOffset = vec3("--line-offset", p.linePointCartesianOffset);
    p.burgersScale = (float)options.number("--bscale", p.burgersScale);
    p.burgersOverrideMagnitude = (float)options.number("--bmag", p.burgersOverrideMagnitude);
    p.mixedCharacterAngleDeg = (float)options.number("--mixed-angle", p.mixedCharacterAngleDeg);
    p.poissonRatio = (float)options.number("--nu", p.poissonRatio);
    p.coreRadius = (float)options.number("--core", p.coreRadius);
    p.cutoffRadius = (float)options.number("--cutoff", p.cutoffRadius);
    p.lineHalfLength = (float)options.number("--line-half", p.lineHalfLength);
    p.cylinderRadius = (float)options.number("--cyl-radius", p.cylinderRadius);
    p.sphereRadius = (float)options.number("--sphere-radius", p.sphereRadius);
    p.ellipsoidRadii = vec3("--ellipsoid", p.ellipsoidRadii);
    if (options.has("--poly2d"))
    {
        // "x1 y1;x2 y2;..."
        std::vector<glm::vec2> points;
        std::stringstream all(options.text("--poly2d"));
        for (std::string point; std::getline(all, point, ';');)
        {
            std::istringstream xy(point);
            float x = 0.0f, y = 0.0f;
            if (xy >> x >> y) points.push_back(glm::vec2(x, y));
        }
        if (points.size() >= 3) p.freeformPoints = points;
    }
    p.dipole = options.has("--dipole");
    {
        const auto offset = options.numbers("--dipole-offset");
        if (offset.size() == 2) p.dipoleOffset = glm::vec2((float)offset[0], (float)offset[1]);
    }
    p.anisotropicElasticity = options.has("--anisotropic");
    p.elasticSymmetry = lower(options.text("--elastic-symmetry", "cubic")) == "hexagonal"
        ? DislocationParams::ElasticSymmetry::Hexagonal : DislocationParams::ElasticSymmetry::Cubic;
    p.elasticC11 = options.number("--elastic-c11", 0.0);
    p.elasticC12 = options.number("--elastic-c12", 0.0);
    p.elasticC44 = options.number("--elastic-c44", 0.0);
    p.elasticC13 = options.number("--elastic-c13", 0.0);
    p.elasticC33 = options.number("--elastic-c33", 0.0);
    p.elasticNoiseAmplitude = options.number("--elastic-noise", p.elasticNoiseAmplitude);
    m_params = p;

    // Source: a file named in the options, else the structure entering the step
    // (both are loaded by drawDialog, which has the element radii for the preview).
    m_source = Structure();
    m_sourceLoaded = false;
    m_useCurrentSceneSource = false;
    m_sourceLabel.clear();
    m_detection = {};
    m_result = {};
    m_statusMsg.clear();
    m_statusIsError = false;
    m_sourcePath = options.text("--input");
    const auto* useInput = m_step.parameters.find("use_input");
    m_stepLoadInput = m_sourcePath.empty() && (!useInput || !useInput->isBool() || useInput->boolean())
                      && !m_step.input.atoms.empty();
    if (!m_sourcePath.empty())
        m_pendingDropPath = m_sourcePath;
    m_openRequested = true;
    return true;
}

std::string DislocationBuilderDialog::stepOptions() const
{
    const DislocationParams& p = m_params;
    const auto ivec3 = [](const glm::ivec3& v) { return std::vector<double>{(double)v.x, (double)v.y, (double)v.z}; };
    const auto vec3 = [](const glm::vec3& v) { return std::vector<double>{v.x, v.y, v.z}; };
    atomforge::pipeline::options::Writer out;
    if (!m_sourcePath.empty()) out.add("--input", m_sourcePath);
    out.add("--character", p.character == DislocationCharacter::Screw ? "screw"
                         : p.character == DislocationCharacter::Mixed ? "mixed" : "edge");
    out.add("--shape", p.shape == DislocationShape::Cylinder ? "cylinder"
                     : p.shape == DislocationShape::Sphere ? "sphere"
                     : p.shape == DislocationShape::Ellipsoid ? "ellipsoid"
                     : p.shape == DislocationShape::Freeform2D ? "freeform" : "halfplane");
    if (!p.autoDirections) out.flag("--manual-vectors");
    out.add("--plane", ivec3(p.planeHkl)).add("--burgers", ivec3(p.burgersUvw)).add("--line", ivec3(p.lineUvw));
    if (p.useFractionalLinePoint) out.add("--line-frac", vec3(p.linePointFractional));
    else out.flag("--line-center");
    out.add("--line-offset", vec3(p.linePointCartesianOffset));
    out.add("--bscale", p.burgersScale);
    if (p.burgersOverrideMagnitude > 0.0f) out.add("--bmag", p.burgersOverrideMagnitude);
    if (p.character == DislocationCharacter::Mixed) out.add("--mixed-angle", p.mixedCharacterAngleDeg);
    out.add("--nu", p.poissonRatio).add("--core", p.coreRadius).add("--cutoff", p.cutoffRadius);
    out.add("--line-half", p.lineHalfLength);

    // Region parameters of the chosen shape.
    if (p.shape == DislocationShape::Cylinder) out.add("--cyl-radius", p.cylinderRadius);
    else if (p.shape == DislocationShape::Sphere) out.add("--sphere-radius", p.sphereRadius);
    else if (p.shape == DislocationShape::Ellipsoid) out.add("--ellipsoid", vec3(p.ellipsoidRadii));
    else if (p.shape == DislocationShape::Freeform2D && p.freeformPoints.size() >= 3)
    {
        std::string polygon;
        for (const glm::vec2& point : p.freeformPoints)
            polygon += (polygon.empty() ? "" : ";") + out.number(point.x) + " " + out.number(point.y);
        out.add("--poly2d", polygon);
    }

    if (p.dipole)
        out.flag("--dipole").add("--dipole-offset", std::vector<double>{p.dipoleOffset.x, p.dipoleOffset.y});
    if (p.anisotropicElasticity)
    {
        const bool hexagonal = p.elasticSymmetry == DislocationParams::ElasticSymmetry::Hexagonal;
        out.flag("--anisotropic").add("--elastic-symmetry", hexagonal ? "hexagonal" : "cubic");
        out.add("--elastic-c11", p.elasticC11).add("--elastic-c12", p.elasticC12).add("--elastic-c44", p.elasticC44);
        if (hexagonal) out.add("--elastic-c13", p.elasticC13).add("--elastic-c33", p.elasticC33);
        out.add("--elastic-noise", p.elasticNoiseAmplitude);
    }
    return out.str();
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
    m_sourcePath = path;

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
    m_sourcePath.clear();
    m_sourceLabel = "scene";

    analyzeSource(scene);
    rebuildPreviewBuffers(m_sourcePreview, scene, radii, shininess);
    autoFitPreviewCamera(m_sourcePreview);

    m_result = {};
    m_outputPreview.dirty = true;
}

void DislocationBuilderDialog::useStepInput(const std::vector<float>& radii,
                                            const std::vector<float>& shininess)
{
    loadFromScene(m_step.input, radii, shininess);  // checks it and shows the preview
    if (m_step.input.atoms.empty() || !m_step.input.hasUnitCell)
        return;
    m_source = m_step.input;  // a copy: the dialog's scene argument is not the step input
    m_useCurrentSceneSource = false;
    m_sourceLabel = "step input";
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
    if (m_stepLoadInput)
    {
        m_stepLoadInput = false;
        useStepInput(elementRadii, elementShininess);
    }

    if (m_openRequested)
    {
        ImGui::OpenPopup("Insert Dislocation");
        m_openRequested = false;
    }

    m_isOpen = ImGui::IsPopupOpen("Insert Dislocation");

    responsive::windowSize(ImVec2(1240.0f, 780.0f), ImGuiCond_FirstUseEver);
    responsive::windowConstraints(ImVec2(960.0f, 600.0f), ImVec2(3200.0f, 3200.0f));

    bool keepOpen = true;
    if (!responsive::beginModal("Insert Dislocation", &keepOpen, ImGuiWindowFlags_NoCollapse))
    {
        m_isOpen = false;
        m_step.finish();
        return;
    }
    m_isOpen = true;

    if (m_step.active())
        ImGui::TextColored(ImVec4(0.45f, 0.75f, 1.0f, 1.0f),
                           "Editing pipeline step %s: Update step writes these settings into the pipeline.", m_step.step.c_str());

    ImGui::TextDisabled("Isotropic-elastic dislocation insertion with automatic FCC/HCP/BCC detection and shape-controlled application region.");
    ImGui::SameLine();
    ImGui::TextColored(validationStateColor(m_result), "%s", validationStateLabel(m_result));
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    const Structure* activeSource = nullptr;
    if (m_sourceLoaded)
        activeSource = m_useCurrentSceneSource ? &structure : &m_source;

    const bool stackPanels = responsive::stacked();
    const float contentHeight = responsive::panelHeight(780);
    const float previewPanelWidth = responsive::previewWidth(680, 440);

    responsive::beginChild("##disloc_previews", ImVec2(previewPanelWidth, contentHeight), false);
    {
        const float splitGap = ImGui::GetStyle().ItemSpacing.y;
        const float upperHeight = (ImGui::GetContentRegionAvail().y - splitGap) * 0.5f;
        const float lowerHeight = ImGui::GetContentRegionAvail().y - upperHeight - splitGap;

        responsive::beginChild("##disloc_source", ImVec2(0.0f, upperHeight), false);
        {
            ImGui::Text("Input Preview");
            ImGui::SameLine();
            ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, responsive::size(10.0f,6.0f));
            if (m_step.active())
            {
                if (responsive::button("Use Step Input"))
                {
                    useStepInput(elementRadii, elementShininess);
                    if (m_sourceLoaded) activeSource = &m_source;
                }
            }
            else if (responsive::button("Use Current Scene"))
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

        responsive::beginChild("##disloc_output", ImVec2(0.0f, lowerHeight), false);
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

    responsive::nextPanel(stackPanels);

    responsive::beginChild("##disloc_controls", ImVec2(0, contentHeight), false);
    {
        ImGui::Text("Dislocation Options");
        ImGui::Separator();

        const bool canGenerate = (activeSource != nullptr) && !m_generationRunning
            && (!m_params.anisotropicElasticity || (m_params.elasticC11 > 0.0 && m_params.elasticC44 > 0.0));
        if (m_generationRunning)
            ImGui::TextDisabled("Generation running in background...");

        if (responsive::button("Detect Lattice", responsive::size(-1.0f,0.0f)))
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
        if (m_params.anisotropicElasticity)
            ImGui::BeginDisabled();
        ImGui::DragFloat("Poisson ratio", &m_params.poissonRatio, 0.005f, 0.05f, 0.49f, "%.3f");
        if (m_params.anisotropicElasticity)
            ImGui::EndDisabled();
        ImGui::DragFloat("Core radius (A)", &m_params.coreRadius, 0.05f, 0.05f, 20.0f, "%.3f");
        ImGui::DragFloat("Cutoff radius (A)", &m_params.cutoffRadius, 0.2f, 0.0f, 100.0f, "%.2f");
        ImGui::DragFloat("Line half length (A)", &m_params.lineHalfLength, 0.5f, 0.1f, 1000000.0f, "%.1f");

        if (m_params.character == DislocationCharacter::Mixed)
            ImGui::SliderFloat("Mixed angle (deg)", &m_params.mixedCharacterAngleDeg, 0.0f, 90.0f, "%.1f");

        ImGui::Separator();
        ImGui::TextDisabled("Anisotropic Elasticity (Stroh formalism)");
        ImGui::Checkbox("Use anisotropic elasticity", &m_params.anisotropicElasticity);
        if (m_params.anisotropicElasticity)
        {
            static const char* kSymmetryLabels[] = {"Cubic", "Hexagonal"};
            int symmetryIndex = (m_params.elasticSymmetry == DislocationParams::ElasticSymmetry::Hexagonal) ? 1 : 0;
            if (ImGui::Combo("Elastic symmetry", &symmetryIndex, kSymmetryLabels, 2))
                m_params.elasticSymmetry = (symmetryIndex == 1)
                    ? DislocationParams::ElasticSymmetry::Hexagonal
                    : DislocationParams::ElasticSymmetry::Cubic;

            float c11 = (float)m_params.elasticC11;
            float c12 = (float)m_params.elasticC12;
            float c44 = (float)m_params.elasticC44;
            if (ImGui::DragFloat("C11 (GPa)", &c11, 1.0f, 0.0f, 2000.0f, "%.1f")) m_params.elasticC11 = c11;
            if (ImGui::DragFloat("C12 (GPa)", &c12, 1.0f, 0.0f, 2000.0f, "%.1f")) m_params.elasticC12 = c12;
            if (ImGui::DragFloat("C44 (GPa)", &c44, 1.0f, 0.0f, 2000.0f, "%.1f")) m_params.elasticC44 = c44;

            if (m_params.elasticSymmetry == DislocationParams::ElasticSymmetry::Hexagonal)
            {
                float c13 = (float)m_params.elasticC13;
                float c33 = (float)m_params.elasticC33;
                if (ImGui::DragFloat("C13 (GPa)", &c13, 1.0f, 0.0f, 2000.0f, "%.1f")) m_params.elasticC13 = c13;
                if (ImGui::DragFloat("C33 (GPa)", &c33, 1.0f, 0.0f, 2000.0f, "%.1f")) m_params.elasticC33 = c33;
                ImGui::TextWrapped("The structure's THIRD cell vector must be the c-axis.");
            }

            float noise = (float)m_params.elasticNoiseAmplitude;
            if (ImGui::DragFloat("Noise amplitude", &noise, 0.0001f, 0.0f, 0.1f, "%.4f")) m_params.elasticNoiseAmplitude = noise;
            ImGui::TextWrapped("Breaks the sextic formalism's degeneracy at isotropic or high-symmetry orientations (e.g. line = plane normal).");

            if (m_params.elasticC11 <= 0.0 || m_params.elasticC44 <= 0.0)
                ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.2f, 1.0f), "C11 and C44 must be positive to generate.");
        }

        ImGui::Separator();
        ImGui::TextDisabled("Dipole");
        ImGui::Checkbox("Insert as dipole (adds an opposite-Burgers-vector partner)", &m_params.dipole);
        if (m_params.dipole)
        {
            ImGui::DragFloat2("Dipole offset (A)", &m_params.dipoleOffset.x, 0.2f, -500.0f, 500.0f, "%.2f");
            ImGui::TextWrapped("Offset of the second (opposite-sign) dislocation, in the local slip-plane axes. A dipole's net Burgers vector is zero, so it is periodicity-compatible unlike a single dislocation.");
        }

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
            if (responsive::button("Add Point"))
            {
                if (m_params.freeformPoints.empty())
                    m_params.freeformPoints.push_back(glm::vec2(0.0f, 0.0f));
                else
                    m_params.freeformPoints.push_back(m_params.freeformPoints.back() + glm::vec2(1.0f, 0.0f));
            }
            ImGui::SameLine();
            if (responsive::button("Remove Last") && !m_params.freeformPoints.empty())
                m_params.freeformPoints.pop_back();

            const float pointsHeight = 112.0f;
            if (responsive::beginChild("##freeform_points", ImVec2(0.0f, pointsHeight), true))
            {
                for (int i = 0; i < (int)m_params.freeformPoints.size(); ++i)
                {
                    ImGui::PushID(i + 24000);
                    ImGui::SetNextItemWidth(responsive::dp(-1.0f));
                    ImGui::DragFloat2("##pt", &m_params.freeformPoints[i].x, 0.05f, -500.0f, 500.0f, "%.2f");
                    ImGui::PopID();
                }
                ImGui::EndChild();
            }
        }

        ImGui::Spacing();
        if (!canGenerate)
            ImGui::BeginDisabled();
        if (responsive::button("Generate Dislocation", responsive::size(-1.0f,0.0f)))
            generateDislocation(*activeSource);
        if (!canGenerate)
            ImGui::EndDisabled();

        // Editing a step: the settings go into the pipeline, the scene is left alone.
        if (m_step.active() && responsive::button("Update step", responsive::size(-1.0f,0.0f)))
        {
            atomforge::pipeline::Json parameters = m_step.parameters;
            parameters["options"] = stepOptions();
            parameters["use_input"] = m_sourcePath.empty();  // the step input, unless a file was loaded
            m_step.commitParameters(parameters);
            m_step.finish();
            ImGui::CloseCurrentPopup();
        }

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

            if (!m_step.active() && responsive::button("Replace Main Scene", responsive::size(-1.0f,0.0f)))
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
    if (responsive::button("Close", ImVec2(buttonWidth, 0.0f)))
    {
        m_step.finish();
        ImGui::CloseCurrentPopup();
    }

    ImGui::EndPopup();
}
