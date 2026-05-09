#pragma once

#include "algorithms/DislocationBuilder.h"
#include "graphics/SceneBuffers.h"
#include "graphics/ShadowMap.h"
#include "io/StructureLoader.h"

#include <GL/glew.h>
#include <glm/glm.hpp>

#include <functional>
#include <future>
#include <string>
#include <vector>

struct Renderer;
struct SphereMesh;
struct CylinderMesh;

struct DislocationBuilderDialog
{
    DislocationBuilderDialog();
    ~DislocationBuilderDialog();

    void initRenderResources(Renderer& renderer);
    void drawMenuItem(bool enabled);

    void drawDialog(Structure& structure,
                    const std::vector<glm::vec3>& elementColors,
                    const std::vector<float>& elementRadii,
                    const std::vector<float>& elementShininess,
                    const std::function<void(Structure&)>& updateBuffers);

    bool isOpen() const { return m_isOpen; }
    void feedDroppedFile(const std::string& path);

private:
    struct PreviewState
    {
        SceneBuffers buffers;
        ShadowMap shadow = {};
        GLuint fbo = 0;
        GLuint colorTex = 0;
        GLuint depthRbo = 0;
        int width = 0;
        int height = 0;
        float yaw = 45.0f;
        float pitch = 35.0f;
        float distance = 10.0f;
        bool dirty = true;
        SphereMesh* sphere = nullptr;
        CylinderMesh* cylinder = nullptr;
        GLuint dislocationLineVAO = 0;
        GLuint dislocationLineVBO = 0;
        glm::vec3 dislocationLineP0 = glm::vec3(0.0f);
        glm::vec3 dislocationLineP1 = glm::vec3(0.0f);
        size_t dislocationLoopPointCount = 0;
        bool hasDislocationLine = false;
    };

    bool tryLoadFile(const std::string& path,
                     const std::vector<float>& radii,
                     const std::vector<float>& shininess);
    void loadFromScene(const Structure& scene,
                       const std::vector<float>& radii,
                       const std::vector<float>& shininess);
    void analyzeSource(const Structure& source);
    void generateDislocation(const Structure& source);
    void pollGenerationResult(const std::vector<float>& radii,
                              const std::vector<float>& shininess);

    void releasePreview(PreviewState& preview);
    void ensurePreviewFBO(PreviewState& preview, int width, int height);
    void rebuildPreviewBuffers(PreviewState& preview,
                               const Structure& structure,
                               const std::vector<float>& radii,
                               const std::vector<float>& shininess,
                               const DislocationResult* dislocationOverlay = nullptr);
    void autoFitPreviewCamera(PreviewState& preview);
    void renderPreviewToFBO(PreviewState& preview, int width, int height);

    bool m_openRequested = false;
    bool m_isOpen = false;
    bool m_glReady = false;
    bool m_sourceLoaded = false;
    bool m_useCurrentSceneSource = false;

    Structure m_source;
    std::string m_sourceLabel;
    std::string m_pendingDropPath;

    DislocationDetectionResult m_detection;
    DislocationParams m_params;
    DislocationResult m_result;

    std::string m_statusMsg;
    bool m_statusIsError = false;

    Renderer* m_renderer = nullptr;
    PreviewState m_sourcePreview;
    PreviewState m_outputPreview;

    bool m_generationRunning = false;
    std::future<DislocationResult> m_generationFuture;
};
