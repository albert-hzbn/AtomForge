#pragma once

#include "algorithms/StackingFaultBuilder.h"
#include "graphics/SceneBuffers.h"
#include "graphics/ShadowMap.h"
#include "io/StructureLoader.h"

#include <GL/glew.h>
#include <glm/glm.hpp>

#include <functional>
#include <string>
#include <vector>

struct Renderer;
struct SphereMesh;
struct CylinderMesh;

struct StackingFaultBuilderDialog
{
    StackingFaultBuilderDialog();
    ~StackingFaultBuilderDialog();

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
            SphereMesh*   sphere   = nullptr;
            CylinderMesh* cylinder = nullptr;
    };

    bool tryLoadFile(const std::string& path,
                     const std::vector<float>& radii,
                     const std::vector<float>& shininess);
    void loadFromScene(const Structure& scene,
                       const std::vector<float>& radii,
                       const std::vector<float>& shininess);
    void analyzeSource();
    void regenerateSequence(const std::vector<float>& radii,
                            const std::vector<float>& shininess);
    bool exportSequence();

    void releasePreview(PreviewState& preview);
    void ensurePreviewFBO(PreviewState& preview, int width, int height);
    void rebuildPreviewBuffers(PreviewState& preview,
                               const Structure& structure,
                               const std::vector<float>& radii,
                               const std::vector<float>& shininess);
    void autoFitPreviewCamera(PreviewState& preview);
    void renderPreviewToFBO(PreviewState& preview, int width, int height);

    bool m_openRequested = false;
    bool m_isOpen = false;
    bool m_glReady = false;
    bool m_sourceLoaded = false;

    Structure m_source;
    std::string m_sourcePath;
    std::string m_sourceLabel;
    std::string m_pendingDropPath;

    StackingFaultDetectionResult m_detection;
    StackingFaultParams m_params;
    StackingFaultResult m_result;
    int m_selectedStructureIndex = 0;

    std::string m_statusMsg;
    bool m_statusIsError = false;

    char m_exportDir[1024];
    char m_exportFolderName[256];
    int m_exportFormatIndex = 0;

    Renderer* m_renderer = nullptr;
    PreviewState m_sourcePreview;
    PreviewState m_outputPreview;
};