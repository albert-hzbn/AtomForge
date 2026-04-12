#pragma once

#include "graphics/SceneBuffers.h"
#include "graphics/ShadowMap.h"
#include "io/StructureLoader.h"

#include <functional>
#include <string>
#include <vector>

#include <GL/glew.h>
#include <glm/glm.hpp>

struct Renderer;
struct SphereMesh;
struct CylinderMesh;

struct SingleCrystalFillDialog
{
    SingleCrystalFillDialog();
    ~SingleCrystalFillDialog();

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
    bool m_openRequested = false;
    bool m_isOpen = false;

    Structure m_reference;
    std::string m_refFileName;
    std::vector<glm::vec3> m_modelVertices;
    std::vector<unsigned int> m_modelIndices;
    glm::vec3 m_modelHalfExtents = glm::vec3(15.0f);

    std::string m_modelName;
    std::vector<std::string> m_pendingDropPaths;
    std::string m_status;

    // GL preview resources
    Renderer*     m_renderer        = nullptr;
    SphereMesh*   m_previewSphere   = nullptr;
    CylinderMesh* m_previewCylinder = nullptr;
    SceneBuffers  m_previewBuffers;
    ShadowMap     m_previewShadow   = {};
    bool          m_glReady         = false;
    bool          m_previewBufDirty = true;

    // Structure preview FBO
    GLuint m_refFBO      = 0;
    GLuint m_refColorTex = 0;
    GLuint m_refDepthRbo = 0;
    int    m_refFboW     = 0;
    int    m_refFboH     = 0;
    float  m_refCamYaw      = 45.0f;
    float  m_refCamPitch    = 35.0f;
    float  m_refCamDistance  = 10.0f;

    // Model preview FBO
    GLuint m_modelFBO      = 0;
    GLuint m_modelColorTex = 0;
    GLuint m_modelDepthRbo = 0;
    int    m_modelFboW     = 0;
    int    m_modelFboH     = 0;
    float  m_modelCamYaw      = 45.0f;
    float  m_modelCamPitch    = 35.0f;
    float  m_modelCamDistance  = 50.0f;

    // Model surface mesh GL objects
    GLuint m_meshProgram = 0;
    GLuint m_meshVAO = 0;
    GLuint m_meshVBO = 0;
    GLuint m_meshIBO = 0;
    int    m_meshTriCount = 0;

    bool loadReferenceStructure(const std::string& path);
    bool loadModelMesh(const std::string& path);

    void ensureRefFBO(int w, int h);
    void ensureModelFBO(int w, int h);
    void rebuildPreviewBuffers(const std::vector<float>& radii,
                               const std::vector<float>& shininess);
    void autoFitRefCamera();
    void renderRefPreviewToFBO(int w, int h);
    void rebuildMeshSurface();
    void autoFitModelCamera();
    void renderModelPreviewToFBO(int w, int h);
};
