#pragma once

#include "ui/FileBrowser.h"
#include "model/Structure.h"

#include <glm/glm.hpp>

#include <string>

class CloudRenderer;
#include <vector>

struct Renderer;
struct ShadowMap;
struct SceneBuffers;

struct ImageExportView
{
    int width = 0;
    int height = 0;
    glm::mat4 projection = glm::mat4(1.0f);
    glm::mat4 view = glm::mat4(1.0f);
    glm::mat4 lightMVP = glm::mat4(1.0f);
    glm::vec3 lightPosition = glm::vec3(0.0f);
    glm::vec3 cameraPosition = glm::vec3(0.0f);
};

// Renders the scene off-screen into RGBA pixels (rows top to bottom), as used
// for raster image export; GIF animation frames are captured this way.
// A large dataset to draw instead of the scene buffers in exported images
// (nullptr for normal structures). The export waits for its atoms to stream in.
void setImageExportCloud(CloudRenderer* cloud);

bool renderSceneToRgba(const ImageExportView& view,
                       const glm::vec4& backgroundColor,
                       bool showBonds,
                       bool showAtoms,
                       bool showBoundingBox,
                       const SceneBuffers& sceneBuffers,
                       Renderer& renderer,
                       const ShadowMap& shadow,
                       std::vector<unsigned char>& rgbaPixels,
                       std::string& errorMessage);

bool exportStructureImage(const ImageExportRequest& request,
                          const ImageExportView& view,
                          const glm::vec4& backgroundColor,
                          bool showBonds,
                          bool showAtoms,
                          bool showBoundingBox,
                          const SceneBuffers& sceneBuffers,
                          Renderer& renderer,
                          const ShadowMap& shadow,
                          const Structure& structure,
                          std::string& errorMessage);
