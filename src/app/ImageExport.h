#pragma once

#include "ui/FileBrowser.h"

#include <glm/glm.hpp>

#include <string>

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

bool exportStructureImage(const ImageExportRequest& request,
                          const ImageExportView& view,
                          const glm::vec4& backgroundColor,
                          bool showBonds,
                          bool showAtoms,
                          bool showBoundingBox,
                          const SceneBuffers& sceneBuffers,
                          Renderer& renderer,
                          const ShadowMap& shadow,
                          std::string& errorMessage);
