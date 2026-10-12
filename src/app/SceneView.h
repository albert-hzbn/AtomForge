#pragma once

#include <glm/glm.hpp>

struct GLFWwindow;
struct Camera;
struct SceneBuffers;

struct FrameView
{
    int framebufferWidth = 0;
    int framebufferHeight = 0;
    int windowWidth = 0;
    int windowHeight = 0;

    glm::mat4 projection = glm::mat4(1.0f);
    // The projection for the whole window (image export), before
    // fitProjectionToArea fits the scene into the part of it left free by
    // docked windows.
    glm::mat4 fullProjection = glm::mat4(1.0f);
    // That part, in window coordinates (x0, y0, x1, y1); all zero: the whole window.
    float area[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    glm::mat4 view = glm::mat4(1.0f);
    glm::mat4 lightMVP = glm::mat4(1.0f);
    glm::vec3 lightPosition = glm::vec3(0.0f);
    glm::vec3 cameraPosition = glm::vec3(0.0f);
};

bool updateViewport(GLFWwindow* window, FrameView& frame);

// Keeps the projection for the whole window in fullProjection and fits the
// scene into frame.area: scaled by the area's share of the height and centred
// in it. Drawing, picking and overlays all use the one projection, so they
// stay in the whole-window coordinates they already use.
void fitProjectionToArea(FrameView& frame);

void buildFrameView(
    Camera& camera,
    const SceneBuffers& sceneBuffers,
    bool useOrthographicView,
    FrameView& frame);
