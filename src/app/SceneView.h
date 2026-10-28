#pragma once

#include "model/DisplaySettings.h"

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
    // Camera distances of the nearest and farthest scene points (bounding sphere), for depth cueing.
    float depthNear = 0.0f, depthFar = 1.0f;
};

struct Renderer;
// The display settings the renderer applies per frame (lighting, depth cueing,
// outlines, cell lines, selection rings); background is the clear colour.
void applyDisplaySettings(Renderer& renderer, const atomforge::display::DisplaySettings& display,
                          const FrameView& frame, const glm::vec3& background);

bool updateViewport(GLFWwindow* window, FrameView& frame);

// Keeps the projection for the whole window in fullProjection and fits the
// scene into frame.area: scaled by the area's share of the height and centred
// in it. Drawing, picking and overlays all use the one projection, so they
// stay in the whole-window coordinates they already use.
void fitProjectionToArea(FrameView& frame);

// The light comes from the direction the display settings give, fixed in
// the scene or turning with the camera.
void buildFrameView(
    Camera& camera,
    const SceneBuffers& sceneBuffers,
    bool useOrthographicView,
    FrameView& frame,
    const atomforge::display::DisplaySettings& display = atomforge::display::DisplaySettings());
