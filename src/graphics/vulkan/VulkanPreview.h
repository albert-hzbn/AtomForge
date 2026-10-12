#pragma once
// What the dialogs' Vulkan previews share: the view and lights of a
// FrameView and the Renderer's settings, and the standard preview of a
// structure (bonds, atoms and the cell box, as the OpenGL previews draw it).

#ifdef ATOMFORGE_VULKAN
#include "app/SceneView.h"
#include "graphics/Renderer.h"
#include "graphics/vulkan/VulkanCanvas.h"

#include "imgui.h"

namespace vkr
{
inline SceneView viewOf(const FrameView& frame)
{
    SceneView view;
    view.projection = frame.projection;
    view.view = frame.view;
    view.lightMVP = frame.lightMVP;
    view.lightPosition = frame.lightPosition;
    view.cameraPosition = frame.cameraPosition;
    return view;
}

// The main view's lighting and materials; previews have no shadows.
inline SceneLighting lightingOf(const Renderer* renderer)
{
    SceneLighting lighting;
    if (renderer)
    {
        lighting.ambient = renderer->lightAmbient;
        lighting.saturation = renderer->lightSaturation;
        lighting.contrast = renderer->lightContrast;
        lighting.specularIntensity = renderer->materialSpecularIntensity;
        lighting.shininessScale = renderer->materialShininessScale;
        lighting.shininessFloor = renderer->materialShininessFloor;
    }
    lighting.shadowStrength = 0.0f;
    return lighting;
}

inline glm::vec3 colorOf(const ImVec4& color)
{
    return glm::vec3(color.x, color.y, color.z);
}

// Starts a preview of frame's size and view; end() it after drawing.
inline bool beginPreview(Canvas& canvas, const Renderer* renderer, const FrameView& frame, const ImVec4& background)
{
    if (!canvas.begin(frame.framebufferWidth, frame.framebufferHeight, colorOf(background))) return false;
    canvas.setView(viewOf(frame), lightingOf(renderer));
    return true;
}

// Bonds, atoms (full spheres) and the cell box.
inline void drawStructure(Canvas& canvas, const SceneBuffers& buffers, bool box = true)
{
    SceneFlags flags;
    flags.spheresOnly = true;
    flags.box = box;
    canvas.drawScene(buffers, flags, glm::vec3(0.85f));
}

inline void drawStandardPreview(Canvas& canvas, const Renderer* renderer, const SceneBuffers& buffers, const FrameView& frame,
                                const ImVec4& background)
{
    if (!beginPreview(canvas, renderer, frame, background)) return;
    drawStructure(canvas, buffers);
    canvas.end();
}
}
#endif
