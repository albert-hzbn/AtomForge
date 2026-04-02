#include "app/SceneView.h"

#include "Camera.h"
#include "graphics/SceneBuffers.h"

#include <GLFW/glfw3.h>

#include <algorithm>
#include <cmath>
#include <glm/gtc/matrix_transform.hpp>

namespace
{
constexpr float kVerticalFovDeg = 45.0f;

float estimateSceneRadius(const SceneBuffers& sceneBuffers)
{
    float maxRadius = 0.0f;

    if (!sceneBuffers.atomPositions.empty())
    {
        for (size_t i = 0; i < sceneBuffers.atomPositions.size(); ++i)
        {
            const float radius = (i < sceneBuffers.atomRadii.size()) ? sceneBuffers.atomRadii[i] : 0.0f;
            const float distance = glm::length(sceneBuffers.atomPositions[i] - sceneBuffers.orbitCenter) + radius;
            maxRadius = std::max(maxRadius, distance);
        }
    }
    else if (!sceneBuffers.boxLines.empty())
    {
        for (const glm::vec3& point : sceneBuffers.boxLines)
        {
            const float distance = glm::length(point - sceneBuffers.orbitCenter);
            maxRadius = std::max(maxRadius, distance);
        }
    }

    return std::max(maxRadius, 1.0f);
}

glm::mat4 makeOrthographicProjection(float aspect,
                                     float cameraDistance,
                                     float sceneRadius,
                                     float verticalFov,
                                     float depthPadding)
{
    const float halfHeight = std::max(0.1f, cameraDistance * std::tan(verticalFov * 0.5f));
    const float halfWidth = halfHeight * aspect;
    const float depthRange = std::max(1000.0f, cameraDistance + sceneRadius + depthPadding);
    return glm::ortho(
        -halfWidth,
        halfWidth,
        -halfHeight,
        halfHeight,
        -depthRange,
        depthRange);
}

glm::mat4 makePerspectiveProjection(float aspect,
                                    float cameraDistance,
                                    float sceneRadius,
                                    float verticalFov,
                                    float depthPadding)
{
    const float nearestSurface = cameraDistance - sceneRadius;
    const float nearClip = (nearestSurface > 0.0f)
        ? std::max(0.01f, nearestSurface * 0.25f)
        : 0.01f;
    const float farClip = std::max(nearClip + 100.0f, cameraDistance + sceneRadius + depthPadding);

    return glm::perspective(verticalFov, aspect, nearClip, farClip);
}

glm::vec3 computeCameraOffset(float distance, float yawDeg, float pitchDeg)
{
    const float yaw = glm::radians(yawDeg);
    const float pitch = glm::radians(pitchDeg);
    return glm::vec3(
        distance * std::cos(pitch) * std::sin(yaw),
        distance * std::sin(pitch),
        distance * std::cos(pitch) * std::cos(yaw));
}

    glm::vec3 computeCameraUp(float yawDeg, float pitchDeg)
    {
        const float yaw = glm::radians(yawDeg);
        const float pitch = glm::radians(pitchDeg);
        return glm::normalize(glm::vec3(
        -std::sin(pitch) * std::sin(yaw),
        std::cos(pitch),
        -std::sin(pitch) * std::cos(yaw)));
    }

glm::vec3 buildLightPosition(const SceneBuffers& sceneBuffers)
{
    return sceneBuffers.orbitCenter + glm::vec3(40.0f, 40.0f, 40.0f);
}
}

bool updateViewport(GLFWwindow* window, FrameView& frame)
{
    glfwGetFramebufferSize(window, &frame.framebufferWidth, &frame.framebufferHeight);
    if (frame.framebufferWidth == 0 || frame.framebufferHeight == 0)
    {
        glfwSwapBuffers(window);
        return false;
    }

    glfwGetWindowSize(window, &frame.windowWidth, &frame.windowHeight);
    if (frame.windowWidth == 0)
        frame.windowWidth = frame.framebufferWidth;
    if (frame.windowHeight == 0)
        frame.windowHeight = frame.framebufferHeight;

    return true;
}

void buildFrameView(
    Camera& camera,
    const SceneBuffers& sceneBuffers,
    bool useOrthographicView,
    FrameView& frame)
{
    const float aspect = (float)frame.framebufferWidth / (float)frame.framebufferHeight;
    const float verticalFov = glm::radians(kVerticalFovDeg);
    const float sceneRadius = estimateSceneRadius(sceneBuffers);
    const float depthPadding = std::max(10.0f, sceneRadius * 0.25f);

    if (useOrthographicView)
    {
        frame.projection = makeOrthographicProjection(
            aspect,
            camera.distance,
            sceneRadius,
            verticalFov,
            depthPadding);
    }
    else
    {
        frame.projection = makePerspectiveProjection(
            aspect,
            camera.distance,
            sceneRadius,
            verticalFov,
            depthPadding);
    }

    const glm::vec3 cameraOffset = computeCameraOffset(camera.distance, camera.yaw, camera.pitch);

    const glm::vec3 target = sceneBuffers.orbitCenter + camera.panOffset;

    frame.cameraPosition = target + cameraOffset;
    // Apply roll: rotate the canonical up vector around the forward (look) axis.
    const glm::vec3 baseUp  = computeCameraUp(camera.yaw, camera.pitch);
    const glm::vec3 forward = glm::normalize(target - frame.cameraPosition);
    const glm::vec3 right   = glm::normalize(glm::cross(forward, baseUp));
    const float rollRad     = glm::radians(camera.roll);
    const glm::vec3 rolledUp = std::cos(rollRad) * baseUp + std::sin(rollRad) * right;
    frame.view = glm::lookAt(
        frame.cameraPosition,
        target,
        rolledUp);

    const glm::mat4 lightProjection = glm::perspective(glm::radians(kVerticalFovDeg), 1.0f, 0.1f, 1000.0f);
    frame.lightPosition = buildLightPosition(sceneBuffers);
    frame.lightMVP = lightProjection * glm::lookAt(
        frame.lightPosition,
        sceneBuffers.orbitCenter,
        glm::vec3(0, 1, 0));
}
