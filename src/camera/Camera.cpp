#include "Camera.h"
#include "imgui.h"

#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>

#include <algorithm>
#include <cmath>

constexpr float Camera::kMinDistance;
constexpr float Camera::kMaxDistance;

Camera* Camera::instance = nullptr;

void Camera::mouseButton(GLFWwindow*,int button,int action,int)
{
    if (ImGui::GetIO().WantCaptureMouse)
        return;

    if(button == GLFW_MOUSE_BUTTON_LEFT)
    {
        if (action == GLFW_PRESS)
        {
            instance->mouseDown = true;
            instance->dragAccum = 0.0f;
        }
        else if (action == GLFW_RELEASE)
        {
            instance->mouseDown = false;
            if (instance->dragAccum < 4.0f)
            {
                instance->pendingClick = true;
                instance->clickX      = instance->lastX;
                instance->clickY      = instance->lastY;
            }
        }
    }

    if (button == GLFW_MOUSE_BUTTON_RIGHT)
    {
        if (action == GLFW_PRESS)
        {
            instance->rightMouseDown = true;
            instance->rightDragAccum = 0.0f;
        }
        else if (action == GLFW_RELEASE)
        {
            instance->rightMouseDown = false;
            if (instance->rightDragAccum < 4.0f)
            {
                instance->pendingRightClick = true;
                instance->rightClickX      = instance->lastX;
                instance->rightClickY      = instance->lastY;
            }
        }
    }

    if (button == GLFW_MOUSE_BUTTON_MIDDLE)
    {
        if (action == GLFW_PRESS)
            instance->middleMouseDown = true;
        else if (action == GLFW_RELEASE)
            instance->middleMouseDown = false;
    }
}

void Camera::cursor(GLFWwindow*,double x,double y)
{
    if (ImGui::GetIO().WantCaptureMouse)
    {
        instance->lastX = x;
        instance->lastY = y;
        return;
    }

    if(instance->mouseDown)
    {
        float dx = x - instance->lastX;
        float dy = y - instance->lastY;

        instance->dragAccum += std::abs(dx) + std::abs(dy);

        if (instance->allowOrbit)
        {
            instance->yaw   -= dx * instance->sensitivity;
            instance->pitch += dy * instance->sensitivity;
        }
    }

    // Right-drag pans unless disabled by the app (e.g. box-select mode).
    if (instance->rightMouseDown)
    {
        const float dx = (float)(x - instance->lastX);
        const float dy = (float)(y - instance->lastY);

        instance->rightDragAccum += std::abs(dx) + std::abs(dy);

        const float yawR   = glm::radians(instance->yaw);
        const float pitchR = glm::radians(instance->pitch);

        const glm::vec3 forward(
            std::cos(pitchR) * std::sin(yawR),
            std::sin(pitchR),
            std::cos(pitchR) * std::cos(yawR));
        const glm::vec3 worldUp(0.0f, 1.0f, 0.0f);
        const glm::vec3 right = glm::normalize(glm::cross(forward, worldUp));
        const glm::vec3 up    = glm::normalize(glm::cross(right, forward));

        if (instance->allowPan)
        {
            const float panScale = 0.01f;
            instance->panOffset += right * (dx * panScale);
            instance->panOffset += up    * (dy * panScale);
        }
    }

    if (instance->middleMouseDown)
    {
        const float dy = (float)(y - instance->lastY);
        // Middle-button drag zoom: proportional to distance for consistent feel.
        float factor = dy * instance->zoomSpeed * 0.004f * instance->distance;
        instance->distance += factor;
        instance->distance = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, instance->distance));
    }

    instance->lastX = x;
    instance->lastY = y;
}

void Camera::scroll(GLFWwindow*, double, double y)
{
    if (ImGui::GetIO().WantCaptureMouse)
        return;

    // Accumulate scroll delta into a velocity bucket.
    // Both discrete mouse-wheel notches (y ≈ ±1) and the many rapid
    // sub-unit events emitted by laptop touchpads land here and are
    // applied together per frame, giving smooth zoom for both devices.
    instance->scrollVelocity += (float)y * instance->zoomSpeed;
}

void Camera::applyScrollVelocity(double now)
{
    // Compute delta-time, guarding against the first-frame spike.
    const double dt = (lastFrameTime > 0.0)
                      ? std::min(now - lastFrameTime, 0.1)
                      : 0.0;
    lastFrameTime = now;

    if (std::abs(scrollVelocity) < 1e-5f)
        return;

    // Exponential decay: time constant ≈ 1/12 s (≈ 83 ms).
    // Over an infinite number of frames the total zoom consumed equals
    // the full accumulated velocity, preserving the same zoom-per-notch
    // as the original instantaneous formula regardless of frame rate.
    constexpr float kDecay = 12.0f;
    const float retain   = std::exp(-kDecay * (float)dt);
    const float consumed = scrollVelocity * (1.0f - retain);
    scrollVelocity      *= retain;
    if (std::abs(scrollVelocity) < 0.001f)
        scrollVelocity = 0.0f;

    // Apply proportional zoom: same 0.15 coefficient as before.
    const float factor = 1.0f - consumed * 0.15f;
    distance *= factor;
    distance  = std::max(kMinDistance, std::min(kMaxDistance, distance));
}