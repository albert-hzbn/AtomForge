#include "Camera.h"
#include "imgui.h"

#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>

#include <algorithm>
#include <cmath>

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

        instance->yaw   -= dx * instance->sensitivity;
        instance->pitch += dy * instance->sensitivity;
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

    // Proportional zoom: each scroll step scales distance by a fixed ratio.
    // This makes zooming feel equally responsive when close up or far away.
    float factor = 1.0f - (float)y * instance->zoomSpeed * 0.15f;
    instance->distance *= factor;

    if(instance->distance < Camera::kMinDistance)
        instance->distance = Camera::kMinDistance;

    if(instance->distance > Camera::kMaxDistance)
        instance->distance = Camera::kMaxDistance;
}