#include "Camera.h"
#include "imgui.h"
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

    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE)
    {
        instance->pendingRightClick = true;
        instance->rightClickX      = instance->lastX;
        instance->rightClickY      = instance->lastY;
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

        // limit pitch so camera never flips
        if(instance->pitch > 89.0f) instance->pitch = 89.0f;
        if(instance->pitch < -89.0f) instance->pitch = -89.0f;
    }

    instance->lastX = x;
    instance->lastY = y;
}

void Camera::scroll(GLFWwindow*, double, double y)
{
    if (ImGui::GetIO().WantCaptureMouse)
        return;

    instance->distance -= y * instance->zoomSpeed;

    if(instance->distance < Camera::kMinDistance)
        instance->distance = Camera::kMinDistance;

    if(instance->distance > Camera::kMaxDistance)
        instance->distance = Camera::kMaxDistance;
}