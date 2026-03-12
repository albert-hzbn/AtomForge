#include "Camera.h"
#include <cmath>

Camera* Camera::instance = nullptr;

void Camera::mouseButton(GLFWwindow*,int button,int action,int)
{
    if(button == GLFW_MOUSE_BUTTON_LEFT)
        instance->mouseDown = action == GLFW_PRESS;
}

void Camera::cursor(GLFWwindow*,double x,double y)
{
    if(instance->mouseDown)
    {
        float dx = x - instance->lastX;
        float dy = y - instance->lastY;

        instance->yaw   += dx * instance->sensitivity;
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
    instance->distance -= y * instance->zoomSpeed;

    if(instance->distance < 2.0f)
        instance->distance = 2.0f;

    if(instance->distance > 50.0f)
        instance->distance = 50.0f;
}