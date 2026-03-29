#include "app/WindowSetup.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>

#include "Camera.h"

GLFWwindow* createMainWindow()
{
    if (!glfwInit())
        return nullptr;

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    GLFWwindow* window = glfwCreateWindow(1280, 800, "AtomForge", nullptr, nullptr);
    if (!window)
    {
        // Fallback for systems without OpenGL 4.6 support.
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
        window = glfwCreateWindow(1280, 800, "AtomForge", nullptr, nullptr);
    }

    if (!window)
    {
        // Final fallback for older systems.
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
#ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
        window = glfwCreateWindow(1280, 800, "AtomForge", nullptr, nullptr);
        if (!window)
        {
            glfwTerminate();
            return nullptr;
        }
    }

    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK)
    {
        glfwDestroyWindow(window);
        glfwTerminate();
        return nullptr;
    }

    if (GLEW_VERSION_4_6)
        std::cout << "[OpenGL] Running with OpenGL 4.6 (compute path available)." << std::endl;
    else if (GLEW_VERSION_4_3)
        std::cout << "[OpenGL] Running with OpenGL 4.3+ (compute path available)." << std::endl;
    else
        std::cout << "[OpenGL] OpenGL 4.3 not available, using legacy instancing path." << std::endl;

    glEnable(GL_DEPTH_TEST);
    return window;
}

void configureCameraCallbacks(GLFWwindow* window, Camera& camera)
{
    Camera::instance = &camera;
    glfwSetMouseButtonCallback(window, Camera::mouseButton);
    glfwSetCursorPosCallback(window, Camera::cursor);
    glfwSetScrollCallback(window, Camera::scroll);
}