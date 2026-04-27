#pragma once

#include <glm/glm.hpp>

struct GLFWwindow;

class Camera
{
public:
    static constexpr float kMinDistance = 2.0f;
    static constexpr float kMaxDistance = 5000.0f;

    // trackball parameters
    float yaw   = 45.0f;
    float pitch = 35.2643897f;
    float roll  = 0.0f;  // degrees, right-hand rotation around forward axis (set by orbital rotation)

    float distance = 10.0f;

    float sensitivity = 0.3f;
    float zoomSpeed   = 0.5f;

    // Accumulated scroll velocity for smooth zoom (supports both mouse wheel
    // and touchpad precision scrolling). Applied and decayed each frame.
    float scrollVelocity = 0.0f;
    double lastFrameTime = 0.0;  // glfwGetTime() at start of previous frame

    // Must be called once per frame (before buildFrameView) to apply
    // accumulated scroll velocity with exponential decay.
    void applyScrollVelocity(double now);

    bool mouseDown = false;
    bool middleMouseDown = false;
    bool rightMouseDown = false;

    glm::vec3 panOffset = glm::vec3(0.0f);
    bool allowPan = true;
    bool allowOrbit = true;

    double lastX = 0;
    double lastY = 0;

    // Accumulated Manhattan drag distance since last left-press.
    float dragAccum = 0.0f;

    // Accumulated Manhattan drag distance since last right-press.
    float rightDragAccum = 0.0f;

    // Set to true for one frame when a short left-click is detected.
    bool   pendingClick = false;
    double clickX       = 0;
    double clickY       = 0;

    // Set to true for one frame when right button is released.
    bool   pendingRightClick = false;
    double rightClickX       = 0;
    double rightClickY       = 0;

    static Camera* instance;

    static void mouseButton(GLFWwindow*,int button,int action,int);
    static void cursor(GLFWwindow*,double x,double y);
    static void scroll(GLFWwindow*,double,double y);

};