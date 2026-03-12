#pragma once
#include <GLFW/glfw3.h>

class Camera
{
public:

    // trackball parameters
    float yaw   = 0.0f;
    float pitch = 20.0f;

    float distance = 10.0f;

    float sensitivity = 0.3f;
    float zoomSpeed   = 0.5f;

    bool mouseDown = false;

    double lastX = 0;
    double lastY = 0;

    static Camera* instance;

    static void mouseButton(GLFWwindow*,int button,int action,int);
    static void cursor(GLFWwindow*,double x,double y);
    static void scroll(GLFWwindow*,double,double y);

};