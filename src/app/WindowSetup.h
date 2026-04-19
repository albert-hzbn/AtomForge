#pragma once

struct GLFWwindow;
struct Camera;
struct SplashScreen;

SplashScreen* createSplashScreen();
void updateSplashScreen(SplashScreen* splash, float progress, const char* status);
void destroySplashScreen(SplashScreen* splash);
GLFWwindow* createMainWindow();
void showMainWindow(GLFWwindow* window);
void configureCameraCallbacks(GLFWwindow* window, Camera& camera);