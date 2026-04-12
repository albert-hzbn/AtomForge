#pragma once

struct GLFWwindow;

// Initialise Dear ImGui, install GLFW/OpenGL3 backends, and apply the
// application style.  Must be called after a valid OpenGL context exists.
void initImGui(GLFWwindow* window);

// Detach backends and destroy the ImGui context.
void shutdownImGui();

// Apply the dark (default) or light theme.  Safe to call at any time.
void applyDarkTheme();
void applyLightTheme();
