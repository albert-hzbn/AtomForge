#pragma once

struct GLFWwindow;

// Initialise Dear ImGui, install GLFW/OpenGL3 backends, and apply the
// application style.  Must be called after a valid OpenGL context exists.
void initImGui(GLFWwindow* window);

// Recompute UI scaling from the current window/display resolution.  Call this
// once per frame before starting a new ImGui frame so HiDPI changes are picked
// up when the window moves between monitors.
void updateImGuiScale(GLFWwindow* window);

// User multiplier relative to automatic monitor scaling. Applied next frame,
// together with widget spacing and on-demand font scaling.
float interfaceScale();
void setInterfaceScale(float multiplier);

// Detach backends and destroy the ImGui context.
void shutdownImGui();

// Apply the dark (default) or light theme.  Safe to call at any time.
void applyDarkTheme();
void applyLightTheme();
