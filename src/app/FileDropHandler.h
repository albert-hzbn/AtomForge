#pragma once

#include <string>
#include <vector>

struct GLFWwindow;
struct EditorState;

// Installs the GLFW drop callback that pushes dropped file paths into
// pendingDrops. The caller is responsible for routing them to the active tab.
void installDropFileCallback(GLFWwindow* window, std::vector<std::string>& pendingDrops);
void processDroppedFiles(EditorState& state);
