#pragma once
// The last frame shown in the main window, read back for screenshots (the UI
// tour, tutorial captures and the step-dialog check), with either renderer.

#include <vector>

namespace graphics
{
// Keeps a copy of every presented frame so lastFrame() can read it (Vulkan;
// OpenGL reads its front buffer). Call once the window exists.
void enableFrameCapture();

// Vertical sync of the main window (glfwSwapInterval with OpenGL).
void setVsync(bool on);

// Waits until the GPU has finished all submitted work (glFinish).
void waitIdle();

// RGB pixels, top row first; false when there is no frame to read.
bool lastFrame(std::vector<unsigned char>& rgb, int& width, int& height);
}
