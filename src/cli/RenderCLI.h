#pragma once

// Headless PNG/JPG/SVG rendering of a structure file (--render mode).
// Requires a working OpenGL context (offscreen GLFW window); prints a clear
// error and returns non-zero if no display/GPU driver is available.
int runRenderCLI(int argc, char* argv[]);
