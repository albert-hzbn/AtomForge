#pragma once

#include <GL/glew.h>

#include <vector>

class CylinderMesh
{
public:
    GLuint vao = 0;
    GLuint vbo = 0;
    int slices = 0;
    int vertexCount = 0;
    // The triangles (x, y, z per vertex), kept for the Vulkan renderer.
    std::vector<float> vertices;

    explicit CylinderMesh(int slices);
};