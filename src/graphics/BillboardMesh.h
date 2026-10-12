#pragma once
#include <GL/glew.h>

#include <vector>

// Billboard quad mesh for imposters (N > 10M atoms)
class BillboardMesh
{
public:
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint ebo = 0;
    // The geometry, kept for the Vulkan renderer (OpenGL buffers are made
    // only when OpenGL draws).
    std::vector<float> vertices;         // x, y, z per vertex
    std::vector<unsigned int> indices;   // triangles

    int indexCount;  // 6 indices for 2 triangles (1 quad)

    BillboardMesh();
    ~BillboardMesh();
};
