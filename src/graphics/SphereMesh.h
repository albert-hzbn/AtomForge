#pragma once
#include <GL/glew.h>

#include <vector>

class SphereMesh
{
public:

    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint ebo = 0;     // Element Buffer Object for indexed rendering
    // The geometry, kept for the Vulkan renderer (OpenGL buffers are made
    // only when OpenGL draws).
    std::vector<float> vertices;         // x, y, z per vertex
    std::vector<unsigned int> indices;   // triangles

    int stacks;
    int slices;

    int vertexCount;
    int indexCount;  // Number of indices (for glDrawElementsInstanced)

    SphereMesh(int stacks, int slices);
    ~SphereMesh();
};