#pragma once
#include <GL/glew.h>

#include <vector>

// Subdivided icosahedron mesh for mid-range instancing (100k–10M atoms).
// Two subdivision levels: 162 vertices, 320 faces (960 indices).
// Much cheaper than a full 40×40 UV sphere but smooth enough to look round.
class LowPolyMesh
{
public:
    GLuint vao = 0;
    GLuint vbo = 0;
    GLuint ebo = 0;
    // The geometry, kept for the Vulkan renderer (OpenGL buffers are made
    // only when OpenGL draws).
    std::vector<float> vertices;         // x, y, z per vertex
    std::vector<unsigned int> indices;   // triangles

    int vertexCount;
    int indexCount;

    LowPolyMesh();
    ~LowPolyMesh();
};
