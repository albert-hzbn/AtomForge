#pragma once
#include <GL/glew.h>

// Subdivided icosahedron mesh for mid-range instancing (100k–10M atoms).
// Two subdivision levels: 162 vertices, 320 faces (960 indices).
// Much cheaper than a full 40×40 UV sphere but smooth enough to look round.
class LowPolyMesh
{
public:
    GLuint vao;
    GLuint vbo;
    GLuint ebo;

    int vertexCount;
    int indexCount;

    LowPolyMesh();
    ~LowPolyMesh();
};
