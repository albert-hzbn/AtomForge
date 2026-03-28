#pragma once
#include <GL/glew.h>

// Billboard quad mesh for imposters (N > 10M atoms)
class BillboardMesh
{
public:
    GLuint vao;
    GLuint vbo;
    GLuint ebo;

    int indexCount;  // 6 indices for 2 triangles (1 quad)

    BillboardMesh();
    ~BillboardMesh();
};
