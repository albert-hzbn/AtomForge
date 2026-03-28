#pragma once
#include <GL/glew.h>

class SphereMesh
{
public:

    GLuint vao;
    GLuint vbo;
    GLuint ebo;     // Element Buffer Object for indexed rendering

    int stacks;
    int slices;

    int vertexCount;
    int indexCount;  // Number of indices (for glDrawElementsInstanced)

    SphereMesh(int stacks, int slices);
    ~SphereMesh();
};