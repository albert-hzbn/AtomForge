#pragma once
#include <GL/glew.h>

class SphereMesh
{
public:

    GLuint vao;
    GLuint vbo;

    int stacks;
    int slices;

    int vertexCount;

    SphereMesh(int stacks,int slices);

    void draw();
};