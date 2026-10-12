#include "CylinderMesh.h"
#include "graphics/GraphicsBackend.h"

#include <cmath>
#include <vector>

namespace
{
void appendVertex(std::vector<float>& vertices, float x, float y, float z)
{
    vertices.push_back(x);
    vertices.push_back(y);
    vertices.push_back(z);
}

void createCylinder(std::vector<float>& vertices, int slices)
{
    const float halfLength = 0.5f;
    const float twoPi = 6.28318530717958647692f;

    for (int index = 0; index < slices; ++index)
    {
        float angle0 = twoPi * (float)index / (float)slices;
        float angle1 = twoPi * (float)(index + 1) / (float)slices;

        float x0 = std::cos(angle0);
        float y0 = std::sin(angle0);
        float x1 = std::cos(angle1);
        float y1 = std::sin(angle1);

        appendVertex(vertices, x0, y0, -halfLength);
        appendVertex(vertices, x0, y0, halfLength);
        appendVertex(vertices, x1, y1, halfLength);

        appendVertex(vertices, x0, y0, -halfLength);
        appendVertex(vertices, x1, y1, halfLength);
        appendVertex(vertices, x1, y1, -halfLength);

        appendVertex(vertices, 0.0f, 0.0f, -halfLength);
        appendVertex(vertices, x1, y1, -halfLength);
        appendVertex(vertices, x0, y0, -halfLength);

        appendVertex(vertices, 0.0f, 0.0f, halfLength);
        appendVertex(vertices, x0, y0, halfLength);
        appendVertex(vertices, x1, y1, halfLength);
    }
}
}

CylinderMesh::CylinderMesh(int sliceCount)
    : slices(sliceCount)
{
    createCylinder(vertices, slices);
    vertexCount = (int)(vertices.size() / 3);
    if (!graphics::usingOpenGL()) return;

    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 vertices.size() * sizeof(float),
                 vertices.data(),
                 GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}