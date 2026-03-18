#include "SphereMesh.h"
#include <cmath>
#include <vector>

static void createSphere(std::vector<float>& v,int stacks,int slices)
{
    constexpr float kPi = 3.14159265358979323846f;

    for(int i=0;i<stacks;i++)
    {
        float lat0 = kPi*(-0.5f + i/(float)stacks);
        float lat1 = kPi*(-0.5f + (i+1)/(float)stacks);

        for(int j=0;j<slices;j++)
        {
            float lng0 = 2*kPi*j/(float)slices;
            float lng1 = 2*kPi*(j+1)/(float)slices;

            float x0 = cos(lng0);
            float y0 = sin(lng0);

            float x1 = cos(lng1);
            float y1 = sin(lng1);

            float z0 = sin(lat0);
            float zr0 = cos(lat0);

            float z1 = sin(lat1);
            float zr1 = cos(lat1);

            /* triangle 1 */

            v.push_back(x0*zr0);
            v.push_back(y0*zr0);
            v.push_back(z0);

            v.push_back(x0*zr1);
            v.push_back(y0*zr1);
            v.push_back(z1);

            v.push_back(x1*zr1);
            v.push_back(y1*zr1);
            v.push_back(z1);

            /* triangle 2 */

            v.push_back(x0*zr0);
            v.push_back(y0*zr0);
            v.push_back(z0);

            v.push_back(x1*zr1);
            v.push_back(y1*zr1);
            v.push_back(z1);

            v.push_back(x1*zr0);
            v.push_back(y1*zr0);
            v.push_back(z0);
        }
    }
}

SphereMesh::SphereMesh(int st,int sl)
{
    stacks = st;
    slices = sl;

    std::vector<float> vertices;

    createSphere(vertices, stacks, slices);

    vertexCount = vertices.size() / 3;

    glGenVertexArrays(1,&vao);
    glGenBuffers(1,&vbo);

    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER,vbo);

    glBufferData(GL_ARRAY_BUFFER,
                 vertices.size()*sizeof(float),
                 vertices.data(),
                 GL_STATIC_DRAW);

    glVertexAttribPointer(0,
                          3,
                          GL_FLOAT,
                          GL_FALSE,
                          3*sizeof(float),
                          (void*)0);

    glEnableVertexAttribArray(0);

    glBindVertexArray(0);
}
