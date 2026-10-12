#version 450
// A full-screen triangle; screen.y is +1 at the top of the image, as in OpenGL.
layout(location = 0) out vec2 screen;
void main()
{
    vec2 p = vec2((gl_VertexIndex<<1)&2, gl_VertexIndex&2);
    screen = p*2.0-1.0;
    gl_Position = vec4(screen.x, -screen.y, 0.0, 1.0);
}
