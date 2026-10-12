#version 450
// A chunk's box for its occlusion query.
layout(push_constant) uniform Box { mat4 mvp; } box;
layout(location = 0) in vec3 position;
void main()
{
    gl_Position = box.mvp * vec4(position, 1.0);
}
