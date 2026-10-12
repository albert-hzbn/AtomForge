#version 450
layout(push_constant) uniform Push { vec4 color; vec4 extra; } push;

layout(location = 0) out vec4 color;
void main()
{
    color = vec4(push.color.rgb, 1.0);
}
