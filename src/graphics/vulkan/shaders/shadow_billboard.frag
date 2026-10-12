#version 450
layout(location = 0) in vec2 quadUv;
void main()
{
    if (dot(quadUv, quadUv) > 1.0) discard;
}
