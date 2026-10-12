#version 450
// The disc of a shadow point sprite.
void main()
{
    vec2 c = gl_PointCoord * 2.0 - 1.0;
    if (dot(c, c) > 1.0) discard;
}
