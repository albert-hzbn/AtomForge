#version 450
// The depth pre-pass of a point-sprite atom: its disc, no colour.
void main()
{
    vec2 c = gl_PointCoord * 2.0 - 1.0;
    if (dot(c, c) > 1.0) discard;
}
