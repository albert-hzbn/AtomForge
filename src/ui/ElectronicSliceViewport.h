#pragma once

#include "electronic/Volume.h"

// Independent planar camera and contour level for the electronic 2D view.
class ElectronicSliceViewport
{
public:
    void draw(const atomforge::electronic::Grid& grid, float low, float high, int palette);
    void invalidate() { m_dirty = true; }
    void reset() { m_zoom = 1; m_angle = 0; m_pan = glm::vec2(0); }

private:
    int m_axis = 2;
    float m_position = .5f;
    float m_level = .1f;
    float m_zoom = 1;
    float m_angle = 0;
    glm::vec2 m_pan{0};
    bool m_dirty = true;
    atomforge::electronic::Grid m_slice;
    std::vector<glm::dvec3> m_contours;
};
