#pragma once

#include "electronic/Volume.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace atomforge::electronic
{
// Clip a convex polygon to the displayed cell.
inline std::vector<glm::dvec3> clipToCell(const Grid& grid, std::vector<glm::dvec3> polygon)
{
    const auto inverse=glm::inverse(grid.cell);
    for (auto& p : polygon) p=inverse*(p-grid.origin);
    for (int axis=0;axis<3;++axis) for (int side=0;side<2;++side)
    {
        std::vector<glm::dvec3> clipped;
        if (polygon.empty()) return {};
        auto previous=polygon.back();
        auto distance=[&](glm::dvec3 p)
        {
            const double d=side ? 1-p[axis] : p[axis];
            return std::abs(d)<1e-12 ? 0.0 : d;
        };
        double before=distance(previous);
        for (auto current : polygon)
        {
            const double after=distance(current);
            if ((before>=0)!=(after>=0))
                clipped.push_back(previous+(current-previous)*(before/(before-after)));
            if (after>=0) clipped.push_back(current);
            previous=current; before=after;
        }
        polygon.clear();
        for (auto p : clipped)
            if (polygon.empty() || glm::length(p-polygon.back())>1e-10) polygon.push_back(p);
        if (polygon.size()>1 && glm::length(polygon.front()-polygon.back())<1e-10) polygon.pop_back();
    }
    for (auto& p : polygon) p=grid.origin+grid.cell*p;
    return polygon;
}

struct SlicePlane
{
    glm::dvec3 origin, u, v;
    std::vector<glm::dvec3> boundary;
};

// Rectangular preview samples are clipped when drawn. Extend finite grids at
// their boundary only for samples outside that polygon, never for analysis.
inline Grid sampleSlicePlane(const Grid& grid, const SlicePlane& plane, int nu, int nv)
{
    grid.validate();
    if (nu<2 || nv<2 || static_cast<double>(nu)*nv>10000000)
        throw std::invalid_argument("Invalid section preview dimensions");
    Grid out;
    out.shape={nu,nv,2}; out.origin=plane.origin;
    out.cell=glm::dmat3(plane.u,plane.v,glm::normalize(glm::cross(plane.u,plane.v)));
    out.unit=grid.unit; out.name="section";
    out.values.resize(static_cast<std::size_t>(nu)*nv*2);
    const auto inverse=glm::inverse(grid.cell);
    for (int y=0;y<nv;++y) for (int x=0;x<nu;++x)
    {
        auto point=out.position(x,y,0);
        if (!grid.periodic) point=grid.origin+grid.cell*glm::clamp(inverse*(point-grid.origin),glm::dvec3(0),glm::dvec3(1));
        out.values[out.index(x,y,0)]=out.values[out.index(x,y,1)]=grid.sample(point);
    }
    out.validate();
    return out;
}

inline bool clipSegmentToCell(const Grid& grid, glm::dvec3& a, glm::dvec3& b)
{
    const auto inverse=glm::inverse(grid.cell);
    const auto start=inverse*(a-grid.origin), direction=inverse*(b-a);
    double enter=0, leave=1;
    for (int axis=0;axis<3;++axis)
    {
        if (std::abs(direction[axis])<1e-14)
        {
            if (start[axis]<0 || start[axis]>1) return false;
        }
        else
        {
            const double x=-start[axis]/direction[axis], y=(1-start[axis])/direction[axis];
            enter=std::max(enter,std::min(x,y)); leave=std::min(leave,std::max(x,y));
        }
    }
    if (leave<=enter) return false;
    const auto delta=b-a;
    b=a+delta*leave; a+=delta*enter;
    return true;
}

// The point is relative to the cell origin in fractional coordinates; the
// normal is Cartesian. The rectangle spans the cell's projection onto the plane.
inline SlicePlane slicePlane(const Grid& grid, glm::dvec3 point, glm::dvec3 normal)
{
    for (int i=0;i<3;++i)
        if (!std::isfinite(point[i]) || !std::isfinite(normal[i]))
            throw std::invalid_argument("Plane point and normal must be finite.");
    if (glm::length(normal)<1e-12) throw std::invalid_argument("Plane normal must be nonzero.");
    normal=glm::normalize(normal);
    const auto ex=glm::normalize(glm::cross(normal,std::abs(normal.z)<.9 ? glm::dvec3(0,0,1) : glm::dvec3(0,1,0)));
    const auto ey=glm::cross(normal,ex);
    const auto center=grid.origin+grid.cell*point;
    glm::dvec2 low(1e300), high(-1e300);
    for (int z=0;z<2;++z) for (int y=0;y<2;++y) for (int x=0;x<2;++x)
    {
        const auto delta=grid.origin+grid.cell*glm::dvec3(x,y,z)-center;
        const glm::dvec2 p(glm::dot(delta,ex),glm::dot(delta,ey));
        low=glm::min(low,p); high=glm::max(high,p);
    }
    SlicePlane plane{center+ex*low.x+ey*low.y,ex*(high.x-low.x),ey*(high.y-low.y),{}};
    plane.boundary=clipToCell(grid,{plane.origin,plane.origin+plane.u,plane.origin+plane.u+plane.v,plane.origin+plane.v});
    glm::dvec3 area(0);
    for (std::size_t i=1;i+1<plane.boundary.size();++i)
        area+=glm::cross(plane.boundary[i]-plane.boundary[0],plane.boundary[i+1]-plane.boundary[0]);
    if (glm::length(area)<1e-12*glm::length(plane.u)*glm::length(plane.v)) plane.boundary.clear();
    return plane;
}
}
