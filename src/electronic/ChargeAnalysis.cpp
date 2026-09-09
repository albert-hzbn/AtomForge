#include "electronic/ChargeAnalysis.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace atomforge::electronic
{
namespace
{
void requireMask(const Grid& mask)
{
    mask.validate();
    if (mask.unit!="mask") throw std::invalid_argument("Create a threshold mask before using Boolean operations");
    for (double v : mask.values)
        if (v!=0 && v!=1) throw std::invalid_argument("Mask values must be zero or one");
}
void requireDensity(const Grid& grid)
{
    grid.validate();
    if (grid.unit!="e/A^3") throw std::invalid_argument("Charge analysis requires electron density in e/A^3");
}
}

Grid densityDifference(const Grid& a, const Grid& b, double weight)
{
    requireDensity(a); requireDensity(b);
    auto out=arithmetic(a,scale(b,weight),"subtract");
    out.name="density_difference";
    return out;
}

Grid thresholdMask(const Grid& grid, double low, double high)
{
    grid.validate();
    if (!std::isfinite(low) || !std::isfinite(high) || low>high)
        throw std::invalid_argument("Mask bounds must be finite and ordered");
    auto out=grid; out.unit="mask"; out.name="threshold_mask";
    for (auto& v : out.values) v=v>=low && v<=high ? 1 : 0;
    return out;
}

Grid booleanMask(const Grid& a, const Grid& b, const std::string& operation)
{
    requireAligned(a,b); requireMask(a); requireMask(b);
    if (operation!="union" && operation!="intersection" && operation!="difference" && operation!="xor")
        throw std::invalid_argument("Unknown Boolean mask operation");
    auto out=a; out.name="mask_"+operation;
    for (std::size_t i=0;i<out.values.size();++i)
    {
        const bool x=a.values[i]!=0, y=b.values[i]!=0;
        out.values[i]=operation=="union" ? x||y : operation=="intersection" ? x&&y : operation=="difference" ? x&&!y : x!=y;
    }
    return out;
}

Grid applyMask(const Grid& grid, const Grid& mask)
{
    requireAligned(grid,mask); requireMask(mask);
    auto out=grid; out.name="masked_"+grid.name;
    for (std::size_t i=0;i<out.values.size();++i) out.values[i]*=mask.values[i];
    return out;
}

std::array<Grid,2> splitDensity(const Grid& grid)
{
    requireDensity(grid);
    std::array<Grid,2> out{grid,grid};
    out[0].name="accumulation"; out[1].name="depletion_magnitude";
    for (std::size_t i=0;i<grid.values.size();++i)
    {
        out[0].values[i]=std::max(0.0,grid.values[i]);
        out[1].values[i]=std::max(0.0,-grid.values[i]);
    }
    return out;
}

std::array<double,3> chargeSummary(const Grid& grid)
{
    const auto parts=splitDensity(grid);
    return {integrate(parts[0]),integrate(parts[1]),integrate(grid)};
}

std::vector<glm::dvec2> cumulativeCharge(const Grid& grid, int axis)
{
    requireDensity(grid);
    const auto profile=planarAverage(grid,axis);
    const int a=(axis+1)%3, b=(axis+2)%3;
    const double area=glm::length(glm::cross(grid.cell[a],grid.cell[b]));
    const double height=std::abs(glm::determinant(grid.cell))/area;
    const double step=height/(grid.periodic ? grid.shape[axis] : grid.shape[axis]-1);
    std::vector<glm::dvec2> out{{0,0}};
    double sum=0;
    // Periodic bins use rectangle quadrature; finite samples use trapezoids,
    // matching integrate() and including the complete cell in the final row.
    if (grid.periodic)
        for (std::size_t i=0;i<profile.size();++i)
        { sum+=profile[i].y*area*step; out.emplace_back((i+1)*step,sum); }
    else
        for (std::size_t i=1;i<profile.size();++i)
        { sum+=(profile[i-1].y+profile[i].y)*.5*area*step; out.emplace_back(i*step,sum); }
    return out;
}
}
