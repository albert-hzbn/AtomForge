#include "electronic/Volume.h"
#include "electronic/DisplayRange.h"
#include "electronic/SlicePlane.h"
#include "electronic/ChargeAnalysis.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace atomforge::electronic;

namespace
{
void close(double actual, double expected, double tolerance = 1e-9)
{
    if (!std::isfinite(actual) || std::abs(actual-expected) > tolerance)
        throw std::runtime_error("Expected " + std::to_string(expected) + ", got " + std::to_string(actual));
}
template<class Work> void rejects(Work work)
{
    try { work(); } catch (const std::exception&) { return; }
    throw std::runtime_error("Invalid input was accepted");
}
Grid constant(bool periodic, int n = 8)
{
    Grid g;
    g.shape = {n,n,n};
    g.periodic = periodic;
    g.unit = "e/A^3";
    g.cell = glm::dmat3({2,0,0},{0.7,3,0},{0.2,0.4,4});
    g.values.assign(static_cast<std::size_t>(n)*n*n,2);
    return g;
}
}

int main()
{
    try
    {
        Grid cube;
        for (bool periodic : {false,true})
        {
            auto density=constant(periodic);
            for (int z=0;z<8;++z) for (int y=0;y<8;++y) for (int x=0;x<8;++x)
                density.values[density.index(x,y,z)]=x<4 ? -2 : 2;
            const auto positive=thresholdMask(density,0,2), negative=thresholdMask(density,-2,0);
            const auto both=booleanMask(positive,negative,"union");
            close(integrate(both),24);
            close(integrate(booleanMask(positive,negative,"intersection")),0);
            close(integrate(booleanMask(positive,negative,"xor")),24);
            close(integrate(booleanMask(both,positive,"difference")),12);
            close(integrate(applyMask(density,positive)),24);
            const auto summary=chargeSummary(density);
            close(summary[0],24); close(summary[1],24); close(summary[2],0);
            close(integrate(densityDifference(constant(periodic),constant(periodic),.5)),24);
            for (int axis=0;axis<3;++axis)
            {
                const auto profile=cumulativeCharge(density,axis);
                close(profile.front().y,0); close(profile.back().y,integrate(density));
                close(cumulativeCharge(constant(periodic),axis).back().y,48);
            }
            rejects([&] { booleanMask(density,positive,"union"); });
            rejects([&] { thresholdMask(density,2,-2); });
            rejects([&] { cumulativeCharge(density,3); });
            auto shifted=positive; shifted.origin.x+=.1;
            rejects([&] { applyMask(density,shifted); });
            auto potential=density; potential.unit="eV";
            rejects([&] { chargeSummary(potential); });
        }
        cube.values.assign(8,1);
        const auto diagonal=slicePlane(cube,{.5,.5,.5},{1,1,1});
        if (diagonal.boundary.size()!=6) throw std::runtime_error("Diagonal cube intersection must be a hexagon");
        for (auto p : diagonal.boundary) close(p.x+p.y+p.z,1.5);
        if (!slicePlane(cube,{2,2,2},{1,1,1}).boundary.empty()) throw std::runtime_error("Outside plane intersects cell");
        if (!slicePlane(cube,{0,0,0},{1,1,1}).boundary.empty()) throw std::runtime_error("Point contact is not a section");
        if (slicePlane(cube,{0,0,0},{0,0,1}).boundary.size()!=4) throw std::runtime_error("Cell face intersection was lost");
        rejects([&] { slicePlane(cube,{.5,.5,.5},{0,0,0}); });
        auto skew=constant(false); skew.origin={3,-2,7};
        const glm::dvec3 normal(2,-1,3), center=skew.origin+skew.cell*glm::dvec3(.5);
        const auto oblique=slicePlane(skew,{.5,.5,.5},normal);
        const auto preview=sampleSlicePlane(skew,oblique,33,33);
        for (double value : preview.values) close(value,2);
        if (oblique.boundary.size()<3) throw std::runtime_error("Skewed cell intersection is empty");
        for (auto p : oblique.boundary)
        {
            close(glm::dot(p-center,normal),0);
            const auto f=glm::inverse(skew.cell)*(p-skew.origin);
            for (int axis=0;axis<3;++axis) if (f[axis]<-1e-9 || f[axis]>1+1e-9) throw std::runtime_error("Section extends outside skewed cell");
        }
        glm::dvec3 start(-1,.5,.5), end(2,.5,.5);
        if (!clipSegmentToCell(cube,start,end)) throw std::runtime_error("Crossing contour was discarded");
        close(start.x,0); close(end.x,1);
        start={-1,2,.5}; end={2,2,.5};
        if (clipSegmentToCell(cube,start,end)) throw std::runtime_error("Outside contour was retained");
        const auto signedRange=displayRange({-3,-1,1,4});
        close(signedRange.suggested,0);
        if (signedRange.contains(-4) || signedRange.contains(5) || !signedRange.contains(-3) || !signedRange.contains(4))
            throw std::runtime_error("Display range validation failed");
        close(displayRange({2,2,2}).suggested,2);
        std::vector<double> skewed(1000,.02); skewed[0]=0; skewed.back()=1000;
        close(displayRange(skewed).suggested,.02);
        const auto sparse=displayRange({0,0,0,0,0,0,0,0,0,1});
        if (!(sparse.suggested>0 && sparse.suggested<1)) throw std::runtime_error("Sparse density estimate is outside usable range");
        for (bool periodic : {false,true})
        {
            const auto g = constant(periodic);
            close(integrate(g),48);
            close(g.sample(g.origin+g.cell*glm::dvec3(0.37,0.62,0.24)),2);
            for (const auto& d : gradient(g)) for (double v : d.values) close(v,0);
            for (double v : laplacian(g).values) close(v,0);
            const auto smoothed = smooth(g,0.3,1);
            close(integrate(smoothed),48);
            const auto basins = voronoiIntegrate(g,{{0,0,0},{1,0,0}});
            close(basins[0].integral+basins[1].integral,48);
            close(basins[0].volume+basins[1].volume,24);
            for (auto row : planarAverage(g,1)) close(row.y,2);
            const auto energy = energyDensity(g);
            const double expected = 0.3*std::pow(3*pi*pi,2.0/3)*std::pow(2*std::pow(bohr,3),5.0/3)*hartree/std::pow(bohr,3);
            close(energy[0].values[0],expected);
            close(energy[1].values[0],-2*expected);
            close(energy[2].values[0],-expected);
        }
        auto polynomial = constant(false);
        for (int z = 0; z < 8; ++z) for (int y = 0; y < 8; ++y) for (int x = 0; x < 8; ++x)
        {
            const auto p = polynomial.position(x,y,z);
            polynomial.values[polynomial.index(x,y,z)] = p.x*p.x+2*p.y*p.y+3*p.z*p.z+p.x*p.y;
        }
        const auto grad = gradient(polynomial);
        const auto lap = laplacian(polynomial);
        for (int z = 0; z < 8; ++z) for (int y = 0; y < 8; ++y) for (int x = 0; x < 8; ++x)
        {
            const auto p = polynomial.position(x,y,z);
            const auto i = polynomial.index(x,y,z);
            close(grad[0].values[i],2*p.x+p.y,1e-10);
            close(grad[1].values[i],4*p.y+p.x,1e-10);
            close(grad[2].values[i],6*p.z,1e-10);
            close(lap.values[i],12,1e-9);
        }
        for (int n : {7,8})
        {
            auto g = constant(true,n);
            g.origin = {0.1,0.2,0.3};
            for (int z = 0; z < n; ++z) for (int y = 0; y < n; ++y) for (int x = 0; x < n; ++x)
                g.values[g.index(x,y,z)] = 2+std::cos(2*pi*x/n);
            std::vector<glm::ivec3> hkl;
            for (int z = -n/2; z <= (n-1)/2; ++z) for (int y = -n/2; y <= (n-1)/2; ++y) for (int x = -n/2; x <= (n-1)/2; ++x) hkl.push_back({x,y,z});
            const auto factors = structureFactors(g,hkl);
            const auto restored = fourierSynthesis(g,factors);
            for (std::size_t i = 0; i < g.values.size(); ++i) close(restored.values[i],g.values[i],1e-10);
            const auto p = patterson(g);
            close(p.values[0],24*4.5);
            close(integrate(p),48*48);
            rejects([&] { structureFactors(g,{{n,0,0}}); });
        }
        auto plane = constant(false);
        for (int z = 0; z < 8; ++z) for (int y = 0; y < 8; ++y) for (int x = 0; x < 8; ++x)
            plane.values[plane.index(x,y,z)] = plane.position(x,y,z).x;
        const auto mesh = isosurface(plane,1.0);
        if (mesh.vertices.empty() || mesh.vertices.size()%3) throw std::runtime_error("Missing planar isosurface");
        for (auto p : mesh.vertices) close(p.x,1);
        const auto slice = section(plane,{0.5,0.2,0.2},{1,0,0},{0,1,0},8,8);
        const auto segments = contours(slice,1);
        if (segments.empty() || segments.size()%2) throw std::runtime_error("Missing contour");
        for (auto p : segments) close(p.x,1);
        // All three lattice-plane orientations, including boundary slices in
        // a skewed cell, must preserve the selected fractional coordinate.
        auto linear = constant(false);
        for (int z=0;z<8;++z) for (int y=0;y<8;++y) for (int x=0;x<8;++x)
            linear.values[linear.index(x,y,z)]=(x+2.0*y+3.0*z)/7;
        for (int axis=0;axis<3;++axis) for (double position : {0.0,.5,1.0})
        {
            const int a=axis==0 ? 1 : 0, b=axis==2 ? 1 : 2;
            const auto cut=section(linear,linear.origin+linear.cell[axis]*position,linear.cell[a],linear.cell[b],9,9);
            for (int y=0;y<9;++y) for (int x=0;x<9;++x)
            {
                glm::dvec3 fraction(0); fraction[axis]=position; fraction[a]=x/8.0; fraction[b]=y/8.0;
                close(cut.values[cut.index(x,y,0)],fraction.x+2*fraction.y+3*fraction.z);
            }
            for (double level : {2.0,3.0}) for (auto point : contours(cut,level))
            {
                const auto fraction=glm::inverse(linear.cell)*(point-linear.origin);
                close(fraction[axis],position);
                close(fraction.x+2*fraction.y+3*fraction.z,level);
            }
        }
        // Conventional NaCl cell, nearest-neighbor distance 1 Angstrom.
        const glm::dmat3 cell(2.0);
        const std::vector<glm::dvec3> sites{{0,0,0},{0,1,1},{1,0,1},{1,1,0},{1,0,0},{0,1,0},{0,0,1},{1,1,1}};
        const std::vector<double> charges{1,1,1,1,-1,-1,-1,-1};
        const auto e = ewald(cell,sites,charges,1.2,8,16);
        const double expected = -1.747564594633182*hartree*bohr;
        close(e.potentials[0],expected,1e-8);
        close(e.energy,4*expected,1e-8);
        const auto converged = ewald(cell,sites,charges,1.0,10,18);
        close(e.energy,converged.energy,1e-8);
        rejects([&] { ewald(cell,{{0,0,0}},{1},1,5,8); });
        auto invalid = constant(true);
        invalid.values[0] = -1;
        rejects([&] { energyDensity(invalid); });
        rejects([&] { patterson(polynomial); });
        rejects([&] { arithmetic(polynomial,constant(true),"subtract"); });
        rejects([&] { polynomial.sample({-1,-1,-1}); });
        rejects([&] { planarAverage(polynomial,3); });
        std::cout << "Electronic analytic regressions passed\n";
    }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
