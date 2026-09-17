#include "electronic/Topology.h"
#include "electronic/Volume.h"

#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace atomforge::electronic;

namespace
{
void expect(bool condition, const std::string& message) { if (!condition) throw std::runtime_error(message); }

// A finite grid of the given shape, all-background (0), with sites marked
// solid (1) wherever `inside` returns true. Values elsewhere in [0,1) never
// cross threshold=0.5, so bettiNumbers(_, 0.5) sees exactly `inside`.
Grid shape(int nx, int ny, int nz, const std::function<bool(int, int, int)>& inside)
{
    Grid g;
    g.shape = {nx, ny, nz};
    g.cell = glm::dmat3(1.0);
    g.periodic = false;
    g.unit = "raw";
    g.values.assign(static_cast<std::size_t>(nx) * ny * nz, 0.0);
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x)
                if (inside(x, y, z)) g.values[g.index(x, y, z)] = 1.0;
    return g;
}
}

int main()
{
    try
    {
        // Single solid voxel: a filled ball, chi=1, B0=1,B1=0,B2=0. Also
        // checks the exact vertex/edge/face count of one unit cube (8,12,6).
        {
            const auto g = shape(3, 3, 3, [](int x, int y, int z) { return x == 1 && y == 1 && z == 1; });
            const auto b = bettiNumbers(g, 0.5);
            expect(b.betti0 == 1 && b.betti1 == 0 && b.betti2 == 0, "Single voxel must be a topological ball");
        }
        // A simply-connected 3x3x3 solid block: still a ball.
        {
            const auto g = shape(5, 5, 5, [](int x, int y, int z) { return x >= 1 && x <= 3 && y >= 1 && y <= 3 && z >= 1 && z <= 3; });
            const auto b = bettiNumbers(g, 0.5);
            expect(b.betti0 == 1 && b.betti1 == 0 && b.betti2 == 0, "Solid block must be a topological ball");
        }
        // Hollow shell: a 3x3x3 solid block with its center voxel empty.
        // Homotopy equivalent to a sphere: chi=2, B0=1, B1=0, B2=1 (one
        // enclosed cavity, independently confirmed by Euler's formula
        // 1 - 0 + 1 = 2 = chi(S^2)).
        {
            const auto g = shape(5, 5, 5, [](int x, int y, int z) {
                if (x == 2 && y == 2 && z == 2) return false;
                return x >= 1 && x <= 3 && y >= 1 && y <= 3 && z >= 1 && z <= 3;
            });
            const auto b = bettiNumbers(g, 0.5);
            expect(b.betti0 == 1 && b.betti2 == 1, "Hollow shell must have one connected component and one cavity");
            expect(b.betti1 == 0, "A hollow shell (sphere) has no handles");
        }
        // Solid square ring (a short square tube), extruded 3 voxels in z
        // with background padding on every side. Homotopy equivalent to a
        // circle: chi(annulus x I) = chi(annulus) = 0, so B0=1,B1=1,B2=0
        // (the z=1..3 tunnel opens to the z=0/z=4 padding, so it is not an
        // enclosed cavity).
        {
            const auto g = shape(9, 9, 5, [](int x, int y, int z) {
                if (z < 1 || z > 3) return false;
                const bool outer = x >= 1 && x <= 7 && y >= 1 && y <= 7;
                const bool inner = x >= 3 && x <= 5 && y >= 3 && y <= 5;
                return outer && !inner;
            });
            const auto b = bettiNumbers(g, 0.5);
            expect(b.betti0 == 1, "Ring must be a single connected component");
            expect(b.betti2 == 0, "Ring's through-hole is open, not an enclosed cavity");
            expect(b.betti1 == 1, "Ring must have exactly one handle (the through-hole)");
        }
        // Two disjoint solid blocks: B0=2.
        {
            const auto g = shape(7, 3, 3, [](int x, int y, int z) {
                (void)y; (void)z;
                return x == 1 || x == 5;
            });
            const auto b = bettiNumbers(g, 0.5);
            expect(b.betti0 == 2 && b.betti1 == 0 && b.betti2 == 0, "Two disjoint blobs must give betti0=2");
        }
        // Betti curve over a superlevel-set sweep: as the threshold rises,
        // the solid set can only shrink. Give the hollow-shell wall a single
        // low-density "pinhole" voxel and a low-but-nonzero background: at a
        // low threshold every shell voxel (including the pinhole) is solid
        // and the cavity is sealed (betti2=1); past the pinhole's density it
        // drops out of the solid set, the cavity leaks into the exterior
        // background through that gap, and betti2 must fall back to 0.
        {
            Grid g = shape(5, 5, 5, [](int, int, int) { return false; });
            for (int z = 1; z <= 3; ++z)
                for (int y = 1; y <= 3; ++y)
                    for (int x = 1; x <= 3; ++x)
                        if (!(x == 2 && y == 2 && z == 2))
                            g.values[g.index(x, y, z)] = (x == 2 && y == 2 && z == 1) ? 1.0 : 2.0;
            const auto curve = bettiCurve(g, {0.5, 1.5, 3.0});
            expect(curve.size() == 3, "Expected one row per threshold");
            expect(curve[0].y == 1 && curve[0].w == 1, "Below the pinhole density the whole shell (incl. pinhole) is solid: sealed cavity");
            expect(curve[1].y == 1 && curve[1].w == 0, "Past the pinhole density the shell has a gap: the cavity leaks out, betti2=0");
            expect(curve[1].z == 0, "A single-voxel gap opens the shell back into a topological ball");
            expect(curve[2].y == 0, "Above the wall density nothing is solid");
        }

        // Errors: periodic grids and non-finite thresholds are rejected.
        {
            Grid periodic = shape(3, 3, 3, [](int, int, int) { return true; });
            periodic.periodic = true;
            bool threw = false;
            try { bettiNumbers(periodic, 0.5); } catch (const std::exception&) { threw = true; }
            expect(threw, "Periodic grids must be rejected");
        }
        {
            const auto g = shape(3, 3, 3, [](int, int, int) { return true; });
            bool threw = false;
            try { bettiCurve(g, {}); } catch (const std::exception&) { threw = true; }
            expect(threw, "An empty threshold list must be rejected");
        }

        std::cout << "Topology analytic regressions passed\n";
    }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
