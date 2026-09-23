// Validates prepareDrag (src/algorithms/DragPrep.h), a native
// reimplementation of BABEL's documented prepareDrag program
// (mainPrepareDrag.f90): linear interpolation between two configurations
// plus a per-atom reaction-coordinate (constraint) direction, for use in an
// external constrained-minimization ("drag") migration-barrier calculation.
#include "algorithms/DragPrep.h"
#include "model/Structure.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace atomforge;

namespace
{
void expectClose(double a, double b, double tol, const char* what)
{
    if (std::abs(a - b) > tol)
        throw std::runtime_error(std::string(what) + ": expected " + std::to_string(b) + ", got " + std::to_string(a));
}
}

int main()
{
    try
    {
        Structure initial;
        initial.hasUnitCell = true;
        initial.cellVectors = {{{10.0, 0, 0}, {0, 10.0, 0}, {0, 0, 10.0}}};
        initial.atoms.push_back({"Cu", 29, 1.0, 2.0, 3.0});
        initial.atoms.push_back({"Cu", 29, 5.0, 5.0, 5.0});
        initial.atoms.push_back({"Cu", 29, 0.2, 5.0, 5.0}); // will migrate "the short way" across x=0 boundary

        Structure final_ = initial;
        final_.atoms[0].x += 2.0; final_.atoms[0].y += 1.0; final_.atoms[0].z -= 0.5;
        final_.atoms[1].x += 0.3;
        final_.atoms[2].x = 9.8; // 0.2 -> 9.8 the "long way" is +9.6, the short way (wrapped) is -0.4

        // --- 1. zeta=0 and zeta=1 must reproduce the endpoints exactly.
        {
            const DragPrepResult r0 = prepareDrag(initial, final_, 0.0, false, true);
            if (!r0.success) throw std::runtime_error("zeta=0 failed: " + r0.message);
            for (size_t i = 0; i < initial.atoms.size(); ++i)
            {
                expectClose(r0.interpolated.atoms[i].x, initial.atoms[i].x, 1e-12, "zeta=0 x");
                expectClose(r0.interpolated.atoms[i].y, initial.atoms[i].y, 1e-12, "zeta=0 y");
                expectClose(r0.interpolated.atoms[i].z, initial.atoms[i].z, 1e-12, "zeta=0 z");
            }

            const DragPrepResult r1 = prepareDrag(initial, final_, 1.0, false, true);
            if (!r1.success) throw std::runtime_error("zeta=1 failed: " + r1.message);
            // Without clipDisplacement, zeta=1 must reproduce `final_` exactly
            // (including atom 2's "long way round" motion).
            for (size_t i = 0; i < initial.atoms.size(); ++i)
            {
                expectClose(r1.interpolated.atoms[i].x, final_.atoms[i].x, 1e-9, "zeta=1 x");
                expectClose(r1.interpolated.atoms[i].y, final_.atoms[i].y, 1e-9, "zeta=1 y");
                expectClose(r1.interpolated.atoms[i].z, final_.atoms[i].z, 1e-9, "zeta=1 z");
            }
        }

        // --- 2. zeta=0.5 must be the exact linear midpoint.
        {
            const DragPrepResult r = prepareDrag(initial, final_, 0.5, false, true);
            if (!r.success) throw std::runtime_error("zeta=0.5 failed: " + r.message);
            expectClose(r.interpolated.atoms[0].x, 2.0, 1e-9, "midpoint atom0 x");
            expectClose(r.interpolated.atoms[0].y, 2.5, 1e-9, "midpoint atom0 y");
            expectClose(r.interpolated.atoms[0].z, 2.75, 1e-9, "midpoint atom0 z");
        }

        // --- 3. clipDisplacement must route atom 2 the SHORT way around the
        // periodic boundary (toward x=-0.2 wrapped to x=9.8, i.e. a net -0.4
        // step) rather than the long way (+9.6), so its interpolated
        // position at zeta=0.5 is close to x=0 (wrapped), not x=5.
        {
            const DragPrepResult r = prepareDrag(initial, final_, 0.5, true, true);
            if (!r.success) throw std::runtime_error("clipped zeta=0.5 failed: " + r.message);
            expectClose(r.constraintDirection[2].x, -0.4, 1e-9, "clipped constraint direction x");
            double wrapped = r.interpolated.atoms[2].x;
            while (wrapped < 0) wrapped += 10.0;
            while (wrapped >= 10.0) wrapped -= 10.0;
            expectClose(wrapped, 0.0, 1e-6, "clipped midpoint (wrapped)");
        }

        // --- 4. Bad inputs must be rejected.
        {
            const DragPrepResult bad = prepareDrag(initial, final_, 1.5, false, true);
            if (bad.success) throw std::runtime_error("zeta outside [0,1] should have been refused");
        }

        std::cout << "Drag prep regressions passed\n";
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
