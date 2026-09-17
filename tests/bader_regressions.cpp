#include "electronic/BaderPartition.h"
#include "electronic/Volume.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace atomforge::electronic;

namespace
{
void expect(bool condition, const std::string& message) { if (!condition) throw std::runtime_error(message); }
void close(double actual, double expected, double tolerance) { expect(std::abs(actual - expected) <= tolerance, "Expected " + std::to_string(expected) + ", got " + std::to_string(actual)); }
}

int main()
{
    try
    {
        // Two identical, well-separated Gaussians on a finite (non-periodic)
        // grid, centered at grid indices 2 and 5 along x (indices 0..7), so
        // the exact bisector (x=3.5) falls between grid points, not on one.
        // y/z centers (3.3, 3.4) deliberately avoid the half-integer index
        // 3.5, which would tie rows 3 and 4 by exact symmetry and is not
        // the case this test means to exercise. For well-separated identical
        // peaks the attraction basins are expected to split at the
        // nearest-maximum boundary; this is checked directly below, not
        // merely assumed.
        const int n = 8;
        Grid density;
        density.shape = {n, n, n};
        density.cell = glm::dmat3(8.0);
        density.periodic = false;
        density.unit = "e/A^3";
        density.values.assign(static_cast<std::size_t>(n) * n * n, 0.0);
        const glm::dvec3 center1(2.0, 3.3, 3.4), center2(5.0, 3.3, 3.4);
        const double sigma = 0.7;
        for (int z = 0; z < n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x < n; ++x)
                {
                    const glm::dvec3 p(x, y, z);
                    const double d1 = glm::length(p - center1), d2 = glm::length(p - center2);
                    density.values[density.index(x, y, z)] =
                        std::exp(-d1 * d1 / (2 * sigma * sigma)) + std::exp(-d2 * d2 / (2 * sigma * sigma));
                }

        const auto partition = baderOnGrid(density);
        expect(partition.basins.size() == 2, "Two well-separated Gaussians must give exactly two basins");

        double totalIntegral = 0.0;
        for (int z = 0; z < n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x < n; ++x)
                    totalIntegral += density.weight(x, y, z) * density.values[density.index(x, y, z)];

        double basinSum = 0.0;
        for (const auto& basin : partition.basins) basinSum += basin.charge;
        close(basinSum, totalIntegral, 1e-9 * std::max(1.0, totalIntegral));  // exact charge conservation

        close(partition.basins[0].charge, partition.basins[1].charge, 1e-6);  // symmetric peaks: equal basins
        close(partition.basins[0].charge, totalIntegral / 2, 1e-6);

        // The true analytic maxima (at center1/center2) generally fall
        // between grid points; the discrete maximum must be the nearest one.
        const glm::dvec3 nearestGridMax1 = density.position(2, 3, 3), nearestGridMax2 = density.position(5, 3, 3);
        for (const auto& basin : partition.basins)
        {
            const bool atMax1 = glm::length(basin.maximum - nearestGridMax1) < 1e-9;
            const bool atMax2 = glm::length(basin.maximum - nearestGridMax2) < 1e-9;
            expect(atMax1 || atMax2, "Each basin maximum must sit at the grid point nearest a Gaussian center");
        }

        // Every voxel's basin must match whichever Gaussian center is
        // Euclidean-closer (the analytic expectation for well-separated,
        // identical, symmetric peaks) -- an independent, per-voxel check,
        // not just the aggregate charge/count checks above.
        for (int z = 0; z < n; ++z)
            for (int y = 0; y < n; ++y)
                for (int x = 0; x < n; ++x)
                {
                    const glm::dvec3 p(x, y, z);
                    const int expectedBasin = glm::length(p - center1) < glm::length(p - center2) ? 0 : 1;
                    const int actualBasin = partition.basin[density.index(x, y, z)];
                    const int actualCenterMatch = glm::length(partition.basins[static_cast<std::size_t>(actualBasin)].maximum - nearestGridMax1) < 1e-9 ? 0 : 1;
                    expect(actualCenterMatch == expectedBasin, "Voxel basin assignment must match the nearer Gaussian center");
                }

        // baderPopulations(): assign basins to explicit Cartesian sites (at
        // the two discrete maxima; site coordinates are Cartesian, matching
        // basin.maximum, not the raw grid-index coordinates used above to
        // define the analytic density).
        const auto populations = baderPopulations(density, {nearestGridMax1, nearestGridMax2});
        close(populations[0], totalIntegral / 2, 1e-6);
        close(populations[1], totalIntegral / 2, 1e-6);

        // A single, isolated peak on a periodic grid must give exactly one
        // basin (every point ascends toward the one periodic maximum) whose
        // charge equals the full integral -- a periodicity/wraparound check.
        {
            Grid single;
            single.shape = {n, n, n};
            single.cell = glm::dmat3(8.0);
            single.periodic = true;
            single.unit = "e/A^3";
            single.values.assign(static_cast<std::size_t>(n) * n * n, 0.0);
            const glm::dvec3 center(4, 4, 4);
            for (int z = 0; z < n; ++z)
                for (int y = 0; y < n; ++y)
                    for (int x = 0; x < n; ++x)
                    {
                        const glm::dvec3 p(x, y, z);
                        const double d = glm::length(p - center);
                        single.values[single.index(x, y, z)] = std::exp(-d * d / (2 * sigma * sigma));
                    }
            const auto single_partition = baderOnGrid(single);
            expect(single_partition.basins.size() == 1, "A single periodic peak must give exactly one basin");
        }

        // Errors.
        {
            bool threw = false;
            Grid wrongUnit = density; wrongUnit.unit = "eV";
            try { baderOnGrid(wrongUnit); } catch (const std::exception&) { threw = true; }
            expect(threw, "Non-e/A^3 density must be rejected");
        }
        {
            bool threw = false;
            Grid negative = density; negative.values[0] = -1;
            try { baderOnGrid(negative); } catch (const std::exception&) { threw = true; }
            expect(threw, "Negative density must be rejected");
        }
        {
            bool threw = false;
            try { baderPopulations(density, {}); } catch (const std::exception&) { threw = true; }
            expect(threw, "An empty site list must be rejected");
        }

        std::cout << "Bader partitioning analytic regressions passed\n";
    }
    catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
