#include "algorithms/SurfaceBuilder.h"
#include "algorithms/CSLComputation.h"

#include <algorithm>
#include <array>
#include <numeric>
#include <cmath>
#include <utility>
#include <vector>

namespace atomforge
{

namespace
{
struct IntVec3 { int n[3]; };

long long detM3(const int m[3][3])
{
    return (long long)m[0][0] * (m[1][1] * (long long)m[2][2] - m[1][2] * (long long)m[2][1])
         - (long long)m[0][1] * (m[1][0] * (long long)m[2][2] - m[1][2] * (long long)m[2][0])
         + (long long)m[0][2] * (m[1][0] * (long long)m[2][1] - m[1][1] * (long long)m[2][0]);
}

// Cartesian vector for an integer combination of the grain's lattice vectors.
void toCartesian(const IntVec3& v, const double cell[3][3], double out[3])
{
    for (int c = 0; c < 3; ++c)
        out[c] = v.n[0] * cell[0][c] + v.n[1] * cell[1][c] + v.n[2] * cell[2][c];
}

double lengthOf(const IntVec3& v, const double cell[3][3])
{
    double c[3];
    toCartesian(v, cell, c);
    return std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
}

// Maps Miller indices from the original cell to a reduced cell sharing its
// orientation (M = reduced * original^-1 must be rational). Returns false
// when M is not rational with small denominators (e.g. a rotated cell).
bool transformMillerIndices(const std::array<std::array<double, 3>, 3>& original,
                            const std::array<std::array<double, 3>, 3>& reduced, int miller[3])
{
    double c[3][3], inv[3][3];
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) c[i][j] = original[i][j];
    const double det = c[0][0] * (c[1][1] * c[2][2] - c[1][2] * c[2][1])
                     - c[0][1] * (c[1][0] * c[2][2] - c[1][2] * c[2][0])
                     + c[0][2] * (c[1][0] * c[2][1] - c[1][1] * c[2][0]);
    if (std::abs(det) < 1e-12) return false;
    inv[0][0] = (c[1][1] * c[2][2] - c[1][2] * c[2][1]) / det;
    inv[0][1] = (c[0][2] * c[2][1] - c[0][1] * c[2][2]) / det;
    inv[0][2] = (c[0][1] * c[1][2] - c[0][2] * c[1][1]) / det;
    inv[1][0] = (c[1][2] * c[2][0] - c[1][0] * c[2][2]) / det;
    inv[1][1] = (c[0][0] * c[2][2] - c[0][2] * c[2][0]) / det;
    inv[1][2] = (c[0][2] * c[1][0] - c[0][0] * c[1][2]) / det;
    inv[2][0] = (c[1][0] * c[2][1] - c[1][1] * c[2][0]) / det;
    inv[2][1] = (c[0][1] * c[2][0] - c[0][0] * c[2][1]) / det;
    inv[2][2] = (c[0][0] * c[1][1] - c[0][1] * c[1][0]) / det;
    double transformed[3] = {0, 0, 0};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            double m = 0;
            for (int k = 0; k < 3; ++k) m += reduced[i][k] * inv[k][j];
            // Reduced-cell coefficients are multiples of 1/2, 1/3, 1/4 or 1/6.
            if (std::abs(m * 12.0 - std::round(m * 12.0)) > 1e-5) return false;
            transformed[i] += m * miller[j];
        }
    long long scaled[3];
    for (int i = 0; i < 3; ++i) scaled[i] = std::llround(transformed[i] * 12.0);
    long long g = 0;
    for (long long v : scaled) g = std::gcd(g, std::llabs(v));
    if (g == 0) return false;
    for (int i = 0; i < 3; ++i) miller[i] = static_cast<int>(scaled[i] / g);
    return true;
}

bool isMultiple(const IntVec3& a, const IntVec3& b)
{
    // True when a and b are parallel (one is an integer/rational multiple
    // of the other), tested via the 2x2 minor cross-products.
    return a.n[0] * b.n[1] - a.n[1] * b.n[0] == 0
        && a.n[0] * b.n[2] - a.n[2] * b.n[0] == 0
        && a.n[1] * b.n[2] - a.n[2] * b.n[1] == 0;
}
}

SurfaceResult buildSurface(const Structure& source, const SurfaceParams& params)
{
    SurfaceResult result;
    if (params.h == 0 && params.k == 0 && params.l == 0)
    {
        result.message = "Miller indices (h,k,l) cannot all be zero.";
        return result;
    }
    if (params.layers < 1)
    {
        result.message = "layers must be at least 1.";
        return result;
    }
    if (params.nmax < 1)
    {
        result.message = "nmax must be at least 1.";
        return result;
    }
    if (!source.hasUnitCell || source.atoms.empty())
    {
        result.message = "Structure needs a unit cell and at least one atom.";
        return result;
    }

    // Miller indices refer to the input cell. They transform like lattice
    // vectors, so a reduced basis P = M C needs (h' k' l') = M (h k l).
    int miller[3] = {params.h, params.k, params.l};
    Structure prepped = source;
    if (params.primitiveInput && reduceToPrimitive(prepped, params.primitiveSymprec))
    {
        if (!transformMillerIndices(source.cellVectors, prepped.cellVectors, miller))
        {
            // The reduction rotated the cell: build from the input cell instead.
            prepped = source;
            miller[0] = params.h; miller[1] = params.k; miller[2] = params.l;
        }
    }

    Grain grain = structureToGrain(prepped);

    // -- 1. Smallest independent pair of integer lattice vectors lying in
    // the (hkl) plane (h*n0 + k*n1 + l*n2 == 0), Gauss-reduced to a short,
    // close-to-orthogonal 2D basis.
    std::vector<IntVec3> inPlane;
    for (int n0 = -params.nmax; n0 <= params.nmax; ++n0)
    for (int n1 = -params.nmax; n1 <= params.nmax; ++n1)
    for (int n2 = -params.nmax; n2 <= params.nmax; ++n2)
    {
        if (n0 == 0 && n1 == 0 && n2 == 0) continue;
        if (miller[0] * n0 + miller[1] * n1 + miller[2] * n2 != 0) continue;
        inPlane.push_back({{n0, n1, n2}});
    }
    if (inPlane.size() < 2)
    {
        result.message = "Could not find two independent in-plane lattice vectors; raise nmax.";
        return result;
    }
    std::sort(inPlane.begin(), inPlane.end(), [&](const IntVec3& a, const IntVec3& b) {
        return lengthOf(a, grain.cell) < lengthOf(b, grain.cell);
    });
    IntVec3 v1 = inPlane.front();
    IntVec3 v2{};
    bool foundV2 = false;
    for (const auto& cand : inPlane)
    {
        if (isMultiple(cand, v1)) continue;
        v2 = cand;
        foundV2 = true;
        break;
    }
    if (!foundV2)
    {
        result.message = "Could not find two independent in-plane lattice vectors; raise nmax.";
        return result;
    }
    // 2D (Gauss/Lagrange) lattice reduction toward a short, near-orthogonal
    // in-plane basis, operating on the integer coefficients directly.
    for (int iter = 0; iter < 32; ++iter)
    {
        if (lengthOf(v2, grain.cell) < lengthOf(v1, grain.cell)) std::swap(v1, v2);
        double c1[3], c2[3];
        toCartesian(v1, grain.cell, c1);
        toCartesian(v2, grain.cell, c2);
        const double dot12 = c1[0] * c2[0] + c1[1] * c2[1] + c1[2] * c2[2];
        const double dot11 = c1[0] * c1[0] + c1[1] * c1[1] + c1[2] * c1[2];
        if (dot11 < 1e-12) break;
        const int m = (int)std::lround(dot12 / dot11);
        if (m == 0) break;
        for (int c = 0; c < 3; ++c) v2.n[c] -= m * v1.n[c];
    }

    // -- 2. Smallest out-of-plane step: minimize how many (hkl) planes it
    // skips (|dot| closest to the ideal single-layer step), then length.
    IntVec3 v3{};
    bool foundV3 = false;
    long long bestDot = 0;
    double bestLen = 0.0;
    for (int n0 = -params.nmax; n0 <= params.nmax; ++n0)
    for (int n1 = -params.nmax; n1 <= params.nmax; ++n1)
    for (int n2 = -params.nmax; n2 <= params.nmax; ++n2)
    {
        const long long dot = (long long)miller[0] * n0 + (long long)miller[1] * n1 + (long long)miller[2] * n2;
        if (dot == 0) continue;
        IntVec3 cand{{n0, n1, n2}};
        int m[3][3] = {{v1.n[0], v1.n[1], v1.n[2]}, {v2.n[0], v2.n[1], v2.n[2]}, {n0, n1, n2}};
        if (detM3(m) == 0) continue;
        const double len = lengthOf(cand, grain.cell);
        const long long absDot = dot < 0 ? -dot : dot;
        if (!foundV3 || absDot < bestDot || (absDot == bestDot && len < bestLen))
        {
            foundV3 = true;
            bestDot = absDot;
            bestLen = len;
            v3 = cand;
        }
    }
    if (!foundV3)
    {
        result.message = "Could not find an out-of-plane lattice vector completing the basis; raise nmax.";
        return result;
    }

    int m[3][3] = {{v1.n[0], v1.n[1], v1.n[2]}, {v2.n[0], v2.n[1], v2.n[2]}, {v3.n[0], v3.n[1], v3.n[2]}};
    Grain slabUnit = makeSupercell(grain, m);
    result.atomsPerLayer = (int)slabUnit.atoms.size();
    // Stack along the true (possibly tilted) lattice translation so each
    // layer keeps its in-plane registry (e.g. ABAB stacking), and only then
    // make the stacking direction normal to the surface.
    Grain slab = makeSupercellDiag(slabUnit, 1, 1, params.layers);
    if (!isOrthogonal(slab.cell))
        slab = setOrthogonalGrain(slab, 2);

    // -- 3. Vacuum, centered along the surface normal (cell direction 2).
    if (params.vacuum > 0.0)
    {
        const double len = cellVecLen(slab.cell[2]);
        if (len > 1e-9)
        {
            const double unit[3] = {slab.cell[2][0] / len, slab.cell[2][1] / len, slab.cell[2][2] / len};
            const double newLen = len + params.vacuum;
            for (int c = 0; c < 3; ++c) slab.cell[2][c] = unit[c] * newLen;
            const double shift[3] = {unit[0] * (params.vacuum * 0.5), unit[1] * (params.vacuum * 0.5), unit[2] * (params.vacuum * 0.5)};
            for (auto& atom : slab.atoms)
            {
                atom.x += shift[0];
                atom.y += shift[1];
                atom.z += shift[2];
            }
        }
    }

    result.structure = grainToStructure(slab);
    result.success = true;
    result.message = "Built a (" + std::to_string(params.h) + " " + std::to_string(params.k) + " "
        + std::to_string(params.l) + ") surface: " + std::to_string(params.layers) + " layers, "
        + std::to_string(result.structure.atoms.size()) + " atoms, " + std::to_string(params.vacuum) + " A vacuum.";
    return result;
}

}
