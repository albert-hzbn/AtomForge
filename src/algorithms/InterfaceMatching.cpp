#include "algorithms/InterfaceMatching.h"

#include "algorithms/InterfaceBuilder.h"
#include "algorithms/SurfaceBuilder.h"

#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <numeric>
#include <tuple>

namespace atomforge::interfaces
{
namespace
{
constexpr double kPi = 3.14159265358979323846;

double cross2(const glm::dvec2& a, const glm::dvec2& b) { return a.x * b.y - a.y * b.x; }

// Rotates every atom and cell vector of a structure by the rotation taking
// a to x and a x b to z (rows of the matrix: the new axes).
void orient(Structure& s)
{
    const glm::dvec3 a(s.cellVectors[0][0], s.cellVectors[0][1], s.cellVectors[0][2]);
    const glm::dvec3 b(s.cellVectors[1][0], s.cellVectors[1][1], s.cellVectors[1][2]);
    const glm::dvec3 x = glm::normalize(a);
    const glm::dvec3 z = glm::normalize(glm::cross(a, b));
    const glm::dvec3 y = glm::cross(z, x);
    const auto rotate = [&](const glm::dvec3& v) { return glm::dvec3(glm::dot(x, v), glm::dot(y, v), glm::dot(z, v)); };
    for (auto& atom : s.atoms)
    {
        const glm::dvec3 p = rotate(glm::dvec3(atom.x, atom.y, atom.z));
        atom.x = p.x;
        atom.y = p.y;
        atom.z = p.z;
    }
    for (auto& row : s.cellVectors)
    {
        const glm::dvec3 p = rotate(glm::dvec3(row[0], row[1], row[2]));
        row = {p.x, p.y, p.z};
    }
    const glm::dvec3 o = rotate(glm::dvec3(s.cellOffset[0], s.cellOffset[1], s.cellOffset[2]));
    s.cellOffset = {o.x, o.y, o.z};
    // The stacking vector points up (c and -c span the same lattice).
    if (s.cellVectors[2][2] < 0)
        for (double& v : s.cellVectors[2]) v = -v;
    // Exact zeros where the rotation leaves rounding noise.
    s.cellVectors[0][1] = s.cellVectors[0][2] = 0.0;
    s.cellVectors[1][2] = 0.0;
}

// The strain taking B's cell (columns v) onto A's (u): the stretch of the
// polar decomposition F = U V^-1 = S R (S symmetric) less the identity, in a
// frame with u1 along x; the rotation R that remains goes to rotationDeg.
// Unlike the plain component ratios it does not mistake a rotation for strain.
bool polarStrain(const glm::dvec2& u1, const glm::dvec2& u2, const glm::dvec2& v1, const glm::dvec2& v2,
                 double strain[3], double& rotationDeg)
{
    const glm::dmat2 U(u1, u2), V(v1, v2);
    if (std::abs(glm::determinant(V)) < 1e-12) return false;
    const glm::dmat2 F = U * glm::inverse(V);
    const glm::dmat2 B = F * glm::transpose(F);
    const double root = std::sqrt(std::max(0.0, glm::determinant(B)));
    const double norm = std::sqrt(B[0][0] + B[1][1] + 2.0 * root);
    if (norm < 1e-12) return false;
    const glm::dmat2 S = (B + root * glm::dmat2(1.0)) / norm;
    const glm::dmat2 R = glm::inverse(S) * F;
    rotationDeg = std::atan2(R[0][1], R[0][0]) * 180.0 / kPi;
    const double phi = std::atan2(u1.y, u1.x), c = std::cos(phi), s = std::sin(phi);
    const glm::dmat2 Q(glm::dvec2(c, -s), glm::dvec2(s, c));  // rotation by -phi
    const glm::dmat2 E = Q * (S - glm::dmat2(1.0)) * glm::transpose(Q);
    strain[0] = E[0][0];
    strain[1] = E[1][1];
    strain[2] = E[1][0];
    return true;
}

long long gcdl(long long a, long long b) { return std::gcd(std::llabs(a), std::llabs(b)); }

// The canonical basis {(p, 0), (q, s)}, 0 <= q < p, of the sublattice of Z^2
// spanned by two integer rows (with a positive determinant).
std::tuple<long long, long long, long long> canonical(const int m[2][2])
{
    const long long det = std::llabs(static_cast<long long>(m[0][0]) * m[1][1] - static_cast<long long>(m[0][1]) * m[1][0]);
    long long s = gcdl(m[0][1], m[1][1]);
    if (s == 0 || det == 0) return {0, 0, 0};
    const long long p = det / s;
    // x m0y + y m1y = s by the extended Euclidean algorithm.
    long long oldR = m[0][1], r = m[1][1], oldX = 1, xx = 0, oldY = 0, yy = 1;
    while (r != 0)
    {
        const long long q = oldR / r;
        std::tie(oldR, r) = std::make_tuple(r, oldR - q * r);
        std::tie(oldX, xx) = std::make_tuple(xx, oldX - q * xx);
        std::tie(oldY, yy) = std::make_tuple(yy, oldY - q * yy);
    }
    if (oldR < 0) { oldX = -oldX; oldY = -oldY; }
    long long q = (oldX * m[0][0] + oldY * m[1][0]) % p;
    if (q < 0) q += p;
    return {p, q, s};
}

struct Candidate
{
    int sa[2];          // A lattice coordinates of u
    int ij[2];          // B lattice coordinates of v
    glm::dvec2 u, v;    // the A vector and the (rotated) B vector it matches
};
}

bool prepareSlab(const Side& side, Structure& slab, std::string& error)
{
    if (!side.crystal.hasUnitCell || side.crystal.atoms.empty())
    {
        error = "needs a structure with a unit cell";
        return false;
    }
    if (side.cleave)
    {
        if (side.h == 0 && side.k == 0 && side.l == 0)
        {
            error = "(0 0 0) is not a plane";
            return false;
        }
        SurfaceParams params;
        params.h = side.h;
        params.k = side.k;
        params.l = side.l;
        params.layers = std::max(1, side.layers);
        params.vacuum = 0.0;
        const SurfaceResult result = buildSurface(side.crystal, params);
        if (!result.success)
        {
            error = result.message;
            return false;
        }
        slab = result.structure;
    }
    else
        slab = repeatLayersZ(side.crystal, std::max(1, side.layers));
    const glm::dvec3 a(slab.cellVectors[0][0], slab.cellVectors[0][1], slab.cellVectors[0][2]);
    const glm::dvec3 b(slab.cellVectors[1][0], slab.cellVectors[1][1], slab.cellVectors[1][2]);
    if (glm::length(glm::cross(a, b)) < 1e-8)
    {
        error = "its a and b cell vectors are parallel";
        return false;
    }
    orient(slab);
    return true;
}

std::vector<Match> findMatches(const Structure& slabA, const Structure& slabB, const SearchParams& params)
{
    std::vector<Match> matches;
    const glm::dvec2 a1(slabA.cellVectors[0][0], slabA.cellVectors[0][1]);
    const glm::dvec2 a2(slabA.cellVectors[1][0], slabA.cellVectors[1][1]);
    const glm::dvec2 b1(slabB.cellVectors[0][0], slabB.cellVectors[0][1]);
    const glm::dvec2 b2(slabB.cellVectors[1][0], slabB.cellVectors[1][1]);
    const double detA = cross2(a1, a2);
    if (std::abs(detA) < 1e-10 || std::abs(cross2(b1, b2)) < 1e-10 || params.angleStep <= 0 || params.maxRepeats < 1)
        return matches;
    const int atomsA = static_cast<int>(slabA.atoms.size()), atomsB = static_cast<int>(slabB.atoms.size());
    const double minAngle = std::clamp(params.minCellAngle, 0.0, 89.0);
    // Per-vector misfit allowed before pairing (the cell's mean strain is checked after).
    const double vectorTolerance = 3.0 * params.maxStrain + 1e-3;

    // Best match per pair of sublattices (A's, B's), over all rotations.
    std::map<std::tuple<long long, long long, long long, long long, long long, long long>, Match> best;
    const double span = params.angleMax - params.angleMin;
    // A rotation by 180 degrees maps every 2D lattice onto itself.
    const int steps = std::max(1, static_cast<int>(std::floor(span / params.angleStep + 1e-9)) + (span < 180.0 - 1e-9 ? 1 : 0));
    std::vector<Candidate> kept;
    for (int step = 0; step < steps; ++step)
    {
        const double theta = params.angleMin + step * params.angleStep;
        const double c = std::cos(theta * kPi / 180.0), s = std::sin(theta * kPi / 180.0);
        const glm::dvec2 c1(c * b1.x - s * b1.y, s * b1.x + c * b1.y);
        const glm::dvec2 c2(c * b2.x - s * b2.y, s * b2.x + c * b2.y);
        kept.clear();
        for (int i = -params.maxRepeats; i <= params.maxRepeats; ++i)
            for (int j = -params.maxRepeats; j <= params.maxRepeats; ++j)
            {
                if (i == 0 && j == 0) continue;
                const glm::dvec2 v = static_cast<double>(i) * c1 + static_cast<double>(j) * c2;
                // A lattice coordinates of v, rounded: the A vector it lies on.
                const double fa = cross2(v, a2) / detA, fb = cross2(a1, v) / detA;
                const int sa0 = static_cast<int>(std::lround(fa)), sa1 = static_cast<int>(std::lround(fb));
                if (sa0 == 0 && sa1 == 0) continue;
                const glm::dvec2 u = static_cast<double>(sa0) * a1 + static_cast<double>(sa1) * a2;
                if (glm::length(v - u) > vectorTolerance * glm::length(u)) continue;
                kept.push_back({{sa0, sa1}, {i, j}, u, v});
            }
        for (std::size_t p = 0; p < kept.size(); ++p)
            for (std::size_t q = p + 1; q < kept.size(); ++q)
            {
                Candidate first = kept[p], second = kept[q];
                const double crossU = cross2(first.u, second.u), crossV = cross2(first.v, second.v);
                if (crossU * crossV <= 0) continue;  // not the same handedness (or collinear)
                if (crossU < 0) std::swap(first, second);
                const double area = std::abs(crossU);
                if (area > params.maxArea || area < 1e-8) continue;
                // A short, nearly orthogonal (Lagrange-reduced) basis of the same cell.
                for (int guard = 0; guard < 64; ++guard)
                {
                    if (glm::dot(second.u, second.u) < glm::dot(first.u, first.u)) std::swap(first, second);
                    const double mu = std::round(glm::dot(first.u, second.u) / glm::dot(first.u, first.u));
                    if (mu == 0) break;
                    const int k = static_cast<int>(mu);
                    second.u -= mu * first.u;
                    second.v -= mu * first.v;
                    second.sa[0] -= k * first.sa[0]; second.sa[1] -= k * first.sa[1];
                    second.ij[0] -= k * first.ij[0]; second.ij[1] -= k * first.ij[1];
                }
                if (cross2(first.u, second.u) < 0)
                {
                    second.u = -second.u; second.v = -second.v;
                    second.sa[0] = -second.sa[0]; second.sa[1] = -second.sa[1];
                    second.ij[0] = -second.ij[0]; second.ij[1] = -second.ij[1];
                }
                const double cosine = glm::dot(first.u, second.u) / (glm::length(first.u) * glm::length(second.u));
                const double angle = std::acos(std::clamp(cosine, -1.0, 1.0)) * 180.0 / kPi;
                if (angle < minAngle || angle > 180.0 - minAngle) continue;
                double e[3] = {}, residual = 0;
                if (!polarStrain(first.u, second.u, first.v, second.v, e, residual)) continue;
                const double e11 = e[0], e22 = e[1], e12 = e[2];
                const double mean = (std::abs(e11) + std::abs(e22) + std::abs(e12)) / 3.0;
                if (mean > params.maxStrain) continue;
                Match match;
                match.a[0][0] = first.sa[0]; match.a[0][1] = first.sa[1];
                match.a[1][0] = second.sa[0]; match.a[1][1] = second.sa[1];
                match.b[0][0] = first.ij[0]; match.b[0][1] = first.ij[1];
                match.b[1][0] = second.ij[0]; match.b[1][1] = second.ij[1];
                match.cellsA = match.a[0][0] * match.a[1][1] - match.a[0][1] * match.a[1][0];
                match.cellsB = match.b[0][0] * match.b[1][1] - match.b[0][1] * match.b[1][0];
                if (match.cellsA <= 0 || match.cellsB <= 0) continue;
                // The rotation B ends up with: the searched one and what the strain leaves.
                match.thetaDeg = std::fmod(theta + residual + 540.0, 360.0) - 180.0;
                match.strain[0] = e11; match.strain[1] = e22; match.strain[2] = e12;
                match.meanStrain = mean;
                match.area = area;
                match.cellAngle = angle;
                match.atoms = atomsA * match.cellsA + atomsB * match.cellsB;
                const auto [pa, qa, sa] = canonical(match.a);
                const auto [pb, qb, sb] = canonical(match.b);
                const auto key = std::make_tuple(pa, qa, sa, pb, qb, sb);
                // The lowest strain; on a tie (a symmetry of the lattices) the smallest rotation.
                const auto found = best.find(key);
                if (found == best.end() || mean < found->second.meanStrain - 1e-12 ||
                    (mean < found->second.meanStrain + 1e-12 && std::abs(match.thetaDeg) < std::abs(found->second.thetaDeg) - 1e-9))
                    best[key] = match;
            }
    }

    // Symmetry-equivalent matches (same cell shape, strain and atoms) once.
    std::vector<Match> candidates;
    for (const auto& entry : best) candidates.push_back(entry.second);
    // Of equivalent ones, the least rotated is kept.
    std::stable_sort(candidates.begin(), candidates.end(), [](const Match& x, const Match& y) { return std::abs(x.thetaDeg) < std::abs(y.thetaDeg); });
    std::map<std::tuple<int, long long, long long, long long, long long>, bool> seen;
    for (const auto& match : candidates)
    {
        const double l1 = std::hypot(match.a[0][0] * a1.x + match.a[0][1] * a2.x, match.a[0][0] * a1.y + match.a[0][1] * a2.y);
        const double l2 = std::hypot(match.a[1][0] * a1.x + match.a[1][1] * a2.x, match.a[1][0] * a1.y + match.a[1][1] * a2.y);
        const auto shape = std::make_tuple(match.atoms, std::llround(std::min(l1, l2) * 1e3), std::llround(std::max(l1, l2) * 1e3),
                                           std::llround(std::min(match.cellAngle, 180.0 - match.cellAngle) * 1e2),
                                           std::llround(match.meanStrain * 1e6));
        if (seen.emplace(shape, true).second) matches.push_back(match);
    }
    std::sort(matches.begin(), matches.end(), [](const Match& x, const Match& y) {
        if (x.atoms != y.atoms) return x.atoms < y.atoms;
        if (x.meanStrain != y.meanStrain) return x.meanStrain < y.meanStrain;
        return x.thetaDeg < y.thetaDeg;
    });
    // A supercell of a smaller match (its cell made of whole copies of that
    // one, at the same rotation and no lower strain) adds nothing: dropped.
    const auto within = [](const int big[2][2], const int small[2][2]) {
        // big = X small with X integer (rows: lattice vectors).
        const long long det = static_cast<long long>(small[0][0]) * small[1][1] - static_cast<long long>(small[0][1]) * small[1][0];
        if (det == 0) return false;
        for (int r = 0; r < 2; ++r)
        {
            const long long x0 = static_cast<long long>(big[r][0]) * small[1][1] - static_cast<long long>(big[r][1]) * small[1][0];
            const long long x1 = static_cast<long long>(big[r][1]) * small[0][0] - static_cast<long long>(big[r][0]) * small[0][1];
            if (x0 % det != 0 || x1 % det != 0) return false;
        }
        return true;
    };
    std::vector<Match> distinct;
    for (const auto& match : matches)
    {
        bool redundant = false;
        for (const auto& smaller : distinct)
        {
            if (smaller.cellsA >= match.cellsA || match.cellsA % smaller.cellsA != 0 || match.cellsB % smaller.cellsB != 0 ||
                match.cellsA / smaller.cellsA != match.cellsB / smaller.cellsB)
                continue;
            if (match.meanStrain < smaller.meanStrain - 1e-9) continue;
            // The same cell repeated, or (at the same strain) a symmetry-equivalent
            // rotation of such a repeat.
            const bool sameRotation = std::abs(std::fmod(match.thetaDeg - smaller.thetaDeg + 540.0, 360.0) - 180.0) <= 0.5;
            if ((sameRotation && within(match.a, smaller.a) && within(match.b, smaller.b)) ||
                std::abs(match.meanStrain - smaller.meanStrain) < 1e-6) { redundant = true; break; }
        }
        if (!redundant) distinct.push_back(match);
    }
    matches.swap(distinct);
    for (auto& match : matches)
        if (std::abs(match.thetaDeg) < 0.005) match.thetaDeg = 0.0;
    if (matches.size() > params.maxMatches) matches.resize(params.maxMatches);
    double lowest = 1e300;
    for (auto& match : matches)
        if (match.meanStrain < lowest - 1e-12)
        {
            match.pareto = true;
            lowest = match.meanStrain;
        }
    return matches;
}

int suggestedMatch(const std::vector<Match>& matches)
{
    int chosen = -1;
    double score = 1e300;
    for (std::size_t i = 0; i < matches.size(); ++i)
    {
        if (!matches[i].pareto) continue;
        const double value = matches[i].meanStrain * 100.0 + matches[i].area / 25.0;
        if (value < score)
        {
            score = value;
            chosen = static_cast<int>(i);
        }
    }
    return chosen;
}

Structure inPlaneSupercell(const Structure& slab, const int m[2][2])
{
    const glm::dvec2 a(slab.cellVectors[0][0], slab.cellVectors[0][1]);
    const glm::dvec2 b(slab.cellVectors[1][0], slab.cellVectors[1][1]);
    const glm::dvec2 A = static_cast<double>(m[0][0]) * a + static_cast<double>(m[0][1]) * b;
    const glm::dvec2 B = static_cast<double>(m[1][0]) * a + static_cast<double>(m[1][1]) * b;
    Structure result = slab;
    result.atoms.clear();
    result.cellVectors[0] = {A.x, A.y, 0.0};
    result.cellVectors[1] = {B.x, B.y, 0.0};
    const glm::dmat2 base(a, b), cell(A, B);
    if (std::abs(glm::determinant(cell)) < 1e-12) return result;
    const glm::dmat2 inverseBase = glm::inverse(base), inverseCell = glm::inverse(cell);
    const int lo0 = std::min({0, m[0][0], m[1][0], m[0][0] + m[1][0]}) - 1, hi0 = std::max({0, m[0][0], m[1][0], m[0][0] + m[1][0]}) + 1;
    const int lo1 = std::min({0, m[0][1], m[1][1], m[0][1] + m[1][1]}) - 1, hi1 = std::max({0, m[0][1], m[1][1], m[0][1] + m[1][1]}) + 1;
    constexpr double eps = 1e-8;
    for (const auto& atom : slab.atoms)
    {
        // The atom's image in the slab's own cell first.
        glm::dvec2 f = inverseBase * glm::dvec2(atom.x, atom.y);
        f -= glm::floor(f + eps);
        const glm::dvec2 home = base * f;
        for (int n0 = lo0; n0 <= hi0; ++n0)
            for (int n1 = lo1; n1 <= hi1; ++n1)
            {
                const glm::dvec2 p = home + static_cast<double>(n0) * a + static_cast<double>(n1) * b;
                const glm::dvec2 g = inverseCell * p;
                if (g.x < -eps || g.x >= 1.0 - eps || g.y < -eps || g.y >= 1.0 - eps) continue;
                AtomSite copy = atom;
                copy.x = p.x;
                copy.y = p.y;
                result.atoms.push_back(copy);
            }
    }
    return result;
}

Structure buildInterface(const Structure& slabA, const Structure& slabB, const Match& match, const Assembly& assembly)
{
    Structure A = inPlaneSupercell(slabA, match.a);
    Structure B = inPlaneSupercell(slabB, match.b);

    // In-plane cells as columns: A's (U) and B's before rotation (V).
    const glm::dmat2 U(glm::dvec2(A.cellVectors[0][0], A.cellVectors[0][1]), glm::dvec2(A.cellVectors[1][0], A.cellVectors[1][1]));
    const glm::dmat2 V(glm::dvec2(B.cellVectors[0][0], B.cellVectors[0][1]), glm::dvec2(B.cellVectors[1][0], B.cellVectors[1][1]));
    const double c = std::cos(match.thetaDeg * kPi / 180.0), s = std::sin(match.thetaDeg * kPi / 180.0);
    const glm::dmat2 R(glm::dvec2(c, s), glm::dvec2(-s, c));
    const glm::dmat2 rotatedV = R * V;
    // The common cell: A's, B's rotated, or their average.
    const glm::dmat2 W = assembly.strainOn == StrainOn::B ? U : assembly.strainOn == StrainOn::A ? rotatedV : 0.5 * (U + rotatedV);
    // Each side's in-plane map onto W (B's includes its rotation).
    const glm::dmat2 mapA = W * glm::inverse(U), mapB = W * glm::inverse(V);
    const auto transform = [](Structure& structure, const glm::dmat2& m) {
        const double F[2][2] = {{m[0][0], m[1][0]}, {m[0][1], m[1][1]}};
        structure = applyTransform2D(structure, F);
    };
    transform(A, mapA);
    transform(B, mapB);

    const auto zRange = [](const Structure& structure) {
        double low = 1e300, high = -1e300;
        for (const auto& atom : structure.atoms)
        {
            low = std::min(low, atom.z);
            high = std::max(high, atom.z);
        }
        return std::make_pair(low, high);
    };
    const auto [lowA, highA] = zRange(A);
    const auto [lowB, highB] = zRange(B);
    const double heightA = highA - lowA, heightB = highB - lowB;
    const double gap = std::max(0.0, assembly.gap), vacuum = std::max(0.0, assembly.vacuum);
    // Periodic: half a gap below A, so the gap B-top to A-bottom (through the
    // boundary) equals the one between A and B. With vacuum: half of it below.
    const double base = vacuum > 0 ? 0.5 * vacuum : 0.5 * gap;
    const double height = heightA + gap + heightB + (vacuum > 0 ? vacuum : gap);
    const glm::dvec2 shift = assembly.shift[0] * W[0] + assembly.shift[1] * W[1];

    Structure result;
    result.hasUnitCell = true;
    result.cellVectors = {{{W[0].x, W[0].y, 0.0}, {W[1].x, W[1].y, 0.0}, {0.0, 0.0, height}}};
    result.cellOffset = {0.0, 0.0, 0.0};
    result.atoms.reserve(A.atoms.size() + B.atoms.size());
    for (auto atom : A.atoms)
    {
        atom.z = atom.z - lowA + base;
        result.atoms.push_back(atom);
    }
    for (auto atom : B.atoms)
    {
        atom.x += shift.x;
        atom.y += shift.y;
        atom.z = atom.z - lowB + base + heightA + gap;
        result.atoms.push_back(atom);
    }
    // In-plane positions wrapped into the cell.
    const glm::dmat2 inverseW = glm::inverse(W);
    for (auto& atom : result.atoms)
    {
        glm::dvec2 f = inverseW * glm::dvec2(atom.x, atom.y);
        f -= glm::floor(f + 1e-9);
        const glm::dvec2 p = W * f;
        atom.x = p.x;
        atom.y = p.y;
    }
    return result;
}
}
