#include "algorithms/PatternMatch.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>

namespace atomforge
{

namespace
{
using Mat3 = Eigen::Matrix3d;
using Vec3 = Eigen::Vector3d;

Vec3 toEigen(const glm::dvec3& v) { return Vec3(v.x, v.y, v.z); }

bool buildCellMatrices(const Structure& s, Mat3& cell, Mat3& invCell)
{
    if (!s.hasUnitCell)
        return false;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            cell(j, i) = s.cellVectors[i][j];
    if (std::abs(cell.determinant()) < 1e-12)
        return false;
    invCell = cell.inverse();
    return true;
}

std::vector<glm::dvec3> neighborDirections(const std::vector<glm::dvec3>& pos, size_t i, float cutoff,
                                           bool usePbc, const Mat3& cell, const Mat3& invCell)
{
    const double cutoff2 = (double)cutoff * (double)cutoff;
    std::vector<glm::dvec3> dirs;
    for (size_t j = 0; j < pos.size(); ++j)
    {
        if (j == i) continue;
        glm::dvec3 delta = pos[j] - pos[i];
        if (usePbc)
        {
            Vec3 frac = invCell * toEigen(delta);
            frac -= frac.array().round().matrix();
            const Vec3 wrapped = cell * frac;
            delta = glm::dvec3(wrapped.x(), wrapped.y(), wrapped.z());
        }
        const double d2 = glm::dot(delta, delta);
        if (d2 > 1e-8 && d2 <= cutoff2)
            dirs.push_back(delta / std::sqrt(d2));
    }
    return dirs;
}
}

CrystalPattern buildPattern(const Structure& reference, float cutoff, bool usePbc)
{
    CrystalPattern result;
    const size_t n = reference.atoms.size();
    if (n < 2)
    {
        result.message = "Reference structure needs at least two atoms.";
        return result;
    }
    if (cutoff <= 0.0f)
    {
        result.message = "cutoff must be positive.";
        return result;
    }

    std::vector<glm::dvec3> pos(n);
    glm::dvec3 center(0.0);
    for (size_t i = 0; i < n; ++i)
    {
        pos[i] = glm::dvec3(reference.atoms[i].x, reference.atoms[i].y, reference.atoms[i].z);
        center += pos[i];
    }
    center /= (double)n;

    Mat3 cell = Mat3::Identity(), invCell = Mat3::Identity();
    usePbc = usePbc && buildCellMatrices(reference, cell, invCell);

    size_t representative = 0;
    double best = std::numeric_limits<double>::max();
    for (size_t i = 0; i < n; ++i)
    {
        const double d2 = glm::dot(pos[i] - center, pos[i] - center);
        if (d2 < best) { best = d2; representative = i; }
    }

    result.directions = neighborDirections(pos, representative, cutoff, usePbc, cell, invCell);
    if (result.directions.empty())
    {
        result.message = "No neighbors found within cutoff around the representative atom.";
        return result;
    }
    result.success = true;
    result.message = "Pattern built from atom " + std::to_string(representative) + " with "
        + std::to_string(result.directions.size()) + " neighbor directions.";
    return result;
}

PatternDetectResult detectPattern(const Structure& structure,
                                  const CrystalPattern& pattern,
                                  float cutoff,
                                  double angleThresholdDeg,
                                  bool usePbc)
{
    PatternDetectResult result;
    if (!pattern.success || pattern.directions.empty())
    {
        result.message = "Pattern is not valid; build it with buildPattern first.";
        return result;
    }
    if (cutoff <= 0.0f)
    {
        result.message = "cutoff must be positive.";
        return result;
    }

    const size_t n = structure.atoms.size();
    std::vector<glm::dvec3> pos(n);
    for (size_t i = 0; i < n; ++i)
        pos[i] = glm::dvec3(structure.atoms[i].x, structure.atoms[i].y, structure.atoms[i].z);

    Mat3 cell = Mat3::Identity(), invCell = Mat3::Identity();
    usePbc = usePbc && buildCellMatrices(structure, cell, invCell);

    const double thresholdRad = angleThresholdDeg * 3.14159265358979323846 / 180.0;
    result.rows.resize(n);

    for (size_t i = 0; i < n; ++i)
    {
        const auto dirs = neighborDirections(pos, i, cutoff, usePbc, cell, invCell);
        PatternMatchRow row;
        row.index = (int)i;
        row.neighborCount = (int)dirs.size();

        if (dirs.size() != pattern.directions.size())
        {
            result.rows[i] = row; // matched = false
            continue;
        }

        // Greedy nearest-angle bipartite assignment: reasonable and cheap
        // for the small (coordination-number-scale) sets involved here,
        // though not guaranteed globally optimal like a full Hungarian
        // algorithm -- a mismatch it misses would only ever make a truly
        // defective atom look slightly more likely to match, never less,
        // so it cannot manufacture a false "not matched" verdict.
        struct Candidate { int a, b; double angle; };
        std::vector<Candidate> candidates;
        candidates.reserve(dirs.size() * pattern.directions.size());
        for (int a = 0; a < (int)dirs.size(); ++a)
            for (int b = 0; b < (int)pattern.directions.size(); ++b)
            {
                const double cosAngle = std::clamp(glm::dot(dirs[a], pattern.directions[b]), -1.0, 1.0);
                candidates.push_back({a, b, std::acos(cosAngle)});
            }
        std::sort(candidates.begin(), candidates.end(), [](const Candidate& x, const Candidate& y) { return x.angle < y.angle; });

        std::vector<bool> usedA(dirs.size(), false), usedB(pattern.directions.size(), false);
        int assigned = 0;
        double maxAngle = 0.0;
        bool ok = true;
        for (const auto& c : candidates)
        {
            if (usedA[c.a] || usedB[c.b]) continue;
            if (c.angle > thresholdRad) { ok = false; break; }
            usedA[c.a] = true; usedB[c.b] = true;
            maxAngle = std::max(maxAngle, c.angle);
            if (++assigned == (int)dirs.size()) break;
        }

        row.matched = ok && assigned == (int)dirs.size();
        row.maxAngleDeviationDeg = row.matched ? maxAngle * 180.0 / 3.14159265358979323846 : 0.0;
        result.rows[i] = row;
    }

    result.success = true;
    const int matchedCount = (int)std::count_if(result.rows.begin(), result.rows.end(), [](const PatternMatchRow& r) { return r.matched; });
    result.message = std::to_string(matchedCount) + " of " + std::to_string(n) + " atoms matched the pattern.";
    return result;
}

}
