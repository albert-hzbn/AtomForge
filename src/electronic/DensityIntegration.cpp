#include "electronic/Volume.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace atomforge::electronic
{
double integrate(const Grid& grid)
{
    grid.validate();
    double sum = 0, correction = 0;
    for (int z = 0; z < grid.shape[2]; ++z)
        for (int y = 0; y < grid.shape[1]; ++y)
            for (int x = 0; x < grid.shape[0]; ++x)
            {
                const double v = grid.values[grid.index(x, y, z)] * grid.weight(x, y, z) - correction;
                const double next = sum + v;
                correction = (next - sum) - v;
                sum = next;
            }
    return sum;
}

double integrateSphere(const Grid& grid, glm::dvec3 center, double radius)
{
    grid.validate();
    if (!std::isfinite(radius) || radius <= 0) throw std::invalid_argument("Sphere radius must be positive");
    for (int i = 0; i < 3; ++i) if (!std::isfinite(center[i])) throw std::invalid_argument("Non-finite sphere center");
    double sum = 0;
    for (int z = 0; z < grid.shape[2]; ++z)
        for (int y = 0; y < grid.shape[1]; ++y)
            for (int x = 0; x < grid.shape[0]; ++x)
            {
                auto d = grid.position(x, y, z) - center;
                if (grid.periodic) d = minimumImage(d, grid.cell);
                if (glm::dot(d, d) <= radius * radius) sum += grid.values[grid.index(x, y, z)] * grid.weight(x, y, z);
            }
    return sum;
}

std::vector<Basin> voronoiIntegrate(const Grid& grid, const std::vector<glm::dvec3>& sites)
{
    grid.validate();
    if (sites.empty()) throw std::invalid_argument("Voronoi integration requires sites");
    for (auto site : sites) for (int a = 0; a < 3; ++a)
        if (!std::isfinite(site[a])) throw std::invalid_argument("Non-finite Voronoi site");
    if (static_cast<double>(sites.size()) * grid.values.size() > 1e9) throw std::invalid_argument("Voronoi workload too large; reduce grid or site count");
    std::vector<Basin> out(sites.size());
    for (int z = 0; z < grid.shape[2]; ++z)
        for (int y = 0; y < grid.shape[1]; ++y)
            for (int x = 0; x < grid.shape[0]; ++x)
            {
                double best = std::numeric_limits<double>::infinity();
                std::vector<std::size_t> nearest;
                const auto p = grid.position(x, y, z);
                for (std::size_t i = 0; i < sites.size(); ++i)
                {
                    auto d = p - sites[i];
                    if (grid.periodic) d = minimumImage(d, grid.cell);
                    const double d2 = glm::dot(d, d);
                    if (d2 < best - 1e-10) { best = d2; nearest = {i}; }
                    else if (std::abs(d2 - best) <= 1e-10) nearest.push_back(i);
                }
                const double w = grid.weight(x, y, z) / nearest.size();
                for (auto i : nearest) { out[i].volume += w; out[i].integral += w * grid.values[grid.index(x, y, z)]; }
            }
    return out;
}

std::vector<Peak> peaks(const Grid& grid, int limit)
{
    grid.validate();
    if (limit < 0) throw std::invalid_argument("Peak limit must be nonnegative");
    std::vector<Peak> out;
    for (int z = 0; z < grid.shape[2]; ++z)
        for (int y = 0; y < grid.shape[1]; ++y)
            for (int x = 0; x < grid.shape[0]; ++x)
            {
                const auto index = grid.index(x, y, z);
                const double value = grid.values[index];
                bool maximum = true, strictlyGreater = false;
                for (int dz = -1; dz <= 1; ++dz)
                    for (int dy = -1; dy <= 1; ++dy)
                        for (int dx = -1; dx <= 1; ++dx)
                        {
                            auto p = glm::ivec3(x + dx, y + dy, z + dz);
                            bool valid = true;
                            for (int a = 0; a < 3; ++a)
                                if (grid.periodic) p[a] = (p[a] + grid.shape[a]) % grid.shape[a];
                                else if (p[a] < 0 || p[a] >= grid.shape[a]) valid = false;
                            if (!valid) continue;
                            const auto j = grid.index(p.x, p.y, p.z);
                            if (grid.values[j] > value || (grid.values[j] == value && j < index)) maximum = false;
                            if (grid.values[j] < value) strictlyGreater = true;
                        }
                if (maximum && strictlyGreater) out.push_back({grid.position(x, y, z), value});
            }
    std::stable_sort(out.begin(), out.end(), [](const Peak& a, const Peak& b) { return a.value > b.value; });
    if (limit && out.size() > static_cast<std::size_t>(limit)) out.resize(limit);
    return out;
}
}
