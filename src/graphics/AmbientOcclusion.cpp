#include "graphics/AmbientOcclusion.h"

#include <algorithm>
#include <cmath>
#include <thread>
#include <unordered_map>

namespace atomforge::display
{
namespace
{
struct CellHash
{
    std::size_t operator()(const glm::ivec3& c) const
    {
        return static_cast<std::size_t>(c.x) * 73856093u ^ static_cast<std::size_t>(c.y) * 19349663u ^
               static_cast<std::size_t>(c.z) * 83492791u;
    }
};

glm::ivec3 cellOf(const glm::vec3& p, float size)
{
    return glm::ivec3(static_cast<int>(std::floor(p.x / size)), static_cast<int>(std::floor(p.y / size)),
                      static_cast<int>(std::floor(p.z / size)));
}
}  // namespace

const std::vector<glm::vec3>& occlusionDirections()
{
    // A Fibonacci lattice of 48 points: even coverage without a preferred axis.
    static const std::vector<glm::vec3> directions = [] {
        constexpr int count = 48;
        const float golden = 3.14159265358979f * (3.0f - std::sqrt(5.0f));
        std::vector<glm::vec3> out;
        for (int i = 0; i < count; ++i) {
            const float y = 1.0f - 2.0f * (static_cast<float>(i) + 0.5f) / count;
            const float r = std::sqrt(std::max(0.0f, 1.0f - y * y));
            const float a = golden * static_cast<float>(i);
            out.emplace_back(r * std::cos(a), y, r * std::sin(a));
        }
        return out;
    }();
    return directions;
}

std::vector<float> ambientOcclusion(const std::vector<glm::vec3>& centres, const std::vector<float>& radii,
                                    float reach, std::size_t maxAtoms)
{
    const std::size_t n = centres.size();
    if (n == 0 || n > maxAtoms || radii.size() != n || !(reach > 0.0f)) return {};
    float largest = 0.0f;
    for (float r : radii) largest = std::max(largest, r);
    // Every sphere that a ray can meet lies within 2 * largest + reach of the centre.
    const float range = 2.0f * largest + reach;
    const float cellSize = std::max(range, 1e-3f);
    std::unordered_map<glm::ivec3, std::vector<std::size_t>, CellHash> grid;
    grid.reserve(n);
    for (std::size_t i = 0; i < n; ++i) grid[cellOf(centres[i], cellSize)].push_back(i);
    const auto& directions = occlusionDirections();

    std::vector<float> occlusion(n, 0.0f);
    const auto work = [&](std::size_t begin, std::size_t end) {
        std::vector<std::size_t> near;
        std::vector<char> blocked(directions.size());
        for (std::size_t i = begin; i < end; ++i) {
            const glm::vec3 c = centres[i];
            const float limit = radii[i] + reach;
            near.clear();
            const glm::ivec3 home = cellOf(c, cellSize);
            for (int dx = -1; dx <= 1; ++dx)
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dz = -1; dz <= 1; ++dz) {
                        const auto it = grid.find(home + glm::ivec3(dx, dy, dz));
                        if (it == grid.end()) continue;
                        for (std::size_t j : it->second)
                            if (j != i && glm::length(centres[j] - c) < limit + radii[j]) near.push_back(j);
                    }
            std::fill(blocked.begin(), blocked.end(), 0);
            for (std::size_t j : near) {
                const glm::vec3 v = centres[j] - c;
                const float distance2 = glm::dot(v, v);
                const float rj2 = radii[j] * radii[j];
                for (std::size_t k = 0; k < directions.size(); ++k) {
                    if (blocked[k]) continue;
                    const float t = glm::dot(v, directions[k]);
                    if (t <= 0.0f) continue;
                    const float miss2 = distance2 - t * t;
                    if (miss2 >= rj2) continue;
                    // The first point of the ray inside sphere j lies within the reach.
                    if (t - std::sqrt(rj2 - miss2) < limit) blocked[k] = 1;
                }
            }
            int count = 0;
            for (char b : blocked) count += b;
            occlusion[i] = static_cast<float>(count) / static_cast<float>(directions.size());
        }
    };
    const std::size_t threads = std::clamp<std::size_t>(std::thread::hardware_concurrency(), 1, 16);
    if (n < 2000 || threads == 1) {
        work(0, n);
        return occlusion;
    }
    std::vector<std::thread> pool;
    const std::size_t chunk = (n + threads - 1) / threads;
    for (std::size_t t = 0; t < threads; ++t) {
        const std::size_t begin = t * chunk, end = std::min(n, begin + chunk);
        if (begin < end) pool.emplace_back(work, begin, end);
    }
    for (auto& thread : pool) thread.join();
    return occlusion;
}
void darkenOccluded(const std::vector<glm::vec3>& centres, const std::vector<float>& radii,
                    std::vector<glm::vec3>& colours, float strength)
{
    if (!(strength > 0.0f) || centres.empty() || radii.size() != centres.size() || colours.size() != centres.size()) return;
    float mean = 0.0f;
    for (float r : radii) mean += r;
    mean /= static_cast<float>(radii.size());
    const std::vector<float> occlusion = ambientOcclusion(centres, radii, 4.0f * mean);
    for (std::size_t i = 0; i < occlusion.size(); ++i)
        colours[i] *= 1.0f - std::min(strength, 1.0f) * occlusion[i];
}
}  // namespace atomforge::display
