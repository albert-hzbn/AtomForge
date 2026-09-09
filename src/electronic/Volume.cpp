#include "electronic/Volume.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace atomforge::electronic
{
void Grid::validate() const
{
    std::size_t count = 1;
    for (int n : shape)
    {
        if (n < 2 || n > 100000 || count > 100000000 / static_cast<std::size_t>(n))
            throw std::invalid_argument("Grid dimensions must be >= 2 and total samples <= 100 million");
        count *= n;
    }
    if (values.size() != count) throw std::invalid_argument("Grid value count does not match dimensions");
    for (int i = 0; i < 3; ++i)
    {
        if (!std::isfinite(origin[i])) throw std::invalid_argument("Non-finite grid origin");
        for (int j = 0; j < 3; ++j)
            if (!std::isfinite(cell[i][j])) throw std::invalid_argument("Non-finite lattice");
    }
    const double determinant = glm::determinant(cell);
    if (!std::isfinite(determinant) || std::abs(determinant) < 1e-12)
        throw std::invalid_argument("Singular grid lattice");
    for (double value : values)
        if (!std::isfinite(value)) throw std::invalid_argument("Non-finite grid value");
}

std::size_t Grid::index(int x, int y, int z) const
{
    return (static_cast<std::size_t>(z) * shape[1] + y) * shape[0] + x;
}

glm::dmat3 Grid::steps() const
{
    glm::dmat3 result = cell;
    for (int i = 0; i < 3; ++i) result[i] /= shape[i] - (periodic ? 0 : 1);
    return result;
}

glm::dvec3 Grid::position(int x, int y, int z) const
{
    return origin + steps() * glm::dvec3(x, y, z);
}

double Grid::weight(int x, int y, int z) const
{
    double w = std::abs(glm::determinant(steps()));
    if (!periodic)
    {
        if (x == 0 || x == shape[0] - 1) w *= 0.5;
        if (y == 0 || y == shape[1] - 1) w *= 0.5;
        if (z == 0 || z == shape[2] - 1) w *= 0.5;
    }
    return w;
}

double Grid::sample(const glm::dvec3& point) const
{
    glm::dvec3 q = glm::inverse(steps()) * (point - origin);
    glm::ivec3 low, high;
    glm::dvec3 t;
    for (int i = 0; i < 3; ++i)
    {
        if (!std::isfinite(q[i])) throw std::invalid_argument("Non-finite sampling point");
        if (periodic) q[i] -= std::floor(q[i] / shape[i]) * shape[i];
        else
        {
            if (q[i] < -1e-8 || q[i] > shape[i] - 1 + 1e-8)
                throw std::out_of_range("Sampling point lies outside finite grid");
            q[i] = std::clamp(q[i], 0.0, static_cast<double>(shape[i] - 1));
        }
        low[i] = std::min(static_cast<int>(std::floor(q[i])), shape[i] - 1);
        high[i] = periodic ? (low[i] + 1) % shape[i] : std::min(low[i] + 1, shape[i] - 1);
        t[i] = q[i] - low[i];
    }
    double value = 0.0;
    for (int z = 0; z < 2; ++z)
        for (int y = 0; y < 2; ++y)
            for (int x = 0; x < 2; ++x)
                value += values[index(x ? high.x : low.x, y ? high.y : low.y, z ? high.z : low.z)]
                    * (x ? t.x : 1 - t.x) * (y ? t.y : 1 - t.y) * (z ? t.z : 1 - t.z);
    return value;
}

void requireAligned(const Grid& a, const Grid& b)
{
    a.validate();
    b.validate();
    if (a.shape != b.shape || a.periodic != b.periodic || glm::length(a.origin - b.origin) > 1e-8)
        throw std::invalid_argument("Grids are not aligned; resample explicitly first");
    for (int i = 0; i < 3; ++i)
        if (glm::length(a.cell[i] - b.cell[i]) > 1e-8)
            throw std::invalid_argument("Grid lattices differ; resample explicitly first");
}

Grid arithmetic(const Grid& a, const Grid& b, const std::string& operation)
{
    requireAligned(a, b);
    if ((operation == "add" || operation == "subtract") && a.unit != b.unit)
        throw std::invalid_argument("Addition and subtraction require matching units");
    Grid out = a;
    out.name = operation;
    if (operation == "multiply") out.unit = "(" + a.unit + ")*(" + b.unit + ")";
    else if (operation == "divide") out.unit = a.unit == b.unit ? "1" : "(" + a.unit + ")/(" + b.unit + ")";
    else if (operation != "add" && operation != "subtract") throw std::invalid_argument("Unknown arithmetic operation");
    for (std::size_t i = 0; i < out.values.size(); ++i)
    {
        if (operation == "add") out.values[i] += b.values[i];
        else if (operation == "subtract") out.values[i] -= b.values[i];
        else if (operation == "multiply") out.values[i] *= b.values[i];
        else
        {
            if (b.values[i] == 0) throw std::invalid_argument("Division by zero in grid arithmetic");
            out.values[i] /= b.values[i];
        }
    }
    out.validate();
    return out;
}

Grid scale(const Grid& a, double factor)
{
    a.validate();
    if (!std::isfinite(factor)) throw std::invalid_argument("Scale must be finite");
    Grid out = a;
    for (double& value : out.values) value *= factor;
    out.validate();
    return out;
}

Grid resample(const Grid& source, const Grid& target)
{
    source.validate();
    target.validate();
    Grid out = target;
    out.name = source.name;
    out.unit = source.unit;
    for (int z = 0; z < out.shape[2]; ++z)
        for (int y = 0; y < out.shape[1]; ++y)
            for (int x = 0; x < out.shape[0]; ++x)
                out.values[out.index(x, y, z)] = source.sample(out.position(x, y, z));
    return out;
}

glm::dvec3 minimumImage(glm::dvec3 delta, const glm::dmat3& cell)
{
    // Enumerate all translations that can improve the initial candidate. Unlike
    // component rounding alone, this also works for non-reduced triclinic cells.
    const auto inverse = glm::inverse(cell);
    const glm::dvec3 f = inverse * delta;
    glm::dvec3 best = cell * (f - glm::round(f));
    const double radius = glm::length(best);
    glm::ivec3 low, high;
    std::size_t count = 1;
    for (int i = 0; i < 3; ++i)
    {
        double bound = glm::length(glm::dvec3(inverse[0][i], inverse[1][i], inverse[2][i])) * radius + 1e-10;
        if (!std::isfinite(f[i]) || bound > 10000 || std::abs(f[i]) > 1e8)
            throw std::invalid_argument("Ill-conditioned lattice or invalid coordinates");
        low[i] = static_cast<int>(std::ceil(f[i] - bound));
        high[i] = static_cast<int>(std::floor(f[i] + bound));
        count *= std::max(1, high[i] - low[i] + 1);
        if (count > 1000000) throw std::invalid_argument("Lattice requires excessive periodic images; reduce the cell");
    }
    for (int z = low.z; z <= high.z; ++z)
        for (int y = low.y; y <= high.y; ++y)
            for (int x = low.x; x <= high.x; ++x)
            {
                const glm::dvec3 candidate = delta - cell * glm::dvec3(x, y, z);
                if (glm::dot(candidate, candidate) < glm::dot(best, best)) best = candidate;
            }
    return best;
}

Grid periodicEndpoints(const Grid& source, double tolerance)
{
    source.validate();
    if (!std::isfinite(tolerance) || tolerance < 0) throw std::invalid_argument("Endpoint tolerance must be nonnegative and finite");
    if (source.periodic) return source;
    Grid out = source;
    for (int& n : out.shape)
    {
        if (n < 3) throw std::invalid_argument("Need >=3 samples to remove periodic endpoint planes");
        --n;
    }
    for (int z = 0; z < source.shape[2]; ++z)
        for (int y = 0; y < source.shape[1]; ++y)
            for (int x = 0; x < source.shape[0]; ++x)
            {
                const double a = source.values[source.index(x,y,z)];
                const double b = source.values[source.index(x%out.shape[0],y%out.shape[1],z%out.shape[2])];
                if (std::abs(a-b) > tolerance * std::max({1.0,std::abs(a),std::abs(b)}))
                    throw std::invalid_argument("Opposite endpoint planes do not match");
            }
    out.periodic = true;
    out.values.resize(static_cast<std::size_t>(out.shape[0])*out.shape[1]*out.shape[2]);
    for (int z = 0; z < out.shape[2]; ++z)
        for (int y = 0; y < out.shape[1]; ++y)
            for (int x = 0; x < out.shape[0]; ++x)
                out.values[out.index(x,y,z)] = source.values[source.index(x,y,z)];
    return out;
}
}
