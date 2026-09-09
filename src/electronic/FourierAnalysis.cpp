#include "electronic/Volume.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace atomforge::electronic
{
namespace
{
using Complex = std::complex<double>;
void fft(std::vector<Complex>& a, double sign)
{
    const std::size_t n = a.size();
    if ((n & (n - 1)) != 0)
    {
        // Non-power-of-two axes use an exact DFT; workload is checked by caller.
        const auto original = a;
        for (std::size_t k = 0; k < n; ++k)
        {
            Complex sum = 0;
            for (std::size_t j = 0; j < n; ++j) sum += original[j] * std::polar(1.0, sign * 2 * pi * k * j / n);
            a[k] = sum;
        }
        return;
    }
    for (std::size_t i = 1, j = 0; i < n; ++i)
    {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (std::size_t length = 2; length <= n; length <<= 1)
    {
        const Complex root = std::polar(1.0, sign * 2 * pi / length);
        for (std::size_t start = 0; start < n; start += length)
        {
            Complex w = 1;
            for (std::size_t j = 0; j < length / 2; ++j)
            {
                const Complex even = a[start + j], odd = a[start + j + length / 2] * w;
                a[start + j] = even + odd;
                a[start + j + length / 2] = even - odd;
                w *= root;
            }
        }
    }
}
void transform(const Grid& grid, std::vector<Complex>& data, double sign)
{
    double work = 0;
    for (int n : grid.shape) work += ((n & (n - 1)) ? n : std::log2(n)) * static_cast<double>(data.size());
    if (work > 2e8) throw std::invalid_argument("Fourier workload too large; resample to power-of-two dimensions");
    for (int axis = 0; axis < 3; ++axis)
    {
        const int b = (axis + 1) % 3, c = (axis + 2) % 3;
        std::vector<Complex> row(grid.shape[axis]);
        for (int j = 0; j < grid.shape[b]; ++j)
            for (int k = 0; k < grid.shape[c]; ++k)
            {
                glm::ivec3 p(0);
                p[b] = j; p[c] = k;
                for (int i = 0; i < grid.shape[axis]; ++i) { p[axis] = i; row[i] = data[grid.index(p.x, p.y, p.z)]; }
                fft(row, sign);
                for (int i = 0; i < grid.shape[axis]; ++i) { p[axis] = i; data[grid.index(p.x, p.y, p.z)] = row[i]; }
            }
    }
}
void requirePeriodic(const Grid& grid)
{
    grid.validate();
    if (!grid.periodic) throw std::invalid_argument("Fourier tools require periodic endpoint-excluded grids; convert explicitly");
}
}

std::vector<Reflection> structureFactors(const Grid& grid, const std::vector<glm::ivec3>& indices)
{
    requirePeriodic(grid);
    std::vector<Complex> data(grid.values.begin(), grid.values.end());
    transform(grid, data, 1);
    const double voxel = std::abs(glm::determinant(grid.cell)) / grid.values.size();
    const auto origin = glm::inverse(grid.cell) * grid.origin;
    std::vector<Reflection> out;
    for (auto h : indices)
    {
        auto q = h;
        for (int a = 0; a < 3; ++a)
        {
            if (h[a] < -grid.shape[a] / 2 || h[a] > (grid.shape[a] - 1) / 2)
                throw std::invalid_argument("Reflection outside unique sampled Nyquist range");
            q[a] = (h[a] + grid.shape[a]) % grid.shape[a];
        }
        out.push_back({h, data[grid.index(q.x, q.y, q.z)] * voxel * std::polar(1.0, 2 * pi * glm::dot(glm::dvec3(h), origin))});
    }
    return out;
}

Grid fourierSynthesis(const Grid& geometry, const std::vector<Reflection>& reflections)
{
    requirePeriodic(geometry);
    std::vector<Complex> data(geometry.values.size(), 0);
    std::vector<bool> occupied(data.size(), false);
    const auto origin = glm::inverse(geometry.cell) * geometry.origin;
    for (const auto& r : reflections)
    {
        auto q = r.hkl;
        for (int a = 0; a < 3; ++a)
        {
            if (q[a] < -geometry.shape[a] / 2 || q[a] > (geometry.shape[a] - 1) / 2)
                throw std::invalid_argument("Reflection outside unique sampled Nyquist range");
            q[a] = (q[a] + geometry.shape[a]) % geometry.shape[a];
        }
        const auto i = geometry.index(q.x, q.y, q.z);
        if (occupied[i]) throw std::invalid_argument("Duplicate reflection");
        if (!std::isfinite(r.value.real()) || !std::isfinite(r.value.imag())) throw std::invalid_argument("Non-finite structure factor");
        occupied[i] = true;
        data[i] = r.value * std::polar(1.0, -2 * pi * glm::dot(glm::dvec3(r.hkl), origin));
    }
    transform(geometry, data, -1);
    Grid out = geometry;
    out.name = "Fourier synthesis";
    const double volume = std::abs(glm::determinant(out.cell));
    double realMax = 0, imaginaryMax = 0;
    for (auto v : data) { realMax = std::max(realMax, std::abs(v.real())); imaginaryMax = std::max(imaginaryMax, std::abs(v.imag())); }
    if (imaginaryMax > 1e-8 * std::max(1.0, realMax))
        throw std::invalid_argument("Structure factors do not produce a real density; supply conjugate/Friedel pairs explicitly");
    for (std::size_t i = 0; i < data.size(); ++i) out.values[i] = data[i].real() / volume;
    out.validate();
    return out;
}

Grid patterson(const Grid& grid)
{
    requirePeriodic(grid);
    std::vector<Complex> data(grid.values.begin(), grid.values.end());
    transform(grid, data, 1);
    for (auto& v : data) v = std::norm(v);
    transform(grid, data, -1);
    Grid out = grid;
    out.origin = glm::dvec3(0);
    out.name = "Patterson autocorrelation";
    out.unit = "(" + grid.unit + ")^2*A^3";
    const double factor = std::abs(glm::determinant(grid.cell)) / (static_cast<double>(data.size()) * data.size());
    for (std::size_t i = 0; i < data.size(); ++i) out.values[i] = data[i].real() * factor;
    out.validate();
    return out;
}

std::vector<Reflection> atomicStructureFactors(const Grid& geometry,
    const std::vector<glm::dvec3>& positions, const std::vector<double>& occupancy,
    const std::vector<double>& bFactors, const std::vector<glm::ivec3>& indices,
    const std::vector<std::complex<double>>& scattering)
{
    requirePeriodic(geometry);
    const auto count = positions.size();
    if (occupancy.size() != count || bFactors.size() != count || scattering.size() != indices.size() * count)
        throw std::invalid_argument("Atomic factors require one occupancy/B per site and one scattering factor per reflection/site");
    const auto inverse = glm::inverse(geometry.cell);
    std::vector<Reflection> out;
    for (std::size_t r = 0; r < indices.size(); ++r)
    {
        const glm::dvec3 h(indices[r]);
        const auto reciprocal = glm::transpose(inverse) * h;
        const double s2 = glm::dot(reciprocal, reciprocal) / 4;
        Complex sum = 0;
        for (std::size_t j = 0; j < count; ++j)
        {
            if (!std::isfinite(occupancy[j]) || occupancy[j] < 0 || occupancy[j] > 1 || !std::isfinite(bFactors[j]) || bFactors[j] < 0)
                throw std::invalid_argument("Occupancy must be 0..1 and B must be nonnegative");
            const auto factor = scattering[r * count + j];
            const auto p = inverse * positions[j];
            for (int a = 0; a < 3; ++a) if (!std::isfinite(p[a])) throw std::invalid_argument("Non-finite atomic position");
            if (!std::isfinite(factor.real()) || !std::isfinite(factor.imag())) throw std::invalid_argument("Non-finite scattering factor");
            sum += occupancy[j] * factor * std::exp(-bFactors[j] * s2) * std::polar(1.0, 2 * pi * glm::dot(h, p));
        }
        out.push_back({indices[r], sum});
    }
    return out;
}
}
