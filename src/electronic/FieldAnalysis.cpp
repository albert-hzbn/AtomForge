#include "electronic/Volume.h"
#include "util/TaskControl.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace atomforge::electronic
{
namespace
{
Grid derivative(const Grid& g, int axis, bool second)
{
    Grid out = g;
    const int n = g.shape[axis];
    if (!g.periodic && n < (second ? 4 : 3))
        throw std::invalid_argument("Finite differences require >=3 samples per axis (>=4 for Laplacian)");
    for (int z = 0; z < g.shape[2]; ++z)
        for (int y = 0; y < g.shape[1]; ++y)
            for (int x = 0; x < g.shape[0]; ++x)
            {
                if (x==0 && y==0) atomforge::taskProgress(double(z)/g.shape[2]);
                glm::ivec3 p(x, y, z);
                const int i = p[axis];
                auto at = [&](int j)
                {
                    auto q = p;
                    q[axis] = g.periodic ? (j % n + n) % n : j;
                    return g.values[g.index(q.x, q.y, q.z)];
                };
                double value;
                if (g.periodic || (i > 0 && i < n - 1))
                    value = second ? at(i + 1) - 2 * at(i) + at(i - 1) : (at(i + 1) - at(i - 1)) * 0.5;
                else if (i == 0)
                    value = second ? 2 * at(0) - 5 * at(1) + 4 * at(2) - at(3) : (-3 * at(0) + 4 * at(1) - at(2)) * 0.5;
                else
                    value = second ? 2 * at(i) - 5 * at(i - 1) + 4 * at(i - 2) - at(i - 3) : (3 * at(i) - 4 * at(i - 1) + at(i - 2)) * 0.5;
                out.values[g.index(x, y, z)] = value;
            }
    return out;
}
}

std::array<Grid, 3> gradient(const Grid& grid)
{
    grid.validate();
    std::array<Grid, 3> result{grid, grid, grid};
    const std::array<Grid, 3> d{derivative(grid, 0, false), derivative(grid, 1, false), derivative(grid, 2, false)};
    const auto transform = glm::transpose(glm::inverse(grid.steps()));
    for (std::size_t i = 0; i < grid.values.size(); ++i)
    {
        const auto v = transform * glm::dvec3(d[0].values[i], d[1].values[i], d[2].values[i]);
        for (int axis = 0; axis < 3; ++axis) result[axis].values[i] = v[axis];
    }
    for (int axis = 0; axis < 3; ++axis)
    {
        result[axis].name = "gradient " + std::to_string(axis);
        result[axis].unit = "(" + grid.unit + ")/A";
        result[axis].validate();
    }
    return result;
}

Grid laplacian(const Grid& grid)
{
    grid.validate();
    Grid result = grid;
    std::fill(result.values.begin(), result.values.end(), 0.0);
    const auto inverse = glm::inverse(grid.steps());
    const auto metric = inverse * glm::transpose(inverse);
    for (int a = 0; a < 3; ++a)
        for (int b = a; b < 3; ++b)
        {
            if (std::abs(metric[a][b]) < 1e-16) continue;
            const auto d = a == b ? derivative(grid, a, true) : derivative(derivative(grid, a, false), b, false);
            for (std::size_t i = 0; i < result.values.size(); ++i)
                result.values[i] += (a == b ? 1 : 2) * metric[a][b] * d.values[i];
        }
    result.name = "Laplacian";
    result.unit = "(" + grid.unit + ")/A^2";
    result.validate();
    return result;
}

std::array<Grid, 3> energyDensity(const Grid& density, double floor)
{
    if (density.unit != "e/A^3") throw std::invalid_argument("Energy conversion requires electron density in e/A^3");
    if (!std::isfinite(floor) || floor <= 0) throw std::invalid_argument("Density floor must be positive and finite");
    const auto grad = gradient(density);
    const auto lap = laplacian(density);
    std::array<Grid, 3> result{density, density, density};
    for (std::size_t i = 0; i < density.values.size(); ++i)
    {
        if (density.values[i] < 0) throw std::invalid_argument("Energy conversion requires nonnegative density, not spin/difference density");
        const double rho = density.values[i] * std::pow(bohr, 3);
        if (density.values[i] <= floor)
        {
            for (auto& f : result) f.values[i] = 0;
            continue;
        }
        double norm = 0;
        for (const auto& f : grad) norm += f.values[i] * f.values[i] * std::pow(bohr, 8);
        const double l = lap.values[i] * std::pow(bohr, 5);
        const double g = 0.3 * std::pow(3 * pi * pi, 2.0 / 3) * std::pow(rho, 5.0 / 3) + norm / (72 * rho) + l / 12;
        const double v = l / 4 - 2 * g;
        const double factor = hartree / std::pow(bohr, 3);
        result[0].values[i] = g * factor;
        result[1].values[i] = v * factor;
        result[2].values[i] = (g + v) * factor;
    }
    const char* names[] = {"kinetic energy density (gradient expansion)", "potential energy density (local virial)", "total energy density (local virial)"};
    for (int i = 0; i < 3; ++i) { result[i].name = names[i]; result[i].unit = "eV/A^3"; result[i].validate(); }
    return result;
}

namespace
{
void validateDensity(const Grid& density, double floor)
{
    if (density.unit != "e/A^3") throw std::invalid_argument("This descriptor requires electron density in e/A^3");
    if (!std::isfinite(floor) || floor <= 0) throw std::invalid_argument("Density floor must be positive and finite");
    for (double v : density.values)
        if (v < 0) throw std::invalid_argument("This descriptor requires nonnegative density, not spin/difference density");
}

// Middle eigenvalue (lambda_2) of a real symmetric 3x3 matrix via the closed-form
// trigonometric solution (Smith, Comm. ACM 4, 168, 1961), avoiding an iterative solver.
double middleEigenvalue(double a11, double a22, double a33, double a12, double a13, double a23)
{
    const double p1 = a12 * a12 + a13 * a13 + a23 * a23;
    if (p1 < 1e-300) { const double lo = std::min({a11, a22, a33}), hi = std::max({a11, a22, a33}); return a11 + a22 + a33 - lo - hi; }
    const double q = (a11 + a22 + a33) / 3;
    const double p = std::sqrt(((a11 - q) * (a11 - q) + (a22 - q) * (a22 - q) + (a33 - q) * (a33 - q) + 2 * p1) / 6);
    const double b11 = (a11 - q) / p, b22 = (a22 - q) / p, b33 = (a33 - q) / p, b12 = a12 / p, b13 = a13 / p, b23 = a23 / p;
    double r = (b11 * (b22 * b33 - b23 * b23) - b12 * (b12 * b33 - b23 * b13) + b13 * (b12 * b23 - b22 * b13)) / 2;
    r = std::clamp(r, -1.0, 1.0);
    const double phi = std::acos(r) / 3;
    const double eig1 = q + 2 * p * std::cos(phi), eig3 = q + 2 * p * std::cos(phi + 2 * pi / 3);
    return 3 * q - eig1 - eig3;
}

struct DensityHessian { std::array<Grid, 3> grad; std::array<Grid, 6> h; };  // grad: d(rho)/dx,dy,dz; h: xx,yy,zz,xy,xz,yz

// Second Cartesian derivatives from differentiating the engine's own metric-aware
// gradient a second time; mixed partials are symmetrized to cancel first-order
// finite-difference asymmetry. This is a finite-difference Hessian, not an
// analytic one: exact for quadratic fields, approximate otherwise.
DensityHessian densityHessian(const Grid& density)
{
    auto grad = gradient(density);
    const auto gx = gradient(grad[0]), gy = gradient(grad[1]), gz = gradient(grad[2]);
    Grid hxy = gx[1], hxz = gx[2], hyz = gy[2];
    for (std::size_t i = 0; i < hxy.values.size(); ++i)
    {
        hxy.values[i] = (gx[1].values[i] + gy[0].values[i]) / 2;
        hxz.values[i] = (gx[2].values[i] + gz[0].values[i]) / 2;
        hyz.values[i] = (gy[2].values[i] + gz[1].values[i]) / 2;
    }
    return {grad, {gx[0], gy[1], gz[2], hxy, hxz, hyz}};
}
}

Grid reducedDensityGradient(const Grid& density, double floor)
{
    validateDensity(density, floor);
    const auto grad = gradient(density);
    Grid result = density;
    const double constant = 2 * std::pow(3 * pi * pi, 1.0 / 3);
    for (std::size_t i = 0; i < density.values.size(); ++i)
    {
        const double rho = density.values[i];
        if (rho <= floor) { result.values[i] = 0; continue; }
        const double gradNorm = std::sqrt(grad[0].values[i] * grad[0].values[i] + grad[1].values[i] * grad[1].values[i]
                                          + grad[2].values[i] * grad[2].values[i]);
        result.values[i] = gradNorm / (constant * std::pow(rho, 4.0 / 3));
    }
    result.name = "reduced density gradient"; result.unit = "dimensionless"; result.validate();
    return result;
}

Grid signedDensity(const Grid& density, double floor)
{
    validateDensity(density, floor);
    const auto hessian = densityHessian(density);
    Grid result = density;
    for (std::size_t i = 0; i < density.values.size(); ++i)
    {
        const double rho = density.values[i];
        if (rho <= floor) { result.values[i] = 0; continue; }
        const double lambda2 = middleEigenvalue(hessian.h[0].values[i], hessian.h[1].values[i], hessian.h[2].values[i],
                                                 hessian.h[3].values[i], hessian.h[4].values[i], hessian.h[5].values[i]);
        result.values[i] = lambda2 == 0 ? 0.0 : std::copysign(rho, lambda2);
    }
    result.name = "sign(lambda_2) * density"; result.unit = "e/A^3"; result.validate();
    return result;
}

Grid dori(const Grid& density, double floor)
{
    validateDensity(density, floor);
    const auto hessian = densityHessian(density);
    Grid result = density;
    for (std::size_t i = 0; i < density.values.size(); ++i)
    {
        const double rho = density.values[i];
        if (rho <= floor) { result.values[i] = 0; continue; }
        const double g[3] = {hessian.grad[0].values[i], hessian.grad[1].values[i], hessian.grad[2].values[i]};
        const double h[3][3] = {
            {hessian.h[0].values[i], hessian.h[3].values[i], hessian.h[4].values[i]},
            {hessian.h[3].values[i], hessian.h[1].values[i], hessian.h[5].values[i]},
            {hessian.h[4].values[i], hessian.h[5].values[i], hessian.h[2].values[i]},
        };
        const double k = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) / (rho * rho);
        double gradGammaSquared = 0;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
            {
                const double jacobian = (h[a][b] * rho - g[a] * g[b]) / (rho * rho);
                gradGammaSquared += jacobian * jacobian;
            }
        const double theta = gradGammaSquared / (k * k * k);
        result.values[i] = theta / (1 + theta);
    }
    result.name = "DORI"; result.unit = "dimensionless"; result.validate();
    return result;
}

Grid smooth(const Grid& grid, double sigma, int radius)
{
    grid.validate();
    if (!std::isfinite(sigma) || sigma <= 0 || radius < 1 || radius > 50)
        throw std::invalid_argument("Smoothing needs sigma > 0 Angstrom and radius 1..50 grid steps");
    Grid out = grid;
    std::vector<std::pair<glm::ivec3, double>> kernel;
    for (int z = -radius; z <= radius; ++z)
        for (int y = -radius; y <= radius; ++y)
            for (int x = -radius; x <= radius; ++x)
            {
                const auto delta = grid.steps() * glm::dvec3(x, y, z);
                kernel.push_back({{x, y, z}, std::exp(-glm::dot(delta, delta) / (2 * sigma * sigma))});
            }
    if (static_cast<double>(kernel.size()) * grid.values.size() > 2e9)
        throw std::invalid_argument("Smoothing workload too large; reduce radius or grid resolution");
    for (int z = 0; z < grid.shape[2]; ++z)
        for (int y = 0; y < grid.shape[1]; ++y)
            for (int x = 0; x < grid.shape[0]; ++x)
            {
                if (x==0 && y==0) atomforge::taskProgress(double(z)/grid.shape[2]);
                double sum = 0, weights = 0;
                for (const auto& k : kernel)
                {
                    auto p = glm::ivec3(x, y, z) + k.first;
                    bool valid = true;
                    for (int a = 0; a < 3; ++a)
                        if (grid.periodic) p[a] = (p[a] % grid.shape[a] + grid.shape[a]) % grid.shape[a];
                        else if (p[a] < 0 || p[a] >= grid.shape[a]) valid = false;
                    if (valid) { sum += k.second * grid.values[grid.index(p.x, p.y, p.z)]; weights += k.second; }
                }
                out.values[grid.index(x, y, z)] = sum / weights;
            }
    out.name = "Gaussian smoothed " + grid.name;
    return out;
}

std::vector<glm::dvec2> lineProfile(const Grid& grid, glm::dvec3 start, glm::dvec3 end, int count)
{
    grid.validate();
    if (count < 2 || count > 1000000) throw std::invalid_argument("Line profile count must be 2..1000000");
    std::vector<glm::dvec2> out;
    const double length = glm::length(end - start);
    for (int i = 0; i < count; ++i)
    {
        const double t = static_cast<double>(i) / (count - 1);
        out.push_back({length * t, grid.sample(start + t * (end - start))});
    }
    return out;
}

std::vector<glm::dvec2> planarAverage(const Grid& grid, int axis)
{
    grid.validate();
    if (axis < 0 || axis > 2) throw std::invalid_argument("Axis must be 0, 1 or 2");
    std::vector<glm::dvec2> out(grid.shape[axis], glm::dvec2(0));
    std::vector<double> weights(out.size(), 0);
    const auto inverse = glm::inverse(grid.cell);
    const double height = 1.0 / glm::length(glm::dvec3(inverse[0][axis], inverse[1][axis], inverse[2][axis]));
    for (int z = 0; z < grid.shape[2]; ++z)
        for (int y = 0; y < grid.shape[1]; ++y)
            for (int x = 0; x < grid.shape[0]; ++x)
            {
                if (x==0 && y==0) atomforge::taskProgress(double(z)/grid.shape[2]);
                const int i = glm::ivec3(x, y, z)[axis];
                const double w = grid.weight(x, y, z);
                out[i].y += grid.values[grid.index(x, y, z)] * w;
                weights[i] += w;
            }
    for (std::size_t i = 0; i < out.size(); ++i)
    {
        out[i].x = height * i / (grid.shape[axis] - (grid.periodic ? 0 : 1));
        out[i].y /= weights[i];
    }
    return out;
}

std::vector<glm::dvec2> macroscopicAverage(const Grid& grid, int axis, int window)
{
    auto out = planarAverage(grid, axis);
    const auto raw = out;
    if (window < 1 || window % 2 == 0 || window > static_cast<int>(out.size()))
        throw std::invalid_argument("Macroscopic window must be odd and no longer than the profile");
    for (int i = 0; i < static_cast<int>(out.size()); ++i)
    {
        double sum = 0;
        int count = 0;
        for (int d = -window / 2; d <= window / 2; ++d)
        {
            int j = i + d;
            if (grid.periodic) j = (j + static_cast<int>(out.size())) % static_cast<int>(out.size());
            if (j >= 0 && j < static_cast<int>(out.size())) { sum += raw[j].y; ++count; }
        }
        out[i].y = sum / count;
    }
    return out;
}

Grid section(const Grid& grid, glm::dvec3 origin, glm::dvec3 u, glm::dvec3 v, int nu, int nv)
{
    grid.validate();
    if (nu < 2 || nv < 2 || static_cast<double>(nu) * nv > 10000000 || glm::length(glm::cross(u, v)) < 1e-12)
        throw std::invalid_argument("Section requires independent spanning vectors and valid dimensions");
    Grid out;
    out.shape = {nu, nv, 2};
    out.origin = origin;
    out.cell = glm::dmat3(u, v, glm::normalize(glm::cross(u, v)));
    out.name = "section";
    out.unit = grid.unit;
    out.values.resize(static_cast<std::size_t>(nu) * nv * 2);
    // Both planes hold the same samples; z=0 is the actual section surface.
    for (int y = 0; y < nv; ++y)
        for (int x = 0; x < nu; ++x)
        {
            const double value = grid.sample(origin + u * (static_cast<double>(x) / (nu - 1)) + v * (static_cast<double>(y) / (nv - 1)));
            out.values[out.index(x, y, 0)] = out.values[out.index(x, y, 1)] = value;
        }
    out.validate();
    return out;
}

std::vector<glm::dvec3> contours(const Grid& slice, double level)
{
    slice.validate();
    if (!std::isfinite(level)) throw std::invalid_argument("Contour level must be finite");
    std::vector<glm::dvec3> out;
    // Consistent triangulation removes the ambiguous marching-square saddle.
    const int triangles[2][3] = {{0, 1, 2}, {0, 2, 3}};
    for (int y = 0; y < slice.shape[1] - 1; ++y)
        for (int x = 0; x < slice.shape[0] - 1; ++x)
        {
            const glm::ivec2 p[4] = {{x, y}, {x + 1, y}, {x + 1, y + 1}, {x, y + 1}};
            for (const auto& tri : triangles)
            {
                std::vector<glm::dvec3> hits;
                for (int e = 0; e < 3; ++e)
                {
                    const auto a = p[tri[e]], b = p[tri[(e + 1) % 3]];
                    const double va = slice.values[slice.index(a.x, a.y, 0)], vb = slice.values[slice.index(b.x, b.y, 0)];
                    if ((va < level) != (vb < level)) hits.push_back(glm::mix(slice.position(a.x, a.y, 0), slice.position(b.x, b.y, 0), (level - va) / (vb - va)));
                }
                if (hits.size() == 2) out.insert(out.end(), hits.begin(), hits.end());
            }
        }
    return out;
}
}
