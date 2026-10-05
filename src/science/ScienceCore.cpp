#include "science/ScienceCore.h"
#include "util/TaskControl.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>

namespace atomforge::science
{
NdArray::NdArray(std::vector<std::size_t> dimensions, double fill) : shape(std::move(dimensions))
{
    std::size_t count = 1;
    for (std::size_t extent : shape) count *= extent;
    values.assign(count, fill);
}

namespace
{
void flatten(const Json& value, std::size_t depth, const std::vector<std::size_t>& shape,
             std::vector<double>& out, const std::string& name)
{
    if (depth < shape.size()) {
        if (!value.isArray() || value.size() != shape[depth])
            throw std::runtime_error(name + " must be a rectangular numeric array");
        for (const auto& item : value.items()) flatten(item, depth + 1, shape, out, name);
        return;
    }
    if (value.isNumber()) out.push_back(value.number());
    else if (value.isBool()) out.push_back(value.boolean() ? 1.0 : 0.0);
    else if (value.isArray()) throw std::runtime_error(name + " must be a rectangular numeric array");
    else throw std::runtime_error(name + " must contain only numbers");
}
}

NdArray toArray(const Json& value, const std::string& name)
{
    NdArray result;
    for (const Json* current = &value; current->isArray(); current = &current->items()[0]) {
        result.shape.push_back(current->size());
        if (!current->size()) break;
    }
    flatten(value, 0, result.shape, result.values, name);
    std::size_t count = 1;
    for (std::size_t extent : result.shape) count *= extent;
    if (count != result.values.size()) throw std::runtime_error(name + " must be a rectangular numeric array");
    return result;
}

namespace
{
Json arrayJson(const NdArray& array, std::size_t axis, std::size_t offset)
{
    if (axis == array.shape.size()) return Json(array.values[offset]);
    std::size_t stride = 1;
    for (std::size_t k = axis + 1; k < array.shape.size(); ++k) stride *= array.shape[k];
    Json result = Json::array();
    result.items().reserve(array.shape[axis]);
    for (std::size_t i = 0; i < array.shape[axis]; ++i) result.push(arrayJson(array, axis + 1, offset + i * stride));
    return result;
}
}

Json toJson(const NdArray& array)
{
    if (array.shape.empty()) return array.values.empty() ? Json() : Json(array.values[0]);
    return arrayJson(array, 0, 0);
}

Json toJson(const std::vector<double>& values)
{
    Json result = Json::array();
    result.items().reserve(values.size());
    for (double value : values) result.push(value);
    return result;
}

Json toJson(const std::vector<int>& values)
{
    Json result = Json::array();
    for (int value : values) result.push(value);
    return result;
}

Json toJson(const std::vector<bool>& values)
{
    Json result = Json::array();
    for (bool value : values) result.push(Json(static_cast<bool>(value)));
    return result;
}

Json toJson(const Vec3& value)
{
    return Json::array({value[0], value[1], value[2]});
}

Json toJson(const Mat3& value)
{
    return Json::array({toJson(value[0]), toJson(value[1]), toJson(value[2])});
}

NdArray finiteArray(const NdArray& array, const std::string& name, int ndim)
{
    bool finite = array.size() > 0;
    for (double value : array.values) finite = finite && std::isfinite(value);
    if (!finite || (ndim >= 0 && array.ndim() != static_cast<std::size_t>(ndim)))
        throw std::runtime_error(name + " must be a nonempty finite array of the required dimension");
    return array;
}

double positive(double value, const std::string& name)
{
    if (!std::isfinite(value) || value <= 0) throw std::runtime_error(name + " must be finite and positive");
    return value;
}

long long integer(double value, const std::string& name, long long minimum)
{
    if (!std::isfinite(value) || value != std::floor(value) || value < static_cast<double>(minimum) || std::abs(value) > 9e15)
        throw std::runtime_error(name + " must be an integer >= " + std::to_string(minimum));
    return static_cast<long long>(value);
}

Mat3 cellMatrix(const NdArray& cell)
{
    const NdArray matrix = finiteArray(cell, "cell", 2);
    if (matrix.shape[0] != 3 || matrix.shape[1] != 3)
        throw std::runtime_error("cell must contain three independent, well-conditioned row vectors");
    Mat3 result{};
    Eigen::Matrix3d eigen;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) eigen(i, j) = result[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = matrix(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
    const Eigen::JacobiSVD<Eigen::Matrix3d> svd(eigen);
    const auto singular = svd.singularValues();
    if (std::abs(eigen.determinant()) < 1e-12 || singular(2) <= 0 || singular(0) / singular(2) > 1e10)
        throw std::runtime_error("cell must contain three independent, well-conditioned row vectors");
    return result;
}

std::vector<Vec3> points(const NdArray& array, const std::string& name)
{
    const NdArray checked = finiteArray(array, name, 2);
    if (checked.shape[1] != 3) throw std::runtime_error(name + " must have shape (atoms, 3)");
    std::vector<Vec3> result(checked.shape[0]);
    for (std::size_t i = 0; i < result.size(); ++i) result[i] = {checked(i, 0), checked(i, 1), checked(i, 2)};
    return result;
}

double determinant(const Mat3& m)
{
    return dot(m[0], cross(m[1], m[2]));
}

Mat3 inverse(const Mat3& m)
{
    const double det = determinant(m);
    if (det == 0.0) throw std::runtime_error("Singular 3x3 matrix");
    Mat3 result{};
    result[0] = {(m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det, (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det, (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det};
    result[1] = {(m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det, (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det, (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det};
    result[2] = {(m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det, (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det, (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det};
    return result;
}

Mat3 multiply(const Mat3& a, const Mat3& b)
{
    Mat3 result{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) result[i][j] += a[i][k] * b[k][j];
    return result;
}

Mat3 transpose(const Mat3& m)
{
    Mat3 result{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) result[i][j] = m[j][i];
    return result;
}

Vec3 rowTimes(const Vec3& v, const Mat3& m)
{
    return {v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0],
            v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1],
            v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2]};
}

Vec3 fractional(const Vec3& cartesian, const Mat3& cell)
{
    return rowTimes(cartesian, inverse(cell));
}

double cellVolume(const Mat3& cell)
{
    return std::abs(determinant(cell));
}

Mat3 identity()
{
    return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
}

MinimumImage::MinimumImage(const Mat3& cell, const Pbc& pbc) : m_pbc(pbc)
{
    m_active = pbc[0] || pbc[1] || pbc[2];
    if (!m_active) return;
    m_basis = cell;
    // Pairwise Gauss reduction of the periodic vectors keeps the lattice and
    // makes a one-shell image search sufficient for skewed cells.
    for (int pass = 0; pass < 100; ++pass) {
        bool changed = false;
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                if (a == b || !pbc[a] || !pbc[b]) continue;
                const double length = dot(m_basis[a], m_basis[a]);
                const double factor = std::round(dot(m_basis[a], m_basis[b]) / length);
                if (factor != 0.0 && dot(m_basis[b], m_basis[b]) > dot(sub(m_basis[b], scale(m_basis[a], factor)), sub(m_basis[b], scale(m_basis[a], factor))) + 1e-12) {
                    m_basis[b] = sub(m_basis[b], scale(m_basis[a], factor));
                    changed = true;
                }
            }
        if (!changed) break;
    }
    m_inverse = inverse(m_basis);
}

Vec3 MinimumImage::operator()(const Vec3& vector) const
{
    if (!m_active) return vector;
    Vec3 f = rowTimes(vector, m_inverse);
    for (int k = 0; k < 3; ++k)
        if (m_pbc[k]) f[k] -= std::round(f[k]);
    const Vec3 base = rowTimes(f, m_basis);
    Vec3 best = base;
    double bestNorm = dot(base, base);
    const int range[3] = {m_pbc[0] ? 1 : 0, m_pbc[1] ? 1 : 0, m_pbc[2] ? 1 : 0};
    for (int i = -range[0]; i <= range[0]; ++i)
        for (int j = -range[1]; j <= range[1]; ++j)
            for (int k = -range[2]; k <= range[2]; ++k) {
                if (!i && !j && !k) continue;
                const Vec3 candidate = add(base, rowTimes({double(i), double(j), double(k)}, m_basis));
                const double candidateNorm = dot(candidate, candidate);
                if (candidateNorm < bestNorm - 1e-12) { best = candidate; bestNorm = candidateNorm; }
            }
    return best;
}

std::vector<Neighbor> neighborList(const std::vector<Vec3>& positions, const Mat3& cell, const Pbc& pbc, double cutoff)
{
    std::vector<Neighbor> result;
    if (positions.empty()) return result;
    const bool periodic = pbc[0] || pbc[1] || pbc[2];
    const Mat3 inv = periodic ? inverse(cell) : identity();
    // Wrap periodic components so that one image shell bound is exact.
    std::vector<Vec3> wrapped(positions.size());
    std::vector<std::array<int, 3>> offsets(positions.size(), {0, 0, 0});
    for (std::size_t i = 0; i < positions.size(); ++i) {
        wrapped[i] = positions[i];
        if (!periodic) continue;
        const Vec3 f = rowTimes(positions[i], inv);
        Vec3 shift{0, 0, 0};
        for (int k = 0; k < 3; ++k)
            if (pbc[k]) { offsets[i][k] = static_cast<int>(std::floor(f[k])); shift[k] = offsets[i][k]; }
        wrapped[i] = sub(positions[i], rowTimes(shift, cell));
    }
    int range[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k)
        if (pbc[k]) {
            const Vec3 reciprocal{inv[0][k], inv[1][k], inv[2][k]};
            range[k] = static_cast<int>(std::ceil(cutoff * norm(reciprocal)));
            if (range[k] > 200) throw std::runtime_error("Neighbor cutoff is too large for the periodic cell");
        }
    Vec3 low = wrapped[0], high = wrapped[0];
    for (const auto& p : wrapped)
        for (int k = 0; k < 3; ++k) { low[k] = std::min(low[k], p[k]); high[k] = std::max(high[k], p[k]); }
    struct Image { int atom; std::array<int, 3> shift; Vec3 position; };
    std::vector<Image> images;
    for (int a = -range[0]; a <= range[0]; ++a)
        for (int b = -range[1]; b <= range[1]; ++b)
            for (int c = -range[2]; c <= range[2]; ++c) {
                const Vec3 translation = rowTimes({double(a), double(b), double(c)}, cell);
                for (std::size_t j = 0; j < wrapped.size(); ++j) {
                    const Vec3 p = add(wrapped[j], translation);
                    bool inside = true;
                    for (int k = 0; k < 3; ++k) inside = inside && p[k] >= low[k] - cutoff && p[k] <= high[k] + cutoff;
                    if (inside) images.push_back({static_cast<int>(j), {a, b, c}, p});
                }
                taskCheckpoint();
            }
    const double binSize = cutoff;
    auto key = [&](const Vec3& p, int dx, int dy, int dz) {
        const long long x = static_cast<long long>(std::floor((p[0] - low[0] + cutoff) / binSize)) + dx;
        const long long y = static_cast<long long>(std::floor((p[1] - low[1] + cutoff) / binSize)) + dy;
        const long long z = static_cast<long long>(std::floor((p[2] - low[2] + cutoff) / binSize)) + dz;
        return (x * 2000003LL + y) * 2000003LL + z;
    };
    std::unordered_map<long long, std::vector<int>> bins;
    for (std::size_t n = 0; n < images.size(); ++n) bins[key(images[n].position, 0, 0, 0)].push_back(static_cast<int>(n));
    const double cutoff2 = cutoff * cutoff;
    for (std::size_t i = 0; i < wrapped.size(); ++i) {
        if (i % 256 == 0) taskCheckpoint();
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    const auto found = bins.find(key(wrapped[i], dx, dy, dz));
                    if (found == bins.end()) continue;
                    for (int index : found->second) {
                        const Image& image = images[static_cast<std::size_t>(index)];
                        if (image.atom == static_cast<int>(i) && !image.shift[0] && !image.shift[1] && !image.shift[2]) continue;
                        const Vec3 d = sub(image.position, wrapped[i]);
                        if (dot(d, d) >= cutoff2) continue;
                        Neighbor neighbor;
                        neighbor.i = static_cast<int>(i);
                        neighbor.j = image.atom;
                        for (int k = 0; k < 3; ++k)
                            neighbor.shift[k] = image.shift[k] + offsets[static_cast<std::size_t>(image.atom)][k] - offsets[i][k];
                        neighbor.vector = d;
                        result.push_back(neighbor);
                    }
                }
    }
    return result;
}

namespace
{
using Complex = std::complex<double>;

void radix2(std::vector<Complex>& data, bool inverse)
{
    const std::size_t n = data.size();
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(data[i], data[j]);
    }
    const double pi = std::acos(-1.0);
    for (std::size_t length = 2; length <= n; length <<= 1) {
        const double angle = 2 * pi / static_cast<double>(length) * (inverse ? 1 : -1);
        for (std::size_t start = 0; start < n; start += length)
            for (std::size_t k = 0; k < length / 2; ++k) {
                const Complex w = std::polar(1.0, angle * static_cast<double>(k));
                const Complex u = data[start + k], v = data[start + k + length / 2] * w;
                data[start + k] = u + v;
                data[start + k + length / 2] = u - v;
            }
    }
}
}

std::vector<std::complex<double>> fft(std::vector<std::complex<double>> data, bool inverse)
{
    const std::size_t n = data.size();
    if (n <= 1) return data;
    if ((n & (n - 1)) == 0) {
        radix2(data, inverse);
        if (inverse) for (auto& value : data) value /= static_cast<double>(n);
        return data;
    }
    // Bluestein chirp-z transform.
    std::size_t m = 1;
    while (m < 2 * n - 1) m <<= 1;
    const double pi = std::acos(-1.0);
    std::vector<Complex> chirp(n);
    for (std::size_t k = 0; k < n; ++k) {
        const double phase = pi * static_cast<double>((k * k) % (2 * n)) / static_cast<double>(n) * (inverse ? 1 : -1);
        chirp[k] = std::polar(1.0, phase);
    }
    std::vector<Complex> a(m), b(m);
    for (std::size_t k = 0; k < n; ++k) a[k] = data[k] * chirp[k];
    b[0] = std::conj(chirp[0]);
    for (std::size_t k = 1; k < n; ++k) b[k] = b[m - k] = std::conj(chirp[k]);
    radix2(a, false);
    radix2(b, false);
    for (std::size_t k = 0; k < m; ++k) a[k] *= b[k];
    radix2(a, true);
    for (std::size_t k = 0; k < n; ++k) {
        data[k] = a[k] / static_cast<double>(m) * chirp[k];
        if (inverse) data[k] /= static_cast<double>(n);
    }
    return data;
}

Pbc parsePbc(const Json& value)
{
    if (!value.isArray() || value.size() != 3) throw std::runtime_error("pbc must contain three booleans");
    Pbc result{};
    for (std::size_t k = 0; k < 3; ++k) {
        if (!value.items()[k].isBool()) throw std::runtime_error("pbc must contain three booleans");
        result[k] = value.items()[k].boolean();
    }
    return result;
}
}
