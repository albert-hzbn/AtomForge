#pragma once

#include "science/Json.h"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <string>
#include <vector>

namespace atomforge::science
{
using Vec3 = std::array<double, 3>;
// Row-vector cell convention: cell[i] is the i-th lattice vector (Angstrom).
using Mat3 = std::array<Vec3, 3>;
using Pbc = std::array<bool, 3>;

// Dense row-major numerical array with an explicit shape.
struct NdArray
{
    std::vector<std::size_t> shape;
    std::vector<double> values;

    NdArray() = default;
    explicit NdArray(std::vector<std::size_t> dimensions, double fill = 0.0);

    std::size_t ndim() const { return shape.size(); }
    std::size_t size() const { return values.size(); }
    double& operator()(std::size_t i) { return values[i]; }
    double operator()(std::size_t i) const { return values[i]; }
    double& operator()(std::size_t i, std::size_t j) { return values[i * shape[1] + j]; }
    double operator()(std::size_t i, std::size_t j) const { return values[i * shape[1] + j]; }
    double& operator()(std::size_t i, std::size_t j, std::size_t k) { return values[(i * shape[1] + j) * shape[2] + k]; }
    double operator()(std::size_t i, std::size_t j, std::size_t k) const { return values[(i * shape[1] + j) * shape[2] + k]; }
};

// Converts a rectangular numeric JSON array (or a scalar, giving shape {}) to an array.
NdArray toArray(const Json& value, const std::string& name);
Json toJson(const NdArray& array);
Json toJson(const std::vector<double>& values);
Json toJson(const std::vector<int>& values);
Json toJson(const std::vector<bool>& values);
Json toJson(const Vec3& value);
Json toJson(const Mat3& value);

// Validation helpers matching the documented tool contracts.
NdArray finiteArray(const NdArray& array, const std::string& name, int ndim = -1);
double positive(double value, const std::string& name);
long long integer(double value, const std::string& name, long long minimum = 1);
Mat3 cellMatrix(const NdArray& cell);
std::vector<Vec3> points(const NdArray& array, const std::string& name);

// Small 3-vector and row-vector-cell algebra.
inline Vec3 add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
inline Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline Vec3 scale(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
inline double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
inline Vec3 cross(const Vec3& a, const Vec3& b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double determinant(const Mat3& m);
Mat3 inverse(const Mat3& m);
Mat3 multiply(const Mat3& a, const Mat3& b);
Mat3 transpose(const Mat3& m);
// Row vector times matrix: v M.
Vec3 rowTimes(const Vec3& v, const Mat3& m);
// Cartesian -> fractional coordinates for a row-vector cell.
Vec3 fractional(const Vec3& cartesian, const Mat3& cell);
double cellVolume(const Mat3& cell);
Mat3 identity();

// General-cell minimum image along periodic axes (Gauss-reduced lattice search).
class MinimumImage
{
public:
    MinimumImage() = default;
    MinimumImage(const Mat3& cell, const Pbc& pbc);
    Vec3 operator()(const Vec3& vector) const;
    bool active() const { return m_active; }

private:
    bool m_active = false;
    Mat3 m_basis{};
    Mat3 m_inverse{};
    Pbc m_pbc{};
};

struct Neighbor
{
    int i = 0;
    int j = 0;
    std::array<int, 3> shift{};
    Vec3 vector{};
};

// All pairs i->j within cutoff, including periodic images (no self pair at zero shift).
// The returned vector is r_j + shift*cell - r_i for the supplied (unwrapped) positions.
std::vector<Neighbor> neighborList(const std::vector<Vec3>& positions, const Mat3& cell, const Pbc& pbc, double cutoff);

// Discrete Fourier transforms of arbitrary length (Bluestein for non powers of two).
std::vector<std::complex<double>> fft(std::vector<std::complex<double>> data, bool inverse = false);

// Parses a JSON boolean triple used for periodic boundaries.
Pbc parsePbc(const Json& value);
}
