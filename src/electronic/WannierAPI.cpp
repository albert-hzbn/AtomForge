#include "electronic/Wannier.h"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#define AF_EXPORT extern "C" __declspec(dllexport)
#else
#define AF_EXPORT extern "C" __attribute__((visibility("default")))
#endif

// Versioned C boundary: exceptions never cross ctypes. Mirrors the handle/guard
// convention in PythonAPI.cpp, kept independent since a WannierHamiltonian is
// not a Grid/Volume.
namespace
{
using namespace atomforge::electronic;
thread_local std::string wannierError;
template<class Work> bool guard(Work work)
{
    wannierError.clear();
    try { work(); return true; }
    catch (const std::exception& e) { wannierError = e.what(); }
    catch (...) { wannierError = "Unknown native Wannier error"; }
    return false;
}
WannierHamiltonian& model(void* h)
{
    if (!h) throw std::invalid_argument("Null Wannier model handle");
    return *static_cast<WannierHamiltonian*>(h);
}
}

AF_EXPORT const char* af_wannier_error() { return wannierError.c_str(); }
AF_EXPORT void af_wannier_free(void* h) { delete static_cast<WannierHamiltonian*>(h); }

AF_EXPORT void* af_wannier_load(const char* path)
{
    void* out = nullptr;
    guard([&]
    {
        if (!path) throw std::invalid_argument("Null path");
        out = new WannierHamiltonian(readWannierHr(path));
    });
    return out;
}

AF_EXPORT int af_wannier_num_wann(void* h) { return h ? static_cast<WannierHamiltonian*>(h)->numWann : 0; }

// matrix: num_wann*num_wann interleaved (real, imag) pairs out, row-major H(k).
AF_EXPORT int af_wannier_hamiltonian(void* h, const double* kpoint, double* matrix)
{
    return guard([&]
    {
        if (!kpoint || !matrix) throw std::invalid_argument("Null buffer");
        const auto values = blochHamiltonian(model(h), {kpoint[0], kpoint[1], kpoint[2]});
        for (std::size_t i = 0; i < values.size(); ++i) { matrix[2 * i] = values[i].real(); matrix[2 * i + 1] = values[i].imag(); }
    });
}

// kpoints: count*3 fractional-reciprocal doubles; energies: count*num_wann doubles out, ascending per k point.
AF_EXPORT int af_wannier_bands(void* h, const double* kpoints, std::size_t count, double* energies)
{
    return guard([&]
    {
        if (!kpoints || !energies) throw std::invalid_argument("Null buffer");
        std::vector<glm::dvec3> points(count);
        for (std::size_t i = 0; i < count; ++i) points[i] = {kpoints[i * 3], kpoints[i * 3 + 1], kpoints[i * 3 + 2]};
        const auto values = wannierBands(model(h), points);
        std::copy(values.begin(), values.end(), energies);
    });
}

// curvature: num_wann doubles out, ascending band order, Hamiltonian-gauge (Kubo) term.
AF_EXPORT int af_wannier_berry_curvature(void* h, const double* kpoint, int planeA, int planeB, double* curvature)
{
    return guard([&]
    {
        if (!kpoint || !curvature) throw std::invalid_argument("Null buffer");
        const auto values = berryCurvature(model(h), {kpoint[0], kpoint[1], kpoint[2]}, planeA, planeB);
        std::copy(values.begin(), values.end(), curvature);
    });
}

AF_EXPORT int af_wannier_chern_number(void* h, int band, int planeA, int planeB, int grid, double fixed, double* result)
{
    return guard([&]
    {
        if (!result) throw std::invalid_argument("Null buffer");
        *result = chernNumber(model(h), band, planeA, planeB, grid, fixed);
    });
}
