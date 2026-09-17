#include "electronic/BaderPartition.h"
#include "electronic/Volume.h"

#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#define AF_EXPORT extern "C" __declspec(dllexport)
#else
#define AF_EXPORT extern "C" __attribute__((visibility("default")))
#endif

// Versioned C boundary: exceptions never cross ctypes. Mirrors the handle/
// guard convention in PythonAPI.cpp/WannierAPI.cpp/LobsterAPI.cpp; takes raw
// grid buffers (as af_electronic_create does) rather than another file's
// private Result handle, so this stays independent of PythonAPI.cpp.
namespace
{
using namespace atomforge::electronic;
struct BaderHandle { BaderPartition partition; glm::dmat3 cell{1.0}; bool periodic = false; };
thread_local std::string baderError;
template<class Work> bool guard(Work work)
{
    baderError.clear();
    try { work(); return true; }
    catch (const std::exception& e) { baderError = e.what(); }
    catch (...) { baderError = "Unknown native Bader error"; }
    return false;
}
BaderHandle& handle(void* h)
{
    if (!h) throw std::invalid_argument("Null Bader partition handle");
    return *static_cast<BaderHandle*>(h);
}
}

AF_EXPORT const char* af_bader_error() { return baderError.c_str(); }
AF_EXPORT void af_bader_free(void* h) { delete static_cast<BaderHandle*>(h); }

AF_EXPORT void* af_bader_partition(const int* shape, const double* cell, const double* origin,
    int periodic, const double* values, std::size_t count, const char* unit)
{
    void* out = nullptr;
    guard([&]
    {
        if (!shape || !cell || !origin || !values || !unit) throw std::invalid_argument("Invalid grid buffers");
        Grid g;
        for (int a = 0; a < 3; ++a)
        {
            g.shape[a] = shape[a]; g.origin[a] = origin[a];
            for (int b = 0; b < 3; ++b) g.cell[a][b] = cell[a * 3 + b];
        }
        g.periodic = periodic != 0;
        g.unit = unit;
        g.values.assign(values, values + count);
        auto result = std::make_unique<BaderHandle>();
        result->cell = g.cell;
        result->periodic = g.periodic;
        result->partition = baderOnGrid(g);
        out = result.release();
    });
    return out;
}
AF_EXPORT int af_bader_num_basins(void* h) { return h ? static_cast<int>(static_cast<BaderHandle*>(h)->partition.basins.size()) : 0; }
// basinIds: count doubles out (one per voxel, x-fastest order, matching the source grid).
AF_EXPORT int af_bader_basin_ids(void* h, double* basinIds)
{
    return guard([&]
    {
        if (!basinIds) throw std::invalid_argument("Null buffer");
        const auto& ids = handle(h).partition.basin;
        for (std::size_t i = 0; i < ids.size(); ++i) basinIds[i] = ids[i];
    });
}
AF_EXPORT int af_bader_basin(void* h, int index, double* charge, double* volume, double* maximum)
{
    return guard([&]
    {
        if (!charge || !volume || !maximum) throw std::invalid_argument("Null buffer");
        const auto& basin = handle(h).partition.basins.at(index);
        *charge = basin.charge; *volume = basin.volume;
        for (int a = 0; a < 3; ++a) maximum[a] = basin.maximum[a];
    });
}
// sites: siteCount*3 Cartesian doubles; populations: siteCount doubles out.
AF_EXPORT int af_bader_populations(void* h, const double* sites, std::size_t siteCount, double* populations)
{
    return guard([&]
    {
        if (!sites || !populations || siteCount == 0) throw std::invalid_argument("Invalid site buffer");
        auto& info = handle(h);
        std::vector<glm::dvec3> points(siteCount);
        for (std::size_t i = 0; i < siteCount; ++i) points[i] = {sites[i * 3], sites[i * 3 + 1], sites[i * 3 + 2]};
        for (std::size_t i = 0; i < siteCount; ++i) populations[i] = 0.0;
        for (const auto& basin : info.partition.basins)
        {
            std::size_t nearest = 0;
            double bestDistance2 = std::numeric_limits<double>::infinity();
            for (std::size_t s = 0; s < points.size(); ++s)
            {
                glm::dvec3 delta = basin.maximum - points[s];
                if (info.periodic) delta = minimumImage(delta, info.cell);
                const double d2 = glm::dot(delta, delta);
                if (d2 < bestDistance2) { bestDistance2 = d2; nearest = s; }
            }
            populations[nearest] += basin.charge;
        }
    });
}
