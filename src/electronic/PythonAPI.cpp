#include "electronic/Volume.h"

#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>

#ifdef _WIN32
#define AF_EXPORT extern "C" __declspec(dllexport)
#else
#define AF_EXPORT extern "C" __attribute__((visibility("default")))
#endif

// Versioned C boundary: exceptions never cross ctypes. Each returned handle is
// independently owned; getters expose borrowed data only until handle deletion.
namespace
{
using namespace atomforge::electronic;
struct Result { Volume volume; std::vector<double> table; int columns = 0; };
thread_local std::string error;
template<class Work> void* guard(Work work)
{
    error.clear();
    try { return work(); }
    catch (const std::exception& e) { error = e.what(); }
    catch (...) { error = "Unknown native post-processing error"; }
    return nullptr;
}
Result& result(void* h)
{
    if (!h) throw std::invalid_argument("Null result handle");
    return *static_cast<Result*>(h);
}
const Grid& grid(void* h)
{
    const auto& r = result(h);
    if (r.volume.fields.empty()) throw std::invalid_argument("Operation requires a grid");
    return r.volume.fields[0];
}
int exactInt(double value)
{
    if (!std::isfinite(value) || value != std::floor(value) || std::abs(value) > 100000000)
        throw std::invalid_argument("Expected a bounded integer parameter");
    return static_cast<int>(value);
}
}

AF_EXPORT int af_electronic_version() { return 1; }
AF_EXPORT const char* af_electronic_error() { return error.c_str(); }
AF_EXPORT void af_electronic_free(void* h) { delete static_cast<Result*>(h); }
AF_EXPORT void* af_electronic_load(const char* path, const char* quantity, const char* coordinateUnit)
{
    return guard([&]() -> void*
    {
        if (!path || !quantity || !coordinateUnit) throw std::invalid_argument("Null load argument");
        auto out = std::make_unique<Result>();
        out->volume = loadVolume(path, quantity, coordinateUnit);
        return out.release();
    });
}
AF_EXPORT void* af_electronic_create(const int* shape, const double* cell, const double* origin,
    int periodic, const double* values, std::size_t count, const char* unit)
{
    return guard([&]() -> void*
    {
        if (!shape || !cell || !origin || !values || !unit || count > 100000000) throw std::invalid_argument("Invalid grid buffers");
        Grid g;
        for (int a = 0; a < 3; ++a)
        {
            g.shape[a] = shape[a]; g.origin[a] = origin[a];
            for (int b = 0; b < 3; ++b) g.cell[a][b] = cell[a * 3 + b];
        }
        g.periodic = periodic != 0;
        g.unit = unit;
        g.values.assign(values, values + count);
        g.validate();
        auto out = std::make_unique<Result>();
        out->volume.fields.push_back(std::move(g));
        return out.release();
    });
}
AF_EXPORT void* af_electronic_select(void* h, int index)
{
    return guard([&]() -> void*
    {
        auto out = std::make_unique<Result>();
        out->volume.sites = result(h).volume.sites;
        out->volume.fields.push_back(result(h).volume.fields.at(index));
        return out.release();
    });
}
AF_EXPORT int af_electronic_fields(void* h) { return h ? static_cast<int>(static_cast<Result*>(h)->volume.fields.size()) : 0; }
AF_EXPORT const char* af_electronic_label(void* h, int field, int unit)
{
    if (!h || field < 0 || field >= af_electronic_fields(h)) return "";
    const auto& g = static_cast<Result*>(h)->volume.fields[field];
    return unit ? g.unit.c_str() : g.name.c_str();
}
AF_EXPORT std::size_t af_electronic_data(void* h, double* destination)
{
    if (!h) return 0;
    const auto& r = *static_cast<Result*>(h);
    const auto& data = r.volume.fields.empty() ? r.table : r.volume.fields[0].values;
    if (destination) std::copy(data.begin(), data.end(), destination);
    return data.size();
}
AF_EXPORT int af_electronic_geometry(void* h, int* shape, double* cell, double* origin)
{
    if (!h || static_cast<Result*>(h)->volume.fields.empty()) return -1;
    const auto& g = static_cast<Result*>(h)->volume.fields[0];
    for (int a = 0; a < 3; ++a)
    {
        if (shape) shape[a] = g.shape[a];
        if (origin) origin[a] = g.origin[a];
        for (int b = 0; b < 3; ++b) if (cell) cell[a * 3 + b] = g.cell[a][b];
    }
    return g.periodic ? 1 : 0;
}
AF_EXPORT int af_electronic_columns(void* h) { return h ? static_cast<Result*>(h)->columns : 0; }
AF_EXPORT std::size_t af_electronic_sites(void* h, double* destination)
{
    if (!h) return 0;
    const auto& sites = static_cast<Result*>(h)->volume.sites;
    if (destination)
        for (std::size_t i = 0; i < sites.size(); ++i)
        {
            destination[i * 4] = sites[i].number;
            for (int a = 0; a < 3; ++a) destination[i * 4 + a + 1] = sites[i].position[a];
        }
    return sites.size();
}
AF_EXPORT void* af_electronic_with_sites(void* h, const double* sites, std::size_t count)
{
    return guard([&]() -> void*
    {
        if (count > 10000000 || (count && !sites)) throw std::invalid_argument("Invalid site buffer");
        auto out = std::make_unique<Result>(result(h));
        out->volume.sites.clear();
        for (std::size_t i = 0; i < count; ++i)
        {
            const int z = exactInt(sites[i * 4]);
            if (z < 0 || z > 118) throw std::invalid_argument("Atomic number outside 0..118");
            const glm::dvec3 p(sites[i*4+1],sites[i*4+2],sites[i*4+3]);
            for (int a = 0; a < 3; ++a) if (!std::isfinite(p[a])) throw std::invalid_argument("Non-finite site");
            out->volume.sites.push_back({z, p});
        }
        return out.release();
    });
}
AF_EXPORT int af_electronic_save(void* h, const char* path, const char* format)
{
    return guard([&]() -> void*
    {
        if (!path || !format) throw std::invalid_argument("Null output argument");
        const auto& r = result(h);
        if (std::string(format) == "obj")
        {
            if (r.columns != 4) throw std::invalid_argument("OBJ export requires an isosurface result");
            Mesh m;
            for (std::size_t i = 0; i < r.table.size(); i += 4) { m.vertices.push_back({r.table[i],r.table[i+1],r.table[i+2]}); m.colors.push_back(r.table[i+3]); }
            saveMesh(m, path);
        }
        else saveVolume(r.volume, path, format);
        return h;
    }) ? 1 : 0;
}

AF_EXPORT void* af_electronic_calculate(void* h, void* other, const char* operation, const double* parameters, std::size_t count)
{
    return guard([&]() -> void*
    {
        if (!operation || (count && !parameters)) throw std::invalid_argument("Null calculation arguments");
        for (std::size_t i = 0; i < count; ++i) if (!std::isfinite(parameters[i])) throw std::invalid_argument("Non-finite calculation parameter");
        const auto& g = grid(h);
        const std::string op(operation);
        auto out = std::make_unique<Result>();
        out->volume.sites = result(h).volume.sites;
        auto need = [&](std::size_t n) { if (count != n) throw std::invalid_argument("Incorrect parameter count for " + op); };
        auto p = [&](std::size_t i) { if (i >= count) throw std::invalid_argument("Missing parameter for " + op); return parameters[i]; };
        auto vec = [&](std::size_t i) { return glm::dvec3(p(i),p(i+1),p(i+2)); };
        auto field = [&](Grid f) { out->volume.fields.push_back(std::move(f)); };
        auto profile = [&](const std::vector<glm::dvec2>& rows)
        { out->columns = 2; for (auto row : rows) { out->table.push_back(row.x); out->table.push_back(row.y); } };
        if (op == "add" || op == "subtract" || op == "multiply" || op == "divide") { need(0); field(arithmetic(g,grid(other),op)); }
        else if (op == "scale") { need(1); field(scale(g,p(0))); }
        else if (op == "resample") { need(0); field(resample(g,grid(other))); }
        else if (op == "periodic") { need(1); field(periodicEndpoints(g,p(0))); }
        else if (op == "smooth") { need(2); field(smooth(g,p(0),exactInt(p(1)))); }
        else if (op == "gradient") { need(0); for (auto& f : gradient(g)) field(std::move(f)); }
        else if (op == "laplacian") { need(0); field(laplacian(g)); }
        else if (op == "energy") { need(1); for (auto& f : energyDensity(g,p(0))) field(std::move(f)); }
        else if (op == "line") { need(7); profile(lineProfile(g,vec(0),vec(3),exactInt(p(6)))); }
        else if (op == "planar") { need(1); profile(planarAverage(g,exactInt(p(0)))); }
        else if (op == "macro") { need(2); profile(macroscopicAverage(g,exactInt(p(0)),exactInt(p(1)))); }
        else if (op == "section") { need(11); field(section(g,vec(0),vec(3),vec(6),exactInt(p(9)),exactInt(p(10)))); }
        else if (op == "contours")
        {
            need(1); out->columns = 3;
            for (auto v : contours(g,p(0))) for (int a = 0; a < 3; ++a) out->table.push_back(v[a]);
        }
        else if (op == "integrate") { need(0); out->columns = 1; out->table = {integrate(g)}; }
        else if (op == "sphere") { need(4); out->columns = 1; out->table = {integrateSphere(g,vec(0),p(3))}; }
        else if (op == "voronoi")
        {
            if (count % 3) throw std::invalid_argument("Expected XYZ sites");
            std::vector<glm::dvec3> sites;
            for (std::size_t i = 0; i < count; i += 3) sites.push_back(vec(i));
            out->columns = 2;
            for (auto b : voronoiIntegrate(g,sites)) { out->table.push_back(b.integral); out->table.push_back(b.volume); }
        }
        else if (op == "peaks")
        {
            need(1); out->columns = 4;
            for (auto peak : peaks(g,exactInt(p(0)))) { for (int a = 0; a < 3; ++a) out->table.push_back(peak.position[a]); out->table.push_back(peak.value); }
        }
        else if (op == "patterson") { need(0); field(patterson(g)); }
        else if (op == "factors")
        {
            if (count % 3) throw std::invalid_argument("Expected hkl triples");
            std::vector<glm::ivec3> indices;
            for (std::size_t i = 0; i < count; i += 3) indices.push_back({exactInt(p(i)),exactInt(p(i+1)),exactInt(p(i+2))});
            out->columns = 5;
            for (const auto& r : structureFactors(g,indices))
                out->table.insert(out->table.end(),{static_cast<double>(r.hkl.x),static_cast<double>(r.hkl.y),static_cast<double>(r.hkl.z),r.value.real(),r.value.imag()});
        }
        else if (op == "synthesis")
        {
            if (count % 5) throw std::invalid_argument("Expected h k l real imag rows");
            std::vector<Reflection> reflections;
            for (std::size_t i = 0; i < count; i += 5) reflections.push_back({{exactInt(p(i)),exactInt(p(i+1)),exactInt(p(i+2))},{p(i+3),p(i+4)}});
            field(fourierSynthesis(g,reflections));
        }
        else if (op == "atomic_factors")
        {
            const auto n = result(h).volume.sites.size();
            if (!n || count < n * 2 || (count - n * 2) % (3 + n * 2)) throw std::invalid_argument("Invalid atomic-factor parameter layout");
            std::vector<glm::dvec3> positions;
            std::vector<double> occupancy, b;
            std::vector<glm::ivec3> indices;
            std::vector<std::complex<double>> factors;
            for (std::size_t i = 0; i < n; ++i) { positions.push_back(result(h).volume.sites[i].position); occupancy.push_back(p(i*2)); b.push_back(p(i*2+1)); }
            for (std::size_t i = n * 2; i < count; i += 3 + n * 2)
            {
                indices.push_back({exactInt(p(i)),exactInt(p(i+1)),exactInt(p(i+2))});
                for (std::size_t j = 0; j < n; ++j) factors.push_back({p(i+3+j*2),p(i+4+j*2)});
            }
            out->columns = 5;
            for (const auto& r : atomicStructureFactors(g,positions,occupancy,b,indices,factors))
                out->table.insert(out->table.end(),{static_cast<double>(r.hkl.x),static_cast<double>(r.hkl.y),static_cast<double>(r.hkl.z),r.value.real(),r.value.imag()});
        }
        else if (op == "ewald")
        {
            need(3 + result(h).volume.sites.size());
            std::vector<glm::dvec3> sites;
            std::vector<double> charges;
            for (std::size_t i = 0; i < result(h).volume.sites.size(); ++i) { sites.push_back(result(h).volume.sites[i].position); charges.push_back(p(i+3)); }
            const auto e = ewald(g.cell,sites,charges,p(0),p(1),p(2));
            out->columns = 1; out->table = {e.energy};
            out->table.insert(out->table.end(),e.potentials.begin(),e.potentials.end());
        }
        else if (op == "isosurface")
        {
            need(1); out->columns = 4;
            const auto mesh = isosurface(g,p(0),other ? &grid(other) : nullptr);
            for (std::size_t i = 0; i < mesh.vertices.size(); ++i)
            { const auto v = mesh.vertices[i]; out->table.insert(out->table.end(),{v.x,v.y,v.z,mesh.colors[i]}); }
        }
        else throw std::invalid_argument("Unknown electronic operation: " + op);
        return out.release();
    });
}
