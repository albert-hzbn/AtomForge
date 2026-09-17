#include "electronic/Lobster.h"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>

#ifdef _WIN32
#define AF_EXPORT extern "C" __declspec(dllexport)
#else
#define AF_EXPORT extern "C" __attribute__((visibility("default")))
#endif

// Versioned C boundary: exceptions never cross ctypes. Mirrors the handle/guard
// convention in PythonAPI.cpp and WannierAPI.cpp.
namespace
{
using namespace atomforge::electronic;
thread_local std::string lobsterError;
template<class Work> bool guard(Work work)
{
    lobsterError.clear();
    try { work(); return true; }
    catch (const std::exception& e) { lobsterError = e.what(); }
    catch (...) { lobsterError = "Unknown native LOBSTER error"; }
    return false;
}
CohpData& cohpcar(void* h)
{
    if (!h) throw std::invalid_argument("Null COHPCAR handle");
    return *static_cast<CohpData*>(h);
}
IcohpList& icohplist(void* h)
{
    if (!h) throw std::invalid_argument("Null ICOHPLIST handle");
    return *static_cast<IcohpList*>(h);
}
}

AF_EXPORT const char* af_lobster_error() { return lobsterError.c_str(); }

AF_EXPORT void* af_cohpcar_load(const char* path)
{
    void* out = nullptr;
    guard([&] { if (!path) throw std::invalid_argument("Null path"); out = new CohpData(readCohpcar(path)); });
    return out;
}
AF_EXPORT void af_cohpcar_free(void* h) { delete static_cast<CohpData*>(h); }
// spinPolarized/numEnergies/numBonds out; returns 1/0. fermiEnergy via af_cohpcar_fermi_energy.
AF_EXPORT int af_cohpcar_info(void* h, int* spinPolarized, std::size_t* numEnergies, int* numBonds)
{
    return guard([&]
    {
        if (!spinPolarized || !numEnergies || !numBonds) throw std::invalid_argument("Null buffer");
        const auto& data = cohpcar(h);
        *spinPolarized = data.spinPolarized ? 1 : 0;
        *numEnergies = data.energies.size();
        *numBonds = static_cast<int>(data.bonds.size());
    });
}
AF_EXPORT double af_cohpcar_fermi_energy(void* h) { return h ? static_cast<CohpData*>(h)->fermiEnergy : 0.0; }
AF_EXPORT int af_cohpcar_energies(void* h, double* buffer)
{
    return guard([&]
    {
        if (!buffer) throw std::invalid_argument("Null buffer");
        const auto& energies = cohpcar(h).energies;
        std::copy(energies.begin(), energies.end(), buffer);
    });
}
AF_EXPORT int af_cohpcar_bond(void* h, int bond, int* atom1, int* atom2, double* length)
{
    return guard([&]
    {
        if (!atom1 || !atom2 || !length) throw std::invalid_argument("Null buffer");
        const auto& info = cohpcar(h).bonds.at(bond);
        *atom1 = info.atom1; *atom2 = info.atom2; *length = info.length;
    });
}
// cohp/icohp: numEnergies doubles out each. spin: 0 (up/only channel) or 1 (down).
AF_EXPORT int af_cohpcar_average(void* h, int spin, double* cohp, double* icohp)
{
    return guard([&]
    {
        if (!cohp || !icohp) throw std::invalid_argument("Null buffer");
        const auto& data = cohpcar(h);
        std::copy(data.averageCohp.at(spin).begin(), data.averageCohp.at(spin).end(), cohp);
        std::copy(data.averageIcohp.at(spin).begin(), data.averageIcohp.at(spin).end(), icohp);
    });
}
AF_EXPORT int af_cohpcar_bond_curve(void* h, int bond, int spin, double* cohp, double* icohp)
{
    return guard([&]
    {
        if (!cohp || !icohp) throw std::invalid_argument("Null buffer");
        const auto& data = cohpcar(h);
        std::copy(data.cohp.at(spin).at(bond).begin(), data.cohp.at(spin).at(bond).end(), cohp);
        std::copy(data.icohp.at(spin).at(bond).begin(), data.icohp.at(spin).at(bond).end(), icohp);
    });
}

AF_EXPORT void* af_icohplist_load(const char* path)
{
    void* out = nullptr;
    guard([&] { if (!path) throw std::invalid_argument("Null path"); out = new IcohpList(readIcohplist(path)); });
    return out;
}
AF_EXPORT void af_icohplist_free(void* h) { delete static_cast<IcohpList*>(h); }
AF_EXPORT int af_icohplist_count(void* h) { return h ? static_cast<int>(static_cast<IcohpList*>(h)->entries.size()) : 0; }
AF_EXPORT int af_icohplist_spin_polarized(void* h) { return h && static_cast<IcohpList*>(h)->spinPolarized ? 1 : 0; }
AF_EXPORT int af_icohplist_entry(void* h, int index, int* atom1, int* atom2, double* length, int* numBonds, double* icohpUp, double* icohpDown)
{
    return guard([&]
    {
        if (!atom1 || !atom2 || !length || !numBonds || !icohpUp || !icohpDown) throw std::invalid_argument("Null buffer");
        const auto& entry = icohplist(h).entries.at(index);
        *atom1 = entry.atom1; *atom2 = entry.atom2; *length = entry.length; *numBonds = entry.numBonds;
        *icohpUp = entry.icohp[0]; *icohpDown = entry.icohp[1];
    });
}
