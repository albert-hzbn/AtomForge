// Reference results of the C++ electronic library (electronic.json and saved
// volumes) for the Rust af-electronic tests: file reading, grid algebra, field
// analysis, Fourier tools, charge analysis, Bader, topology, LOBSTER and
// Wannier90 on the fixtures of make_electronic_fixtures.py.
// Usage (from the repository root): golden_electronic <fixture dir> <output dir>.

#include "electronic/BaderPartition.h"
#include "electronic/ChargeAnalysis.h"
#include "electronic/Lobster.h"
#include "electronic/Topology.h"
#include "electronic/Volume.h"
#include "electronic/Wannier.h"
#include "science/Json.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>

using namespace atomforge::electronic;
using atomforge::science::Json;

namespace
{
Json vec3(glm::dvec3 v) { return Json::array({v.x, v.y, v.z}); }
Json numbers(const std::vector<double>& values)
{
    Json out = Json::array();
    for (double v : values) out.push(v);
    return out;
}
Json gridJson(const Grid& g)
{
    Json o = Json::object();
    o["shape"] = Json::array({g.shape[0], g.shape[1], g.shape[2]});
    Json cell = Json::array();
    for (int c = 0; c < 3; ++c) for (int r = 0; r < 3; ++r) cell.push(g.cell[c][r]);
    o["cell_columns"] = cell;
    o["origin"] = vec3(g.origin);
    o["periodic"] = g.periodic;
    o["name"] = g.name;
    o["unit"] = g.unit;
    o["values"] = numbers(g.values);
    return o;
}
Json volumeJson(const Volume& v)
{
    Json o = Json::object();
    Json sites = Json::array();
    for (const auto& s : v.sites) sites.push(Json::array({s.number, s.position.x, s.position.y, s.position.z}));
    o["sites"] = sites;
    Json fields = Json::array();
    for (const auto& f : v.fields) fields.push(gridJson(f));
    o["fields"] = fields;
    return o;
}
Json vec2s(const std::vector<glm::dvec2>& list)
{
    Json out = Json::array();
    for (auto p : list) out.push(Json::array({p.x, p.y}));
    return out;
}
Json vec3s(const std::vector<glm::dvec3>& list)
{
    Json out = Json::array();
    for (auto p : list) out.push(vec3(p));
    return out;
}
Json reflections(const std::vector<Reflection>& list)
{
    Json out = Json::array();
    for (const auto& r : list) out.push(Json::array({r.hkl.x, r.hkl.y, r.hkl.z, r.value.real(), r.value.imag()}));
    return out;
}

Json cases = Json::array();
void attempt(const std::string& name, const std::function<Json()>& work)
{
    Json entry = Json::object();
    entry["name"] = name;
    try { entry["result"] = work(); }
    catch (const std::exception& e) { entry["error"] = std::string(e.what()); }
    cases.push(entry);
}
}

int main(int argc, char** argv)
{
    if (argc < 3) { std::fprintf(stderr, "usage: golden_electronic <fixture dir> <output dir>\n"); return 1; }
    const std::string in = std::string(argv[1]) + "/", out = std::string(argv[2]) + "/";
    std::filesystem::create_directories(out + "electronic_saved");

    // Reading.
    const char* loads[][3] = {
        {"density.cube", "auto", "bohr"}, {"density.cube", "density", "bohr"}, {"density.cube", "raw", "angstrom"},
        {"negative_axis.cube", "density", "bohr"}, {"two_channel.cube", "potential", "bohr"}, {"orbitals.cube", "auto", "bohr"},
        {"truncated.cube", "auto", "bohr"}, {"extra.cube", "auto", "bohr"}, {"D_exponent.cube", "auto", "bohr"},
        {"density.cube", "elf", "bohr"}, {"density.cube", "charge", "bohr"}, {"density.cube", "auto", "nm"},
        {"CHGCAR", "auto", "bohr"}, {"CHGCAR", "raw", "bohr"}, {"LOCPOT", "auto", "bohr"}, {"ELFCAR", "auto", "bohr"},
        {"CHGCAR_four", "auto", "bohr"}, {"CHGCAR_three", "auto", "bohr"}, {"CHGCAR_traj", "auto", "bohr"}, {"grid.raw", "auto", "bohr"},
        {"grid.xsf", "auto", "bohr"}, {"grid.xsf", "density", "bohr"}, {"nogrid.xsf", "auto", "bohr"}, {"anim.xsf", "auto", "bohr"},
        {"missing.cube", "auto", "bohr"},
    };
    for (const auto& l : loads)
        attempt(std::string("load ") + l[0] + " " + l[1] + " " + l[2], [&] { return volumeJson(loadVolume(in + l[0], l[1], l[2])); });

    const Volume cubeVolume = loadVolume(in + "density.cube", "density", "bohr");
    const Volume chgcar = loadVolume(in + "CHGCAR", "auto", "bohr");
    const Volume locpot = loadVolume(in + "LOCPOT", "auto", "bohr");
    const Volume xsf = loadVolume(in + "grid.xsf", "auto", "bohr");
    const Grid& D = cubeVolume.fields[0];
    const Grid& P = chgcar.fields[0];
    const Grid& S = chgcar.fields[1];
    const Grid& L = locpot.fields[0];
    std::vector<glm::dvec3> cubeSites, chgSites;
    for (const auto& s : cubeVolume.sites) cubeSites.push_back(s.position);
    for (const auto& s : chgcar.sites) chgSites.push_back(s.position);

    // Writing (files compared byte for byte by the Rust test).
    const char* saves[][3] = {{"cube", "cube", "density_cube"}, {"cube", "xsf", "density_xsf"}, {"chg", "vasp", "chgcar_vasp"},
                              {"chg", "cube", "chgcar_cube"}, {"chg", "xsf", "chgcar_xsf"}, {"loc", "cube", "locpot_cube"},
                              {"cube", "vasp", "error_vasp"}, {"chg", "png", "error_png"}};
    for (const auto& s : saves)
        attempt(std::string("save ") + s[2], [&] {
            const Volume& v = std::string(s[0]) == "cube" ? cubeVolume : std::string(s[0]) == "chg" ? chgcar : locpot;
            saveVolume(v, out + "electronic_saved/" + s[2], s[1]);
            return Json(true);
        });

    // Grid algebra.
    for (const char* op : {"add", "subtract", "multiply", "divide", "power"})
        attempt(std::string("arithmetic ") + op, [&] { return gridJson(arithmetic(P, S, op)); });
    attempt("arithmetic divide P", [&] { return gridJson(arithmetic(P, scale(P, 0.25), "divide")); });
    attempt("arithmetic units", [&] { return gridJson(arithmetic(P, L, "add")); });
    attempt("arithmetic aligned", [&] { return gridJson(arithmetic(D, P, "multiply")); });
    attempt("scale", [&] { return gridJson(scale(P, 2.5)); });
    attempt("scale nan", [&] { return gridJson(scale(P, std::nan(""))); });
    attempt("resample finite", [&] {
        Grid target;
        target.shape = {5, 4, 3};
        target.cell = glm::dmat3(1.5, 0, 0, 0.1, 1.2, 0, 0, 0.05, 0.9);
        target.origin = D.origin + glm::dvec3(0.4, 0.5, 0.3);
        target.values.assign(60, 0.0);
        return gridJson(resample(D, target));
    });
    attempt("resample periodic", [&] {
        Grid target = P;
        target.shape = {5, 6, 4};
        target.values.assign(120, 0.0);
        return gridJson(resample(P, target));
    });
    attempt("resample outside", [&] {
        Grid target = D;
        target.origin += glm::dvec3(50, 0, 0);
        return gridJson(resample(D, target));
    });
    attempt("periodic endpoints", [&] {
        Grid finite = P;
        finite.periodic = false;
        finite.shape = {P.shape[0] + 1, P.shape[1] + 1, P.shape[2] + 1};
        finite.values.assign(static_cast<std::size_t>(finite.shape[0]) * finite.shape[1] * finite.shape[2], 0.0);
        for (int z = 0; z < finite.shape[2]; ++z)
            for (int y = 0; y < finite.shape[1]; ++y)
                for (int x = 0; x < finite.shape[0]; ++x)
                    finite.values[finite.index(x, y, z)] = P.values[P.index(x % P.shape[0], y % P.shape[1], z % P.shape[2])];
        return gridJson(periodicEndpoints(finite));
    });
    attempt("periodic endpoints mismatch", [&] { return gridJson(periodicEndpoints(D)); });
    attempt("minimum images", [&] {
        Json list = Json::array();
        const glm::dmat3 skew(4, 0, 0, 3.5, 1.2, 0, 0.7, -0.4, 2.0);
        for (glm::dvec3 d : {glm::dvec3(3.9, 0.2, -0.1), glm::dvec3(-7.1, 2.5, 1.9), glm::dvec3(0.1, 0.1, 0.1), glm::dvec3(12.3, -8.8, 4.4)})
            list.push(vec3(minimumImage(d, skew)));
        return list;
    });

    // Integration.
    attempt("integrate D", [&] { return Json(integrate(D)); });
    attempt("integrate P", [&] { return Json(integrate(P)); });
    attempt("integrate sphere D", [&] { return Json(integrateSphere(D, cubeSites[0], 0.9)); });
    attempt("integrate sphere P", [&] { return Json(integrateSphere(P, chgSites[1], 1.3)); });
    attempt("integrate sphere bad", [&] { return Json(integrateSphere(P, chgSites[1], -1)); });
    attempt("voronoi D", [&] {
        Json list = Json::array();
        for (const auto& b : voronoiIntegrate(D, cubeSites)) list.push(Json::array({b.integral, b.volume}));
        return list;
    });
    attempt("voronoi P", [&] {
        Json list = Json::array();
        for (const auto& b : voronoiIntegrate(P, chgSites)) list.push(Json::array({b.integral, b.volume}));
        return list;
    });
    for (int limit : {0, 3})
        for (int which = 0; which < 2; ++which)
            attempt("peaks " + std::to_string(which) + " " + std::to_string(limit), [&] {
                Json list = Json::array();
                for (const auto& p : peaks(which ? P : D, limit)) list.push(Json::array({p.position.x, p.position.y, p.position.z, p.value}));
                return list;
            });

    // Field analysis.
    for (int which = 0; which < 2; ++which)
    {
        const Grid& g = which ? P : D;
        const std::string tag = which ? " P" : " D";
        attempt("gradient" + tag, [&] { Json l = Json::array(); for (const auto& c : gradient(g)) l.push(gridJson(c)); return l; });
        attempt("laplacian" + tag, [&] { return gridJson(laplacian(g)); });
        attempt("energy density" + tag, [&] { Json l = Json::array(); for (const auto& c : energyDensity(g, 1e-6)) l.push(gridJson(c)); return l; });
        attempt("rdg" + tag, [&] { return gridJson(reducedDensityGradient(g, 1e-6)); });
        attempt("signed density" + tag, [&] { return gridJson(signedDensity(g, 1e-6)); });
        attempt("dori" + tag, [&] { return gridJson(dori(g, 1e-6)); });
        attempt("smooth" + tag, [&] { return gridJson(smooth(g, 0.6, 2)); });
        for (int axis = 0; axis < 3; ++axis)
            attempt("planar average" + tag + std::to_string(axis), [&] { return vec2s(planarAverage(g, axis)); });
        attempt("macroscopic" + tag, [&] { return vec2s(macroscopicAverage(g, 2, 3)); });
    }
    attempt("energy density S", [&] { Json l = Json::array(); for (const auto& c : energyDensity(S, 1e-6)) l.push(gridJson(c)); return l; });
    attempt("rdg units", [&] { return gridJson(reducedDensityGradient(L, 1e-6)); });
    attempt("smooth bad", [&] { return gridJson(smooth(P, 0.5, 60)); });
    attempt("line profile D", [&] { return vec2s(lineProfile(D, D.position(1, 1, 1), D.position(9, 8, 6), 17)); });
    attempt("line profile P", [&] { return vec2s(lineProfile(P, glm::dvec3(-1, -2, 0.5), glm::dvec3(7, 9, 4), 25)); });
    attempt("macroscopic bad", [&] { return vec2s(macroscopicAverage(P, 2, 2)); });
    attempt("section and contours", [&] {
        const Grid slice = section(D, D.position(1, 1, 3), glm::dvec3(2.0, 0.1, 0), glm::dvec3(0.05, 1.5, 0.2), 14, 11);
        Json o = Json::object();
        o["section"] = gridJson(slice);
        o["contours"] = vec3s(contours(slice, 0.05));
        return o;
    });
    attempt("section P", [&] { return gridJson(section(P, glm::dvec3(0.3, 0.2, 0.1), glm::dvec3(8, 1, 0), glm::dvec3(-1, 7, 2), 9, 7)); });
    attempt("section bad", [&] { return gridJson(section(P, glm::dvec3(0), glm::dvec3(1, 0, 0), glm::dvec3(2, 0, 0), 9, 7)); });

    // Fourier.
    std::vector<glm::ivec3> all;
    for (int h = -P.shape[0] / 2; h <= (P.shape[0] - 1) / 2; ++h)
        for (int k = -P.shape[1] / 2; k <= (P.shape[1] - 1) / 2; ++k)
            for (int l = -P.shape[2] / 2; l <= (P.shape[2] - 1) / 2; ++l) all.push_back({h, k, l});
    attempt("structure factors", [&] { return reflections(structureFactors(P, {{0, 0, 0}, {1, 0, 0}, {-1, 2, 1}, {3, -3, 2}})); });
    attempt("structure factors range", [&] { return reflections(structureFactors(P, {{4, 0, 0}})); });
    attempt("structure factors finite", [&] { return reflections(structureFactors(D, {{0, 0, 0}})); });
    attempt("fourier synthesis", [&] { return gridJson(fourierSynthesis(P, structureFactors(P, all))); });
    attempt("fourier synthesis partial", [&] { return gridJson(fourierSynthesis(P, structureFactors(P, {{1, 0, 0}}))); });
    attempt("fourier synthesis duplicate", [&] { return gridJson(fourierSynthesis(P, structureFactors(P, {{1, 0, 0}, {1, 0, 0}}))); });
    attempt("patterson", [&] { return gridJson(patterson(P)); });
    attempt("patterson odd", [&] {
        Grid odd = P;
        odd.shape = {5, 6, 3};
        odd.values.assign(90, 0.0);
        for (std::size_t i = 0; i < odd.values.size(); ++i) odd.values[i] = std::sin(0.37 * static_cast<double>(i)) + 1.0;
        return gridJson(patterson(odd));
    });
    attempt("atomic structure factors", [&] {
        std::vector<std::complex<double>> f;
        const std::vector<glm::ivec3> hkl = {{1, 0, 0}, {1, 1, 0}, {2, -1, 1}};
        for (std::size_t r = 0; r < hkl.size(); ++r)
            for (std::size_t j = 0; j < chgSites.size(); ++j) f.push_back({29.0 - 3.0 * r - 10.0 * j, 0.5 * r});
        return reflections(atomicStructureFactors(P, chgSites, {1.0, 0.8}, {0.5, 1.2}, hkl, f));
    });
    attempt("ewald", [&] {
        const glm::dmat3 cell(5.64, 0, 0, 0, 5.64, 0, 0, 0, 5.64);
        const std::vector<glm::dvec3> sites = {{0, 0, 0}, {2.82, 2.82, 0}, {2.82, 0, 2.82}, {0, 2.82, 2.82}, {2.82, 0, 0}, {0, 2.82, 0}, {0, 0, 2.82}, {2.82, 2.82, 2.82}};
        const auto r = ewald(cell, sites, {1, 1, 1, 1, -1, -1, -1, -1}, 0.6, 9.0, 4.0);
        Json o = Json::object();
        o["potentials"] = numbers(r.potentials);
        o["energy"] = r.energy;
        return o;
    });
    attempt("ewald skew", [&] {
        const glm::dmat3 cell(4.0, 0, 0, 1.3, 3.7, 0, 0.6, -0.5, 4.4);
        const auto r = ewald(cell, {{0.1, 0.2, 0.3}, {2.4, 1.9, 2.0}, {1.0, 3.0, 3.5}}, {2, -1.5, -0.5}, 0.8, 7.0, 5.0);
        Json o = Json::object();
        o["potentials"] = numbers(r.potentials);
        o["energy"] = r.energy;
        return o;
    });
    attempt("ewald charged", [&] { return Json(ewald(glm::dmat3(3.0), {{0, 0, 0}}, {1.0}, 0.5, 5, 5).energy); });

    // Charge analysis.
    const Grid P2 = scale(P, 0.9);
    attempt("density difference", [&] { return gridJson(densityDifference(P, P2, 0.5)); });
    attempt("density difference units", [&] { return gridJson(densityDifference(P, L)); });
    const Grid maskA = thresholdMask(P, 0.5, 3.0), maskB = thresholdMask(S, -0.1, 0.1);
    attempt("threshold mask", [&] { return gridJson(maskA); });
    for (const char* op : {"union", "intersection", "difference", "xor", "nand"})
        attempt(std::string("boolean ") + op, [&] { return gridJson(booleanMask(maskA, maskB, op)); });
    attempt("boolean unmasked", [&] { return gridJson(booleanMask(P, maskB, "union")); });
    attempt("apply mask", [&] { return gridJson(applyMask(P, maskA)); });
    attempt("split", [&] { Json l = Json::array(); for (const auto& g : splitDensity(S)) l.push(gridJson(g)); return l; });
    attempt("charge summary", [&] { const auto s = chargeSummary(S); return Json::array({s[0], s[1], s[2]}); });
    attempt("cumulative P", [&] { return vec2s(cumulativeCharge(P, 2)); });
    attempt("cumulative D", [&] { return vec2s(cumulativeCharge(D, 0)); });

    // Bader.
    for (int which = 0; which < 2; ++which)
        attempt(std::string("bader ") + (which ? "P" : "D"), [&] {
            const auto partition = baderOnGrid(which ? P : D);
            Json o = Json::object();
            Json ids = Json::array();
            for (int id : partition.basin) ids.push(id);
            o["basin"] = ids;
            Json basins = Json::array();
            for (const auto& b : partition.basins) basins.push(Json::array({b.charge, b.volume, b.maximum.x, b.maximum.y, b.maximum.z}));
            o["basins"] = basins;
            o["populations"] = numbers(baderPopulations(which ? P : D, which ? chgSites : cubeSites));
            return o;
        });
    attempt("bader negative", [&] { return Json(static_cast<int>(baderOnGrid(S).basins.size())); });

    // Topology.
    attempt("betti curve", [&] {
        Json list = Json::array();
        for (auto row : bettiCurve(D, {0.001, 0.01, 0.05, 0.2, 1.0, 5.0})) list.push(Json::array({row.x, row.y, row.z, row.w}));
        return list;
    });
    attempt("betti shell", [&] {
        Grid g;
        g.shape = {9, 9, 9};
        g.values.assign(729, 0.0);
        for (int z = 0; z < 9; ++z)
            for (int y = 0; y < 9; ++y)
                for (int x = 0; x < 9; ++x)
                {
                    const double r = std::sqrt((x - 4.0) * (x - 4.0) + (y - 4.0) * (y - 4.0) + (z - 4.0) * (z - 4.0));
                    g.values[g.index(x, y, z)] = (r > 1.5 && r < 3.5) ? 1.0 : 0.0;
                }
        const auto b = bettiNumbers(g, 0.5);
        return Json::array({b.betti0, b.betti1, b.betti2});
    });
    attempt("betti periodic", [&] { return Json(bettiNumbers(P, 1.0).betti0); });

    // LOBSTER.
    for (const char* name : {"COHPCAR.lobster", "COHPCAR.spin", "COHPCAR.noheaders", "COHPCAR.badheader", "COHPCAR.nolabel", "COHPCAR.short", "missing"})
        attempt(std::string("cohpcar ") + name, [&] {
            const auto c = readCohpcar(in + name);
            Json o = Json::object();
            o["spin"] = c.spinPolarized;
            o["fermi"] = c.fermiEnergy;
            o["energies"] = numbers(c.energies);
            Json bonds = Json::array();
            for (const auto& b : c.bonds) bonds.push(Json::array({b.atom1, b.atom2, b.length}));
            o["bonds"] = bonds;
            Json curves = Json::array();
            for (std::size_t s = 0; s < c.cohp.size(); ++s)
            {
                curves.push(numbers(c.averageCohp[s]));
                curves.push(numbers(c.averageIcohp[s]));
                for (std::size_t b = 0; b < c.cohp[s].size(); ++b) { curves.push(numbers(c.cohp[s][b])); curves.push(numbers(c.icohp[s][b])); }
            }
            o["curves"] = curves;
            return o;
        });
    for (const char* name : {"ICOHPLIST.lobster", "ICOHPLIST.spin", "ICOHPLIST.old", "ICOHPLIST.bad", "missing"})
        attempt(std::string("icohplist ") + name, [&] {
            const auto l = readIcohplist(in + name);
            Json o = Json::object();
            o["spin"] = l.spinPolarized;
            Json entries = Json::array();
            for (const auto& e : l.entries) entries.push(Json::array({e.atom1, e.atom2, e.length, e.numBonds, e.icohp[0], e.icohp[1]}));
            o["entries"] = entries;
            return o;
        });

    // Wannier90.
    for (const char* name : {"short_hr.dat", "missing_hr.dat"})
        attempt(std::string("wannier ") + name, [&] { return Json(readWannierHr(in + name).numWann); });
    const auto hr = readWannierHr(in + "haldane_hr.dat");
    const auto nonHermitian = readWannierHr(in + "nonhermitian_hr.dat");
    attempt("wannier read", [&] {
        Json o = Json::object();
        o["num_wann"] = hr.numWann;
        Json vectors = Json::array();
        for (std::size_t i = 0; i < hr.vectors.size(); ++i)
        {
            Json row = Json::array({hr.vectors[i].x, hr.vectors[i].y, hr.vectors[i].z, hr.degeneracies[i]});
            for (auto v : hr.hamiltonians[i]) { row.push(v.real()); row.push(v.imag()); }
            vectors.push(row);
        }
        o["blocks"] = vectors;
        return o;
    });
    attempt("bloch", [&] {
        Json list = Json::array();
        for (auto v : blochHamiltonian(hr, glm::dvec3(0.13, 0.27, 0.0))) { list.push(v.real()); list.push(v.imag()); }
        return list;
    });
    attempt("bands", [&] { return numbers(wannierBands(hr, {{0, 0, 0}, {0.5, 0, 0}, {1.0 / 3, 1.0 / 3, 0}, {0.2, 0.7, 0.1}})); });
    attempt("bands non-Hermitian", [&] { return numbers(wannierBands(nonHermitian, {{0.1, 0.2, 0.3}})); });
    attempt("berry", [&] { return numbers(berryCurvature(hr, glm::dvec3(0.21, 0.37, 0), 0, 1)); });
    attempt("berry plane", [&] { return numbers(berryCurvature(hr, glm::dvec3(0.21, 0.37, 0), 1, 1)); });
    for (int band : {0, 1})
        attempt("chern " + std::to_string(band), [&] { return Json(chernNumber(hr, band, 0, 1, 12, 0.0)); });
    attempt("chern band", [&] { return Json(chernNumber(hr, 2, 0, 1, 12, 0.0)); });

    std::ofstream(out + "electronic.json") << cases.dump();
    std::printf("wrote %selectronic.json (%zu cases)\n", out.c_str(), cases.size());
    return 0;
}
