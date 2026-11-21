// Reference outputs of the C++ algorithms (algorithms.json) for the Rust
// af-algorithms tests. Each section mirrors a Rust test; inputs are written
// alongside so both sides start from the same structures.
// Usage (from the repository root): golden_algorithms <output dir>.

#include "algorithms/AmorphousBuilder.h"
#include "algorithms/AngularDistributionAnalysis.h"
#include "algorithms/SubstitutionalSolidSolutionBuilder.h"
#include "algorithms/BulkCrystalBuilder.h"
#include "algorithms/CSLComputation.h"
#include "algorithms/CellSculptorAlgo.h"
#include "algorithms/DislocationBuilder.h"
#include "algorithms/DislocationFit.h"
#include "algorithms/NyeTensor.h"
#include "algorithms/VitekMap.h"
#include "algorithms/DragPrep.h"
#include "algorithms/ElasticConstants.h"
#include "algorithms/NanostructureTools.h"
#include "algorithms/PatternMatch.h"
#include "algorithms/SQSBuilder.h"
#include "algorithms/StackingFaultBuilder.h"
#include "algorithms/SurfaceBuilder.h"
#include "algorithms/InterfaceBuilder.h"
#include "algorithms/InterfaceMatching.h"
#include "algorithms/PolyCrystalBuilder.h"
#include "algorithms/NanoCrystalBuilder.h"
#include "ElementData.h"
#include "algorithms/StrainTool.h"
#include "algorithms/VacancyBuilder.h"
#include "algorithms/MeshLoader.h"
#include "algorithms/ShapeModel.h"
#include "algorithms/CommonNeighbourAnalysis.h"
#include "algorithms/RadialDistributionAnalysis.h"
#include "algorithms/ShortRangeOrderAnalysis.h"
#include "algorithms/VoronoiComputation.h"
#include "algorithms/InterstitialVoidAnalysis.h"
#include "science/Json.h"
#include "science/ScienceData.h"
#include "util/ElementData.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <string>
#include <vector>

using atomforge::science::Json;

namespace
{
template <typename F>
Json attempt(F&& f)
{
    try {
        return f();
    } catch (const std::exception& e) {
        Json error = Json::object();
        error["error"] = std::string(e.what());
        return error;
    }
}

Json vec3(const std::array<double, 3>& v) { return Json::array({v[0], v[1], v[2]}); }
Json fvec3(const glm::vec3& v) { return Json::array({double(v.x), double(v.y), double(v.z)}); }

Json structureDump(const Structure& s)
{
    Json o = Json::object();
    Json atoms = Json::array();
    for (const auto& a : s.atoms) {
        Json atom = Json::object();
        atom["symbol"] = a.symbol;
        atom["z"] = a.atomicNumber;
        atom["position"] = Json::array({a.x, a.y, a.z});
        atom["colour"] = Json::array({double(a.r), double(a.g), double(a.b)});
        atom["fixed"] = int(a.fixed);
        atom["magmom"] = a.magmom;
        atom["charge"] = a.charge;
        atom["velocity"] = Json::array({a.velocity[0], a.velocity[1], a.velocity[2]});
        atoms.push(atom);
    }
    o["atoms"] = atoms;
    o["has_cell"] = s.hasUnitCell;
    o["cell"] = Json::array({vec3(s.cellVectors[0]), vec3(s.cellVectors[1]), vec3(s.cellVectors[2])});
    o["pbc"] = Json::array({s.pbc[0], s.pbc[1], s.pbc[2]});
    o["offset"] = vec3(s.cellOffset);
    Json colours = Json::array();
    for (const auto& c : s.grainColors) colours.push(Json::array({double(c[0]), double(c[1]), double(c[2])}));
    o["grain_colors"] = colours;
    Json regions = Json::array();
    for (int r : s.grainRegionIds) regions.push(r);
    o["grain_region_ids"] = regions;
    Json property = Json::array(), frozen = Json::array();
    for (double v : s.atomProperty) property.push(std::isnan(v) ? Json("nan") : Json(v));
    for (double v : s.frozenProperty) frozen.push(std::isnan(v) ? Json("nan") : Json(v));
    o["atom_property"] = property;
    o["frozen_property"] = frozen;
    o["loop_points"] = int(s.dislocationLoopPoints.size());
    o["detection_done"] = s.dislocationDetectionDone;
    return o;
}

AtomSite atom(const std::string& symbol, double x, double y, double z)
{
    AtomSite a;
    a.symbol = symbol;
    a.atomicNumber = atomforge::science::atomicNumber(symbol);
    a.x = x; a.y = y; a.z = z;
    getDefaultElementColor(a.atomicNumber, a.r, a.g, a.b);
    return a;
}

// A lattice of `basis` (fractional) repeated n times in a cell of rows `cell`.
Structure crystal(const std::array<std::array<double, 3>, 3>& cell, const std::vector<std::pair<std::string, std::array<double, 3>>>& basis,
                  int na, int nb, int nc)
{
    Structure s;
    s.hasUnitCell = true;
    for (int r = 0; r < 3; ++r)
        for (int k = 0; k < 3; ++k) s.cellVectors[r][k] = cell[r][k] * (r == 0 ? na : r == 1 ? nb : nc);
    for (int i = 0; i < na; ++i)
        for (int j = 0; j < nb; ++j)
            for (int k = 0; k < nc; ++k)
                for (const auto& [symbol, f] : basis) {
                    const double u = i + f[0], v = j + f[1], w = k + f[2];
                    s.atoms.push_back(atom(symbol, u * cell[0][0] + v * cell[1][0] + w * cell[2][0], u * cell[0][1] + v * cell[1][1] + w * cell[2][1],
                                           u * cell[0][2] + v * cell[1][2] + w * cell[2][2]));
                }
    return s;
}

std::vector<std::pair<std::string, Structure>> structures()
{
    std::vector<std::pair<std::string, Structure>> out;
    const double a = 3.615;
    const std::array<std::array<double, 3>, 3> cubic{{{a, 0, 0}, {0, a, 0}, {0, 0, a}}};
    const std::vector<std::pair<std::string, std::array<double, 3>>> fcc{{"Cu", {0, 0, 0}}, {"Cu", {0, .5, .5}}, {"Cu", {.5, 0, .5}}, {"Cu", {.5, .5, 0}}};
    out.push_back({"fcc 2x2x2", crystal(cubic, fcc, 2, 2, 2)});
    out.push_back({"fcc cell", crystal(cubic, fcc, 1, 1, 1)});
    const double b = 2.87;
    out.push_back({"bcc 3x3x3", crystal({{{b, 0, 0}, {0, b, 0}, {0, 0, b}}}, {{"Fe", {0, 0, 0}}, {"Fe", {.5, .5, .5}}}, 3, 3, 3)});
    const double h = 2.95, c = 4.68;
    out.push_back({"hcp 3x3x2", crystal({{{h, 0, 0}, {-0.5 * h, std::sqrt(3.0) / 2 * h, 0}, {0, 0, c}}},
                                        {{"Ti", {1.0 / 3, 2.0 / 3, 0.25}}, {"Ti", {2.0 / 3, 1.0 / 3, 0.75}}}, 3, 3, 2)});
    const double d = 5.431;
    out.push_back({"diamond", crystal({{{d, 0, 0}, {0, d, 0}, {0, 0, d}}},
                                      {{"Si", {0, 0, 0}}, {"Si", {0, .5, .5}}, {"Si", {.5, 0, .5}}, {"Si", {.5, .5, 0}},
                                       {"Si", {.25, .25, .25}}, {"Si", {.25, .75, .75}}, {"Si", {.75, .25, .75}}, {"Si", {.75, .75, .25}}}, 1, 1, 1)});
    const double n = 5.64;
    out.push_back({"rocksalt", crystal({{{n, 0, 0}, {0, n, 0}, {0, 0, n}}},
                                       {{"Na", {0, 0, 0}}, {"Na", {0, .5, .5}}, {"Na", {.5, 0, .5}}, {"Na", {.5, .5, 0}},
                                        {"Cl", {.5, 0, 0}}, {"Cl", {0, .5, 0}}, {"Cl", {0, 0, .5}}, {"Cl", {.5, .5, .5}}}, 1, 1, 1)});
    // A random triclinic cell with an offset.
    Structure random;
    random.hasUnitCell = true;
    random.cellVectors = {{{7.1, 0, 0}, {1.2, 6.8, 0}, {-0.9, 0.7, 7.4}}};
    random.cellOffset = {0.3, -0.2, 0.5};
    std::mt19937 rng(3);
    std::uniform_real_distribution<double> u(0, 1);
    const char* symbols[] = {"Al", "Ni", "O", "Fe"};
    for (int i = 0; i < 24; ++i) {
        const double f[3] = {u(rng), u(rng), u(rng)};
        double p[3];
        for (int k = 0; k < 3; ++k) p[k] = random.cellOffset[k] + f[0] * random.cellVectors[0][k] + f[1] * random.cellVectors[1][k] + f[2] * random.cellVectors[2][k];
        random.atoms.push_back(atom(symbols[i % 4], p[0], p[1], p[2]));
    }
    out.push_back({"random triclinic", random});
    // A cluster without a cell.
    Structure cluster = crystal(cubic, fcc, 2, 2, 2);
    cluster.hasUnitCell = false;
    out.push_back({"fcc cluster", cluster});
    // Thermally displaced fcc and a random ternary fcc alloy.
    Structure warm = crystal(cubic, fcc, 3, 3, 3);
    std::normal_distribution<double> noise(0.0, 0.08);
    for (auto& atom : warm.atoms) { atom.x += noise(rng); atom.y += noise(rng); atom.z += noise(rng); }
    out.push_back({"warm fcc", warm});
    Structure alloy = crystal(cubic, fcc, 3, 3, 3);
    std::uniform_int_distribution<int> pick(0, 2);
    const char* alloySymbols[] = {"Ni", "Cu", "Al"};
    for (auto& atom : alloy.atoms) {
        const std::string symbol = alloySymbols[pick(rng)];
        atom = ::atom(symbol, atom.x, atom.y, atom.z);
    }
    out.push_back({"fcc alloy", alloy});
    return out;
}

Json voidsJson(const InterstitialVoidDetectionResult& r)
{
    Json o = Json::object();
    o["success"] = r.success;
    o["message"] = r.message;
    Json regions = Json::array();
    for (const auto& v : r.regions) {
        Json g = Json::object();
        g["position"] = fvec3(v.position);
        g["clearance"] = double(v.clearance);
        g["coordination"] = v.coordination;
        g["volume"] = double(v.volume);
        g["volume_bin"] = v.volumeBin;
        g["kind"] = int(v.kind);
        Json neighbours = Json::array();
        for (int i : v.neighborAtomIndices) neighbours.push(i);
        g["neighbors"] = neighbours;
        Json faces = Json::array();
        for (const auto& f : v.faces) faces.push(Json::array({f[0], f[1], f[2]}));
        g["faces"] = faces;
        Json vertices = Json::array();
        for (const auto& f : v.faceVertices) vertices.push(Json::array({fvec3(f[0]), fvec3(f[1]), fvec3(f[2])}));
        g["face_vertices"] = vertices;
        regions.push(g);
    }
    o["regions"] = regions;
    return o;
}

// Positions, symbols and cell of a built structure, compactly.
Json builtJson(const Structure& s)
{
    Json o = Json::object();
    std::string symbols;
    Json positions = Json::array(), colours = Json::array();
    for (const auto& a : s.atoms) {
        symbols += (symbols.empty() ? "" : " ") + a.symbol;
        positions.push(a.x);
        positions.push(a.y);
        positions.push(a.z);
        colours.push(double(a.r));
        colours.push(double(a.g));
        colours.push(double(a.b));
    }
    o["symbols"] = symbols;
    o["positions"] = positions;
    o["colours"] = colours;
    o["has_cell"] = s.hasUnitCell;
    o["cell"] = Json::array({vec3(s.cellVectors[0]), vec3(s.cellVectors[1]), vec3(s.cellVectors[2])});
    o["offset"] = vec3(s.cellOffset);
    return o;
}

Json bulkCase(CrystalSystem system, int number, LatticeParameters params, const std::vector<std::pair<std::string, std::array<double, 3>>>& basis,
              const std::vector<glm::vec3>& colours, int origin)
{
    std::vector<AtomSite> atoms;
    for (const auto& [symbol, f] : basis) {
        AtomSite a;
        a.symbol = symbol;
        a.atomicNumber = atomforge::science::atomicNumber(symbol);
        a.x = f[0]; a.y = f[1]; a.z = f[2];
        atoms.push_back(a);
    }
    Structure s;
    const BulkBuildResult r = buildBulkCrystal(s, system, number, params, atoms, colours, origin);
    Json o = Json::object();
    o["system"] = int(system);
    o["number"] = number;
    o["params"] = Json::array({params.a, params.b, params.c, params.alpha, params.beta, params.gamma});
    Json basisJson = Json::array();
    for (const auto& [symbol, f] : basis) basisJson.push(Json::array({symbol, f[0], f[1], f[2]}));
    o["basis"] = basisJson;
    o["colour_table"] = int(colours.size());
    o["origin"] = origin;
    o["success"] = r.success;
    o["message"] = r.message;
    o["generated"] = r.generatedAtoms;
    o["result_number"] = r.spaceGroupNumber;
    o["symbol"] = r.spaceGroupSymbol;
    o["built"] = builtJson(s);
    return o;
}

Json cnaJson(const atomforge::analysis::CnaResult& r)
{
    Json o = Json::object();
    o["valid"] = r.valid;
    o["message"] = r.message;
    o["atom_count"] = r.atomCount;
    o["pair_count"] = r.pairCount;
    o["pbc_used"] = r.pbcUsed;
    Json signatures = Json::array();
    for (const auto& [sig, count] : r.signatureCounts) signatures.push(Json::array({sig.common, sig.bonds, sig.chain, count}));
    o["signatures"] = signatures;
    Json environments = Json::array();
    for (const auto& [name, count] : r.environmentCounts) environments.push(Json::array({name, count}));
    o["environments"] = environments;
    Json rows = Json::array();
    for (const auto& row : r.atomRows)
        rows.push(Json::array({row.index, row.atomicNumber, row.symbol, row.coordination, row.dominantSignature.common, row.dominantSignature.bonds,
                               row.dominantSignature.chain, row.dominantSignatureCount, row.environment}));
    o["rows"] = rows;
    return o;
}

Json rdfJson(const atomforge::analysis::RdfResult& r)
{
    Json o = Json::object();
    o["valid"] = r.valid;
    o["message"] = r.message;
    o["flags"] = Json::array({r.pbcUsed, r.normalized, r.hasFirstPeak, r.hasFirstMinimum, r.hasDistortionMetrics, r.distortionWindowAuto, r.hasFirstPeakFwhm});
    o["counts"] = Json::array({r.atomCount, r.refCount, r.targetCount, r.binCount});
    o["values"] = Json::array({double(r.volume), double(r.density), double(r.rMin), double(r.rMax), double(r.binWidth), double(r.firstPeakR),
                               double(r.firstPeakValue), double(r.firstMinimumR), double(r.firstMinimumValue), double(r.distortionWindowMin),
                               double(r.distortionWindowMax), double(r.shellMeanDistance), double(r.shellStdDev), double(r.shellRelativeDistortion),
                               double(r.shellCoordination), double(r.firstPeakFwhm)});
    Json bins = Json::array();
    for (const auto& b : r.bins) {
        bins.push(double(b.rCenter));
        bins.push(double(b.g));
        bins.push(double(b.rawCount));
        bins.push(double(b.cumulative));
    }
    o["bins"] = bins;
    return o;
}

Json adfJson(const AdfResult& r)
{
    Json o = Json::object();
    o["valid"] = r.valid;
    o["message"] = r.message;
    o["values"] = Json::array({double(r.rCutoff), r.binCount, double(r.binWidth), r.pbcUsed, r.normalized, r.nAtoms, r.nCentreAtoms, double(r.nTriplets)});
    Json elements = Json::array();
    for (const auto& [symbol, count] : r.elementCounts) elements.push(Json::array({symbol, count}));
    o["elements"] = elements;
    Json bins = Json::array();
    for (const auto& b : r.bins) {
        bins.push(double(b.angleDeg));
        bins.push(double(b.count));
        bins.push(double(b.value));
    }
    o["bins"] = bins;
    Json peaks = Json::array();
    for (const auto& p : r.peaks) peaks.push(Json::array({double(p.angleDeg), double(p.value), p.label}));
    o["peaks"] = peaks;
    Json stats = Json::array();
    for (const auto& [symbol, st] : r.coordStats) stats.push(Json::array({symbol, double(st.mean), double(st.stddev), st.count}));
    o["stats"] = stats;
    return o;
}

Json sroJson(const SroAnalysisResult& r)
{
    Json o = Json::object();
    o["success"] = r.success;
    o["message"] = r.message;
    const auto common = [](int atoms, int shells, const std::vector<double>& distances, const std::map<std::string, double>& composition, const std::string& message) {
        Json c = Json::object();
        c["atoms"] = atoms;
        c["shells"] = shells;
        c["distances"] = atomforge::science::toJson(distances);
        Json comp = Json::array();
        for (const auto& [symbol, x] : composition) comp.push(Json::array({symbol, x}));
        c["composition"] = comp;
        c["message"] = message;
        return c;
    };
    Json wc = common(r.warrenCowley.nAtoms, r.warrenCowley.nShells, r.warrenCowley.shellDistances, r.warrenCowley.composition, r.warrenCowley.message);
    Json entries = Json::array();
    for (const auto& e : r.warrenCowley.entries)
        entries.push(Json::array({e.shell, e.shellDistance, e.centerElement, e.neighborElement, e.nij, e.zi, e.pij, e.cj, e.alpha}));
    wc["entries"] = entries;
    o["warren_cowley"] = wc;
    Json rc = common(r.raoCurtin.nAtoms, r.raoCurtin.nShells, r.raoCurtin.shellDistances, r.raoCurtin.composition, r.raoCurtin.message);
    rc["coordination"] = atomforge::science::toJson(r.raoCurtin.shellCoordination);
    Json rcEntries = Json::array();
    for (const auto& e : r.raoCurtin.entries)
        rcEntries.push(Json::array({e.shell, e.shellDistance, e.elem_i, e.elem_j, e.n_ij, e.n_ji, e.z_shell, e.p_ij, e.ci, e.cj, e.alpha}));
    rc["entries"] = rcEntries;
    o["rao_curtin"] = rc;
    return o;
}

// The CSL bicrystal as the CLI (BuildDefectCLI.cpp) composes it.
Json cslBicrystal(const Structure& input, const int axis[3], int sigmaTarget, int sigmaMax, int planeIdx, int ucA, int ucB, float vacuum, float gap,
                  float overlapD, bool conventional)
{
    Json o = Json::object();
    std::vector<SigmaCandidate> candidates = computeGBInfo(axis, sigmaMax);
    if (candidates.empty()) { o["error"] = "no candidates"; return o; }
    int selIdx = 0;
    if (sigmaTarget > 0) {
        bool found = false;
        for (int i = 0; i < (int)candidates.size(); ++i)
            if (candidates[i].sigma == sigmaTarget) { selIdx = i; found = true; break; }
        if (!found) { o["error"] = "sigma not found"; return o; }
    }
    const SigmaCandidate& sel = candidates[selIdx];
    const int direction = planeIdx;
    Grain inputGrain = structureToGrain(input);
    int cslT[3][3];
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) cslT[i][j] = sel.csl[j][i];
    Grain grainA = makeSupercell(inputGrain, cslT);
    if (!isOrthogonal(grainA.cell)) grainA = setOrthogonalGrain(grainA, direction);
    int scaleB[3] = {1, 1, 1};
    scaleB[direction] = ucB;
    Grain tempA = makeSupercellDiag(grainA, scaleB[0], scaleB[1], scaleB[2]);
    Grain grainB = getBFromA(tempA);
    { int ones[3][3] = {{1,0,0},{0,1,0},{0,0,1}}; grainB = makeSupercell(grainB, ones); }
    int scaleA[3] = {1, 1, 1};
    scaleA[direction] = ucA;
    grainA = makeSupercellDiag(grainA, scaleA[0], scaleA[1], scaleA[2]);
    Grain gb = stackGrains(grainA, grainB, direction, vacuum, gap);
    int removed = removeOverlaps(gb.atoms, overlapD);
    Structure structure = grainToStructure(gb);
    {
        double inv[3][3];
        invertCell(gb.cell, inv);
        std::vector<double> layers;
        for (const auto& a : structure.atoms) {
            double cart[3] = {a.x, a.y, a.z};
            double frac[3];
            cartToFrac(cart, inv, frac);
            layers.push_back(wrapFrac(frac[direction]));
        }
        std::sort(layers.begin(), layers.end());
        double minSpacing = 1.0;
        for (size_t i = 1; i < layers.size(); i++) { double d = layers[i] - layers[i - 1]; if (d > 1e-8 && d < minSpacing) minSpacing = d; }
        structure.pbcBoundaryTol = static_cast<float>(minSpacing * 0.5 + 1e-6);
    }
    if (!conventional) reduceToPrimitiveGB(structure, direction);
    structure.grainRegionIds.assign(structure.atoms.size(), 0);
    {
        double cell[3][3];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) cell[i][j] = structure.cellVectors[i][j];
        double inv[3][3];
        invertCell(cell, inv);
        for (size_t i = 0; i < structure.atoms.size(); ++i) {
            double cart[3] = {structure.atoms[i].x, structure.atoms[i].y, structure.atoms[i].z};
            double frac[3];
            cartToFrac(cart, inv, frac);
            structure.grainRegionIds[i] = (wrapFrac(frac[direction]) >= 0.5) ? 1 : 0;
        }
    }
    o["sigma"] = sel.sigma;
    o["plane"] = Json::array({sel.plane[planeIdx][0], sel.plane[planeIdx][1], sel.plane[planeIdx][2]});
    o["removed"] = removed;
    o["tol"] = double(structure.pbcBoundaryTol);
    o["structure"] = structureDump(structure);
    return o;
}
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::fprintf(stderr, "usage: golden_algorithms <output dir>\n");
        return 1;
    }
    Json root = Json::object();
    const auto inputs = structures();
    Json inputList = Json::array();
    for (const auto& [name, s] : inputs) {
        Json item = Json::object();
        item["name"] = name;
        item["structure"] = structureDump(s);
        inputList.push(item);
    }
    root["inputs"] = inputList;

    Json voids = Json::array();
    struct Case { const char* name; InterstitialVoidDetectionParams params; };
    InterstitialVoidDetectionParams coarse;
    coarse.gridResolution = 9;
    coarse.minClearance = 0.3f;
    InterstitialVoidDetectionParams capped;
    capped.maxVoids = 5;
    capped.minSeparation = 1.0f;
    capped.firstNeighborShellFactor = 1.4f;
    capped.maxNeighborCount = 8;
    for (const Case& c : {Case{"default", {}}, Case{"coarse", coarse}, Case{"capped", capped}})
        for (std::size_t i = 0; i < inputs.size(); ++i) {
            Json item = Json::object();
            item["input"] = int(i);
            item["case"] = c.name;
            item["params"] = Json::array({c.params.gridResolution, c.params.maxVoids, double(c.params.minClearance), double(c.params.minSeparation),
                                          double(c.params.firstNeighborShellFactor), c.params.maxNeighborCount});
            item["result"] = voidsJson(detectInterstitialVoidRegions(inputs[i].second, c.params));
            voids.push(item);
        }
    root["voids"] = voids;

    // Bulk crystals: every space group, both origin choices, error paths.
    {
        std::vector<glm::vec3> colours;
        for (int z = 0; z <= 118; ++z) {
            float r, g, b;
            getDefaultElementColor(z, r, g, b);
            colours.push_back({r * 0.5f, g, b * 0.25f});
        }
        const std::vector<std::pair<std::string, std::array<double, 3>>> basis{
            {"Fe", {0.1, 0.2, 0.3}}, {"O", {0, 0, 0}}, {"Na", {0.5, 0.5, 0.25}}, {"H", {0.3333333, 0.6666667, 0.25}}, {"Cl", {-0.25, 1.125, 0.5}}};
        Json bulk = Json::array();
        Json labels = Json::array();
        for (int number = 0; number <= 231; ++number) {
            Json l = Json::object();
            l["label"] = spaceGroupLabel(number);
            l["two_origins"] = hasTwoOrigins(number);
            labels.push(l);
        }
        for (int number = 1; number <= 230; ++number) {
            CrystalSystem system = CrystalSystem::Triclinic;
            for (int r = 0; r < kNumSpaceGroupRanges; ++r)
                if (number >= kSpaceGroupRanges[r].first && number <= kSpaceGroupRanges[r].last) system = kSpaceGroupRanges[r].system;
            LatticeParameters params{4.1, 4.3, 4.7, 95, 100, 105};
            applySystemConstraints(system, params);
            bulk.push(bulkCase(system, number, params, basis, number % 3 ? colours : std::vector<glm::vec3>{}, 1));
            if (hasTwoOrigins(number)) bulk.push(bulkCase(system, number, params, basis, colours, 2));
        }
        const LatticeParameters cubic{3.6, 3.6, 3.6, 90, 90, 90};
        bulk.push(bulkCase(CrystalSystem::Cubic, 225, cubic, {{"Cu", {0, 0, 0}}}, colours, 3));
        bulk.push(bulkCase(CrystalSystem::Cubic, 225, cubic, {{"Cu", {0, 0, 0}}}, colours, 2));
        bulk.push(bulkCase(CrystalSystem::Cubic, 225, cubic, {}, colours, 1));
        bulk.push(bulkCase(CrystalSystem::Cubic, 0, cubic, {{"Cu", {0, 0, 0}}}, colours, 1));
        bulk.push(bulkCase(CrystalSystem::Cubic, 231, cubic, {{"Cu", {0, 0, 0}}}, colours, 1));
        bulk.push(bulkCase(CrystalSystem::Cubic, 225, {3.6, 3.7, 3.6, 90, 90, 90}, {{"Cu", {0, 0, 0}}}, colours, 1));
        bulk.push(bulkCase(CrystalSystem::Hexagonal, 194, {3.6, 3.6, 5.8, 90, 90, 90}, {{"Mg", {1.0 / 3, 2.0 / 3, 0.25}}}, colours, 1));
        bulk.push(bulkCase(CrystalSystem::Triclinic, 1, {3.6, -1, 3.6, 90, 90, 90}, {{"Cu", {0, 0, 0}}}, colours, 1));
        bulk.push(bulkCase(CrystalSystem::Triclinic, 1, {3.6, 3.6, 3.6, 0, 90, 90}, {{"Cu", {0, 0, 0}}}, colours, 1));
        bulk.push(bulkCase(CrystalSystem::Triclinic, 1, {3.6, 3.6, 3.6, 60, 60, 120}, {{"Cu", {0, 0, 0}}}, colours, 1));
        bulk.push(bulkCase(CrystalSystem::Monoclinic, 14, {3.6, 3.6, 3.6, 80, 90, 90}, {{"Cu", {0, 0, 0}}}, colours, 1));
        bulk.push(bulkCase(CrystalSystem::Orthorhombic, 62, {3.6, 3.6, 3.6, 90, 91, 90}, {{"Cu", {0, 0, 0}}}, colours, 1));
        bulk.push(bulkCase(CrystalSystem::Tetragonal, 139, {3.6, 3.5, 3.6, 90, 90, 90}, {{"Cu", {0, 0, 0}}}, colours, 1));
        bulk.push(bulkCase(CrystalSystem::Trigonal, 166, {3.6, 3.6, 3.6, 90, 90, 90}, {{"Cu", {0, 0, 0}}}, colours, 1));
        root["bulk"] = bulk;
        root["space_group_labels"] = labels;
    }

    // Structure analyses on every input.
    {
        using namespace atomforge::analysis;
        Json cna = Json::array(), rdf = Json::array(), adf = Json::array(), sro = Json::array();
        for (std::size_t i = 0; i < inputs.size(); ++i) {
            const Structure& s = inputs[i].second;
            for (int c = 0; c < 3; ++c) {
                CnaParams p;
                if (c == 1) p.usePbc = false;
                if (c == 2) p.cutoffRadius = 3.0f;
                Json item = Json::object();
                item["input"] = int(i);
                item["params"] = Json::array({double(p.cutoffScale), p.usePbc, double(p.cutoffRadius)});
                item["result"] = cnaJson(computeCna(s, p));
                cna.push(item);
            }
            for (int c = 0; c < 4; ++c) {
                RdfParams p;
                if (c == 1) { p.usePbc = false; p.smoothingPasses = 2; p.rMax = 6.0f; p.binCount = 120; }
                if (c == 2) { p.normalize = false; p.autoDistortionWindow = false; p.refAtomicNumberFilter = 28; p.targetAtomicNumberFilter = 13; p.rMin = 1.0f; }
                if (c == 3) { p.rMax = 12.0f; p.binCount = 300; p.enableDistortionAnalysis = false; }
                Json item = Json::object();
                item["input"] = int(i);
                item["params"] = Json::array({p.refAtomicNumberFilter, p.targetAtomicNumberFilter, p.usePbc, p.normalize, double(p.rMin), double(p.rMax), p.binCount,
                                              p.smoothingPasses, p.enableDistortionAnalysis, p.autoDistortionWindow, double(p.manualDistortionMin), double(p.manualDistortionMax)});
                item["result"] = rdfJson(computeRdf(s, p));
                rdf.push(item);
            }
            for (int c = 0; c < 3; ++c) {
                AdfParams p;
                if (c == 1) { p.usePbc = false; p.rCutoff = 2.9f; p.binCount = 90; p.smoothPasses = 0; p.normalize = false; }
                if (c == 2) { p.centreMode = AdfCentreMode::ByPair; p.centreSymbol = "Ni"; p.neighSymbol1 = "Cu"; p.neighSymbol2 = ""; p.rCutoff = 3.0f; }
                Json item = Json::object();
                item["input"] = int(i);
                item["params"] = Json::array({double(p.rCutoff), p.binCount, p.usePbc, int(p.centreMode), p.centreSymbol, p.neighSymbol1, p.neighSymbol2, p.normalize, p.smoothPasses});
                item["result"] = adfJson(computeADF(s, p));
                adf.push(item);
            }
            for (const auto& [shells, tolerance] : std::vector<std::pair<int, double>>{{3, 0.1}, {6, 0.25}}) {
                Json item = Json::object();
                item["input"] = int(i);
                item["params"] = Json::array({shells, tolerance});
                item["result"] = sroJson(analyzeSRO(s, shells, tolerance));
                sro.push(item);
            }
        }
        root["cna"] = cna;
        root["rdf"] = rdf;
        root["adf"] = adf;
        root["sro"] = sro;
    }

    // Voronoi cells (also with a boundary tolerance that replicates cells).
    {
        Json voronoi = Json::array();
        for (std::size_t i = 0; i < inputs.size(); ++i)
            for (float tol : {0.0f, 0.3f}) {
                Structure s = inputs[i].second;
                s.pbcBoundaryTol = tol;
                const VoronoiDiagram d = computeVoronoi(s);
                Json cells = Json::array();
                for (const auto& cell : d.cells) {
                    Json faces = Json::array();
                    for (const auto& face : cell.faces) {
                        Json vertices = Json::array();
                        for (const auto& v : face.vertices) {
                            vertices.push(double(v.x));
                            vertices.push(double(v.y));
                            vertices.push(double(v.z));
                        }
                        faces.push(vertices);
                    }
                    cells.push(faces);
                }
                Json item = Json::object();
                item["input"] = int(i);
                item["tolerance"] = double(tol);
                item["cells"] = cells;
                voronoi.push(item);
            }
        root["voronoi"] = voronoi;
    }

    // Substitutional solid solutions and amorphous packings.
    {
        Json sss = Json::array();
        const std::vector<std::pair<std::vector<std::pair<int, float>>, unsigned>> compositions = {
            {{{28, 0.5f}, {29, 0.5f}}, 12345u},
            {{{28, 0.2f}, {29, 0.3f}, {13, 0.5f}}, 7u},
            {{{26, 1.0f}, {24, 1.0f}, {28, 1.0f}}, 99u},
            {{{26, 0.333f}, {24, 0.333f}, {28, 0.333f}, {25, 0.001f}}, 3u},
            {{{29, 0.0f}}, 1u},
            {{{0, 0.5f}}, 1u},
            {{{29, -0.5f}}, 1u},
            {{}, 1u},
        };
        for (std::size_t i : {std::size_t(0), std::size_t(8), std::size_t(9)})
            for (const auto& [composition, seed] : compositions) {
                SSSParams params;
                for (const auto& [z, f] : composition) params.composition.push_back({z, f});
                params.seed = seed;
                const SSSResult r = buildSubstitutionalSolidSolution(inputs[i].second, params);
                Json item = Json::object();
                item["input"] = int(i);
                Json comp = Json::array();
                for (const auto& [z, f] : composition) comp.push(Json::array({z, double(f)}));
                item["composition"] = comp;
                item["seed"] = double(seed);
                item["success"] = r.success;
                item["message"] = r.message;
                Json counts = Json::array();
                for (const auto& [z, c] : r.elementCounts) counts.push(Json::array({z, c}));
                item["counts"] = counts;
                item["output"] = structureDump(r.output);
                sss.push(item);
            }
        root["solid_solution"] = sss;

        std::vector<float> radii = makeLiteratureCovalentRadii();
        std::vector<glm::vec3> colours;
        for (int z = 0; z <= 118; ++z) {
            float r, g, b;
            getDefaultElementColor(z, r, g, b);
            colours.push_back({r, g * 0.5f, b});
        }
        Json amorphous = Json::array();
        struct Case { std::vector<std::pair<int, int>> elements; int mode; float a, b, c, density, scale, tolerance; unsigned seed; int attempts; bool periodic; std::vector<AmorphousPairDist> pairs; };
        const std::vector<Case> cases = {
            {{{14, 60}}, 1, 0, 0, 0, 2.33f, 1.0f, 0.75f, 42, 1000, true, {}},
            {{{14, 30}, {8, 60}}, 1, 0, 0, 0, 2.2f, 1.1f, 0.75f, 7, 1000, true, {}},
            {{{26, 40}, {5, 10}}, 0, 9.5f, 10.0f, 10.5f, 0, 1.0f, 0.8f, 3, 200, false, {{26, 5, 1.6f}}},
            {{{6, 200}}, 0, 6.0f, 6.0f, 6.0f, 0, 1.0f, 0.9f, 5, 50, true, {}},
            {{{1, 10}, {8, 5}}, 0, 5.0f, 5.0f, 5.0f, 0, 1.0f, 0.75f, 11, 100, true, {{1, 1, 1.5f}, {8, 1, 0.9f}}},
            {{}, 1, 0, 0, 0, 2.0f, 1.0f, 0.75f, 1, 10, true, {}},
            {{{200, 5}}, 1, 0, 0, 0, 2.0f, 1.0f, 0.75f, 1, 10, true, {}},
            {{{14, -5}}, 1, 0, 0, 0, 2.0f, 1.0f, 0.75f, 1, 10, true, {}},
            {{{14, 0}}, 1, 0, 0, 0, 2.0f, 1.0f, 0.75f, 1, 10, true, {}},
            {{{14, 5}}, 0, 0.05f, 5.0f, 5.0f, 0, 1.0f, 0.75f, 1, 10, true, {}},
        };
        for (const Case& c : cases) {
            AmorphousParams params;
            for (const auto& [z, n] : c.elements) params.elements.push_back({z, n});
            params.boxMode = c.mode == 0 ? AmorphousBoxMode::Manual : AmorphousBoxMode::AutoFromDensity;
            params.boxA = c.a; params.boxB = c.b; params.boxC = c.c;
            params.targetDensity = c.density;
            params.cellScaleFactor = c.scale;
            params.covalentTolerance = c.tolerance;
            params.seed = c.seed;
            params.maxAttempts = c.attempts;
            params.periodic = c.periodic;
            params.pairDistances = c.pairs;
            const AmorphousResult r = buildAmorphousStructure(params, radii, colours);
            Json item = Json::object();
            Json elements = Json::array();
            for (const auto& [z, n] : c.elements) elements.push(Json::array({z, n}));
            item["elements"] = elements;
            item["values"] = Json::array({c.mode, double(c.a), double(c.b), double(c.c), double(c.density), double(c.scale), double(c.tolerance), double(c.seed),
                                          c.attempts, c.periodic});
            Json pairs = Json::array();
            for (const auto& pd : c.pairs) pairs.push(Json::array({pd.z1, pd.z2, double(pd.minDist)}));
            item["pairs"] = pairs;
            item["success"] = r.success;
            item["message"] = r.message;
            item["numbers"] = Json::array({r.requestedAtoms, r.placedAtoms, double(r.actualDensity), double(r.boxA), double(r.boxB), double(r.boxC)});
            item["output"] = structureDump(r.output);
            amorphous.push(item);
        }
        root["amorphous"] = amorphous;
    }

    // Shapes: scripts, distances, materials, bounds, volumes, errors.
    {
        using namespace atomforge::shape;
        std::vector<std::string> scripts;
        for (const auto& [title, script] : templates()) scripts.push_back(script);
        for (const char* extra : {
                 "box(size = 10, center = false); cylinder(h = 8, r1 = 3, r2 = 1, center = false);",
                 "cube([4, 6, 8]); capsule(r = 2, h = 6);",
                 "scale([1, 2, 0.5]) sphere(d = 12);",
                 "mirror([1, 1, 0]) translate([5, 0, 0]) cone(4, 1, 6);",
                 "rotate(a = 30, v = [1, 1, 0]) box(size = [10, 4, 2]);",
                 "offset(r = -1.5) octahedron(d = 8); grow(2) prism(n = 5, apothem = 4, h = 3, rotation = 15);",
                 "linear_array(n = 3, step = [6, 0, 0]) sphere(2); polar_array(n = 4, angle = 90) translate([8, 0, 0]) sphere(1.5);",
                 "smooth_intersection(k = 1) { sphere(6); box(9); } smooth_difference(k = 0.5) { box(8); sphere(4.5); }",
                 "polyhedron(hkl = [1, 1, 0], d = 7); facets(hkl = [[1, 0, 0], [0, 0, 1]], d = [6, 9], family = false);",
                 "A = [1, 2, 3] * 2 - [0.5, 0.5, 0.5]; translate(A / 2) sphere(r = 2 ^ 2 / 2 + sqrt(4) - min(1, 2, -1) + max(3) + atan2(1, 1) / 45);",
                 "x = cos(60) + sin(30) + tan(45) + asin(1) / 90 + acos(0) / 90 + atan(1) / 45 + abs(-1) + floor(1.7) + ceil(0.2) + round(0.5) + exp(0) + ln(1) + log(1) + pow(2, 3) + pi / PI + true; sphere(r = x);",
                 "hollow(t = 1) cylinder(d = 10, h = 4); intersect() { halfspace(hkl = [0, 0, 1], d = 1); sphere(5); } subtract() { box(6); slab(hkl = [1, 1, 1], thickness = 2); }",
                 "region(\"none\") sphere(3); region(\"Cu\") box(4); sphere(6); material(\"Ag\") translate([2, 0, 0]) sphere(2);",
                 "grid_array(n = [2, 3, 1], step = [5, 5, 0]) region(\"Fe\") box(2); box(30);",
                 "/* block\n comment */ # hash comment\n group() { move([1, 2, 3]) ellipsoid(2); torus(4, 1); blend(1) { sphere(2); } }",
                 "plane(normal = [0, 1, 0], d = 2); extrude(points = [[0, 0], [4, 0], [4, 3], [2, 5], [0, 3]], h = 2);",
                 "sphere(r = .5e1); cube(size = 0x4);",
             })
            scripts.push_back(extra);
        for (const char* bad : {
                 "", "region(\"Au\") sphere(3);", "sphere(r = -1);", "sphere(r = 1", "sphere(1, 2);", "sphere(radius = 2);", "foo(1);",
                 "box(size = [1, 2, 3, 4]);", "x = y + 1; sphere(x);", "sphere(r = 'a');", "translate([1, 2]) ;", "union() { sphere(1);",
                 "sphere(r = 1 / 0);", "sphere(r = [1, 2] * [3, 4]);", "sphere(r = [1, 2] + 3);", "sphere(r = 3 / [1]);", "@", "'open",
                 "cylinder(h = 4, r1 = 0, r2 = 0);", "prism(n = 2);", "polygon(points = [[0, 0], [1, 1]]);", "facets(hkl = [1, 1]);",
                 "facets(hkl = [1, 1, 1]);", "facets(hkl = [[1, 0, 0], [0, 1, 0]], d = [1]);", "rotate() sphere(1);", "scale([0, 1, 1]) sphere(1);",
                 "mirror([0, 0, 0]) sphere(1);", "linear_array(n = 1.5) sphere(1);", "grid_array(n = [0, 1, 1]) sphere(1);", "region() sphere(1);",
                 "sphere(r = sqrt(1, 2));", "sphere(r = min());", "sphere(r = bogus(1));", "sphere(r = (1);", "sphere(r = [1, 2);", "3 = 4;",
                 "capsule(h = -1);", "ellipsoid(r = [1, 0, 1]);", "polar_array(axis = [0, 0, 0]) sphere(1);", "facets(hkl = [0, 0, 0], d = 1);",
                 "facets(hkl = [1, 0, 0], d = -1);", "halfspace(normal = [0, 0, 0]);", "slab(normal = [0, 0, 0]);", "box(size = [1, -1, 1]);",
             })
            scripts.push_back(bad);
        const glm::dmat3 hexCell(glm::dvec3(4.0, 0, 0), glm::dvec3(-2.0, 3.4641016151377544, 0), glm::dvec3(0, 0, 6.5));
        Json shapes = Json::array();
        for (const auto& script : scripts)
            for (int withCell = 0; withCell < 2; ++withCell) {
                Json item = Json::object();
                item["script"] = script;
                item["cell"] = bool(withCell);
                item["result"] = attempt([&] {
                    const Model m = Model::parse(script, withCell ? &hexCell : nullptr);
                    Json r = Json::object();
                    r["primitives"] = m.primitives();
                    const Bounds b = m.bounds();
                    r["bounds"] = Json::array({b.low.x, b.low.y, b.low.z, b.high.x, b.high.y, b.high.z, b.finite, b.empty});
                    Json regions = Json::array();
                    for (const auto& e : m.regionElements()) regions.push(e);
                    r["regions"] = regions;
                    r["volume"] = m.volume(2.5);
                    Json samples = Json::array();
                    std::mt19937 points(17);
                    std::uniform_real_distribution<double> u(-25, 25);
                    for (int k = 0; k < 60; ++k) {
                        const glm::dvec3 p(u(points), u(points), u(points) * 0.6);
                        const std::string* material = m.material(p);
                        samples.push(Json::array({p.x, p.y, p.z, m.distance(p), material ? Json(*material) : Json()}));
                    }
                    r["samples"] = samples;
                    return r;
                });
                shapes.push(item);
            }
        root["shapes"] = shapes;
    }

    // Cell sculptor: supercells, slabs and cells from three slabs.
    {
        Json sculptor = Json::array();
        const std::vector<std::vector<std::array<int, 7>>> slabSets = {
            {{1, 1, 1, 1, 0, 2, 0}},
            {{1, 0, 0, 1, -1, 2, 0}, {0, 1, 0, 1, 0, 1, 0}, {0, 0, 1, 1, 0, 1, 0}},
            {{1, 1, 0, 1, 0, 1, 0}, {-1, 1, 0, 1, 0, 1, 0}, {0, 0, 1, 1, 0, 2, 0}},
            {{1, 0, 0, 0, 0, 0, 1}, {0, 1, 1, 0, 0, 0, 1}},
            {{1, 1, 1, 1, 0, 1, 0}, {1, -1, 0, 0, 0, 0, 1}, {1, 1, -2, 1, 0, 1, 0}},
            {{0, 0, 1, 1, 1, 1, 0}, {1, 0, 0, 1, 0, 1, 0}, {0, 1, 0, 1, 0, 1, 0}},
        };
        for (std::size_t i : {std::size_t(1), std::size_t(2), std::size_t(3), std::size_t(5), std::size_t(6), std::size_t(7), std::size_t(9)})
            for (std::size_t set = 0; set < slabSets.size(); ++set) {
                const Structure& source = inputs[i].second;
                std::vector<CellSlabPlane> slabs;
                for (const auto& v : slabSets[set]) {
                    CellSlabPlane slab;
                    slab.h = v[0]; slab.k = v[1]; slab.l = v[2];
                    slab.usePeriodic = v[3] != 0;
                    slab.startPlane = v[4]; slab.nPeriods = v[5];
                    if (v[6]) { slab.d1 = -3.5f; slab.d2 = 4.25f; }
                    slabs.push_back(slab);
                }
                const Structure sc = cscBuildSupercell(source, 3, 2 + int(set % 2), 3);
                Json item = Json::object();
                item["input"] = int(i);
                item["set"] = int(set);
                Json per = Json::array();
                for (const auto& slab : slabs) {
                    const glm::vec3 n = cscNormal(slab, source);
                    const glm::vec3 ref = cscSlabRef(slab, source, cscCentroid(sc));
                    per.push(Json::array({double(n.x), double(n.y), double(n.z), double(cscDhkl(slab, source)), double(cscStructuralDhkl(slab, source)),
                                          double(cscAtomLayerSpacing(slab, source)), double(cscD1(slab, source)), double(cscD2(slab, source)),
                                          double(ref.x), double(ref.y), double(ref.z)}));
                }
                item["slabs"] = per;
                item["bounded"] = cscIsBounded(slabs, source);
                item["supercell"] = structureDump(sc);
                item["cut"] = structureDump(cscApplySlabs(sc, slabs, source));
                std::vector<glm::vec3> poly = {{-5, -5, 0}, {5, -5, 0}, {5, 5, 0}, {-5, 5, 0}};
                Json clipped = Json::array();
                for (const auto& v : cscClipHalfspace(poly, glm::normalize(glm::vec3(1, 0.5f, 0)), 1.25f)) clipped.push(fvec3(v));
                item["clip"] = clipped;
                sculptor.push(item);
            }
        root["sculptor"] = sculptor;
    }

    // Mesh files.
    {
        Json meshes = Json::array();
        for (const char* name : {"quad.obj", "bad_index.obj", "empty.obj", "missing.obj", "ascii.stl", "binary.stl", "broken.stl", "missing.stl"}) {
            const std::string path = std::string("rust/tests/fixtures/mesh/") + name;
            std::vector<glm::vec3> vertices;
            std::vector<unsigned> indices;
            std::string error;
            const bool obj = std::string(name).find(".obj") != std::string::npos;
            const bool ok = obj ? parseObjMesh(path, vertices, indices, error) : parseStlMesh(path, vertices, indices, error);
            Json item = Json::object();
            item["name"] = std::string(name);
            item["ok"] = ok;
            item["error"] = error;
            Json v = Json::array(), ix = Json::array();
            for (const auto& p : vertices) v.push(fvec3(p));
            for (unsigned k : indices) ix.push(double(k));
            item["vertices"] = v;
            item["indices"] = ix;
            meshes.push(item);
        }
        root["meshes"] = meshes;
    }

    // Small tools: elastic tensors, strain, vacancies, SQS, drag, nanowires, core-shell, patterns.
    {
        Json small = Json::object();
        const auto tensorJson = [](const ElasticTensor& t) {
            Json values = Json::array();
            for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) values.push(t.c[i][j][k][l]);
            return values;
        };
        Json elastic = Json::array();
        const ElasticTensor cubic = makeCubicElasticTensor(168.4, 121.4, 75.4);
        const ElasticTensor hex = makeHexagonalElasticTensor(162.4, 92.0, 69.0, 180.7, 46.7);
        double voigt[6][6] = {};
        for (int a = 0; a < 6; ++a) for (int b = a; b < 6; ++b) voigt[a][b] = 10.0 + a * 7.5 - b * 1.25 + (a == b ? 50 : 0);
        const ElasticTensor general = makeVoigtElasticTensor(voigt);
        const ElasticTensor iso = makeCubicElasticTensor(200, 100, 50);
        const glm::dmat3 axes(glm::normalize(glm::dvec3(1, -1, 0)), glm::normalize(glm::dvec3(1, 1, -2)), glm::normalize(glm::dvec3(1, 1, 1)));
        for (const ElasticTensor* t : {&cubic, &hex, &general, &iso}) {
            Json item = Json::object();
            item["tensor"] = tensorJson(*t);
            item["rotated"] = tensorJson(rotateElasticTensor(*t, axes));
            item["isotropic"] = Json::array({isNearlyIsotropic(*t), isNearlyIsotropic(*t, 0.3), isNearlyIsotropic(rotateElasticTensor(*t, axes))});
            elastic.push(item);
        }
        small["elastic"] = elastic;

        Json strain = Json::array();
        for (std::size_t i : {std::size_t(1), std::size_t(3), std::size_t(6), std::size_t(7)}) {
            for (const auto& e : std::vector<std::array<double, 6>>{{0.01, 0, 0, 0, 0, 0}, {-0.02, 0.03, 0.01, 0.005, -0.004, 0.002}}) {
                const auto r = atomforge::applyStrain(inputs[i].second, atomforge::engineeringStrain(e[0], e[1], e[2], e[3], e[4], e[5]));
                Json item = Json::object();
                item["input"] = int(i);
                item["e"] = Json::array({e[0], e[1], e[2], e[3], e[4], e[5]});
                item["success"] = r.success;
                item["message"] = r.message;
                item["structure"] = structureDump(r.structure);
                strain.push(item);
            }
        }
        small["strain"] = strain;

        Json vacancies = Json::array();
        struct VCase { std::size_t input; const char* element; double percent; int count; double separation; unsigned seed; };
        for (const VCase& c : {VCase{0, "", 25, 0, 0, 1}, VCase{0, "Cu", 0, 5, 3.0, 7}, VCase{9, "Ni", 30, 0, 4.0, 3}, VCase{9, "", 0, 20, 2.0, 11},
                               VCase{7, "", 10, 0, 3.0, 5}, VCase{7, "", 0, 0, 0, 1}, VCase{9, "Xe", 10, 0, 0, 1}, VCase{1, "", 0, 9, 0, 1}, VCase{3, "", 200, 0, 0, 1}}) {
            atomforge::VacancyParams params;
            params.element = c.element;
            params.targetPercentage = c.percent;
            params.targetCount = c.count;
            params.minSeparation = c.separation;
            params.seed = c.seed;
            const auto r = atomforge::buildVacancies(inputs[c.input].second, params);
            Json item = Json::object();
            item["values"] = Json::array({int(c.input), std::string(c.element), c.percent, c.count, c.separation, double(c.seed)});
            item["success"] = r.success;
            item["message"] = r.message;
            item["counts"] = Json::array({r.eligibleCount, r.requestedCount, r.removedCount});
            item["structure"] = structureDump(r.structure);
            vacancies.push(item);
        }
        small["vacancies"] = vacancies;

        Json sqs = Json::array();
        struct SCase { std::size_t input; std::vector<std::pair<std::string, double>> composition; int shells, steps; double t0, t1; unsigned seed; };
        for (const SCase& c : {SCase{0, {{"Cu", 0.5}, {"Ni", 0.5}}, 2, 200, 1.0, 0.02, 1}, SCase{8, {{"Fe", 0.34}, {"Cr", 0.33}, {"Ni", 0.33}}, 3, 300, 0.5, 0.01, 4},
                               SCase{2, {{"Fe", 0.75}, {"Al", 0.25}}, 2, 150, 0, 0, 2}, SCase{0, {{"Cu", 1.0}}, 2, 10, 1, 0, 1}, SCase{0, {{"Cu", 0.5}, {"Qq", 0.5}}, 2, 10, 1, 0, 1},
                               SCase{0, {{"Cu", 0.5}, {"Ni", 0.6}}, 2, 10, 1, 0, 1}, SCase{0, {{"Cu", 0.5}, {"Ni", 0.5}}, 2, -1, 1, 0, 1}, SCase{0, {{"Cu", 0.5}, {"Ni", 0.5}}, 1, 1, 1, 0, 1}}) {
            atomforge::SQSParams params;
            for (const auto& [symbol, f] : c.composition) params.composition[symbol] = f;
            params.shells = c.shells; params.steps = c.steps; params.startTemperature = c.t0; params.endTemperature = c.t1; params.seed = c.seed;
            const auto r = atomforge::buildSQS(inputs[c.input].second, params);
            Json item = Json::object();
            Json comp = Json::array();
            for (const auto& [symbol, f] : c.composition) comp.push(Json::array({symbol, f}));
            item["composition"] = comp;
            item["values"] = Json::array({int(c.input), c.shells, c.steps, c.t0, c.t1, double(c.seed)});
            item["success"] = r.success;
            item["message"] = r.message;
            item["numbers"] = Json::array({r.initialObjective, r.finalObjective, r.acceptedSwaps});
            item["structure"] = structureDump(r.structure);
            sqs.push(item);
        }
        small["sqs"] = sqs;

        Json drag = Json::array();
        Structure moved = inputs[7].second;
        for (std::size_t a = 0; a < moved.atoms.size(); ++a) { moved.atoms[a].x += 0.3 * a - 2.0; moved.atoms[a].z -= 6.0 + 0.05 * a; }
        Structure fewer = moved;
        fewer.atoms.pop_back();
        for (const auto& [zeta, clip, pbc, which] : std::vector<std::tuple<double, bool, bool, int>>{{0.5, true, true, 0}, {0.25, false, true, 0}, {1.0, true, false, 0}, {1.5, true, true, 0}, {0.5, true, true, 1}, {0.5, true, true, 2}}) {
            const Structure& final_ = which == 1 ? fewer : moved;
            const Structure& initial = which == 2 ? inputs[3].second /* no atoms? */ : inputs[7].second;
            const auto r = atomforge::prepareDrag(initial, which == 2 ? initial : final_, zeta, clip, pbc);
            Json item = Json::object();
            item["values"] = Json::array({zeta, clip, pbc, which});
            item["success"] = r.success;
            item["message"] = r.message;
            Json dirs = Json::array();
            for (const auto& d : r.constraintDirection) dirs.push(Json::array({d.x, d.y, d.z}));
            item["directions"] = dirs;
            item["interpolated"] = structureDump(r.interpolated);
            drag.push(item);
        }
        small["drag"] = drag;
        small["drag_final"] = structureDump(moved);

        Json wires = Json::array();
        for (std::size_t i : {std::size_t(1), std::size_t(3), std::size_t(6)})
            for (const auto& [axis, radius, sides, vacuum, repeats] : std::vector<std::tuple<int, double, int, double, int>>{{2, 6.0, 0, 5.0, 1}, {0, 8.0, 6, 4.0, 2}, {1, 5.0, 4, 0.0, 3}, {3, 5, 0, 1, 1}, {2, -1, 0, 1, 1}, {2, 0.01, 0, 1, 1}}) {
                atomforge::NanowireParams params;
                params.axis = axis; params.radius = radius; params.sides = sides; params.vacuum = vacuum; params.axisRepeats = repeats;
                const auto r = atomforge::buildNanowire(inputs[i].second, params);
                Json item = Json::object();
                item["values"] = Json::array({int(i), axis, radius, sides, vacuum, repeats});
                item["success"] = r.success;
                item["message"] = r.message;
                item["count"] = r.atomCount;
                item["structure"] = structureDump(r.structure);
                wires.push(item);
            }
        small["wires"] = wires;

        Json shells = Json::array();
        for (const auto& [core, shell, radius, centroid] : std::vector<std::tuple<std::string, std::string, double, bool>>{{"Au", "Ag", 4.0, true}, {"Pt", "Pd", 3.0, false}, {"", "Ag", 4.0, true}, {"Xx", "Ag", 4.0, true}}) {
            atomforge::CoreShellParams params;
            params.coreElement = core; params.shellElement = shell; params.coreRadius = radius; params.useCentroid = centroid;
            params.center[0] = 1.0; params.center[1] = 2.0; params.center[2] = 0.5;
            const auto r = atomforge::applyCoreShell(inputs[7].second, params);
            Json item = Json::object();
            item["values"] = Json::array({core, shell, radius, centroid});
            item["success"] = r.success;
            item["message"] = r.message;
            item["counts"] = Json::array({r.coreCount, r.shellCount});
            item["structure"] = structureDump(r.structure);
            shells.push(item);
        }
        small["core_shell"] = shells;

        Json patterns = Json::array();
        for (const auto& [ref, target, cutoff, threshold, pbc] : std::vector<std::tuple<int, int, float, double, bool>>{{0, 0, 2.7f, 5.0, true}, {0, 8, 2.7f, 12.0, true}, {0, 8, 2.7f, 4.0, true}, {2, 2, 2.6f, 5.0, true}, {0, 7, 2.7f, 10.0, false}, {3, 3, 3.0f, 5.0, true}, {0, 0, -1.0f, 5.0, true}}) {
            const auto pattern = atomforge::buildPattern(inputs[ref].second, cutoff, pbc);
            const auto r = atomforge::detectPattern(inputs[target].second, pattern, cutoff, threshold, pbc);
            Json item = Json::object();
            item["values"] = Json::array({ref, target, double(cutoff), threshold, pbc});
            item["pattern"] = Json::array({pattern.success, pattern.message});
            Json dirs = Json::array();
            for (const auto& d : pattern.directions) dirs.push(Json::array({d.x, d.y, d.z}));
            item["directions"] = dirs;
            item["success"] = r.success;
            item["message"] = r.message;
            Json rows = Json::array();
            for (const auto& row : r.rows) rows.push(Json::array({row.index, row.matched, row.neighborCount, row.maxAngleDeviationDeg}));
            item["rows"] = rows;
            patterns.push(item);
        }
        small["patterns"] = patterns;
        root["small"] = small;
    }

    // Stacking-fault sequences.
    {
        Json faults = Json::array();
        struct FCase { std::size_t input; int plane; int mode; int layers; float interval; float maximum; bool pbc; };
        const std::vector<FCase> cases = {
            {1, 0, 0, 9, 0.25f, 2.0f, true}, {0, 0, 1, 6, 0.5f, 1.5f, true}, {2, 0, 0, 9, 0.25f, 1.0f, true}, {2, 6, 1, 8, 0.5f, 1.0f, true},
            {3, 0, 0, 9, 0.5f, 2.0f, true}, {3, 3, 0, 6, 0.5f, 1.0f, true}, {3, 4, 1, 6, 0.5f, 1.0f, true}, {4, 0, 0, 9, 0.5f, 2.0f, true},
            {5, 0, 0, 9, 0.5f, 2.0f, true}, {6, 0, 0, 9, 0.5f, 2.0f, true}, {7, 0, 0, 9, 0.5f, 2.0f, true}, {8, 0, 0, 9, 1.0f, 2.0f, false},
            {1, 0, 0, 1, 0.5f, 2.0f, true}, {1, 0, 0, 9, 0.0f, 2.0f, true}, {1, 0, 0, 9, 0.5f, -1.0f, true},
        };
        for (const FCase& c : cases) {
            StackingFaultParams params;
            params.plane = StackingFaultPlane(c.plane);
            params.cellMode = StackingFaultCellMode(c.mode);
            params.layerCount = c.layers;
            params.interval = c.interval;
            params.maxDisplacementFactor = c.maximum;
            params.usePbcForDetection = c.pbc;
            const StackingFaultDetectionResult d = detectStackingFaultFamily(inputs[c.input].second, c.pbc);
            const StackingFaultResult r = buildStackingFaultSequence(inputs[c.input].second, params);
            Json item = Json::object();
            item["values"] = Json::array({int(c.input), c.plane, c.mode, c.layers, double(c.interval), double(c.maximum), c.pbc});
            item["detection"] = Json::array({d.success, d.message, int(d.family), d.fccCount, d.hcpCount, d.bccCount, d.recognizedCount});
            item["success"] = r.success;
            item["message"] = r.message;
            item["numbers"] = Json::array({int(r.family), r.detectedPhase, r.detectedLayerCount, r.shiftedAtomCount, double(r.partialDisplacement)});
            Json sequence = Json::array();
            for (const auto& step : r.sequence) {
                Json e = Json::object();
                e["label"] = step.label;
                e["factor"] = double(step.displacementFactor);
                e["structure"] = builtJson(step.structure);
                sequence.push(e);
            }
            item["sequence"] = sequence;
            faults.push(item);
        }
        root["stacking_faults"] = faults;
    }

    // Dislocations: characters, shapes, dipoles, manual directions, anisotropic fields.
    {
        Json dislocations = Json::array();
        struct DCase { std::size_t input; int character; int shape; bool autoDirections; int aniso; bool dipole; bool fractional; float override_; };
        const std::vector<DCase> cases = {
            {0, 0, 0, true, 0, false, true, -1}, {0, 1, 0, true, 0, false, true, -1}, {0, 2, 1, true, 0, false, true, -1}, {0, 0, 2, true, 0, true, false, -1},
            {0, 0, 3, true, 0, false, true, 2.5f}, {0, 2, 4, true, 0, false, true, -1}, {2, 0, 0, true, 0, false, true, -1}, {2, 1, 1, true, 0, false, true, -1},
            {3, 0, 0, true, 0, false, true, -1}, {3, 1, 0, true, 0, true, true, -1}, {8, 0, 0, true, 0, false, true, -1}, {0, 0, 0, false, 0, false, true, -1},
            {0, 0, 0, true, 1, false, true, -1}, {0, 1, 1, true, 1, false, true, -1}, {0, 2, 0, true, 1, true, true, -1}, {3, 0, 0, true, 2, false, true, -1},
            {2, 2, 0, true, 3, false, true, -1}, {7, 0, 0, true, 0, false, true, -1}, {6, 0, 0, true, 0, false, true, -1}, {0, 0, 0, false, 0, false, true, -2},
        };
        for (const DCase& c : cases) {
            DislocationParams params;
            params.character = DislocationCharacter(c.character);
            params.shape = DislocationShape(c.shape);
            params.autoDirections = c.autoDirections;
            if (!c.autoDirections) {
                params.planeHkl = glm::ivec3(1, 1, 1);
                params.burgersUvw = glm::ivec3(1, 0, -1);
                params.lineUvw = c.override_ == -2 ? glm::ivec3(1, 0, -1) : glm::ivec3(1, -2, 1);
            }
            params.useFractionalLinePoint = c.fractional;
            params.linePointCartesianOffset = glm::vec3(0.1f, -0.2f, 0.05f);
            params.burgersOverrideMagnitude = c.override_ > 0 ? c.override_ : -1.0f;
            params.dipole = c.dipole;
            params.dipoleOffset = glm::vec2(4.0f, 1.0f);
            params.cutoffRadius = 12.0f;
            params.mixedCharacterAngleDeg = 30.0f;
            params.cylinderRadius = 6.0f;
            params.sphereRadius = 7.0f;
            params.ellipsoidRadii = glm::vec3(8.0f, 5.0f, 6.0f);
            if (c.aniso) {
                params.anisotropicElasticity = true;
                if (c.aniso == 1) { params.elasticSymmetry = DislocationParams::ElasticSymmetry::Cubic; params.elasticC11 = 168.4; params.elasticC12 = 121.4; params.elasticC44 = 75.4; }
                if (c.aniso == 2) { params.elasticSymmetry = DislocationParams::ElasticSymmetry::Hexagonal; params.elasticC11 = 162.4; params.elasticC12 = 92.0; params.elasticC13 = 69.0; params.elasticC33 = 180.7; params.elasticC44 = 46.7; }
                if (c.aniso == 3) {
                    params.elasticSymmetry = DislocationParams::ElasticSymmetry::General;
                    for (int a = 0; a < 6; ++a) for (int b = a; b < 6; ++b) params.elasticVoigt[a][b] = (a == b ? 200.0 + 10 * a : 0) + (a < 3 && b < 3 && a != b ? 120.0 : 0) + (a == 0 && b == 4 ? 5.0 : 0);
                }
            }
            const DislocationResult r = buildDislocation(inputs[c.input].second, params);
            Json item = Json::object();
            item["values"] = Json::array({int(c.input), c.character, c.shape, c.autoDirections, c.aniso, c.dipole, c.fractional, double(c.override_)});
            item["success"] = r.success;
            item["message"] = r.message;
            item["numbers"] = Json::array({int(r.family), r.shiftedAtomCount, double(r.burgersMagnitude)});
            item["frame"] = Json::array({fvec3(r.linePoint), fvec3(r.lineDirection), fvec3(r.burgersDirection), fvec3(r.e1), fvec3(r.e2)});
            Json loop = Json::array();
            for (const auto& p : r.loopPoints) loop.push(fvec3(p));
            item["loop"] = loop;
            const auto& v = r.validation;
            item["validation"] = Json::array({v.passed, v.finiteCoordinates, v.latticeFamilyPreserved, v.minDistanceSafe, int(v.familyBefore), int(v.familyAfter),
                                              double(v.minInteratomicDistance), double(v.rmsDisplacement), double(v.maxDisplacement), v.message});
            item["output"] = builtJson(r.output);
            dislocations.push(item);
        }
        root["dislocations"] = dislocations;
    }

    // Vitek maps, Nye tensors and fits on inserted dislocations.
    {
        Json analyses = Json::array();
        struct NCase { std::size_t input; int character; };
        for (const NCase& c : {NCase{0, 0}, NCase{2, 1}, NCase{3, 0}}) {
            DislocationParams params;
            params.character = DislocationCharacter(c.character);
            params.cutoffRadius = 12.0f;
            const DislocationResult d = buildDislocation(inputs[c.input].second, params);
            const Structure& reference = inputs[c.input].second;
            Json item = Json::object();
            item["input"] = int(c.input);
            item["deformed"] = structureDump(d.output);
            Json runs = Json::array();
            for (const auto& [cutoff, pbc, burgers] : std::vector<std::tuple<float, bool, double>>{{0.0f, true, 2.556}, {3.0f, false, 0.0}, {2.9f, true, 2.556}}) {
                Json run = Json::object();
                run["values"] = Json::array({double(cutoff), pbc, burgers});
                const auto v = atomforge::computeVitekMap(reference, d.output, glm::dvec3(d.lineDirection), burgers, cutoff, pbc);
                Json pairs = Json::array();
                for (const auto& pr : v.pairs)
                    pairs.push(Json::array({pr.indexI, pr.indexJ, pr.positionI.x, pr.positionI.y, pr.positionI.z, pr.positionJ.x, pr.positionJ.y, pr.positionJ.z,
                                            pr.screwComponent, pr.edgeComponent.x, pr.edgeComponent.y}));
                run["vitek"] = Json::array({v.success, v.message, pairs});
                const auto nye = atomforge::computeNyeTensor(reference, d.output, cutoff > 0 ? cutoff : 3.2f, pbc);
                Json alpha = Json::array();
                for (std::size_t a = 0; a < nye.alpha.size(); ++a) {
                    Json m = Json::array();
                    for (int r = 0; r < 3; ++r) for (int k = 0; k < 3; ++k) m.push(nye.alpha[a][r][k]);
                    m.push(nye.alphaNorm[a]);
                    alpha.push(m);
                }
                run["nye"] = Json::array({nye.success, nye.message, alpha});
                std::vector<glm::dvec3> positions;
                for (const auto& atom : reference.atoms) positions.push_back(glm::dvec3(atom.x, atom.y, atom.z));
                const auto fit = atomforge::fitDislocationFromNye(positions, nye, glm::dvec3(d.lineDirection), 6.25);
                run["fit"] = Json::array({fit.success, fit.message, fit.linePosition.x, fit.linePosition.y, fit.linePosition.z,
                                          fit.burgersVector.x, fit.burgersVector.y, fit.burgersVector.z, fit.totalWeight});
                runs.push(run);
            }
            item["line"] = fvec3(d.lineDirection);
            item["runs"] = runs;
            analyses.push(item);
        }
        root["dislocation_analyses"] = analyses;
    }

    // CSL grain boundaries: candidate tables and bicrystals.
    {
        Json tables = Json::array();
        for (const auto& a : std::vector<std::array<int, 3>>{{0, 0, 1}, {1, 1, 0}, {1, 1, 1}, {1, 1, 2}, {2, 1, 0}, {3, 1, 0}, {2, 2, 0}, {0, 0, 0}}) {
            const int axis[3] = {a[0], a[1], a[2]};
            Json list = Json::array();
            for (const auto& c : computeGBInfo(axis, 51)) {
                Json csl = Json::array();
                for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) csl.push(c.csl[i][j]);
                Json types = Json::array();
                for (int k = 0; k < 3; ++k) types.push(classifyBoundaryType(axis, c.plane[k].data()));
                list.push(Json::array({c.sigma, c.m, c.n, double(c.thetaDeg), csl, types}));
            }
            Json item = Json::object();
            item["axis"] = Json::array({a[0], a[1], a[2]});
            item["candidates"] = list;
            tables.push(item);
        }
        root["csl_tables"] = tables;
        Json bicrystals = Json::array();
        struct BCase { std::size_t input; std::array<int, 3> axis; int sigma; int plane; int ucA, ucB; float vacuum, gap, overlap; bool conventional; };
        for (const BCase& c : {BCase{1, {0, 0, 1}, 5, 0, 1, 1, 0, 0, 0, false}, BCase{1, {0, 0, 1}, 5, 1, 2, 1, 2.0f, 0.5f, 0.8f, true},
                               BCase{1, {1, 1, 0}, 3, 0, 1, 1, 0, 0, 0.5f, false}, BCase{1, {1, 1, 0}, 3, 2, 1, 2, 0, 0, 0, true},
                               BCase{1, {1, 1, 1}, 7, 0, 1, 1, 0, 0, 1.0f, false}, BCase{1, {1, 1, 1}, 3, 2, 1, 1, 3.0f, 0, 0, false},
                               BCase{4, {0, 0, 1}, 5, 0, 1, 1, 0, 0, 1.0f, false}, BCase{5, {1, 1, 0}, 3, 0, 1, 1, 0, 0, 0.5f, false},
                               BCase{1, {0, 0, 1}, 13, 2, 1, 1, 0, 0, 0, false}, BCase{1, {0, 0, 1}, 4, 0, 1, 1, 0, 0, 0, false}}) {
            const int axis[3] = {c.axis[0], c.axis[1], c.axis[2]};
            Json item = Json::object();
            item["values"] = Json::array({int(c.input), c.axis[0], c.axis[1], c.axis[2], c.sigma, c.plane, c.ucA, c.ucB, double(c.vacuum), double(c.gap), double(c.overlap), c.conventional});
            item["result"] = cslBicrystal(inputs[c.input].second, axis, c.sigma, 51, c.plane, c.ucA, c.ucB, c.vacuum, c.gap, c.overlap, c.conventional);
            bicrystals.push(item);
        }
        root["csl_bicrystals"] = bicrystals;
    }

    // Primitive reduction and surface slabs.
    {
        Json primitives = Json::array();
        for (std::size_t i = 0; i < inputs.size(); ++i) {
            Structure s = inputs[i].second;
            const bool ok = reduceToPrimitive(s, 1e-3);
            Json item = Json::object();
            item["input"] = int(i);
            item["ok"] = ok;
            item["structure"] = structureDump(s);
            primitives.push(item);
        }
        root["primitives"] = primitives;
        Json surfaces = Json::array();
        struct SCase { std::size_t input; int h, k, l, layers; double vacuum; int nmax; bool primitive; };
        for (const SCase& c : {SCase{1, 1, 1, 1, 4, 15.0, 8, true}, SCase{1, 1, 1, 0, 3, 10.0, 8, true}, SCase{1, 1, 0, 0, 5, 0.0, 8, false},
                               SCase{1, 2, 1, 1, 2, 12.5, 6, true}, SCase{1, 3, 2, 1, 2, 8.0, 6, false}, SCase{3, 0, 0, 1, 3, 10.0, 6, true},
                               SCase{3, 1, 0, 0, 2, 10.0, 6, false}, SCase{4, 1, 1, 1, 2, 10.0, 6, true}, SCase{5, 1, 0, 0, 2, 10.0, 6, true},
                               SCase{5, 1, 1, 1, 2, 10.0, 6, true}, SCase{2, 1, 1, 0, 2, 10.0, 4, true}, SCase{6, 1, 0, 1, 1, 5.0, 3, false},
                               SCase{1, 0, 0, 0, 2, 10.0, 6, true}, SCase{1, 1, 1, 1, 0, 10.0, 6, true}, SCase{7, 1, 1, 1, 2, 10.0, 6, true}}) {
            atomforge::SurfaceParams params;
            params.h = c.h; params.k = c.k; params.l = c.l; params.layers = c.layers; params.vacuum = c.vacuum; params.nmax = c.nmax; params.primitiveInput = c.primitive;
            const auto r = atomforge::buildSurface(inputs[c.input].second, params);
            Json item = Json::object();
            item["values"] = Json::array({int(c.input), c.h, c.k, c.l, c.layers, c.vacuum, c.nmax, c.primitive});
            item["success"] = r.success;
            item["message"] = r.message;
            item["per_layer"] = r.atomsPerLayer;
            item["structure"] = structureDump(r.structure);
            surfaces.push(item);
        }
        root["surfaces"] = surfaces;
    }

    // Interfaces: slabs, lattice matches and assembled cells.
    {
        namespace ifc = atomforge::interfaces;
        const auto matchJson = [](const ifc::Match& m) {
            return Json::array({m.a[0][0], m.a[0][1], m.a[1][0], m.a[1][1], m.b[0][0], m.b[0][1], m.b[1][0], m.b[1][1], m.thetaDeg,
                                m.strain[0], m.strain[1], m.strain[2], m.meanStrain, m.area, m.cellAngle, m.cellsA, m.cellsB, m.atoms, m.pareto});
        };
        Json cases = Json::array();
        struct ICase { std::size_t ia; bool ca; int ha, ka, la, layersA; std::size_t ib; bool cb; int hb, kb, lb, layersB; int maxRepeats; double maxStrain, maxArea, minAngle, angleMin, angleMax, angleStep; };
        for (const ICase& c : {ICase{1, true, 1, 1, 1, 3, 3, true, 0, 0, 1, 2, 4, 0.03, 80.0, 20.0, 0.0, 180.0, 1.0},
                               ICase{1, true, 1, 0, 0, 2, 2, true, 1, 1, 0, 2, 4, 0.04, 60.0, 20.0, 0.0, 90.0, 2.0},
                               ICase{5, true, 1, 0, 0, 2, 4, true, 1, 0, 0, 2, 3, 0.05, 80.0, 30.0, 0.0, 180.0, 5.0},
                               ICase{6, false, 0, 0, 1, 2, 1, true, 1, 1, 0, 2, 3, 0.05, 60.0, 20.0, 10.0, 50.0, 2.5},
                               ICase{0, true, 1, 1, 1, 2, 8, true, 1, 1, 1, 1, 3, 0.02, 50.0, 20.0, 0.0, 180.0, 3.0},
                               ICase{1, true, 1, 1, 1, 2, 9, true, 1, 1, 1, 1, 6, 0.05, 150.0, 20.0, 0.0, 180.0, 1.0},
                               ICase{3, true, 0, 0, 1, 2, 2, true, 1, 1, 0, 1, 6, 0.04, 120.0, 25.0, 0.0, 180.0, 0.5},
                               ICase{7, true, 1, 1, 1, 2, 1, true, 1, 1, 1, 2, 3, 0.03, 50.0, 20.0, 0.0, 180.0, 1.0},
                               ICase{1, true, 0, 0, 0, 2, 1, true, 1, 1, 1, 2, 3, 0.03, 50.0, 20.0, 0.0, 180.0, 1.0}}) {
            Json item = Json::object();
            item["values"] = Json::array({int(c.ia), c.ca, c.ha, c.ka, c.la, c.layersA, int(c.ib), c.cb, c.hb, c.kb, c.lb, c.layersB,
                                          c.maxRepeats, c.maxStrain, c.maxArea, c.minAngle, c.angleMin, c.angleMax, c.angleStep});
            ifc::Side sa, sb;
            sa.crystal = inputs[c.ia].second; sa.cleave = c.ca; sa.h = c.ha; sa.k = c.ka; sa.l = c.la; sa.layers = c.layersA;
            sb.crystal = inputs[c.ib].second; sb.cleave = c.cb; sb.h = c.hb; sb.k = c.kb; sb.l = c.lb; sb.layers = c.layersB;
            Structure slabA, slabB;
            std::string errA, errB;
            const bool okA = ifc::prepareSlab(sa, slabA, errA), okB = ifc::prepareSlab(sb, slabB, errB);
            item["slab_a"] = okA ? structureDump(slabA) : Json(errA);
            item["slab_b"] = okB ? structureDump(slabB) : Json(errB);
            if (okA && okB) {
                ifc::SearchParams sp;
                sp.maxRepeats = c.maxRepeats; sp.maxStrain = c.maxStrain; sp.maxArea = c.maxArea; sp.minCellAngle = c.minAngle;
                sp.angleMin = c.angleMin; sp.angleMax = c.angleMax; sp.angleStep = c.angleStep;
                const auto matches = ifc::findMatches(slabA, slabB, sp);
                Json list = Json::array();
                for (const auto& m : matches) list.push(matchJson(m));
                item["matches"] = list;
                const int chosen = ifc::suggestedMatch(matches);
                item["suggested"] = chosen;
                Json built = Json::array();
                if (!matches.empty()) {
                    const ifc::Match& m = matches[chosen >= 0 ? chosen : 0];
                    ifc::Assembly asm1;
                    built.push(structureDump(ifc::buildInterface(slabA, slabB, m, asm1)));
                    ifc::Assembly asm2; asm2.strainOn = ifc::StrainOn::A; asm2.gap = 2.5; asm2.vacuum = 10.0; asm2.shift[0] = 0.25; asm2.shift[1] = 0.5;
                    built.push(structureDump(ifc::buildInterface(slabA, slabB, m, asm2)));
                    ifc::Assembly asm3; asm3.strainOn = ifc::StrainOn::Both; asm3.gap = -1.0; asm3.shift[0] = -0.3;
                    built.push(structureDump(ifc::buildInterface(slabA, slabB, matches.back(), asm3)));
                }
                item["built"] = built;
            }
            cases.push(item);
        }
        root["interfaces"] = cases;
        // The structure helpers on their own.
        Json helpers = Json::array();
        helpers.push(structureDump(repeatLayersZ(inputs[3].second, 3)));
        helpers.push(structureDump(repeatLayersZ(inputs[1].second, 1)));
        const double F[2][2] = {{1.02, 0.1}, {-0.05, 0.97}};
        helpers.push(structureDump(applyTransform2D(inputs[6].second, F)));
        helpers.push(structureDump(repeatInterfaceXY(inputs[4].second, 2, 3)));
        helpers.push(structureDump(repeatInterfaceXY(inputs[4].second, 1, 1)));
        const int m1[2][2] = {{2, 1}, {-1, 1}};
        helpers.push(structureDump(ifc::inPlaneSupercell(inputs[6].second, m1)));
        const int m2[2][2] = {{1, 1}, {1, 1}};
        helpers.push(structureDump(ifc::inPlaneSupercell(inputs[1].second, m2)));
        root["interface_helpers"] = helpers;
    }

    // Polycrystals.
    {
        const std::vector<glm::vec3> colors = makeDefaultElementColors();
        Json cases = Json::array();
        struct PCase { std::size_t input; float sx, sy, sz; int grains, seed, mode; int specCount; };
        for (const PCase& c : {PCase{1, 20.0f, 18.0f, 15.0f, 4, 42, 0, 0}, PCase{3, 16.0f, 16.0f, 16.0f, 3, 7, 0, 0}, PCase{5, 17.3f, 12.0f, 11.1f, 5, 123, 1, 3},
                               PCase{4, 14.0f, 14.0f, 14.0f, 3, 5, 2, 3}, PCase{1, 7.22f, 7.22f, 7.22f, 1, 1, 1, 0}, PCase{6, 15.0f, 12.0f, 10.0f, 6, 99, 0, 0},
                               PCase{9, 12.0f, 12.0f, 12.0f, 2, 3, 2, 1}, PCase{7, 10.0f, 10.0f, 10.0f, 2, 1, 0, 0}, PCase{1, 10.0f, 10.0f, 10.0f, 0, 1, 0, 0},
                               PCase{1, 10.0f, -1.0f, 10.0f, 2, 1, 0, 0}, PCase{0, 2000.0f, 2000.0f, 2000.0f, 2, 1, 0, 0}}) {
            PolyParams params;
            params.sizeX = c.sx; params.sizeY = c.sy; params.sizeZ = c.sz; params.numGrains = c.grains; params.seed = c.seed;
            params.orientationMode = static_cast<GrainOrientationMode>(c.mode);
            const int indices[3] = {0, 7, 2};
            for (int k = 0; k < c.specCount; ++k) {
                GrainOrientation o;
                o.grainIndex = indices[k];
                o.phi1 = 10.0f + 37.5f * k; o.Phi = 54.7356f - 11.0f * k; o.phi2 = 45.0f + 90.0f * k;
                params.specifiedOrientations.push_back(o);
            }
            Structure out;
            const PolyBuildResult r = buildPolycrystal(out, inputs[c.input].second, params, colors);
            Json item = Json::object();
            item["values"] = Json::array({int(c.input), double(c.sx), double(c.sy), double(c.sz), c.grains, c.seed, c.mode, c.specCount});
            item["result"] = Json::array({r.success, r.message, r.inputAtoms, r.outputAtoms, r.numGrains});
            item["structure"] = structureDump(out);
            cases.push(item);
        }
        root["polycrystals"] = cases;
    }

    // Nanocrystals: shapes, orientations, a mesh fill and Wulff polyhedra.
    {
        const std::vector<glm::vec3> colors = makeDefaultElementColors();
        const std::vector<glm::vec3> cube = {{-6, -6, -6}, {6, -6, -6}, {6, 6, -6}, {-6, 6, -6}, {-6, -6, 6}, {6, -6, 6}, {6, 6, 6}, {-6, 6, 6}};
        const std::vector<unsigned int> faces = {0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6, 3, 0, 4, 3, 4, 7, 0, 1, 99};
        const auto v3 = [](const glm::vec3& v) { return Json::array({double(v.x), double(v.y), double(v.z)}); };
        const std::vector<std::vector<double>> table = {
            {1.0, 0.0, 0.0, 8.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 5.0, 1.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 20.0, 1.0},
            {1.0, 0.0, 1.0, 9.0, 7.0, 5.0, 2.0, 1.0, 10.0, 20.0, 30.0, 1.0, 0.0, 0.0, 1.0, 5.0, 1.0, 3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 20.0, 1.0},
            {2.0, 0.0, 2.0, 6.0, 5.0, 4.0, 2.0, 2.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 5.0, 1.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 20.0, 1.0},
            {3.0, 0.0, 3.0, 6.0, 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 5.0, 1.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 20.0, 1.0},
            {4.0, 0.0, 4.0, 9.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 5.0, 1.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 20.0, 1.0},
            {5.0, 0.0, 5.0, 10.0, 6.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 5.0, 0.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 20.0, 1.0},
            {1.0, 0.0, 6.0, 8.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 3.0, 1.0, 2.5, 0.0, 1.0, 2.0, 3.0, 0.0, 20.0, 1.0},
            {7.0, 0.0, 0.0, 5.0, 0.0, 0.0, 2.0, 1.0, 15.0, 0.0, 45.0, 1.0, 0.0, 0.0, 1.0, 5.0, 1.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 20.0, 1.0},
            {1.0, 0.0, 7.0, 5.0, 5.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 5.0, 1.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 20.0, 1.2},
            {1.0, 0.0, 7.0, 5.0, 5.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 5.0, 1.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 20.0, 0.0},
            {6.0, 0.0, 3.0, 5.0, 9.0, 0.0, 1.0, 1.0, -15.0, 40.0, 75.0, 1.0, 0.0, 0.0, 1.0, 5.0, 1.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 20.0, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 5.0, 1.0, 5.0, 1.0, 0.0, 0.0, 0.0, 1.0, 12.0, 1.0},
            {3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0, 5.0, 10.0, 15.0, 1.0, 0.0, 0.0, 1.0, 5.0, 1.0, 5.0, 1.0, 0.0, 0.0, 0.0, 2.0, 10.0, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 5.0, 0.0, 5.0, 1.0, 0.0, 0.0, 0.0, 3.0, 9.0, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 5.0, 1.0, 5.0, 1.0, 0.0, 0.0, 0.0, 4.0, 9.0, 1.0},
            {7.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 5.0, 1.0, 5.0, 1.0, 0.0, 0.0, 0.0, 1.0, 9.0, 1.0},
            {1.0, 0.0, 0.0, 8.0, 0.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 5.0, 1.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 20.0, 1.0},
            {9.0, 0.0, 0.0, 500.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 5.0, 1.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 20.0, 1.0},
            {4.0, 0.0, 0.0, 0.1, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 5.0, 1.0, 5.0, 0.0, 0.3, 0.3, 0.3, 0.0, 20.0, 1.0},
            {2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0, 1.0, 4.0, 1.0, 0.0, 0.0, 0.0, 1.0, 7.0, 1.0}
        };
        Json cases = Json::array();
        for (const auto& c : table) {
            NanoParams p;
            p.generationMode = static_cast<NanoGenerationMode>(int(c[1]));
            p.shape = static_cast<NanoShape>(int(c[2]));
            const float s1 = float(c[3]), s2 = float(c[4]), s3 = float(c[5]);
            p.sphereRadius = s1; p.ellipRx = s1; p.ellipRy = s2; p.ellipRz = s3; p.boxHx = s1; p.boxHy = s2; p.boxHz = s3;
            p.cylRadius = s1; p.cylHeight = s2; p.octRadius = s1; p.truncOctRadius = s1; p.truncOctTrunc = s2; p.cuboRadius = s1;
            p.modelHx = s1; p.modelHy = s2; p.modelHz = s3;
            p.cylAxis = int(c[6]);
            p.applyCrystalOrientation = c[7] != 0; p.useMillerOrientation = c[7] == 2;
            p.orientXDeg = float(c[8]); p.orientYDeg = float(c[9]); p.orientZDeg = float(c[10]);
            p.millerH = int(c[11]); p.millerK = int(c[12]); p.millerL = int(c[13]);
            p.autoReplicate = c[14] != 0; p.repA = int(c[15]); p.repB = int(c[15]) + 1; p.repC = int(c[15]);
            p.setOutputCell = c[16] != 0; p.vacuumPadding = float(c[17]);
            p.autoCenterFromAtoms = c[18] != 0; p.cx = float(c[19]); p.cy = float(c[20]); p.cz = float(c[21]);
            const int set = int(c[22]);
            if (set == 1) p.wulffPlanes = {{1, 0, 0, 1.0f}, {1, 1, 1, 0.8f}};
            if (set == 2) p.wulffPlanes = {{0, 0, 1, 1.0f}, {1, 0, 0, 1.1f}, {1, 0, 1, 1.2f}};
            if (set == 3) p.wulffPlanes = {{1, 0, 0, 1.0f}, {2, 2, 0, 1.2f}, {1, 1, 1, 0.9f}};
            if (set == 4) p.wulffPlanes = {{1, 0, 0, 1.0f}, {1, 1, 1, 0.0f}};
            p.wulffMaxRadius = float(c[23]);
            p.modelScale = float(c[24]);
            Json values = Json::array();
            for (double v : c) values.push(v);
            Json item = Json::object();
            item["values"] = values;
            item["bounding"] = Json::array({double(computeBoundingRadius(p)), double(computeShapeHalfExtents(p).hx), double(computeShapeHalfExtents(p).hy), double(computeShapeHalfExtents(p).hz)});
            if (p.generationMode == NanoGenerationMode::WulffConstruction) {
                const WulffPreview w = computeWulffPreview(inputs[std::size_t(c[0])].second, p);
                Json planes = Json::array();
                for (const auto& pl : w.planes) planes.push(Json::array({pl.familyIndex, pl.h, pl.k, pl.l, double(pl.distance), v3(pl.normal)}));
                Json wf = Json::array();
                for (const auto& f : w.faces) {
                    Json verts = Json::array();
                    for (const auto& v : f.vertices) verts.push(v3(v));
                    wf.push(Json::array({f.familyIndex, f.h, f.k, f.l, double(f.distance), v3(f.normal), verts}));
                }
                item["preview"] = Json::array({w.success, w.message, v3(w.center), v3(w.minPoint), v3(w.maxPoint), double(w.maxPlaneDistance), double(w.boundingRadius), planes, wf});
            }
            Structure out;
            const bool mesh = p.shape == NanoShape::MeshModel && p.modelScale > 0;
            const NanoBuildResult r = buildNanocrystal(out, inputs[std::size_t(c[0])].second, p, colors, mesh ? cube : std::vector<glm::vec3>{}, mesh ? faces : std::vector<unsigned int>{});
            item["result"] = Json::array({r.success, r.message, int(r.mode), r.inputAtoms, r.outputAtoms, int(r.shape), double(r.estimatedDiameter), r.repA, r.repB, r.repC, r.tilingUsed, r.repClamped, r.wulffFamilyCount, r.wulffFaceCount});
            item["structure"] = structureDump(out);
            cases.push(item);
        }
        root["nanocrystals"] = cases;
        Json labels = Json::array();
        for (int i = 0; i < kNumShapes; ++i) labels.push(std::string(shapeLabel(static_cast<NanoShape>(i))));
        root["nano_labels"] = labels;
    }

    std::ofstream(std::string(argv[1]) + "/algorithms.json") << root.dump(1);
    std::printf("wrote %s/algorithms.json\n", argv[1]);
    return 0;
}
