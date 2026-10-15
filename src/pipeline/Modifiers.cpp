// The modifier library. Each entry describes its parameters (shown by the
// desktop pipeline panel and accepted by the text syntax) and applies itself
// to the data flowing through the pipeline.
#include "algorithms/InterstitialVoidAnalysis.h"
#include "pipeline/Expression.h"
#include "pipeline/Pipeline.h"
#include "science/ScienceCatalog.h"
#include "science/ScienceCore.h"
#include "science/ScienceData.h"
#include "science/ScienceTools.h"
#include "util/ElementData.h"

#include <algorithm>
#include <deque>
#include <filesystem>
#include <cmath>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>

namespace atomforge::pipeline
{
namespace
{
using science::Mat3;
using science::Vec3;

// ---------------------------------------------------------------- parameter access
const ModifierType& self(const char* id) { return *findModifierType(id); }
double number(const Json& p, const char* id, const char* name) { return parameter(p, self(id), name).number(); }
bool flag(const Json& p, const char* id, const char* name) { return parameter(p, self(id), name).boolean(); }
std::string text(const Json& p, const char* id, const char* name) { return parameter(p, self(id), name).string(); }
Vec3 vector(const Json& p, const char* id, const char* name)
{
    const Json v = parameter(p, self(id), name);
    return {v.items()[0].number(), v.items()[1].number(), v.items()[2].number()};
}

Mat3 cellOf(const Structure& s) { Mat3 m{}; for (int r = 0; r < 3; ++r) m[r] = s.cellVectors[r]; return m; }
Vec3 positionOf(const AtomSite& a) { return {a.x, a.y, a.z}; }
void place(AtomSite& a, const Vec3& x) { a.x = x[0]; a.y = x[1]; a.z = x[2]; }

void requireCell(const PipelineData& data, const char* what)
{
    if (!data.structure.hasUnitCell) throw std::runtime_error(std::string(what) + " needs a structure with a unit cell");
}

// Keeps atoms by index (in the given order), carrying aligned per-atom metadata.
void rebuild(PipelineData& data, const std::vector<std::size_t>& keep)
{
    Structure& s = data.structure;
    const std::size_t n = s.atoms.size();
    Structure out = s;
    out.atoms.clear();
    out.grainColors.clear();
    out.grainRegionIds.clear();
    out.atomProperty.clear();
    std::vector<char> selected;
    for (std::size_t i : keep) {
        out.atoms.push_back(s.atoms[i]);
        if (s.grainColors.size() == n) out.grainColors.push_back(s.grainColors[i]);
        if (s.grainRegionIds.size() == n) out.grainRegionIds.push_back(s.grainRegionIds[i]);
        if (s.atomProperty.size() == n) out.atomProperty.push_back(s.atomProperty[i]);
        selected.push_back(i < data.selected.size() ? data.selected[i] : 0);
    }
    out.dislocationLoopPoints.clear();
    out.dislocationDetectionDone = false;
    data.structure = std::move(out);
    data.selected = std::move(selected);
}

// Indices the operation applies to: the selection, or every atom.
std::vector<std::size_t> targets(const PipelineData& data, bool selectedOnly, const char* what)
{
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < data.structure.atoms.size(); ++i)
        if (!selectedOnly || data.selected[i]) out.push_back(i);
    if (selectedOnly && out.empty()) throw std::runtime_error(std::string(what) + ": no atoms are selected; add a selection modifier above");
    return out;
}

// Combines a new selection with the current one.
void select(PipelineData& data, const std::vector<char>& mask, const std::string& mode)
{
    for (std::size_t i = 0; i < mask.size(); ++i) {
        char& s = data.selected[i];
        if (mode == "replace") s = mask[i];
        else if (mode == "add") s = s || mask[i];
        else if (mode == "subtract") s = s && !mask[i];
        else s = s && mask[i];  // intersect
    }
    data.notes.push_back(std::to_string(data.selectedCount()) + " of " + std::to_string(mask.size()) + " atoms selected");
}

Vec3 minimumImage(const Structure& s, Vec3 d)
{
    if (!s.hasUnitCell) return d;
    const Mat3 cell = cellOf(s);
    Vec3 f = science::fractional(d, cell);
    for (double& v : f) v -= std::round(v);
    return science::rowTimes(f, cell);
}

const char* kModes = "replace|add|subtract|intersect";
#define MODE {"mode", "Combine with the current selection", "choice", "replace", kModes}

// An element symbol in any case ("cu", "CU" -> "Cu"); throws for an unknown one.
std::string standardElement(const std::string& name)
{
    std::string lower = name;
    for (char& c : lower) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    for (int z = 1; z <= 118; ++z) {
        std::string symbol = elementSymbol(z);
        for (char& c : symbol) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (symbol == lower) return elementSymbol(z);
    }
    throw std::runtime_error("Unknown element " + name);
}

// ---------------------------------------------------------------- selection
void selectElement(const Json& p, PipelineData& data)
{
    std::string list = text(p, "select-element", "elements");
    std::replace(list.begin(), list.end(), ',', ' ');
    std::istringstream in(list);
    std::vector<std::string> elements;
    for (std::string e; in >> e;) elements.push_back(standardElement(e));
    if (elements.empty()) throw std::runtime_error("Name at least one element");
    std::vector<char> mask(data.structure.atoms.size(), 0);
    for (std::size_t i = 0; i < mask.size(); ++i)
        mask[i] = std::find(elements.begin(), elements.end(), data.structure.atoms[i].symbol) != elements.end();
    select(data, mask, text(p, "select-element", "mode"));
}

// The value of an expression for every atom.
std::vector<double> evaluateForAtoms(const Expression& expression, const PipelineData& data)
{
    const Structure& s = data.structure;
    const Mat3 inv = s.hasUnitCell ? science::inverse(cellOf(s)) : Mat3{};
    std::vector<double> values(s.atoms.size());
    AtomVariables v;
    v.count = static_cast<double>(s.atoms.size());
    for (std::size_t i = 0; i < s.atoms.size(); ++i) {
        const AtomSite& a = s.atoms[i];
        v.x = a.x; v.y = a.y; v.z = a.z;
        if (s.hasUnitCell) {
            const Vec3 f = science::rowTimes(positionOf(a), inv);
            v.fx = f[0]; v.fy = f[1]; v.fz = f[2];
        }
        v.index = static_cast<double>(i);
        v.property = s.atomProperty.size() == s.atoms.size() ? s.atomProperty[i] : std::numeric_limits<double>::quiet_NaN();
        v.selected = data.selected[i] ? 1 : 0;
        v.atomicNumber = a.atomicNumber;
        v.fixed = static_cast<double>((a.fixed & 1) + ((a.fixed >> 1) & 1) + ((a.fixed >> 2) & 1));
        v.element = a.symbol;
        values[i] = expression.evaluate(v);
    }
    return values;
}

void selectExpression(const Json& p, PipelineData& data)
{
    const Expression expression(text(p, "select-expression", "expression"));
    const auto values = evaluateForAtoms(expression, data);
    std::vector<char> mask(values.size(), 0);
    for (std::size_t i = 0; i < values.size(); ++i) mask[i] = values[i] != 0;
    select(data, mask, text(p, "select-expression", "mode"));
}

// Pairs of atoms closer than a distance: both atoms, or only the later one
// of each pair (so deleting the selection keeps one atom per overlap).
void selectOverlapping(const Json& p, PipelineData& data)
{
    const double distance = number(p, "select-overlapping", "distance");
    if (!(distance > 0)) throw std::runtime_error("select-overlapping needs a positive distance");
    const bool both = text(p, "select-overlapping", "atoms") == "both";
    const Structure& s = data.structure;
    std::vector<Vec3> positions;
    for (const auto& a : s.atoms) positions.push_back(positionOf(a));
    const science::Pbc pbc{s.hasUnitCell, s.hasUnitCell, s.hasUnitCell};
    std::vector<char> mask(s.atoms.size(), 0);
    std::size_t pairs = 0;
    for (const auto& n : science::neighborList(positions, s.hasUnitCell ? cellOf(s) : science::identity(), pbc, distance)) {
        if (n.i >= n.j) continue;  // each pair once; skips an atom's own periodic images
        ++pairs;
        mask[static_cast<std::size_t>(n.j)] = 1;
        if (both) mask[static_cast<std::size_t>(n.i)] = 1;
    }
    data.notes.push_back(std::to_string(pairs) + " overlapping pairs");
    select(data, mask, text(p, "select-overlapping", "mode"));
}

void selectSlab(const Json& p, PipelineData& data)
{
    const std::string axis = text(p, "select-slab", "axis");
    const double low = number(p, "select-slab", "minimum"), high = number(p, "select-slab", "maximum");
    if (!(low <= high)) throw std::runtime_error("The slab minimum must not exceed its maximum");
    const bool fractional = axis == "a" || axis == "b" || axis == "c";
    if (fractional) requireCell(data, "A fractional slab");
    const int k = (axis == "x" || axis == "a") ? 0 : (axis == "y" || axis == "b") ? 1 : 2;
    const Mat3 inv = fractional ? science::inverse(cellOf(data.structure)) : Mat3{};
    std::vector<char> mask(data.structure.atoms.size(), 0);
    for (std::size_t i = 0; i < mask.size(); ++i) {
        const Vec3 x = positionOf(data.structure.atoms[i]);
        double value = fractional ? science::rowTimes(x, inv)[static_cast<std::size_t>(k)] : x[static_cast<std::size_t>(k)];
        if (fractional) value -= std::floor(value);
        mask[i] = value >= low && value <= high;
    }
    select(data, mask, text(p, "select-slab", "mode"));
}

void selectSphere(const Json& p, PipelineData& data)
{
    const Vec3 centre = vector(p, "select-sphere", "center");
    const double radius = number(p, "select-sphere", "radius");
    if (!(radius > 0)) throw std::runtime_error("The sphere radius must be positive");
    std::vector<char> mask(data.structure.atoms.size(), 0);
    for (std::size_t i = 0; i < mask.size(); ++i)
        mask[i] = science::norm(minimumImage(data.structure, science::sub(positionOf(data.structure.atoms[i]), centre))) <= radius;
    select(data, mask, text(p, "select-sphere", "mode"));
}

void selectProperty(const Json& p, PipelineData& data)
{
    const Structure& s = data.structure;
    if (s.atomProperty.size() != s.atoms.size() || s.atoms.empty())
        throw std::runtime_error("No per-atom property: add compute-property above (or colour atoms by a result first)");
    const double low = number(p, "select-property", "minimum"), high = number(p, "select-property", "maximum");
    std::vector<char> mask(s.atoms.size(), 0);
    for (std::size_t i = 0; i < mask.size(); ++i) mask[i] = std::isfinite(s.atomProperty[i]) && s.atomProperty[i] >= low && s.atomProperty[i] <= high;
    select(data, mask, text(p, "select-property", "mode"));
}

void selectRandom(const Json& p, PipelineData& data)
{
    const double fraction = number(p, "select-random", "fraction");
    if (!(fraction >= 0 && fraction <= 1)) throw std::runtime_error("The fraction must be between 0 and 1");
    // Exactly round(fraction * N) atoms, reproducible for a seed.
    const std::size_t n = data.structure.atoms.size();
    std::vector<std::size_t> order(n);
    for (std::size_t i = 0; i < n; ++i) order[i] = i;
    std::mt19937_64 random(static_cast<std::uint64_t>(number(p, "select-random", "seed")));
    std::shuffle(order.begin(), order.end(), random);
    std::vector<char> mask(n, 0);
    for (std::size_t k = 0; k < static_cast<std::size_t>(std::llround(fraction * static_cast<double>(n))); ++k) mask[order[k]] = 1;
    select(data, mask, text(p, "select-random", "mode"));
}

void expandSelection(const Json& p, PipelineData& data)
{
    const double cutoff = number(p, "expand-selection", "cutoff");
    const long long steps = static_cast<long long>(number(p, "expand-selection", "steps"));
    if (!(cutoff > 0) || steps < 1) throw std::runtime_error("expand-selection needs a positive cutoff and at least one step");
    const Structure& s = data.structure;
    std::vector<Vec3> positions;
    for (const auto& a : s.atoms) positions.push_back(positionOf(a));
    const science::Pbc pbc{s.hasUnitCell, s.hasUnitCell, s.hasUnitCell};
    const auto neighbours = science::neighborList(positions, s.hasUnitCell ? cellOf(s) : science::identity(), pbc, cutoff);
    for (long long step = 0; step < steps; ++step) {
        std::vector<char> grown = data.selected;
        for (const auto& n : neighbours)
            if (data.selected[static_cast<std::size_t>(n.i)]) grown[static_cast<std::size_t>(n.j)] = 1;
        data.selected = std::move(grown);
    }
    data.notes.push_back(std::to_string(data.selectedCount()) + " atoms selected after expansion");
}

void invertSelection(const Json&, PipelineData& data)
{
    for (char& s : data.selected) s = !s;
    data.notes.push_back(std::to_string(data.selectedCount()) + " atoms selected");
}

void clearSelection(const Json&, PipelineData& data) { data.clearSelection(); }

// ---------------------------------------------------------------- modification
void deleteSelected(const Json&, PipelineData& data)
{
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < data.structure.atoms.size(); ++i)
        if (!data.selected[i]) keep.push_back(i);
    const std::size_t removed = data.structure.atoms.size() - keep.size();
    rebuild(data, keep);
    data.notes.push_back(std::to_string(removed) + " atoms deleted");
}

// Selective dynamics: hold coordinates of the target atoms fixed in a
// relaxation (written to POSCAR as F flags), or free them again.
void fixAtoms(const Json& p, PipelineData& data)
{
    const std::string axes = text(p, "fix-atoms", "axes");
    std::uint8_t mask = 0;
    for (char c : axes) mask = static_cast<std::uint8_t>(mask | (c == 'x' ? 1 : c == 'y' ? 2 : c == 'z' ? 4 : 0));
    const auto indices = targets(data, text(p, "fix-atoms", "target") == "selected", "fix-atoms");
    for (std::size_t i : indices) data.structure.atoms[i].fixed = static_cast<std::uint8_t>(data.structure.atoms[i].fixed | mask);
    data.notes.push_back(std::to_string(indices.size()) + " atoms fixed along " + axes);
}

void freeAtoms(const Json& p, PipelineData& data)
{
    const auto indices = targets(data, text(p, "free-atoms", "target") == "selected", "free-atoms");
    for (std::size_t i : indices) data.structure.atoms[i].fixed = 0;
    data.notes.push_back(std::to_string(indices.size()) + " atoms free to move");
}

void assignElement(const Json& p, PipelineData& data)
{
    const std::string symbol = standardElement(text(p, "assign-element", "element"));
    const int z = science::atomicNumber(symbol);
    float r = 1, g = 1, b = 1;
    getDefaultElementColor(z, r, g, b);
    const auto indices = targets(data, text(p, "assign-element", "target") == "selected", "assign-element");
    for (std::size_t i : indices) {
        AtomSite& a = data.structure.atoms[i];
        a.symbol = symbol; a.atomicNumber = z; a.r = r; a.g = g; a.b = b;
    }
    data.notes.push_back(std::to_string(indices.size()) + " atoms changed to " + symbol);
}

void displace(const Json& p, PipelineData& data)
{
    const Vec3 shift = vector(p, "displace", "vector");
    for (std::size_t i : targets(data, text(p, "displace", "target") == "selected", "displace"))
        place(data.structure.atoms[i], science::add(positionOf(data.structure.atoms[i]), shift));
}

void randomDisplacement(const Json& p, PipelineData& data)
{
    const double amplitude = number(p, "random-displacement", "amplitude");
    if (!(amplitude >= 0)) throw std::runtime_error("The amplitude must not be negative");
    std::mt19937_64 random(static_cast<std::uint64_t>(number(p, "random-displacement", "seed")));
    std::normal_distribution<double> gaussian(0.0, amplitude);
    for (std::size_t i : targets(data, text(p, "random-displacement", "target") == "selected", "random-displacement"))
        place(data.structure.atoms[i], science::add(positionOf(data.structure.atoms[i]), {gaussian(random), gaussian(random), gaussian(random)}));
}

void slice(const Json& p, PipelineData& data)
{
    Vec3 normal = vector(p, "slice", "normal");
    const double length = science::norm(normal);
    if (!(length > 0)) throw std::runtime_error("The slice normal must not be zero");
    normal = science::scale(normal, 1 / length);
    const double distance = number(p, "slice", "distance"), width = number(p, "slice", "width");
    const bool invert = flag(p, "slice", "invert");
    // Without a width the plane cuts away the half-space above it; with a width a slab is cut.
    std::vector<char> inside(data.structure.atoms.size(), 0);
    for (std::size_t i = 0; i < inside.size(); ++i) {
        const double d = science::dot(positionOf(data.structure.atoms[i]), normal) - distance;
        const bool hit = width > 0 ? std::abs(d) <= 0.5 * width : d > 0;
        inside[i] = hit != invert;
    }
    if (text(p, "slice", "action") == "select") { select(data, inside, "replace"); return; }
    std::vector<std::size_t> keep;
    for (std::size_t i = 0; i < inside.size(); ++i) if (!inside[i]) keep.push_back(i);
    const std::size_t removed = inside.size() - keep.size();
    rebuild(data, keep);
    data.notes.push_back(std::to_string(removed) + " atoms removed by the slice");
}

void addAtom(const Json& p, PipelineData& data)
{
    const std::string symbol = standardElement(text(p, "add-atom", "element"));
    const int z = science::atomicNumber(symbol);
    Vec3 position = vector(p, "add-atom", "position");
    if (flag(p, "add-atom", "fractional")) {
        requireCell(data, "A fractional position");
        position = science::add(science::rowTimes(position, cellOf(data.structure)), data.structure.cellOffset);
    }
    AtomSite atom;
    atom.symbol = symbol;
    atom.atomicNumber = z;
    place(atom, position);
    getDefaultElementColor(z, atom.r, atom.g, atom.b);
    Structure& s = data.structure;
    const std::size_t n = s.atoms.size();
    // New atoms extend the per-atom metadata (no grain colour, region -1, no property).
    if (s.grainColors.size() == n) s.grainColors.push_back({0.6f, 0.6f, 0.6f});
    if (s.grainRegionIds.size() == n) s.grainRegionIds.push_back(-1);
    if (s.atomProperty.size() == n && n) s.atomProperty.push_back(std::numeric_limits<double>::quiet_NaN());
    s.atoms.push_back(atom);
    data.selected.push_back(flag(p, "add-atom", "select") ? 1 : 0);
}

void insertInterstitials(const Json& p, PipelineData& data)
{
    const std::string symbol = standardElement(text(p, "insert-interstitials", "element"));
    const int z = science::atomicNumber(symbol);
    InterstitialVoidDetectionParams params;
    params.gridResolution = static_cast<int>(number(p, "insert-interstitials", "resolution"));
    params.minClearance = static_cast<float>(number(p, "insert-interstitials", "clearance"));
    params.minSeparation = static_cast<float>(number(p, "insert-interstitials", "separation"));
    const auto detection = detectInterstitialVoidRegions(data.structure, params);
    if (!detection.success) throw std::runtime_error(detection.message.empty() ? "No interstitial sites were found" : detection.message);
    // Candidate sites of the requested kind (the same voids as Edit > Add Interstitial Atoms).
    const std::string kind = text(p, "insert-interstitials", "kind");
    std::vector<const InterstitialVoidRegion*> sites;
    for (const auto& region : detection.regions) {
        const bool wanted = kind == "any" || (kind == "tetrahedral" && region.kind == InterstitialVoidKind::Tetrahedral) ||
                            (kind == "octahedral" && region.kind == InterstitialVoidKind::Octahedral) ||
                            (kind == "irregular" && region.kind == InterstitialVoidKind::Irregular);
        if (wanted) sites.push_back(&region);
    }
    if (sites.empty()) throw std::runtime_error("No " + kind + " interstitial sites in this structure");
    const long long count = static_cast<long long>(number(p, "insert-interstitials", "count"));
    if (count < 0) throw std::runtime_error("The count must not be negative");
    if (text(p, "insert-interstitials", "order") == "random") {
        std::mt19937_64 random(static_cast<std::uint64_t>(number(p, "insert-interstitials", "seed")));
        std::shuffle(sites.begin(), sites.end(), random);
    } else {
        std::stable_sort(sites.begin(), sites.end(), [](const auto* a, const auto* b) { return a->clearance > b->clearance; });
    }
    const std::size_t n = count == 0 ? sites.size() : std::min(sites.size(), static_cast<std::size_t>(count));
    for (std::size_t k = 0; k < n; ++k) {
        AtomSite atom;
        atom.symbol = symbol;
        atom.atomicNumber = z;
        place(atom, {sites[k]->position.x, sites[k]->position.y, sites[k]->position.z});
        getDefaultElementColor(z, atom.r, atom.g, atom.b);
        Structure& s = data.structure;
        const std::size_t before = s.atoms.size();
        if (s.grainColors.size() == before) s.grainColors.push_back({0.6f, 0.6f, 0.6f});
        if (s.grainRegionIds.size() == before) s.grainRegionIds.push_back(-1);
        if (s.atomProperty.size() == before && before) s.atomProperty.push_back(std::numeric_limits<double>::quiet_NaN());
        s.atoms.push_back(atom);
        data.selected.push_back(1);
    }
    data.notes.push_back(std::to_string(n) + " of " + std::to_string(sites.size()) + " " + kind + " sites filled with " + symbol);
}

void merge(const Json& p, PipelineData& data)
{
    const std::string file = text(p, "merge", "file");
    if (file.empty()) throw std::runtime_error("Choose the structure file to merge");
    const auto frames = science::readFrames(std::filesystem::u8path(file));
    const Vec3 offset = vector(p, "merge", "offset");
    Structure& s = data.structure;
    const std::size_t n = s.atoms.size();
    const bool colours = s.grainColors.size() == n, regions = s.grainRegionIds.size() == n, property = s.atomProperty.size() == n && n;
    const Structure& other = frames.front().structure;
    for (AtomSite atom : other.atoms) {
        place(atom, science::add(positionOf(atom), offset));
        s.atoms.push_back(atom);
        if (colours) s.grainColors.push_back({0.6f, 0.6f, 0.6f});
        if (regions) s.grainRegionIds.push_back(-1);
        if (property) s.atomProperty.push_back(std::numeric_limits<double>::quiet_NaN());
        data.selected.push_back(flag(p, "merge", "select") ? 1 : 0);
    }
    if (!s.hasUnitCell && other.hasUnitCell && flag(p, "merge", "use_cell")) { s.hasUnitCell = true; s.cellVectors = other.cellVectors; }
    data.notes.push_back(std::to_string(other.atoms.size()) + " atoms merged from " + std::filesystem::u8path(file).filename().u8string());
}

// ---------------------------------------------------------------- cell
void setCell(const Json& p, PipelineData& data)
{
    Mat3 cell{};
    cell[0] = vector(p, "set-cell", "a");
    cell[1] = vector(p, "set-cell", "b");
    cell[2] = vector(p, "set-cell", "c");
    if (std::abs(science::determinant(cell)) < 1e-9) throw std::runtime_error("The cell vectors are coplanar");
    Structure& s = data.structure;
    if (flag(p, "set-cell", "scale_atoms")) {
        requireCell(data, "Scaling atoms with the cell");
        const Mat3 inv = science::inverse(cellOf(s));
        for (auto& a : s.atoms) place(a, science::rowTimes(science::rowTimes(science::sub(positionOf(a), s.cellOffset), inv), cell));
        s.cellOffset = {0, 0, 0};
    }
    for (int r = 0; r < 3; ++r) s.cellVectors[r] = cell[r];
    s.hasUnitCell = true;
}

void replicate(const Json& p, PipelineData& data)
{
    requireCell(data, "replicate");
    const Vec3 counts = vector(p, "replicate", "counts");
    long long n[3];
    for (int k = 0; k < 3; ++k) {
        n[k] = std::llround(counts[static_cast<std::size_t>(k)]);
        if (n[k] < 1 || n[k] > 100) throw std::runtime_error("Replication counts must be between 1 and 100");
    }
    if (static_cast<double>(n[0] * n[1] * n[2]) * static_cast<double>(data.structure.atoms.size()) > 50e6) throw std::runtime_error("The replicated structure would exceed 50 million atoms");
    Structure& s = data.structure;
    const std::size_t count = s.atoms.size();
    const bool colours = s.grainColors.size() == count, regions = s.grainRegionIds.size() == count, property = s.atomProperty.size() == count;
    Structure out = s;
    out.atoms.clear(); out.grainColors.clear(); out.grainRegionIds.clear(); out.atomProperty.clear();
    std::vector<char> selected;
    for (long long i = 0; i < n[0]; ++i)
        for (long long j = 0; j < n[1]; ++j)
            for (long long k = 0; k < n[2]; ++k) {
                const Vec3 shift = science::add(science::add(science::scale(s.cellVectors[0], static_cast<double>(i)), science::scale(s.cellVectors[1], static_cast<double>(j))),
                                                science::scale(s.cellVectors[2], static_cast<double>(k)));
                for (std::size_t a = 0; a < count; ++a) {
                    AtomSite atom = s.atoms[a];
                    place(atom, science::add(positionOf(atom), shift));
                    out.atoms.push_back(atom);
                    if (colours) out.grainColors.push_back(s.grainColors[a]);
                    if (regions) out.grainRegionIds.push_back(s.grainRegionIds[a]);
                    if (property) out.atomProperty.push_back(s.atomProperty[a]);
                    selected.push_back(data.selected[a]);
                }
            }
    for (int r = 0; r < 3; ++r) out.cellVectors[r] = science::scale(s.cellVectors[r], static_cast<double>(n[r]));
    out.dislocationLoopPoints.clear();
    out.dislocationDetectionDone = false;
    data.structure = std::move(out);
    data.selected = std::move(selected);
}

// Integer supercell (Edit > Transform Structure): the new cell vectors are the
// rows of M as combinations of a, b and c; the atoms fill the new cell.
void supercell(const Json& p, PipelineData& data)
{
    requireCell(data, "supercell");
    std::string matrix = text(p, "supercell", "matrix");
    std::replace(matrix.begin(), matrix.end(), ',', ' ');
    std::istringstream in(matrix);
    Mat3 m{};
    long long integers[3][3];
    for (int r = 0; r < 3; ++r)
        for (int k = 0; k < 3; ++k) {
            double value = 0;
            if (!(in >> value) || std::abs(value - std::round(value)) > 1e-9)
                throw std::runtime_error("The supercell matrix needs nine integers, row by row");
            integers[r][k] = std::llround(value);
            m[static_cast<std::size_t>(r)][static_cast<std::size_t>(k)] = static_cast<double>(integers[r][k]);
        }
    const double det = science::determinant(m);
    if (std::abs(det) < 0.5) throw std::runtime_error("The supercell matrix is singular");
    Structure& s = data.structure;
    const std::size_t count = s.atoms.size();
    if (std::abs(det) * static_cast<double>(count) > 50e6) throw std::runtime_error("The supercell would exceed 50 million atoms");
    const Mat3 cell = cellOf(s), inverseCell = science::inverse(cell), inverseM = science::inverse(m);
    const Vec3 origin = {s.cellOffset[0], s.cellOffset[1], s.cellOffset[2]};
    // The old-cell translations that reach the new cell: its corners, f = u M with u in {0,1}^3.
    long long low[3], high[3];
    for (int k = 0; k < 3; ++k) {
        low[k] = 0; high[k] = 0;
        for (int r = 0; r < 3; ++r) (integers[r][k] < 0 ? low[k] : high[k]) += integers[r][k];
    }
    const bool colours = s.grainColors.size() == count, regions = s.grainRegionIds.size() == count, property = s.atomProperty.size() == count;
    Structure out = s;
    out.atoms.clear(); out.grainColors.clear(); out.grainRegionIds.clear(); out.atomProperty.clear();
    std::vector<char> selected;
    constexpr double eps = 1e-7;
    for (std::size_t a = 0; a < count; ++a) {
        Vec3 f = science::rowTimes(science::sub(positionOf(s.atoms[a]), origin), inverseCell);
        for (double& x : f) x -= std::floor(x);
        for (long long i = low[0]; i <= high[0]; ++i)
            for (long long j = low[1]; j <= high[1]; ++j)
                for (long long k = low[2]; k <= high[2]; ++k) {
                    const Vec3 shifted = {f[0] + static_cast<double>(i), f[1] + static_cast<double>(j), f[2] + static_cast<double>(k)};
                    const Vec3 g = science::rowTimes(shifted, inverseM);
                    if (g[0] < -eps || g[0] >= 1 - eps || g[1] < -eps || g[1] >= 1 - eps || g[2] < -eps || g[2] >= 1 - eps) continue;
                    AtomSite atom = s.atoms[a];
                    place(atom, science::add(origin, science::rowTimes(shifted, cell)));
                    out.atoms.push_back(atom);
                    if (colours) out.grainColors.push_back(s.grainColors[a]);
                    if (regions) out.grainRegionIds.push_back(s.grainRegionIds[a]);
                    if (property) out.atomProperty.push_back(s.atomProperty[a]);
                    selected.push_back(data.selected[a]);
                }
    }
    const Mat3 newCell = science::multiply(m, cell);
    for (int r = 0; r < 3; ++r) out.cellVectors[r] = newCell[static_cast<std::size_t>(r)];
    out.dislocationLoopPoints.clear();
    out.dislocationDetectionDone = false;
    data.notes.push_back("Supercell: " + std::to_string(out.atoms.size()) + " atoms (" + std::to_string(std::llround(std::abs(det))) + " cells)");
    data.structure = std::move(out);
    data.selected = std::move(selected);
}

// ---------------------------------------------------------------- scientific tools
// Analyze-menu tools that produce a structure, as steps: the request is the
// step's parameters (JSON-valued ones as text) plus the step's structure.
struct ScienceStep
{
    const char* id;       // modifier id
    const char* tool;     // catalog tool
    const char* title;
    const char* frame;    // which output frame: "last", or a symmetry cell chosen by the "cell" parameter
};

const std::vector<ScienceStep>& scienceSteps()
{
    static const std::vector<ScienceStep> steps = {
        {"relax", "relax", "Relax structure (interatomic potential)", "last"},
        {"nvt-dynamics", "nvt", "NVT molecular dynamics (final frame)", "last"},
        {"npt-dynamics", "npt", "NPT molecular dynamics (final frame)", "last"},
        {"standardize-cell", "symmetry", "Standardize cell (symmetry)", "cell"},
    };
    return steps;
}

bool scalarKind(const std::string& kind) { return kind == "float" || kind == "int" || kind == "bool"; }

void runScienceStep(const ScienceStep& step, const Json& p, PipelineData& data)
{
    const ::ScienceToolDef* tool = ::findScienceTool(step.tool);
    const ModifierType& type = self(step.id);
    Json request = Json::object();
    for (const auto& definition : tool->parameters) {
        const std::string kind = definition.kind;
        if (kind == "structure") { request[definition.name] = science::structureJson(data.structure); continue; }
        const Json value = parameter(p, type, definition.name);
        if (scalarKind(kind)) { request[definition.name] = value; continue; }
        const std::string text = value.isString() ? value.string() : value.dump();
        if (text.find_first_not_of(" \t") == std::string::npos) continue;  // optional and empty
        try { request[definition.name] = Json::parse(text); }
        catch (const std::exception& error) { throw std::runtime_error(std::string(definition.label) + ": " + error.what()); }
    }
    const auto output = science::runTool(step.tool, request);
    if (output.frames.empty()) throw std::runtime_error(std::string(tool->title) + " produced no structure");
    std::size_t index = output.frames.size() - 1;
    if (std::string(step.frame) == "cell") {
        const std::string wanted = parameter(p, type, "cell").string();
        const Json* order = output.result.find("frame_order");
        for (std::size_t i = 0; order && i < order->items().size(); ++i)
            if (order->items()[i].string() == wanted) index = i;
    }
    Structure result = output.frames[index];
    // Same atoms (relaxation, dynamics): keep the selection and per-atom metadata.
    if (result.atoms.size() == data.structure.atoms.size()) {
        Structure& s = data.structure;
        for (std::size_t i = 0; i < s.atoms.size(); ++i) place(s.atoms[i], positionOf(result.atoms[i]));
        s.cellVectors = result.cellVectors;
        s.hasUnitCell = result.hasUnitCell;
    } else {
        data.structure = std::move(result);
        data.clearSelection();
    }
    data.notes.push_back(std::string(tool->title) + ": " + std::to_string(data.structure.atoms.size()) + " atoms");
}

// The steps' types, with the tool's own parameters (structure inputs come from the pipeline).
std::vector<ModifierType> scienceStepTypes()
{
    static std::deque<std::string> texts;
    std::vector<ModifierType> types;
    for (const ScienceStep& step : scienceSteps()) {
        const ::ScienceToolDef* tool = ::findScienceTool(step.tool);
        if (!tool) continue;
        std::vector<ModifierParameter> parameters;
        for (const auto& definition : tool->parameters) {
            const std::string kind = definition.kind;
            if (kind == "structure") continue;
            // JSON-valued inputs (potentials, lists) are edited as their JSON text.
            parameters.push_back({definition.name, definition.label, scalarKind(kind) ? definition.kind : "string", definition.value, ""});
        }
        if (std::string(step.frame) == "cell")
            parameters.push_back({"cell", "Cell to keep", "choice", "primitive", "primitive|conventional|symmetrized"});
        const std::string help = tool->help;
        texts.push_back("Runs " + std::string(tool->title) + " (Analyze > " + tool->category + ") on the structure and keeps its " +
                        (std::string(step.frame) == "cell" ? "chosen standardized cell" : "final structure") + ". " + help.substr(0, help.find('\n')));
        const char* helpText = texts.back().c_str();
        const ScienceStep* captured = &step;
        types.push_back({step.id, step.title, "Simulation", helpText, parameters,
                         [captured](const Json& p, PipelineData& data) { runScienceStep(*captured, p, data); }});
    }
    return types;
}

void transform(const Json& p, PipelineData& data)
{
    Mat3 m{};
    m[0] = vector(p, "transform", "row1");
    m[1] = vector(p, "transform", "row2");
    m[2] = vector(p, "transform", "row3");
    if (std::abs(science::determinant(m)) < 1e-12) throw std::runtime_error("The transformation matrix is singular");
    const bool selectedOnly = text(p, "transform", "target") == "selected";
    // Positions and cell vectors are row vectors: x' = x M.
    for (std::size_t i : targets(data, selectedOnly, "transform"))
        place(data.structure.atoms[i], science::rowTimes(positionOf(data.structure.atoms[i]), m));
    if (flag(p, "transform", "transform_cell") && !selectedOnly && data.structure.hasUnitCell)
        for (auto& row : data.structure.cellVectors) row = science::rowTimes(row, m);
}

void strain(const Json& p, PipelineData& data)
{
    const Vec3 normal = vector(p, "strain", "normal"), shear = vector(p, "strain", "shear");
    // Engineering-Voigt shears (yz, xz, xy) are twice the tensor components.
    Mat3 f = science::identity();
    f[0][0] += normal[0]; f[1][1] += normal[1]; f[2][2] += normal[2];
    f[1][2] = f[2][1] = 0.5 * shear[0];
    f[0][2] = f[2][0] = 0.5 * shear[1];
    f[0][1] = f[1][0] = 0.5 * shear[2];
    for (auto& a : data.structure.atoms) place(a, science::rowTimes(positionOf(a), f));
    if (data.structure.hasUnitCell)
        for (auto& row : data.structure.cellVectors) row = science::rowTimes(row, f);
}

void wrap(const Json&, PipelineData& data)
{
    requireCell(data, "wrap");
    const Mat3 cell = cellOf(data.structure), inv = science::inverse(cell);
    for (auto& a : data.structure.atoms) {
        Vec3 f = science::rowTimes(science::sub(positionOf(a), data.structure.cellOffset), inv);
        for (double& v : f) v -= std::floor(v);
        place(a, science::add(science::rowTimes(f, cell), data.structure.cellOffset));
    }
}

void center(const Json&, PipelineData& data)
{
    requireCell(data, "center");
    if (data.structure.atoms.empty()) return;
    const Mat3 cell = cellOf(data.structure), inv = science::inverse(cell);
    Vec3 low{HUGE_VAL, HUGE_VAL, HUGE_VAL}, high{-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
    for (const auto& a : data.structure.atoms) {
        const Vec3 f = science::rowTimes(positionOf(a), inv);
        for (int k = 0; k < 3; ++k) { low[k] = std::min(low[k], f[k]); high[k] = std::max(high[k], f[k]); }
    }
    Vec3 shift{};
    for (int k = 0; k < 3; ++k) shift[k] = 0.5 - 0.5 * (low[k] + high[k]);
    const Vec3 cartesian = science::add(science::rowTimes(shift, cell), data.structure.cellOffset);
    for (auto& a : data.structure.atoms) place(a, science::add(positionOf(a), cartesian));
}

void addVacuum(const Json& p, PipelineData& data)
{
    requireCell(data, "add-vacuum");
    const std::string axis = text(p, "add-vacuum", "axis");
    const int k = axis == "a" ? 0 : axis == "b" ? 1 : 2;
    const double thickness = number(p, "add-vacuum", "thickness");
    if (thickness < 0) throw std::runtime_error("The vacuum thickness must not be negative");
    auto& row = data.structure.cellVectors[k];
    const double length = science::norm(row);
    const Vec3 direction = science::scale(row, 1 / length);
    row = science::scale(row, (length + thickness) / length);
    if (flag(p, "add-vacuum", "center"))
        for (auto& a : data.structure.atoms) place(a, science::add(positionOf(a), science::scale(direction, 0.5 * thickness)));
}

// Positions, cell and periodicity of the structure as a tool request.
Json atomsRequest(const Structure& s)
{
    Json positions = Json::array();
    for (const auto& a : s.atoms) positions.push(Json::array({a.x, a.y, a.z}));
    Json request = Json::object();
    request["positions"] = positions;
    if (s.hasUnitCell) {
        Json cell = Json::array();
        for (const auto& row : s.cellVectors) cell.push(Json::array({row[0], row[1], row[2]}));
        request["cell"] = cell;
        request["pbc"] = Json::array({true, true, true});
    }
    return request;
}

void setProperty(PipelineData& data, std::vector<double> values, const std::string& name)
{
    data.structure.atomProperty = std::move(values);
    data.structure.atomPropertyName = name;
    std::size_t valid = 0;
    for (double v : data.structure.atomProperty) valid += std::isfinite(v);
    data.notes.push_back(name + " computed for " + std::to_string(valid) + " of " + std::to_string(data.structure.atomProperty.size()) + " atoms");
}

std::vector<double> numbersOf(const Json& list)
{
    std::vector<double> out;
    for (const auto& item : list.items()) out.push_back(item.isNumber() ? item.number() : std::numeric_limits<double>::quiet_NaN());
    return out;
}

// Elastic strain relative to an ideal lattice (polyhedral template matching).
void elasticStrainStep(const Json& p, PipelineData& data)
{
    Json request = atomsRequest(data.structure);
    request["lattice"] = text(p, "elastic-strain", "lattice");
    request["lattice_constant_A"] = number(p, "elastic-strain", "a");
    request["c_over_a"] = number(p, "elastic-strain", "c_over_a");
    const Json result = science::runTool("elastic-strain", request).result;
    const std::string component = text(p, "elastic-strain", "component");
    if (component == "volumetric") { setProperty(data, numbersOf(result.at("volumetric_strain")), "Volumetric elastic strain"); return; }
    if (component == "shear") { setProperty(data, numbersOf(result.at("shear_strain")), "Von Mises elastic shear strain"); return; }
    static const char* names[] = {"xx", "yy", "zz", "yz", "xz", "xy"};
    const std::size_t c = static_cast<std::size_t>(std::find(names, names + 6, component) - names);
    std::vector<double> values;
    for (const auto& row : result.at("strain_voigt").items()) values.push_back(row.items()[c].isNumber() ? row.items()[c].number() : std::numeric_limits<double>::quiet_NaN());
    setProperty(data, values, "Elastic strain " + component);
}

// Displacements from a reference structure file (same atoms, same order).
void displacementStep(const Json& p, PipelineData& data)
{
    const std::string file = text(p, "displacements", "reference");
    if (file.empty()) throw std::runtime_error("Choose the reference structure file");
    const Structure reference = science::readFrames(std::filesystem::u8path(file)).front().structure;
    const Structure& s = data.structure;
    if (reference.atoms.size() != s.atoms.size())
        throw std::runtime_error("The reference has " + std::to_string(reference.atoms.size()) + " atoms, the structure " + std::to_string(s.atoms.size()));
    Json request = Json::object();
    Json current = atomsRequest(s).at("positions");
    request["reference"] = atomsRequest(reference).at("positions");
    request["current"] = current;
    const std::string mapping = text(p, "displacements", "affine_mapping");
    if (reference.hasUnitCell) {
        request["reference_cell"] = atomsRequest(reference).at("cell");
        request["pbc"] = Json::array({true, true, true});
        if (s.hasUnitCell) request["current_cell"] = atomsRequest(s).at("cell");
    }
    if (mapping != "off") request["affine_mapping"] = mapping;
    const Json result = science::runTool("displacements", request).result;
    const std::string component = text(p, "displacements", "component");
    if (component == "magnitude") { setProperty(data, numbersOf(result.at("magnitude_A")), "Displacement magnitude (A)"); return; }
    const std::size_t c = component == "x" ? 0 : component == "y" ? 1 : 2;
    std::vector<double> values;
    for (const auto& row : result.at("displacement_A").items()) values.push_back(row.items()[c].number());
    setProperty(data, values, "Displacement " + component + " (A)");
}

// Statistics of the per-atom property (OVITO's Reduce property), as notes.
void reduceProperty(const Json& p, PipelineData& data)
{
    const Structure& s = data.structure;
    if (s.atomProperty.size() != s.atoms.size() || s.atoms.empty()) throw std::runtime_error("reduce-property needs a per-atom property; add compute-property above");
    const bool selectedOnly = text(p, "reduce-property", "target") == "selected";
    std::vector<double> values;
    for (std::size_t i = 0; i < s.atoms.size(); ++i)
        if ((!selectedOnly || data.selected[i]) && std::isfinite(s.atomProperty[i])) values.push_back(s.atomProperty[i]);
    if (values.empty()) throw std::runtime_error("No atom has a value");
    Json request = Json::object();
    request["values"] = science::toJson(values);
    request["bins"] = number(p, "reduce-property", "bins");
    const Json r = science::runTool("property-histogram", request).result;
    char line[256];
    std::snprintf(line, sizeof line, "%s: n=%.0f sum=%.6g mean=%.6g std=%.6g min=%.6g median=%.6g max=%.6g",
                  s.atomPropertyName.empty() ? "property" : s.atomPropertyName.c_str(), r.at("count").number(), r.at("sum").number(),
                  r.at("mean").number(), r.at("std").number(), r.at("min").number(), r.at("median").number(), r.at("max").number());
    data.notes.push_back(line);
    std::string histogram = "histogram:";
    const auto centers = numbersOf(r.at("centers")), counts = numbersOf(r.at("counts"));
    for (std::size_t b = 0; b < centers.size(); ++b) {
        std::snprintf(line, sizeof line, " %.4g:%.0f", centers[b], counts[b]);
        histogram += line;
    }
    data.notes.push_back(histogram);
}

// Spatial binning of the per-atom property (or atom counts); each atom takes
// the value of its bin, so the profile or map can colour the structure.
void binProperty(const Json& p, PipelineData& data)
{
    const Structure& s = data.structure;
    requireCell(data, "bin-property");
    Json request = atomsRequest(s);
    const bool hasProperty = s.atomProperty.size() == s.atoms.size() && !s.atoms.empty();
    if (hasProperty) request["values"] = science::toJson(s.atomProperty);
    const Vec3 bins = vector(p, "bin-property", "bins");
    request["bins"] = Json::array({bins[0], bins[1], bins[2]});
    request["operation"] = text(p, "bin-property", "operation");
    const Json result = science::runTool("spatial-binning", request).result;
    const auto grid = numbersOf(result.at("values"));
    const Mat3 cell = cellOf(s);
    std::vector<double> perAtom(s.atoms.size());
    const long long n[3] = {std::llround(bins[0]), std::llround(bins[1]), std::llround(bins[2])};
    for (std::size_t i = 0; i < s.atoms.size(); ++i) {
        Vec3 f = science::fractional(positionOf(s.atoms[i]), cell);
        long long index[3];
        for (int d = 0; d < 3; ++d) {
            f[static_cast<std::size_t>(d)] -= std::floor(f[static_cast<std::size_t>(d)]);
            index[d] = std::min(n[d] - 1, static_cast<long long>(f[static_cast<std::size_t>(d)] * static_cast<double>(n[d])));
        }
        perAtom[i] = grid[static_cast<std::size_t>((index[0] * n[1] + index[1]) * n[2] + index[2])];
    }
    const std::string operation = text(p, "bin-property", "operation");
    setProperty(data, perAtom, "Binned " + operation + (hasProperty ? " of " + s.atomPropertyName : std::string(" (atoms)")));
}

// An orthogonal cell just enclosing the atoms (plus padding), with the atoms
// shifted to start at the padding.
void shrinkWrap(const Json& p, PipelineData& data)
{
    Structure& s = data.structure;
    if (s.atoms.empty()) throw std::runtime_error("shrink-wrap needs atoms");
    const double padding = number(p, "shrink-wrap", "padding");
    if (padding < 0) throw std::runtime_error("padding must be nonnegative");
    Vec3 low = positionOf(s.atoms[0]), high = low;
    for (const auto& a : s.atoms)
        for (int d = 0; d < 3; ++d) {
            const double x = positionOf(a)[static_cast<std::size_t>(d)];
            low[static_cast<std::size_t>(d)] = std::min(low[static_cast<std::size_t>(d)], x);
            high[static_cast<std::size_t>(d)] = std::max(high[static_cast<std::size_t>(d)], x);
        }
    for (auto& a : s.atoms) place(a, {a.x - low[0] + padding, a.y - low[1] + padding, a.z - low[2] + padding});
    s.hasUnitCell = true;
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            s.cellVectors[r][c] = r == c ? std::max(high[static_cast<std::size_t>(r)] - low[static_cast<std::size_t>(r)] + 2 * padding, 1e-3) : 0.0;
    s.cellOffset = {0.0, 0.0, 0.0};
    data.notes.push_back("Cell " + std::to_string(s.cellVectors[0][0]) + " x " + std::to_string(s.cellVectors[1][1]) + " x " +
                         std::to_string(s.cellVectors[2][2]) + " Angstrom");
}

// Colours the selected (or all) atoms; the colour stays with the atoms.
void assignColor(const Json& p, PipelineData& data)
{
    const Vec3 colour = vector(p, "assign-color", "color");
    for (double c : colour)
        if (!(c >= 0 && c <= 1)) throw std::runtime_error("Colour components are between 0 and 1 (red green blue)");
    const auto indices = targets(data, text(p, "assign-color", "target") == "selected", "assign-color");
    Structure& s = data.structure;
    for (std::size_t i : indices) {
        s.atoms[i].r = static_cast<float>(colour[0]);
        s.atoms[i].g = static_cast<float>(colour[1]);
        s.atoms[i].b = static_cast<float>(colour[2]);
    }
    if (s.grainColors.size() == s.atoms.size())
        for (std::size_t i : indices) s.grainColors[i] = {static_cast<float>(colour[0]), static_cast<float>(colour[1]), static_cast<float>(colour[2])};
    data.notes.push_back(std::to_string(indices.size()) + " atoms coloured");
}

// ---------------------------------------------------------------- analysis
void computeExpression(const Json& p, PipelineData& data)
{
    const std::string formula = text(p, "compute-expression", "expression");
    const auto values = evaluateForAtoms(Expression(formula), data);
    const bool selectedOnly = text(p, "compute-expression", "target") == "selected";
    Structure& s = data.structure;
    std::vector<double> property = s.atomProperty.size() == s.atoms.size() ? s.atomProperty
                                   : std::vector<double>(s.atoms.size(), std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < values.size(); ++i)
        if (!selectedOnly || data.selected[i]) property[i] = values[i];
    s.atomProperty = std::move(property);
    const std::string name = text(p, "compute-expression", "name");
    s.atomPropertyName = name.empty() ? formula : name;
    data.notes.push_back(s.atomPropertyName + " computed for " + std::to_string(selectedOnly ? data.selectedCount() : values.size()) + " atoms");
}

// Grains from local lattice orientations; atoms are coloured by grain and the
// grain id becomes the per-atom property.
void segmentGrains(const Json& p, PipelineData& data)
{
    Structure& s = data.structure;
    Json request = atomsRequest(s);
    request["merge_threshold_deg"] = number(p, "segment-grains", "threshold");
    request["min_grain_size"] = number(p, "segment-grains", "min_size");
    request["adopt_orphans"] = flag(p, "segment-grains", "adopt");
    const Json result = science::runTool("grain-segmentation", request).result;
    const Json& ids = result.at("grain_id");
    s.atomProperty.clear();
    s.grainColors.clear();
    for (const auto& id : ids.items()) {
        const double g = id.number();
        s.atomProperty.push_back(g);
        // Golden-ratio hues give neighbouring grain ids distinct colours; grey for none.
        if (g <= 0) { s.grainColors.push_back({0.6f, 0.6f, 0.6f}); continue; }
        const double hue = std::fmod(g * 0.6180339887, 1.0) * 6;
        const int sector = static_cast<int>(hue);
        const float f = static_cast<float>(hue - sector), v = 0.9f, lo = 0.25f;
        const float rise = lo + (v - lo) * f, fall = v - (v - lo) * f;
        const std::array<std::array<float, 3>, 6> colours = {{{v, rise, lo}, {fall, v, lo}, {lo, v, rise}, {lo, fall, v}, {rise, lo, v}, {v, lo, fall}}};
        s.grainColors.push_back(colours[static_cast<std::size_t>(sector % 6)]);
    }
    s.atomPropertyName = "Grain id (0 none)";
    data.notes.push_back(std::to_string(static_cast<long long>(result.at("grain_count").number())) + " grains");
}

void computeProperty(const Json& p, PipelineData& data)
{
    const std::string analysis = text(p, "compute-property", "analysis");
    const double cutoff = number(p, "compute-property", "cutoff");
    science::AtomProperty property;
    if (analysis == "structure-type") property = science::analysePerAtom(data.structure, "structure-type", cutoff);
    else if (analysis == "diamond") property = science::analysePerAtom(data.structure, "diamond-structure", cutoff);
    else if (analysis == "centrosymmetry") property = science::analysePerAtom(data.structure, "centrosymmetry", cutoff);
    else if (analysis == "coordination") property = science::analysePerAtom(data.structure, "bond-order", cutoff, "Coordination");
    else if (analysis == "ptm") property = science::analysePerAtom(data.structure, "ptm", cutoff);
    else if (analysis == "cna") property = science::analysePerAtom(data.structure, "cna", cutoff);
    else if (analysis == "voronoi-volume") property = science::analysePerAtom(data.structure, "voronoi", cutoff);
    else if (analysis == "voronoi-coordination") property = science::analysePerAtom(data.structure, "voronoi", cutoff, "Voronoi coordination");
    else if (analysis == "orientation") property = science::analysePerAtom(data.structure, "ptm", cutoff, "Orientation angle (degrees)");
    else if (analysis == "planar-faults") property = science::analysePerAtom(data.structure, "fcc-planar-faults", cutoff);
    else if (analysis == "chill-plus") {
        Structure oxygens = data.structure;
        std::vector<std::size_t> index;
        oxygens.atoms.clear();
        for (std::size_t i = 0; i < data.structure.atoms.size(); ++i)
            if (data.structure.atoms[i].symbol == "O") { oxygens.atoms.push_back(data.structure.atoms[i]); index.push_back(i); }
        if (oxygens.atoms.empty() || oxygens.atoms.size() == data.structure.atoms.size()) property = science::analysePerAtom(data.structure, "chill-plus", cutoff);
        else {
            const auto onOxygens = science::analysePerAtom(oxygens, "chill-plus", cutoff);
            property.name = onOxygens.name;
            property.values.assign(data.structure.atoms.size(), std::numeric_limits<double>::quiet_NaN());
            for (std::size_t k = 0; k < index.size(); ++k) property.values[index[k]] = onOxygens.values[k];
        }
    }
    else if (analysis == "local-entropy") property = science::analysePerAtom(data.structure, "local-entropy", cutoff);
    else if (analysis == "rings") property = science::analysePerAtom(data.structure, "find-rings", cutoff);
    else property = science::analysePerAtom(data.structure, "bond-order", cutoff);
    data.structure.atomProperty = property.values;
    data.structure.atomPropertyName = property.name;
    data.notes.push_back(property.name + " computed for " + std::to_string(property.values.size()) + " atoms");
}
}

namespace
{
std::vector<ModifierType>& registry()
{
    static std::vector<ModifierType> types = [] {
        std::vector<ModifierType> builtIn = {
        // Selection
        {"select-element", "Select by element", "Selection", "Selects atoms of the listed elements.",
         {{"elements", "Elements (e.g. Cu O)", "elements", "Cu", ""}, MODE}, selectElement},
        {"select-expression", "Select by expression", "Selection",
         "Selects atoms for which a condition holds, e.g. fz > 0.5 && element == O.\nVariables: x y z (Angstrom), fx fy fz (fractional), index, count, element,\nZ, property (per-atom property), selected. Operators + - * / ^, < <= > >= == !=,\n&& || ! (and or not), functions abs sqrt exp log sin cos tan floor ceil min max.",
         {{"expression", "Condition", "expression", "x > 0", ""}, MODE}, selectExpression},
        {"select-slab", "Select a slab", "Selection", "Selects atoms between two coordinates along x, y or z (Angstrom) or along a, b or c (fractional, wrapped into [0, 1)).",
         {{"axis", "Axis", "choice", "z", "x|y|z|a|b|c"}, {"minimum", "From", "float", "0", ""}, {"maximum", "To", "float", "5", ""}, MODE}, selectSlab},
        {"select-sphere", "Select a sphere", "Selection", "Selects atoms within a radius of a point (minimum-image distance in periodic cells).",
         {{"center", "Centre (Angstrom)", "vector", "0 0 0", ""}, {"radius", "Radius (Angstrom)", "float", "3", ""}, MODE}, selectSphere},
        {"select-property", "Select by property range", "Selection", "Selects atoms whose per-atom property lies in a range (needs compute-property or a coloured result).",
         {{"minimum", "Minimum", "float", "0", ""}, {"maximum", "Maximum", "float", "1", ""}, MODE}, selectProperty},
        {"select-random", "Select a random fraction", "Selection", "Selects round(fraction x N) atoms at random, reproducibly for a seed.",
         {{"fraction", "Fraction", "float", "0.1", ""}, {"seed", "Seed", "int", "1", ""}, MODE}, selectRandom},
        {"expand-selection", "Expand selection to neighbours", "Selection", "Adds the neighbours (within the cutoff) of selected atoms, repeated for the number of steps.",
         {{"cutoff", "Neighbour cutoff (Angstrom)", "float", "3", ""}, {"steps", "Steps", "int", "1", ""}}, expandSelection},
        {"select-overlapping", "Select overlapping atoms", "Selection",
         "Selects atoms closer than a distance to another atom (minimum image in periodic cells): the later atom of each pair, or both.",
         {{"distance", "Distance (Angstrom)", "float", "1", ""}, {"atoms", "Atoms to select", "choice", "one", "one|both"}, MODE}, selectOverlapping},
        {"invert-selection", "Invert selection", "Selection", "Selects the unselected atoms and deselects the rest.", {}, invertSelection},
        {"clear-selection", "Clear selection", "Selection", "Deselects all atoms.", {}, clearSelection},
        // Modification
        {"delete-selected", "Delete selected atoms", "Modification", "Removes the selected atoms.", {}, deleteSelected},
        {"assign-element", "Assign element", "Modification", "Changes the element of the selected (or all) atoms, e.g. to create substitutions.",
         {{"element", "Element", "element", "Cu", ""}, {"target", "Atoms", "choice", "selected", "selected|all"}}, assignElement},
        {"fix-atoms", "Fix atoms", "Modification",
         "Holds coordinates of the selected (or all) atoms fixed in a relaxation; POSCAR files keep them as selective dynamics (F).",
         {{"axes", "Fixed coordinates", "choice", "xyz", "xyz|xy|xz|yz|x|y|z"}, {"target", "Atoms", "choice", "selected", "selected|all"}}, fixAtoms},
        {"free-atoms", "Free atoms", "Modification", "Lets the selected (or all) atoms move along every coordinate again.",
         {{"target", "Atoms", "choice", "selected", "selected|all"}}, freeAtoms},
        {"displace", "Displace", "Modification", "Moves atoms by a fixed vector.",
         {{"vector", "Displacement (Angstrom)", "vector", "0 0 0", ""}, {"target", "Atoms", "choice", "all", "all|selected"}}, displace},
        {"random-displacement", "Random displacement", "Modification", "Adds Gaussian displacements (standard deviation per component), reproducibly for a seed.",
         {{"amplitude", "Standard deviation (Angstrom)", "float", "0.05", ""}, {"seed", "Seed", "int", "1", ""}, {"target", "Atoms", "choice", "all", "all|selected"}}, randomDisplacement},
        {"add-atom", "Add atom", "Modification", "Adds one atom at a Cartesian (Angstrom) or fractional position.",
         {{"element", "Element", "element", "Cu", ""}, {"position", "Position", "vector", "0 0 0", ""}, {"fractional", "Fractional position", "bool", "false", ""},
          {"select", "Select the new atom", "bool", "true", ""}}, addAtom},
        {"insert-interstitials", "Insert interstitial atoms", "Modification", "Fills interstitial voids (found as in Edit > Add Interstitial Atoms) with atoms: the largest voids first or at random, all of them (count 0) or a number; the new atoms are selected.",
         {{"element", "Element", "element", "H", ""}, {"kind", "Void kind", "choice", "any", "any|tetrahedral|octahedral|irregular"},
          {"count", "Atoms to insert (0 = every site)", "int", "1", ""}, {"order", "Choose sites", "choice", "largest", "largest|random"},
          {"seed", "Seed (random order)", "int", "1", ""}, {"resolution", "Grid resolution", "int", "14", ""},
          {"clearance", "Minimum clearance (Angstrom)", "float", "0.6", ""}, {"separation", "Minimum separation (Angstrom)", "float", "0.5", ""}}, insertInterstitials},
        {"merge", "Merge a structure file", "Modification", "Appends the atoms of another structure file (POSCAR, extended XYZ, JSON or LAMMPS dump), shifted by an offset.",
         {{"file", "Structure file", "string", "", ""}, {"offset", "Offset (Angstrom)", "vector", "0 0 0", ""}, {"select", "Select the merged atoms", "bool", "true", ""},
          {"use_cell", "Take its cell when this structure has none", "bool", "true", ""}}, merge},
        {"slice", "Slice", "Modification", "Cuts with a plane n.r = distance: removes the half-space above it, or with a width the slab |n.r - distance| <= width/2. Invert keeps that region instead; action select only selects it.",
         {{"normal", "Plane normal", "vector", "0 0 1", ""}, {"distance", "Distance from the origin (Angstrom)", "float", "0", ""},
          {"width", "Slab width (0 = half-space)", "float", "0", ""}, {"invert", "Invert", "bool", "false", ""}, {"action", "Action", "choice", "delete", "delete|select"}}, slice},
        // Cell
        {"replicate", "Replicate", "Cell", "Repeats the cell along a, b and c.",
         {{"counts", "Copies along a, b, c", "ivec3", "2 2 2", ""}}, replicate},
        {"supercell", "Transform structure (supercell)", "Cell", "Builds the supercell whose cell vectors are integer combinations of a, b and c (rows of the matrix, as in Edit > Transform Structure); the atoms fill the new cell.",
         {{"matrix", "Matrix rows (nine integers)", "string", "2 0 0 0 2 0 0 0 2", ""}}, supercell},
        {"transform", "Affine transformation", "Cell", "Applies x' = x M to atoms (row vectors) and, for all atoms, optionally to the cell.",
         {{"row1", "Matrix row 1", "vector", "1 0 0", ""}, {"row2", "Matrix row 2", "vector", "0 1 0", ""}, {"row3", "Matrix row 3", "vector", "0 0 1", ""},
          {"transform_cell", "Transform the cell too", "bool", "true", ""}, {"target", "Atoms", "choice", "all", "all|selected"}}, transform},
        {"strain", "Strain", "Cell", "Homogeneous strain of atoms and cell: normal components xx yy zz and engineering shears yz xz xy.",
         {{"normal", "Normal strain xx, yy, zz", "vector", "0 0 0", ""}, {"shear", "Engineering shear yz, xz, xy", "vector", "0 0 0", ""}}, strain},
        {"set-cell", "Set cell vectors", "Cell", "Replaces the lattice vectors; with scale atoms, atoms keep their fractional coordinates.",
         {{"a", "Vector a (Angstrom)", "vector", "1 0 0", ""}, {"b", "Vector b", "vector", "0 1 0", ""}, {"c", "Vector c", "vector", "0 0 1", ""},
          {"scale_atoms", "Scale atoms with the cell", "bool", "true", ""}}, setCell},
        {"wrap", "Wrap into cell", "Cell", "Maps atoms outside the cell back into it.", {}, wrap},
        {"assign-color", "Assign colour", "Modification", "Colours the selected (or all) atoms with a red green blue colour (components 0 to 1).",
         {{"color", "Colour (red green blue, 0-1)", "vector", "1 0 0", ""}, {"target", "Atoms", "choice", "selected", "selected|all"}}, assignColor},
        {"shrink-wrap", "Shrink-wrap cell", "Cell", "Replaces the cell by the orthogonal box that just encloses the atoms, plus padding on every side.",
         {{"padding", "Padding (Angstrom)", "float", "0", ""}}, shrinkWrap},
        {"center", "Center in cell", "Cell", "Shifts the atoms so their extent (in fractional coordinates) is centred in the cell.", {}, center},
        {"add-vacuum", "Add vacuum", "Cell", "Lengthens one cell vector by a thickness, e.g. to turn a bulk cell into a slab; optionally centres the atoms.",
         {{"axis", "Cell vector", "choice", "c", "a|b|c"}, {"thickness", "Vacuum (Angstrom)", "float", "10", ""}, {"center", "Centre the atoms", "bool", "true", ""}}, addVacuum},
        // Analysis
        {"compute-property", "Compute per-atom property", "Analysis", "Computes a per-atom property for colouring and select-property: Ackland-Jones structure type, adaptive CNA, diamond structure, Voronoi volume or coordination, centrosymmetry, coordination,\nSteinhardt q6, polyhedral template matching (structure or orientation angle), fcc planar faults, Chill+ ice types, local entropy or rings\nthrough each atom (bond cutoff). PTM, planar faults, Chill+ and entropy choose their own neighbour range.",
         {{"analysis", "Property", "choice", "structure-type", "structure-type|cna|diamond|centrosymmetry|coordination|q6|ptm|orientation|voronoi-volume|voronoi-coordination|planar-faults|chill-plus|local-entropy|rings"},
          {"cutoff", "Neighbour cutoff (Angstrom)", "float", "3", ""}}, computeProperty},
        {"compute-expression", "Compute property from a formula", "Analysis",
         "Sets the per-atom property to a formula of the atom variables, e.g. sqrt(x^2 + y^2) or property * 2 (see select-expression for the variables);\nwith target selected, only selected atoms change.",
         {{"expression", "Formula", "expression", "z", ""}, {"name", "Property name (empty: the formula)", "string", "", ""},
          {"target", "Atoms", "choice", "all", "all|selected"}}, computeExpression},
        {"elastic-strain", "Elastic strain", "Analysis",
         "Sets the per-atom property to the elastic strain relative to the ideal lattice (polyhedral template matching): volumetric (trace/3), von Mises\nshear or one Green-Lagrange component. Atoms of other structures get no value.",
         {{"lattice", "Lattice", "choice", "fcc", "fcc|hcp|bcc|sc|cubic-diamond|hex-diamond"}, {"a", "Unstrained lattice constant (Angstrom)", "float", "3.615", ""},
          {"c_over_a", "Unstrained c/a (hexagonal)", "float", "1.633", ""},
          {"component", "Property", "choice", "volumetric", "volumetric|shear|xx|yy|zz|yz|xz|xy"}}, elasticStrainStep},
        {"displacements", "Displacement vectors", "Analysis",
         "Sets the per-atom property to the displacement from a reference structure file (same atoms in the same order; minimum image in periodic\ncells; optionally after mapping the cell deformation).",
         {{"reference", "Reference structure file", "string", "", ""}, {"component", "Property", "choice", "magnitude", "magnitude|x|y|z"},
          {"affine_mapping", "Affine mapping", "choice", "off", "off|to_reference|to_current"}}, displacementStep},
        {"bin-property", "Spatial binning", "Analysis",
         "Averages (or sums, counts, ...) the per-atom property over bins along a, b and c (1 = not binned) and gives each atom its bin's value;\nwithout a property the atoms are counted. bins 1 1 20 gives a profile along c.",
         {{"bins", "Bins along a, b, c", "ivec3", "1 1 20", ""}, {"operation", "Operation", "choice", "mean", "mean|sum|density|min|max|count"}}, binProperty},
        {"reduce-property", "Property statistics", "Analysis",
         "Reports the count, sum, mean, standard deviation, minimum, median and maximum and a histogram of the per-atom property (all or selected atoms).",
         {{"target", "Atoms", "choice", "all", "all|selected"}, {"bins", "Histogram bins", "int", "10", ""}}, reduceProperty},
        {"segment-grains", "Segment grains", "Analysis",
         "Finds the grains of a polycrystal from local lattice orientations (polyhedral template matching), colours atoms by grain and sets the\nper-atom property to the grain id (0 for none).",
         {{"threshold", "Misorientation threshold (degrees)", "float", "4", ""}, {"min_size", "Minimum grain size (atoms)", "int", "100", ""},
          {"adopt", "Assign boundary atoms to grains", "bool", "true", ""}}, segmentGrains},
        };
        for (auto& type : scienceStepTypes()) builtIn.push_back(std::move(type));
        builtIn.reserve(256);  // registered types keep their addresses
        return builtIn;
    }();
    return types;
}
}

const std::vector<ModifierType>& modifierTypes() { return registry(); }

void registerModifierType(ModifierType type)
{
    auto& types = registry();
    for (auto& existing : types)
        if (std::string(existing.id) == type.id) { existing = std::move(type); return; }
    if (types.size() == types.capacity()) throw std::runtime_error("Too many modifier types");
    types.push_back(std::move(type));
}
#undef MODE
}
