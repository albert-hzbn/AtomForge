// Analysis, property and colouring steps of the pipeline (in the spirit of
// OVITO's modifiers): per-atom analyses through the scientific tools,
// statistics, binning and selections derived from them.
#include "pipeline/ModifierSupport.h"
#include "pipeline/PipelineSteps.h"
#include "science/AtomProperties.h"
#include "science/ScienceData.h"
#include "science/ScienceTools.h"
#include "util/ElementData.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <stdexcept>

namespace atomforge::pipeline
{
namespace
{
using namespace support;

// Positions, cell and periodic axes of the structure as a tool request.
Json atomsRequest(const Structure& s) { return science::atomsRequest(s, science::periodicAxes(s)); }

// Atoms on the surface of the solid: those a probe sphere touches.
void selectSurface(const Json& p, PipelineData& data)
{
    Json request = atomsRequest(data.structure);
    request["probe_radius_A"] = number(p, "select-surface", "probe_radius");
    request["export_mesh"] = false;
    const Json result = science::runTool("surface-mesh", request).result;
    std::vector<char> mask;
    for (const auto& flag : result.at("surface_atom").items()) mask.push_back(flag.boolean() ? 1 : 0);
    char note[160];
    std::snprintf(note, sizeof note, "Surface area %.6g A^2, solid volume %.6g A^3", result.at("surface_area_A2").number(), result.at("solid_volume_A3").number());
    data.notes.push_back(note);
    select(data, mask, text(p, "select-surface", "mode"));
}

// Atoms inside LAMMPS regions (region commands, ';' between them).
void selectRegion(const Json& p, PipelineData& data)
{
    Json request = atomsRequest(data.structure);
    request["regions"] = text(p, "select-region", "regions");
    request["lattice_scale"] = number(p, "select-region", "lattice_scale");
    request["export_mesh"] = false;
    const Json result = science::runTool("lammps-regions", request).result;
    std::vector<char> mask;
    for (const auto& flag : result.at("inside_any").items()) mask.push_back(flag.boolean() ? 1 : 0);
    select(data, mask, text(p, "select-region", "mode"));
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
    const science::Pbc pbc = science::periodicAxes(s);
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
        request["pbc"] = atomsRequest(reference).at("pbc");
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

// Darkens atom colours by ambient occlusion (covalent radii).
void ambientOcclusionStep(const Json& p, PipelineData& data)
{
    Structure& s = data.structure;
    if (s.atoms.empty()) return;
    static const std::vector<float> covalent = makeLiteratureCovalentRadii();
    Json positions = Json::array();
    std::vector<double> radii;
    for (const auto& a : s.atoms) {
        positions.push(Json::array({a.x, a.y, a.z}));
        const std::size_t z = static_cast<std::size_t>(std::max(0, a.atomicNumber));
        radii.push_back(z < covalent.size() && covalent[z] > 0 ? covalent[z] : 1.0);
    }
    Json request = Json::object();
    request["positions"] = positions;
    request["radii"] = science::toJson(radii);
    request["directions"] = number(p, "ambient-occlusion", "directions");
    const Json result = science::runTool("ambient-occlusion", request).result;
    const double intensity = number(p, "ambient-occlusion", "intensity");
    if (!(intensity >= 0 && intensity <= 1)) throw std::runtime_error("intensity is between 0 and 1");
    const auto& brightness = result.at("brightness").items();
    for (std::size_t i = 0; i < s.atoms.size(); ++i) {
        const float factor = static_cast<float>(1 - intensity + intensity * brightness[i].number());
        s.atoms[i].r *= factor;
        s.atoms[i].g *= factor;
        s.atoms[i].b *= factor;
        if (s.grainColors.size() == s.atoms.size()) for (float& c : s.grainColors[i]) c *= factor;
    }
    std::vector<double> values;
    for (const auto& b : brightness) values.push_back(b.number());
    s.atomProperty = std::move(values);
    s.atomPropertyName = "Ambient brightness";
    data.notes.push_back("Ambient occlusion from " + std::to_string(static_cast<long long>(number(p, "ambient-occlusion", "directions"))) + " directions");
}

// Keeps a copy of the per-atom property for later comparison.
void freezeProperty(const Json&, PipelineData& data)
{
    Structure& s = data.structure;
    if (s.atomProperty.size() != s.atoms.size() || s.atoms.empty()) throw std::runtime_error("freeze-property needs a per-atom property; add compute-property above");
    s.frozenProperty = s.atomProperty;
    s.frozenPropertyName = s.atomPropertyName;
    data.notes.push_back("Kept " + (s.atomPropertyName.empty() ? std::string("the property") : s.atomPropertyName) + " as frozen");
}

void clearProperty(const Json&, PipelineData& data)
{
    data.structure.atomProperty.clear();
    data.structure.atomPropertyName.clear();
}

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

// The analyses of compute-property: the choice, the scientific tool, the
// property it reports (empty: its first one) and, for analyses of one species
// (Chill+ looks at water oxygens), that element. A new analysis needs one row.
struct PropertyChoice
{
    const char* choice;
    const char* tool;
    const char* property;
    const char* element;
};

const std::vector<PropertyChoice>& propertyChoices()
{
    namespace name = science::propertyName;
    static const std::vector<PropertyChoice> choices = {
        {"structure-type", "structure-type", "", ""},
        {"cna", "cna", "", ""},
        {"ptm", "ptm", "", ""},
        {"orientation", "ptm", name::orientationAngle, ""},
        {"ptm-ordering", "ptm", name::ptmOrdering, ""},
        {"vorotop", "vorotop", "", ""},
        {"diamond", "diamond-structure", "", ""},
        {"centrosymmetry", "centrosymmetry", "", ""},
        {"coordination", "bond-order", name::coordination, ""},
        {"q6", "bond-order", "", ""},
        {"bond-valence", "bond-orders", "", ""},
        {"voronoi-volume", "voronoi", "", ""},
        {"voronoi-coordination", "voronoi", name::voronoiCoordination, ""},
        {"surface", "surface-mesh", "", ""},
        {"planar-faults", "fcc-planar-faults", "", ""},
        {"chill-plus", "chill-plus", "", "O"},
        {"local-entropy", "local-entropy", "", ""},
        {"rings", "find-rings", "", ""},
    };
    return choices;
}

// The choices as the parameter's option list ("structure-type|cna|...").
const char* propertyChoiceOptions()
{
    static const std::string options = [] {
        std::string list;
        for (const auto& c : propertyChoices()) list += std::string(list.empty() ? "" : "|") + c.choice;
        return list;
    }();
    return options.c_str();
}

science::AtomProperty analyse(const PropertyChoice& choice, const Structure& structure, double cutoff)
{
    // Only the atoms of the element, when the structure has others too; the rest get no value.
    std::vector<std::size_t> members;
    for (std::size_t i = 0; *choice.element && i < structure.atoms.size(); ++i)
        if (structure.atoms[i].symbol == choice.element) members.push_back(i);
    if (members.empty() || members.size() == structure.atoms.size())
        return science::analysePerAtom(structure, choice.tool, cutoff, choice.property);
    Structure subset = structure;
    subset.atoms.clear();
    for (std::size_t i : members) subset.atoms.push_back(structure.atoms[i]);
    const science::AtomProperty onSubset = science::analysePerAtom(subset, choice.tool, cutoff, choice.property);
    science::AtomProperty property{onSubset.name, std::vector<double>(structure.atoms.size(), std::numeric_limits<double>::quiet_NaN())};
    for (std::size_t k = 0; k < members.size(); ++k) property.values[members[k]] = onSubset.values[k];
    return property;
}

void computeProperty(const Json& p, PipelineData& data)
{
    const std::string analysis = text(p, "compute-property", "analysis");
    const auto& choices = propertyChoices();
    const auto found = std::find_if(choices.begin(), choices.end(), [&](const PropertyChoice& c) { return analysis == c.choice; });
    if (found == choices.end()) throw std::runtime_error("Unknown per-atom property " + analysis + "; choose from " + propertyChoiceOptions());
    const science::AtomProperty property = analyse(*found, data.structure, number(p, "compute-property", "cutoff"));
    data.structure.atomProperty = property.values;
    data.structure.atomPropertyName = property.name;
    data.notes.push_back(property.name + " computed for " + std::to_string(property.values.size()) + " atoms");
}
}

std::vector<ModifierType> analysisStepTypes()
{
    return {
        {"select-overlapping", "Select overlapping atoms", "Selection",
         "Selects atoms closer than a distance to another atom (minimum image in periodic cells): the later atom of each pair, or both.",
         {{"distance", "Distance (Angstrom)", "float", "1", ""}, {"atoms", "Atoms to select", "choice", "one", "one|both"}, kModeParameter}, selectOverlapping},
        {"select-surface", "Select surface atoms", "Selection",
         "Selects the atoms on the surface of the solid: those a probe sphere of the given radius touches (as in Analysis > Surface Mesh); reports the surface area and solid volume.",
         {{"probe_radius", "Probe radius (Angstrom)", "float", "4", ""}, kModeParameter}, selectSurface},
        {"select-region", "Select a LAMMPS region", "Selection",
         "Selects the atoms inside LAMMPS regions given by their region commands, separated by ';' (e.g. region top block INF INF INF INF 10 INF units box);\nunits lattice multiplies coordinates by the lattice spacing.",
         {{"regions", "Region commands", "string", "region r1 sphere 0 0 0 5 units box", ""}, {"lattice_scale", "Lattice spacing (units lattice)", "float", "1", ""}, kModeParameter}, selectRegion},
        {"assign-color", "Assign colour", "Modification", "Colours the selected (or all) atoms with a red green blue colour (components 0 to 1).",
         {{"color", "Colour (red green blue, 0-1)", "vector", "1 0 0", ""}, {"target", "Atoms", "choice", "selected", "selected|all"}}, assignColor},
        {"ambient-occlusion", "Ambient occlusion", "Modification",
         "Darkens the atoms by how hidden they are from all directions (covalent radii), which brings out pores and grooves; the brightness becomes the per-atom property.",
         {{"intensity", "Intensity (0-1)", "float", "0.7", ""}, {"directions", "Directions", "int", "40", ""}}, ambientOcclusionStep},
        {"clear-property", "Remove per-atom property", "Analysis", "Removes the per-atom property (and its colouring).", {}, clearProperty},
        {"freeze-property", "Freeze per-atom property", "Analysis",
         "Keeps a copy of the per-atom property; later expressions compare with it as the variable frozen, e.g. compute-property cna | freeze-property |\nrelax ... | compute-property cna | select-expression \"property != frozen\" selects the atoms whose structure changed.", {}, freezeProperty},
        {"compute-property", "Compute per-atom property", "Analysis", "Computes a per-atom property for colouring and select-property: structure types (Ackland-Jones, adaptive CNA, template matching with its\norientation angle and chemical ordering, VoroTop, diamond), centrosymmetry, coordination, Steinhardt q6, bond valence, Voronoi volume or\ncoordination, surface atoms, fcc planar faults, Chill+ ice types (on the oxygens), local entropy or rings through each atom. The cutoff is\nused by analyses that need one; the others choose their own neighbour range.",
         {{"analysis", "Property", "choice", "structure-type", propertyChoiceOptions()},
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
}
}
