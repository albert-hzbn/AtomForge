#include "science/Symmetry.h"
#include "util/ElementData.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

#ifdef ATOMS_ENABLE_SPGLIB
#include <spglib.h>
#endif

namespace atomforge::science
{
namespace
{
const char* crystalSystem(int number)
{
    if (number <= 2) return "triclinic";
    if (number <= 15) return "monoclinic";
    if (number <= 74) return "orthorhombic";
    if (number <= 142) return "tetragonal";
    if (number <= 167) return "trigonal";
    if (number <= 194) return "hexagonal";
    return "cubic";
}

double wrap(double x) { return x - std::round(x); }

Structure makeStructure(const Mat3& cell, const std::vector<Vec3>& fractions, const std::vector<int>& numbers)
{
    Structure structure;
    structure.hasUnitCell = true;
    for (int r = 0; r < 3; ++r) structure.cellVectors[r] = cell[r];
    for (std::size_t i = 0; i < fractions.size(); ++i) {
        AtomSite atom;
        atom.atomicNumber = numbers[i];
        atom.symbol = elementSymbol(numbers[i]);
        const Vec3 position = rowTimes(fractions[i], cell);
        atom.x = position[0]; atom.y = position[1]; atom.z = position[2];
        getDefaultElementColor(atom.atomicNumber, atom.r, atom.g, atom.b);
        structure.atoms.push_back(atom);
    }
    return structure;
}

Json latticeParameters(const Mat3& cell)
{
    const double a = norm(cell[0]), b = norm(cell[1]), c = norm(cell[2]);
    const double degrees = 180.0 / std::acos(-1.0);
    Json result = Json::object();
    result["a_A"] = a;
    result["b_A"] = b;
    result["c_A"] = c;
    result["alpha_deg"] = std::acos(std::clamp(dot(cell[1], cell[2]) / (b * c), -1.0, 1.0)) * degrees;
    result["beta_deg"] = std::acos(std::clamp(dot(cell[0], cell[2]) / (a * c), -1.0, 1.0)) * degrees;
    result["gamma_deg"] = std::acos(std::clamp(dot(cell[0], cell[1]) / (a * b), -1.0, 1.0)) * degrees;
    result["volume_A3"] = cellVolume(cell);
    return result;
}

#ifdef ATOMS_ENABLE_SPGLIB
// Lattice rows (Angstrom) -> spglib's column convention.
void toSpglib(const Mat3& cell, double lattice[3][3])
{
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) lattice[r][c] = cell[c][r];
}

Mat3 fromSpglib(const double lattice[3][3])
{
    Mat3 cell{};
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) cell[r][c] = lattice[c][r];
    return cell;
}

// spglib standardization (conventional or primitive, idealized).
Structure standardized(const Mat3& cell, std::vector<Vec3> fractions, std::vector<int> numbers, bool primitive, double symprec)
{
    const std::size_t count = fractions.size();
    // The standardized conventional cell can hold up to four times the input atoms.
    fractions.resize(count * 4);
    numbers.resize(count * 4);
    double lattice[3][3];
    toSpglib(cell, lattice);
    const int n = spg_standardize_cell(lattice, reinterpret_cast<double (*)[3]>(fractions.data()), numbers.data(),
                                       static_cast<int>(count), primitive ? 1 : 0, 0, symprec);
    if (n <= 0) throw std::runtime_error("spglib could not standardize the cell");
    fractions.resize(static_cast<std::size_t>(n));
    numbers.resize(static_cast<std::size_t>(n));
    for (auto& f : fractions) for (double& x : f) x -= std::floor(x);
    return makeStructure(fromSpglib(lattice), fractions, numbers);
}
#endif
}

ToolOutput symmetryAnalysis(const StructureInput& input, const Parameters& p)
{
    const Structure& structure = input.structure;
    if (!structure.hasUnitCell || !(input.pbc[0] && input.pbc[1] && input.pbc[2]) || structure.atoms.empty())
        throw std::runtime_error("Symmetry analysis needs a nonempty periodic structure");
    const double symprec = positive(p.number("symprec_A", 1e-3), "symprec_A");
    const double angleTolerance = p.number("angle_tolerance_deg", -1.0);
    const std::string outputCell = p.has("output_cell") ? p.json("output_cell").string() : "conventional";
    if (outputCell != "conventional" && outputCell != "primitive" && outputCell != "symmetrized")
        throw std::runtime_error("output_cell must be conventional, primitive or symmetrized");
#ifndef ATOMS_ENABLE_SPGLIB
    (void)symprec; (void)angleTolerance;
    throw std::runtime_error("This AtomForge build has no spglib symmetry support; rebuild with spglib to analyse symmetry");
#else
    Mat3 cell{};
    for (int r = 0; r < 3; ++r) cell[r] = structure.cellVectors[r];
    const Mat3 inv = inverse(cell);
    std::vector<Vec3> fractions;
    std::vector<int> numbers;
    for (const auto& atom : structure.atoms) {
        fractions.push_back(rowTimes({atom.x, atom.y, atom.z}, inv));
        numbers.push_back(atom.atomicNumber);
    }
    double lattice[3][3];
    toSpglib(cell, lattice);
    SpglibDataset* dataset = spgat_get_dataset(lattice, reinterpret_cast<double (*)[3]>(fractions.data()), numbers.data(),
                                               static_cast<int>(numbers.size()), symprec, angleTolerance);
    if (!dataset) throw std::runtime_error("spglib could not detect the symmetry; try a larger symprec_A");

    Json result = Json::object();
    const int number = dataset->spacegroup_number;
    result["spacegroup_number"] = number;
    result["international_symbol"] = std::string(dataset->international_symbol);
    result["hall_symbol"] = std::string(dataset->hall_symbol);
    result["hall_number"] = dataset->hall_number;
    result["point_group"] = std::string(dataset->pointgroup_symbol);
    result["crystal_system"] = std::string(crystalSystem(number));
    result["operation_count"] = dataset->n_operations;
    result["symprec_A"] = symprec;

    // Per-atom Wyckoff letters, site symmetries and symmetry-equivalent orbits.
    Json wyckoff = Json::array(), site = Json::array(), orbit = Json::array();
    std::vector<int> orbits(numbers.size()), letters(numbers.size());
    std::map<int, int> orbitIndex;
    for (std::size_t i = 0; i < numbers.size(); ++i) {
        letters[i] = dataset->wyckoffs[i];
        wyckoff.push(std::string(1, static_cast<char>('a' + dataset->wyckoffs[i])));
        site.push(std::string(dataset->site_symmetry_symbols[i]));
        const int representative = dataset->equivalent_atoms[i];
        const auto found = orbitIndex.emplace(representative, static_cast<int>(orbitIndex.size())).first;
        orbits[i] = found->second;
        orbit.push(found->second);
    }
    result["wyckoff_letters"] = wyckoff;
    result["site_symmetry"] = site;
    result["orbit_id"] = orbit;
    result["wyckoff_index"] = toJson(letters);
    result["orbit_count"] = static_cast<int>(orbitIndex.size());
    // One row per orbit: element, Wyckoff position, multiplicity in the input cell, site symmetry.
    Json sites = Json::array();
    for (const auto& [representative, index] : orbitIndex) {
        int multiplicity = 0;
        for (int value : orbits) multiplicity += value == index;
        Json row = Json::object();
        row["element"] = std::string(elementSymbol(numbers[static_cast<std::size_t>(representative)]));
        row["wyckoff"] = std::string(1, static_cast<char>('a' + dataset->wyckoffs[representative]));
        row["site_symmetry"] = std::string(dataset->site_symmetry_symbols[representative]);
        row["multiplicity"] = multiplicity;
        row["representative_atom"] = representative;
        sites.push(row);
    }
    result["sites"] = sites;
    Json operations = Json::array();
    for (int k = 0; k < dataset->n_operations; ++k) {
        Json rotation = Json::array();
        for (int r = 0; r < 3; ++r)
            rotation.push(Json::array({dataset->rotations[k][r][0], dataset->rotations[k][r][1], dataset->rotations[k][r][2]}));
        Json op = Json::object();
        op["rotation"] = rotation;
        op["translation"] = Json::array({dataset->translations[k][0], dataset->translations[k][1], dataset->translations[k][2]});
        operations.push(op);
    }
    result["operations"] = operations;

    // Symmetrized positions in the input cell: average each atom's images under all operations.
    std::vector<Vec3> sum(numbers.size(), Vec3{0, 0, 0});
    double largestShift = 0;
    for (int k = 0; k < dataset->n_operations; ++k) {
        for (std::size_t i = 0; i < numbers.size(); ++i) {
            Vec3 image{};
            for (int r = 0; r < 3; ++r)
                image[static_cast<std::size_t>(r)] = dataset->rotations[k][r][0] * fractions[i][0] + dataset->rotations[k][r][1] * fractions[i][1] +
                                                     dataset->rotations[k][r][2] * fractions[i][2] + dataset->translations[k][r];
            std::size_t best = 0;
            double bestDistance = HUGE_VAL;
            for (std::size_t j = 0; j < numbers.size(); ++j) {
                if (numbers[j] != numbers[i]) continue;
                Vec3 delta{wrap(image[0] - fractions[j][0]), wrap(image[1] - fractions[j][1]), wrap(image[2] - fractions[j][2])};
                const double distance = norm(rowTimes(delta, cell));
                if (distance < bestDistance) { bestDistance = distance; best = j; }
            }
            Vec3 delta{wrap(image[0] - fractions[best][0]), wrap(image[1] - fractions[best][1]), wrap(image[2] - fractions[best][2])};
            sum[best] = add(sum[best], add(fractions[best], delta));
        }
    }
    std::vector<Vec3> symmetrized(numbers.size());
    for (std::size_t i = 0; i < numbers.size(); ++i) {
        symmetrized[i] = scale(sum[i], 1.0 / dataset->n_operations);
        largestShift = std::max(largestShift, norm(rowTimes(sub(symmetrized[i], fractions[i]), cell)));
    }
    result["max_symmetrization_shift_A"] = largestShift;
    spg_free_dataset(dataset);

    ToolOutput output;
    Structure conventional = standardized(cell, fractions, numbers, false, symprec);
    Structure primitive = standardized(cell, fractions, numbers, true, symprec);
    Structure symmetrizedStructure = makeStructure(cell, symmetrized, numbers);
    Mat3 conventionalCell{}, primitiveCell{};
    for (int r = 0; r < 3; ++r) { conventionalCell[r] = conventional.cellVectors[r]; primitiveCell[r] = primitive.cellVectors[r]; }
    result["input_cell"] = latticeParameters(cell);
    result["conventional_cell"] = latticeParameters(conventionalCell);
    result["primitive_cell"] = latticeParameters(primitiveCell);
    result["conventional_atom_count"] = static_cast<int>(conventional.atoms.size());
    result["primitive_atom_count"] = static_cast<int>(primitive.atoms.size());
    result["conventional_structure"] = structureJson(conventional);
    result["primitive_structure"] = structureJson(primitive);
    result["output_cell"] = outputCell;
    std::vector<std::pair<std::string, Structure>> cells = {
        {"conventional", conventional}, {"primitive", primitive}, {"symmetrized", symmetrizedStructure}};
    std::stable_partition(cells.begin(), cells.end(), [&](const auto& entry) { return entry.first != outputCell; });
    Json order = Json::array();
    for (auto& [name, frame] : cells) { order.push(name); output.frames.push_back(frame); }
    result["frame_order"] = order;
    output.result = result;
    return output;
#endif
}
}
