#include "science/ReciprocalPath.h"
#include "science/HpkotData.h"
#include "util/ElementData.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>

#ifdef ATOMS_ENABLE_SPGLIB
#include <spglib.h>
#endif

namespace atomforge::science
{
namespace
{
const double kTwoPi = 2 * std::acos(-1.0);

// Arithmetic over the named cell and k-vector parameters of the HPKOT tables.
class Expression
{
public:
    Expression(const std::string& text, const std::map<std::string, double>& variables) : m_text(text), m_variables(variables) {}

    double evaluate()
    {
        const double value = sum();
        skip();
        if (m_position != m_text.size()) throw std::runtime_error("Invalid k-path expression: " + m_text);
        return value;
    }

private:
    void skip() { while (m_position < m_text.size() && m_text[m_position] == ' ') ++m_position; }
    double sum()
    {
        double value = product();
        while (true) {
            skip();
            if (m_position < m_text.size() && (m_text[m_position] == '+' || m_text[m_position] == '-')) {
                const char op = m_text[m_position++];
                const double rhs = product();
                value = op == '+' ? value + rhs : value - rhs;
            } else return value;
        }
    }
    double product()
    {
        double value = factor();
        while (true) {
            skip();
            if (m_position < m_text.size() && (m_text[m_position] == '*' || m_text[m_position] == '/')) {
                const char op = m_text[m_position++];
                const double rhs = factor();
                value = op == '*' ? value * rhs : value / rhs;
            } else return value;
        }
    }
    double factor()
    {
        skip();
        if (m_position >= m_text.size()) throw std::runtime_error("Invalid k-path expression: " + m_text);
        const char c = m_text[m_position];
        if (c == '-' || c == '+') { ++m_position; const double value = factor(); return c == '-' ? -value : value; }
        if (c == '(') {
            ++m_position;
            const double value = sum();
            skip();
            if (m_position >= m_text.size() || m_text[m_position] != ')') throw std::runtime_error("Invalid k-path expression: " + m_text);
            ++m_position;
            return value;
        }
        if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
            std::size_t used = 0;
            const double value = std::stod(m_text.substr(m_position), &used);
            m_position += used;
            return value;
        }
        const std::size_t start = m_position;
        while (m_position < m_text.size() && (std::isalnum(static_cast<unsigned char>(m_text[m_position])) || m_text[m_position] == '_')) ++m_position;
        const std::string name = m_text.substr(start, m_position - start);
        const auto found = m_variables.find(name);
        if (name.empty() || found == m_variables.end()) throw std::runtime_error("Unknown k-path symbol '" + name + "' in " + m_text);
        return found->second;
    }

    const std::string& m_text;
    const std::map<std::string, double>& m_variables;
    std::size_t m_position = 0;
};

struct CellParameters { double a, b, c, cosAlpha, cosBeta, cosGamma; };

CellParameters parameters(const Mat3& cell)
{
    const double a = norm(cell[0]), b = norm(cell[1]), c = norm(cell[2]);
    return {a, b, c, dot(cell[1], cell[2]) / b / c, dot(cell[0], cell[2]) / a / c, dot(cell[0], cell[1]) / a / b};
}

// Rows are reciprocal vectors with real . reciprocal^T = 2 pi I.
Mat3 reciprocalRows(const Mat3& real)
{
    const Mat3 inv = inverse(real);
    Mat3 result{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) result[i][j] = kTwoPi * inv[j][i];
    return result;
}

Mat3 integerMatrix(std::initializer_list<int> values)
{
    Mat3 result{};
    int k = 0;
    for (int value : values) { result[k / 3][k % 3] = value; ++k; }
    return result;
}

// HPKOT Table 3: conventional (a,b,c) -> primitive (a,b,c) P, and P^-1.
void primitiveMatrices(const std::string& bravais, Mat3& p, Mat3& inverseP)
{
    if (bravais == "cF" || bravais == "oF") {
        p = integerMatrix({0, 1, 1, 1, 0, 1, 1, 1, 0}); inverseP = integerMatrix({-1, 1, 1, 1, -1, 1, 1, 1, -1});
        for (auto& row : p) row = scale(row, 0.5);
    } else if (bravais == "cI" || bravais == "tI" || bravais == "oI") {
        p = integerMatrix({-1, 1, 1, 1, -1, 1, 1, 1, -1}); inverseP = integerMatrix({0, 1, 1, 1, 0, 1, 1, 1, 0});
        for (auto& row : p) row = scale(row, 0.5);
    } else if (bravais == "hR") {
        p = integerMatrix({2, -1, -1, 1, 1, -2, 1, 1, 1}); inverseP = integerMatrix({1, 0, 1, -1, 1, 1, 0, -1, 1});
        for (auto& row : p) row = scale(row, 1.0 / 3.0);
    } else if (bravais == "oC") {
        p = integerMatrix({1, 1, 0, -1, 1, 0, 0, 0, 2}); inverseP = integerMatrix({1, -1, 0, 1, 1, 0, 0, 0, 1});
        for (auto& row : p) row = scale(row, 0.5);
    } else if (bravais == "oA") {
        p = integerMatrix({0, 0, 2, 1, 1, 0, -1, 1, 0}); inverseP = integerMatrix({0, 1, -1, 0, 1, 1, 1, 0, 0});
        for (auto& row : p) row = scale(row, 0.5);
    } else if (bravais == "mC") {
        p = integerMatrix({1, -1, 0, 1, 1, 0, 0, 0, 2}); inverseP = integerMatrix({1, 1, 0, -1, 1, 0, 0, 0, 1});
        for (auto& row : p) row = scale(row, 0.5);
    } else {
        p = identity(); inverseP = identity();
    }
}

bool centrosymmetric(const std::string& pointGroup)
{
    static const char* groups[] = {"-1", "2/m", "mmm", "4/m", "4/mmm", "-3", "-3m", "6/m", "6/mmm", "m-3", "m-3m"};
    for (const char* group : groups)
        if (pointGroup == group) return true;
    return false;
}
}

ToolOutput reciprocalPath(const StructureInput& input, double spacing, double symprec, bool timeReversal)
{
    positive(spacing, "spacing_inv_A");
    positive(symprec, "symprec_A");
    const Structure& structure = input.structure;
    if (!structure.hasUnitCell || !(input.pbc[0] && input.pbc[1] && input.pbc[2]) || structure.atoms.empty())
        throw std::runtime_error("Reciprocal paths require a nonempty periodic structure");
#ifndef ATOMS_ENABLE_SPGLIB
    (void)timeReversal;
    throw std::runtime_error("This AtomForge build has no spglib symmetry support; rebuild with spglib to compute symmetry paths");
#else
    Mat3 cell{};
    for (int r = 0; r < 3; ++r) cell[r] = structure.cellVectors[r];
    cellMatrix([&] { NdArray m({3, 3}); for (std::size_t r = 0; r < 3; ++r) for (std::size_t c = 0; c < 3; ++c) m(r, c) = cell[r][c]; return m; }());
    const Mat3 inv = inverse(cell);
    double lattice[3][3];
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) lattice[r][c] = cell[c][r];  // spglib: column vectors
    std::vector<std::array<double, 3>> fractions;
    std::vector<int> types;
    for (const auto& atom : structure.atoms) {
        fractions.push_back(rowTimes({atom.x, atom.y, atom.z}, inv));
        types.push_back(atom.atomicNumber);
    }
    SpglibDataset* dataset = spg_get_dataset(lattice, reinterpret_cast<double (*)[3]>(fractions.data()), types.data(),
                                             static_cast<int>(types.size()), symprec);
    if (!dataset) throw std::runtime_error("Spglib could not detect the symmetry of the system");
    Mat3 conventional{};
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) conventional[r][c] = dataset->std_lattice[c][r];
    std::vector<Vec3> conventionalPositions;
    std::vector<int> conventionalTypes;
    for (int i = 0; i < dataset->n_std_atoms; ++i) {
        conventionalPositions.push_back({dataset->std_positions[i][0], dataset->std_positions[i][1], dataset->std_positions[i][2]});
        conventionalTypes.push_back(dataset->std_types[i]);
    }
    const int number = dataset->spacegroup_number;
    const std::string international = dataset->international_symbol;
    const std::string pointGroup = dataset->pointgroup_symbol;
    spg_free_dataset(dataset);

    const char family = number <= 2 ? 'a' : number <= 15 ? 'm' : number <= 74 ? 'o' : number <= 142 ? 't' : number <= 194 ? 'h' : 'c';
    const std::string bravais = std::string(1, family) + international.substr(0, 1);
    const bool hasInversion = centrosymmetric(pointGroup);
    CellParameters cp = parameters(conventional);
    double a = cp.a, b = cp.b, c = cp.c, cosAlpha = cp.cosAlpha, cosBeta = cp.cosBeta, cosGamma = cp.cosGamma;
    std::string extended;
    if (bravais == "cP") extended = number <= 206 ? "cP1" : "cP2";
    else if (bravais == "cF") extended = number <= 206 ? "cF1" : "cF2";
    else if (bravais == "cI") extended = "cI1";
    else if (bravais == "tP") extended = "tP1";
    else if (bravais == "tI") extended = c <= a ? "tI1" : "tI2";
    else if (bravais == "oP") extended = "oP1";
    else if (bravais == "oF") {
        if (1 / (a * a) > 1 / (b * b) + 1 / (c * c)) extended = "oF1";
        else if (1 / (c * c) > 1 / (a * a) + 1 / (b * b)) extended = "oF2";
        else extended = "oF3";
    } else if (bravais == "oI") {
        // The longest conventional vector selects oI1 (c), oI2 (a) or oI3 (b).
        std::vector<std::pair<double, int>> sorted = {{c, 1}, {b, 3}, {a, 2}};
        std::sort(sorted.begin(), sorted.end());
        extended = "oI" + std::to_string(sorted.back().second);
    } else if (bravais == "oC") extended = a <= b ? "oC1" : "oC2";
    else if (bravais == "oA") extended = b <= c ? "oA1" : "oA2";
    else if (bravais == "hP") {
        static const int hp1[] = {143, 144, 145, 146, 147, 148, 149, 151, 153, 157, 159, 160, 161, 162, 163};
        extended = std::find(std::begin(hp1), std::end(hp1), number) != std::end(hp1) ? "hP1" : "hP2";
    } else if (bravais == "hR") extended = std::sqrt(3.0) * a <= std::sqrt(2.0) * c ? "hR1" : "hR2";
    else if (bravais == "mP") extended = "mP1";
    else if (bravais == "mC") {
        const double sinBeta2 = 1.0 - cosBeta * cosBeta;
        if (b < a * std::sqrt(sinBeta2)) extended = "mC1";
        else extended = -a * cosBeta / c + a * a * sinBeta2 / (b * b) <= 1.0 ? "mC2" : "mC3";
    } else if (bravais == "aP") {
        // Reciprocal-space Niggli reduction, then the HPKOT all-acute/all-obtuse choice.
        Mat3 reciprocal = reciprocalRows(conventional);
        double columns[3][3];
        for (int r = 0; r < 3; ++r)
            for (int k = 0; k < 3; ++k) columns[r][k] = reciprocal[k][r];
        if (!spg_niggli_reduce(columns, 1e-5)) throw std::runtime_error("Niggli reduction failed for the triclinic cell");
        for (int r = 0; r < 3; ++r)
            for (int k = 0; k < 3; ++k) reciprocal[k][r] = columns[r][k];
        const Mat3 real2 = reciprocalRows(reciprocal);
        const auto k2 = parameters(reciprocal);
        const double conditions[3] = {std::abs(k2.b * k2.c * k2.cosAlpha), std::abs(k2.c * k2.a * k2.cosBeta), std::abs(k2.a * k2.b * k2.cosGamma)};
        const Mat3 m2Choices[3] = {integerMatrix({0, 0, 1, 1, 0, 0, 0, 1, 0}), integerMatrix({0, 1, 0, 0, 0, 1, 1, 0, 0}), identity()};
        const int smallest = static_cast<int>(std::min_element(conditions, conditions + 3) - conditions);
        const Mat3 real3 = multiply(transpose(m2Choices[smallest]), real2);
        const auto k3 = parameters(reciprocalRows(real3));
        const bool pa = k3.cosAlpha > 0, pb = k3.cosBeta > 0, pg = k3.cosGamma > 0;
        Mat3 m3 = identity();
        if ((pa && !pb && !pg) || (!pa && pb && pg)) m3 = integerMatrix({1, 0, 0, 0, -1, 0, 0, 0, -1});
        else if ((!pa && pb && !pg) || (pa && !pb && pg)) m3 = integerMatrix({-1, 0, 0, 0, 1, 0, 0, 0, -1});
        else if ((!pa && !pb && pg) || (pa && pb && !pg)) m3 = integerMatrix({-1, 0, 0, 0, -1, 0, 0, 0, 1});
        const Mat3 final = multiply(transpose(m3), real3);
        const auto kf = parameters(reciprocalRows(final));
        if (kf.cosAlpha <= 0 && kf.cosBeta <= 0 && kf.cosGamma <= 0) extended = "aP2";
        else if (kf.cosAlpha >= 0 && kf.cosBeta >= 0 && kf.cosGamma >= 0) extended = "aP3";
        else throw std::runtime_error("Unexpected aP triclinic lattice: neither all-obtuse nor all-acute");
        const Mat3 finalInverse = inverse(final);
        for (auto& position : conventionalPositions) position = rowTimes(rowTimes(position, conventional), finalInverse);
        conventional = final;
        cp = parameters(final);
        a = cp.a; b = cp.b; c = cp.c; cosAlpha = cp.cosAlpha; cosBeta = cp.cosBeta; cosGamma = cp.cosGamma;
    } else throw std::runtime_error("Unknown Bravais lattice " + bravais + " for space group " + std::to_string(number));

    // Primitive cell in the HPKOT convention (positions not wrapped).
    Mat3 p{}, inverseP{};
    primitiveMatrices(bravais, p, inverseP);
    const Mat3 primitive = multiply(transpose(p), conventional);
    std::vector<Vec3> mapped;
    for (const auto& position : conventionalPositions) mapped.push_back(rowTimes(position, transpose(inverseP)));
    const int ratio = static_cast<int>(std::lround(determinant(inverseP)));
    std::vector<int> group(mapped.size(), -1);
    std::vector<std::size_t> representatives;
    for (std::size_t i = 0; i < mapped.size(); ++i) {
        if (group[i] >= 0) continue;
        int members = 0;
        for (std::size_t j = i; j < mapped.size(); ++j) {
            bool same = true;
            for (int k = 0; k < 3; ++k) {
                const double d = mapped[j][k] - mapped[i][k];
                same = same && std::abs(d - std::round(d)) < 1e-6;
            }
            if (!same) continue;
            if (group[j] >= 0 || conventionalTypes[j] != conventionalTypes[i])
                throw std::runtime_error("Problem creating the primitive cell: overlapping atoms of different type");
            group[j] = static_cast<int>(representatives.size());
            ++members;
        }
        if (members != ratio) throw std::runtime_error("Problem creating the primitive cell from the conventional cell");
        representatives.push_back(i);
    }
    Structure primitiveStructure;
    primitiveStructure.hasUnitCell = true;
    for (int r = 0; r < 3; ++r) primitiveStructure.cellVectors[r] = primitive[r];
    for (std::size_t index : representatives) {
        AtomSite atom;
        atom.atomicNumber = conventionalTypes[index];
        atom.symbol = elementSymbol(atom.atomicNumber);
        const Vec3 position = rowTimes(mapped[index], primitive);
        atom.x = position[0]; atom.y = position[1]; atom.z = position[2];
        getDefaultElementColor(atom.atomicNumber, atom.r, atom.g, atom.b);
        primitiveStructure.atoms.push_back(atom);
    }

    const hpkot::Lattice* table = nullptr;
    for (const auto& candidate : hpkot::lattices)
        if (extended == candidate.name) table = &candidate;
    if (!table) throw std::runtime_error("Missing HPKOT data for " + extended);
    std::map<std::string, double> variables = {{"a", a}, {"b", b}, {"c", c}, {"cosalpha", cosAlpha}, {"cosbeta", cosBeta},
                                               {"cosgamma", cosGamma}, {"sinbeta", std::sqrt(1.0 - cosBeta * cosBeta)}};
    for (int i = 0; i < table->parameterCount; ++i)
        variables[table->parameters[i].name] = Expression(table->parameters[i].expression, variables).evaluate();
    std::vector<std::pair<std::string, Vec3>> points;
    for (int i = 0; i < table->pointCount; ++i) {
        const auto& point = table->points[i];
        points.push_back({point.label, {Expression(point.x, variables).evaluate(), Expression(point.y, variables).evaluate(),
                                        Expression(point.z, variables).evaluate()}});
    }
    std::vector<std::pair<std::string, std::string>> path;
    for (int i = 0; i < table->segmentCount; ++i) path.push_back({table->path[i].from, table->path[i].to});
    // Without inversion or time reversal, k and -k are inequivalent: add the -k path.
    const bool augmented = !hasInversion && !timeReversal;
    if (augmented) {
        const auto original = points;
        for (const auto& [label, coordinates] : original)
            if (label != "GAMMA") points.push_back({label + "'", scale(coordinates, -1)});
        const auto originalPath = path;
        for (const auto& [from, to] : originalPath)
            path.push_back({from == "GAMMA" ? from : from + "'", to == "GAMMA" ? to : to + "'"});
    }
    auto coordinate = [&](const std::string& label) {
        for (const auto& point : points)
            if (point.first == label) return point.second;
        throw std::runtime_error("Undefined k-point " + label);
    };
    const Mat3 reciprocal = reciprocalRows(primitive);
    std::vector<Vec3> kFractional;
    std::vector<std::string> labels;
    std::vector<double> distance;
    Json segments = Json::array();
    double previous = 0;
    for (const auto& [from, to] : path) {
        const Vec3 start = coordinate(from), stop = coordinate(to);
        const double length = norm(sub(rowTimes(stop, reciprocal), rowTimes(start, reciprocal)));
        const int samples = std::max(2, static_cast<int>(length / spacing));
        long long segmentStart = static_cast<long long>(labels.size());
        for (int i = 0; i < samples; ++i) {
            if (i == 0 && !labels.empty() && labels.back() == from) { --segmentStart; continue; }
            const double fraction = static_cast<double>(i) / (samples - 1);
            kFractional.push_back(add(start, scale(sub(stop, start), fraction)));
            labels.push_back(i == 0 ? from : (i == samples - 1 ? to : ""));
            distance.push_back(previous + length * fraction);
        }
        previous += length;
        segments.push(Json::array({Json(segmentStart), Json(static_cast<long long>(labels.size()))}));
    }
    ToolOutput output;
    output.frames.push_back(primitiveStructure);
    Json kFractionalJson = Json::array(), kCartesian = Json::array(), labelJson = Json::array(), pathJson = Json::array();
    for (const auto& k : kFractional) {
        kFractionalJson.push(toJson(k));
        kCartesian.push(toJson(rowTimes(k, reciprocal)));
    }
    for (const auto& label : labels) labelJson.push(label);
    for (const auto& [from, to] : path) pathJson.push(Json::array({Json(from), Json(to)}));
    Json special = Json::object();
    for (const auto& [label, coordinates] : points) special[label] = toJson(coordinates);
    Json result = Json::object();
    result["primitive_structure"] = structureJson(primitiveStructure);
    result["spacegroup_number"] = number;
    result["spacegroup_symbol"] = international;
    result["bravais_lattice_extended"] = extended;
    result["has_inversion_symmetry"] = hasInversion;
    result["augmented_path"] = augmented;
    result["path"] = pathJson;
    result["special_points_fractional"] = special;
    result["kpoints_fractional"] = kFractionalJson;
    result["kpoints_inv_A"] = kCartesian;
    result["distance_inv_A"] = toJson(distance);
    result["labels"] = labelJson;
    result["segments"] = segments;
    result["reciprocal_lattice_inv_A"] = toJson(reciprocal);
    output.result = result;
    return output;
#endif
}

ToolOutput dftInputs(const StructureInput& structure, int pointsPerSegment, double symprecA, bool timeReversal)
{
    if (pointsPerSegment < 2 || pointsPerSegment > 1000) throw std::runtime_error("points_per_segment must be between 2 and 1000");
    ToolOutput path = reciprocalPath(structure, 0.05, symprecA, timeReversal);
    const Structure& primitive = path.frames.front();
    const Json& special = path.result.at("special_points_fractional");
    std::vector<std::pair<std::string, std::string>> segments;
    for (const auto& segment : path.result.at("path").items())
        segments.push_back({segment.items()[0].string(), segment.items()[1].string()});
    auto coordinate = [&](const std::string& label) -> const Json& { return special.at(label); };
    auto qeLabel = [](std::string label) { return label == "GAMMA" ? std::string("G") : label; };
    std::ostringstream kpoints, poscar, qe;
    kpoints << std::setprecision(10) << std::fixed;
    poscar << std::setprecision(12) << std::fixed;
    qe << std::setprecision(10) << std::fixed;
    const std::string title = "HPKOT path, space group " + path.result.at("spacegroup_symbol").string() + " (" +
                              std::to_string(static_cast<int>(path.result.at("spacegroup_number").number())) + ")";
    kpoints << title << "\n" << pointsPerSegment << "\nLine-mode\nReciprocal\n";
    for (std::size_t s = 0; s < segments.size(); ++s) {
        if (s) kpoints << "\n";
        for (const auto& label : {segments[s].first, segments[s].second}) {
            const auto& k = coordinate(label).items();
            kpoints << k[0].number() << ' ' << k[1].number() << ' ' << k[2].number() << " ! " << (label == "GAMMA" ? std::string("\\Gamma") : label) << "\n";
        }
    }
    // POSCAR (VASP 5), species grouped in order of first appearance.
    std::vector<std::string> species;
    for (const auto& atom : primitive.atoms)
        if (std::find(species.begin(), species.end(), atom.symbol) == species.end()) species.push_back(atom.symbol);
    poscar << title << "\n1.0\n";
    for (const auto& row : primitive.cellVectors) poscar << "  " << row[0] << ' ' << row[1] << ' ' << row[2] << "\n";
    for (const auto& s : species) poscar << ' ' << s;
    poscar << "\n";
    for (const auto& s : species) poscar << ' ' << std::count_if(primitive.atoms.begin(), primitive.atoms.end(), [&](const AtomSite& a) { return a.symbol == s; });
    poscar << "\nDirect\n";
    Mat3 cell{};
    for (int r = 0; r < 3; ++r) cell[r] = primitive.cellVectors[r];
    for (const auto& s : species)
        for (const auto& atom : primitive.atoms)
            if (atom.symbol == s) {
                const Vec3 f = fractional({atom.x, atom.y, atom.z}, cell);
                poscar << "  " << f[0] << ' ' << f[1] << ' ' << f[2] << "\n";
            }
    // Quantum ESPRESSO: a point weight is the number of steps to the next point; 1 jumps.
    std::vector<std::pair<std::string, int>> sequence;
    for (std::size_t s = 0; s < segments.size(); ++s) {
        if (sequence.empty() || sequence.back().first != segments[s].first) {
            if (!sequence.empty()) sequence.back().second = 1;
            sequence.push_back({segments[s].first, pointsPerSegment});
        } else sequence.back().second = pointsPerSegment;
        sequence.push_back({segments[s].second, 1});
    }
    qe << "CELL_PARAMETERS angstrom\n";
    for (const auto& row : primitive.cellVectors) qe << "  " << row[0] << ' ' << row[1] << ' ' << row[2] << "\n";
    qe << "ATOMIC_SPECIES\n";
    for (const auto& s : species) qe << "  " << s << ' ' << atomicMass(s) << ' ' << s << ".UPF\n";
    qe << "ATOMIC_POSITIONS crystal\n";
    for (const auto& atom : primitive.atoms) {
        const Vec3 f = fractional({atom.x, atom.y, atom.z}, cell);
        qe << "  " << atom.symbol << ' ' << f[0] << ' ' << f[1] << ' ' << f[2] << "\n";
    }
    qe << "K_POINTS crystal_b\n" << sequence.size() << "\n";
    for (const auto& [label, weight] : sequence) {
        const auto& k = coordinate(label).items();
        qe << "  " << k[0].number() << ' ' << k[1].number() << ' ' << k[2].number() << ' ' << weight << " ! " << qeLabel(label) << "\n";
    }
    Json files = Json::object();
    files["KPOINTS"] = kpoints.str();
    files["POSCAR"] = poscar.str();
    files["qe_band_cards.in"] = qe.str();
    Json result = path.result;
    result["points_per_segment"] = pointsPerSegment;
    result["files"] = files;
    path.result = result;
    return path;
}
}
