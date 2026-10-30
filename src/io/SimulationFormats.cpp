#include "io/SimulationFormats.h"

#include "util/ElementData.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace atomforge::io
{
namespace
{
using Vec3 = std::array<double, 3>;
constexpr double kBohr = 0.529177210903;  // Angstrom

std::string lower(std::string s)
{
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string canonical(const std::string& format)
{
    const std::string f = lower(format);
    if (f == "pwscf" || f == "pwi" || f == "in") return "pwi";
    return f;
}

int numberOf(const std::string& symbol)
{
    for (int z = 1; z <= 118; ++z)
        if (symbol == elementSymbol(z)) return z;
    return 0;
}

// An element symbol from a label such as "Fe1", "fe" or "Fe_up".
std::string elementFrom(const std::string& label)
{
    std::string letters;
    for (char c : label) {
        if (!std::isalpha(static_cast<unsigned char>(c))) break;
        letters += c;
        if (letters.size() == 2) break;
    }
    if (letters.empty()) return "";
    letters[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(letters[0])));
    if (letters.size() == 2) {
        letters[1] = static_cast<char>(std::tolower(static_cast<unsigned char>(letters[1])));
        if (numberOf(letters)) return letters;
        letters.resize(1);
    }
    return numberOf(letters) ? letters : "";
}

AtomSite atom(const std::string& label, const Vec3& position)
{
    AtomSite a;
    a.symbol = elementFrom(label);
    if (a.symbol.empty()) {
        // An atomic number (XSF).
        char* end = nullptr;
        const long z = std::strtol(label.c_str(), &end, 10);
        if (*end == 0 && z >= 1 && z <= 118) a.symbol = elementSymbol(static_cast<int>(z));
    }
    if (a.symbol.empty()) throw std::runtime_error("Unknown element '" + label + "'");
    a.atomicNumber = numberOf(a.symbol);
    a.x = position[0];
    a.y = position[1];
    a.z = position[2];
    getDefaultElementColor(a.atomicNumber, a.r, a.g, a.b);
    return a;
}

std::vector<std::string> species(const Structure& s)
{
    std::vector<std::string> out;
    for (const auto& a : s.atoms)
        if (std::find(out.begin(), out.end(), a.symbol) == out.end()) out.push_back(a.symbol);
    return out;
}

Vec3 fromFractional(const Structure& s, const Vec3& f)
{
    Vec3 x{};
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) x[c] += f[r] * s.cellVectors[r][c];
    return x;
}

double length(const Vec3& v) { return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }

void number(std::ostream& out, double value) { out << std::setw(18) << std::fixed << std::setprecision(10) << value; }

void vectorLine(std::ostream& out, const Vec3& v, const char* indent = "  ")
{
    out << indent;
    for (double x : v) number(out, x);
    out << '\n';
}

// ---------------------------------------------------------------- writers
void writePwi(const Structure& s, std::ostream& out)
{
    const auto kinds = species(s);
    out << "&CONTROL\n  calculation = 'scf'\n  prefix = 'atomforge'\n  pseudo_dir = './'\n  outdir = './out'\n/\n"
        << "&SYSTEM\n  ibrav = 0\n  nat = " << s.atoms.size() << "\n  ntyp = " << kinds.size() << "\n  ecutwfc = 40\n  ecutrho = 320\n/\n"
        << "&ELECTRONS\n/\n";
    bool fixed = false;
    for (const auto& a : s.atoms) fixed = fixed || a.fixed;
    if (fixed) out << "&IONS\n/\n";
    out << "ATOMIC_SPECIES\n";
    for (const auto& k : kinds) out << "  " << std::left << std::setw(3) << k << std::right << ' ' << std::fixed << std::setprecision(4)
                                    << elementAtomicMass(numberOf(k)) << ' ' << k << ".UPF\n";
    out << "CELL_PARAMETERS angstrom\n";
    for (const auto& row : s.cellVectors) vectorLine(out, {row[0], row[1], row[2]});
    out << "ATOMIC_POSITIONS angstrom\n";
    for (const auto& a : s.atoms) {
        out << "  " << std::left << std::setw(3) << a.symbol << std::right;
        number(out, a.x); number(out, a.y); number(out, a.z);
        // if_pos: 0 keeps a coordinate fixed.
        if (fixed) out << ' ' << ((a.fixed & 1) ? 0 : 1) << ' ' << ((a.fixed & 2) ? 0 : 1) << ' ' << ((a.fixed & 4) ? 0 : 1);
        out << '\n';
    }
    // About 30 k-points per inverse Angstrom of reciprocal length along periodic axes.
    out << "K_POINTS automatic\n ";
    for (int k = 0; k < 3; ++k) {
        const double len = length({s.cellVectors[k][0], s.cellVectors[k][1], s.cellVectors[k][2]});
        out << ' ' << (s.pbc[static_cast<std::size_t>(k)] ? std::max(1, static_cast<int>(std::ceil(30.0 / len))) : 1);
    }
    out << " 0 0 0\n";
}

void writeCell(const Structure& s, std::ostream& out)
{
    out << "%BLOCK LATTICE_CART\nANG\n";
    for (const auto& row : s.cellVectors) vectorLine(out, {row[0], row[1], row[2]});
    out << "%ENDBLOCK LATTICE_CART\n\n%BLOCK POSITIONS_ABS\nANG\n";
    for (const auto& a : s.atoms) {
        out << "  " << std::left << std::setw(3) << a.symbol << std::right;
        number(out, a.x); number(out, a.y); number(out, a.z);
        out << '\n';
    }
    out << "%ENDBLOCK POSITIONS_ABS\n";
    // Fixed coordinates as linear constraints: constraint, species, index within the species, direction.
    std::vector<std::string> lines;
    std::map<std::string, int> counter;
    int constraint = 0;
    for (const auto& a : s.atoms) {
        const int index = ++counter[a.symbol];
        for (int k = 0; k < 3; ++k)
            if (a.fixed & (1 << k)) {
                std::ostringstream line;
                line << "  " << ++constraint << ' ' << a.symbol << ' ' << index << ' ' << (k == 0) << ' ' << (k == 1) << ' ' << (k == 2);
                lines.push_back(line.str());
            }
    }
    if (!lines.empty()) {
        out << "\n%BLOCK IONIC_CONSTRAINTS\n";
        for (const auto& l : lines) out << l << '\n';
        out << "%ENDBLOCK IONIC_CONSTRAINTS\n";
    }
}

void writeGen(const Structure& s, std::ostream& out)
{
    const auto kinds = species(s);
    out << s.atoms.size() << ' ' << (s.hasUnitCell ? 'S' : 'C') << '\n';
    for (const auto& k : kinds) out << k << ' ';
    out << '\n';
    for (std::size_t i = 0; i < s.atoms.size(); ++i) {
        const auto& a = s.atoms[i];
        out << std::setw(6) << i + 1 << ' ' << std::setw(3)
            << (std::find(kinds.begin(), kinds.end(), a.symbol) - kinds.begin()) + 1;
        number(out, a.x); number(out, a.y); number(out, a.z);
        out << '\n';
    }
    if (s.hasUnitCell) {
        vectorLine(out, {s.cellOffset[0], s.cellOffset[1], s.cellOffset[2]});
        for (const auto& row : s.cellVectors) vectorLine(out, {row[0], row[1], row[2]});
    }
}

void writeXsf(const Structure& s, std::ostream& out)
{
    int periodic = 0;
    for (int k = 0; k < 3; ++k) periodic += s.hasUnitCell && s.pbc[static_cast<std::size_t>(k)];
    if (!s.hasUnitCell) {
        out << "ATOMS\n";
    } else {
        // XCrySDen's dimensionality keywords; a slab is periodic along a and b, a polymer along a.
        const char* kind = periodic == 3 ? "CRYSTAL" : periodic == 2 && s.pbc[0] && s.pbc[1] ? "SLAB"
                         : periodic == 1 && s.pbc[0] ? "POLYMER" : periodic == 0 ? "MOLECULE" : "CRYSTAL";
        out << kind << "\nPRIMVEC\n";
        for (const auto& row : s.cellVectors) vectorLine(out, {row[0], row[1], row[2]});
        out << "PRIMCOORD\n" << s.atoms.size() << " 1\n";
    }
    for (const auto& a : s.atoms) {
        out << std::setw(4) << a.atomicNumber;
        number(out, a.x); number(out, a.y); number(out, a.z);
        out << '\n';
    }
}

void writeCp2k(const Structure& s, std::ostream& out)
{
    std::string periodic;
    for (int k = 0; k < 3; ++k)
        if (s.pbc[static_cast<std::size_t>(k)]) periodic += "XYZ"[k];
    if (periodic.empty()) periodic = "NONE";
    out << "&SUBSYS\n  &CELL\n";
    for (int k = 0; k < 3; ++k) {
        out << "    " << "ABC"[k];
        for (double x : s.cellVectors[k]) number(out, x);
        out << '\n';
    }
    out << "    PERIODIC " << periodic << "\n  &END CELL\n  &COORD\n";
    for (const auto& a : s.atoms) {
        out << "    " << std::left << std::setw(3) << a.symbol << std::right;
        number(out, a.x); number(out, a.y); number(out, a.z);
        out << '\n';
    }
    out << "  &END COORD\n";
    for (const auto& k : species(s))
        out << "  &KIND " << k << "\n    ELEMENT " << k << "\n    BASIS_SET DZVP-MOLOPT-SR-GTH\n    POTENTIAL GTH-PBE\n  &END KIND\n";
    out << "&END SUBSYS\n";
}

// ---------------------------------------------------------------- readers
std::vector<std::string> linesOf(const std::string& text)
{
    std::vector<std::string> lines;
    std::istringstream in(text);
    for (std::string line; std::getline(in, line);) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(line);
    }
    return lines;
}

std::vector<std::string> tokens(const std::string& line)
{
    std::istringstream in(line);
    std::vector<std::string> out;
    for (std::string w; in >> w;) out.push_back(w);
    return out;
}

double toNumber(const std::string& w)
{
    std::string t = w;
    std::replace(t.begin(), t.end(), 'd', 'e');
    std::replace(t.begin(), t.end(), 'D', 'e');
    std::size_t used = 0;
    const double v = std::stod(t, &used);
    if (used != t.size()) throw std::runtime_error("'" + w + "' is not a number");
    return v;
}

Vec3 vector3(const std::vector<std::string>& t, std::size_t from)
{
    if (t.size() < from + 3) throw std::runtime_error("Expected three numbers");
    return {toNumber(t[from]), toNumber(t[from + 1]), toNumber(t[from + 2])};
}

void setCell(Structure& s, const std::array<Vec3, 3>& rows)
{
    s.hasUnitCell = true;
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) s.cellVectors[r][c] = rows[r][c];
}

void readXsf(const std::string& text, Structure& s)
{
    const auto lines = linesOf(text);
    std::array<Vec3, 3> rows{};
    bool haveCell = false;
    std::array<bool, 3> pbc{true, true, true};
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const auto t = tokens(lines[i]);
        if (t.empty() || t[0][0] == '#') continue;
        const std::string key = t[0];
        if (key == "SLAB") pbc = {true, true, false};
        else if (key == "POLYMER") pbc = {true, false, false};
        else if (key == "MOLECULE") pbc = {false, false, false};
        else if (key == "PRIMVEC") {
            for (int r = 0; r < 3; ++r) rows[r] = vector3(tokens(lines.at(i + 1 + r)), 0);
            haveCell = true;
            i += 3;
        } else if (key == "PRIMCOORD") {
            const auto header = tokens(lines.at(i + 1));
            const long n = std::stol(header.at(0));
            for (long k = 0; k < n; ++k) {
                const auto a = tokens(lines.at(i + 2 + k));
                s.atoms.push_back(atom(a.at(0), vector3(a, 1)));
            }
            i += 1 + n;
        } else if (key == "ATOMS") {
            for (std::size_t k = i + 1; k < lines.size(); ++k) {
                const auto a = tokens(lines[k]);
                if (a.size() < 4 || (std::isalpha(static_cast<unsigned char>(a[0][0])) && elementFrom(a[0]).empty())) break;
                try { s.atoms.push_back(atom(a[0], vector3(a, 1))); } catch (const std::exception&) { break; }
                i = k;
            }
        }
    }
    if (haveCell) {
        setCell(s, rows);
        s.pbc = pbc;
    }
}

void readGen(const std::string& text, Structure& s)
{
    std::vector<std::string> lines;
    for (const auto& line : linesOf(text)) {
        const auto t = tokens(line);
        if (!t.empty() && t[0][0] != '#') lines.push_back(line);
    }
    const auto header = tokens(lines.at(0));
    const long n = std::stol(header.at(0));
    const char type = static_cast<char>(std::toupper(static_cast<unsigned char>(header.at(1)[0])));
    const auto kinds = tokens(lines.at(1));
    std::vector<std::pair<std::string, Vec3>> raw;
    for (long k = 0; k < n; ++k) {
        const auto a = tokens(lines.at(2 + k));
        const long kind = std::stol(a.at(1));
        if (kind < 1 || kind > static_cast<long>(kinds.size())) throw std::runtime_error("Species index out of range");
        raw.push_back({kinds[static_cast<std::size_t>(kind - 1)], vector3(a, 2)});
    }
    if (type == 'S' || type == 'F') {
        const Vec3 origin = vector3(tokens(lines.at(2 + n)), 0);
        std::array<Vec3, 3> rows{};
        for (int r = 0; r < 3; ++r) rows[r] = vector3(tokens(lines.at(3 + n + r)), 0);
        setCell(s, rows);
        s.cellOffset = origin;
    }
    for (const auto& [symbol, position] : raw)
        s.atoms.push_back(atom(symbol, type == 'F' ? fromFractional(s, position) : position));
}

// The content of a %BLOCK name ... %ENDBLOCK name (case-insensitive).
std::vector<std::vector<std::string>> block(const std::vector<std::string>& lines, const std::string& name, bool& found)
{
    std::vector<std::vector<std::string>> out;
    found = false;
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const auto t = tokens(lines[i]);
        if (t.size() >= 2 && lower(t[0]) == "%block" && lower(t[1]) == name) {
            found = true;
            for (std::size_t k = i + 1; k < lines.size(); ++k) {
                const auto row = tokens(lines[k]);
                if (!row.empty() && lower(row[0]) == "%endblock") return out;
                if (!row.empty() && row[0][0] != '!' && row[0][0] != '#') out.push_back(row);
            }
        }
    }
    return out;
}

double unitScale(std::vector<std::vector<std::string>>& rows)
{
    if (!rows.empty() && rows.front().size() == 1 && std::isalpha(static_cast<unsigned char>(rows.front()[0][0]))) {
        const std::string unit = lower(rows.front()[0]);
        rows.erase(rows.begin());
        if (unit == "bohr" || unit == "a0") return kBohr;
        if (unit == "nm") return 10.0;
    }
    return 1.0;
}

void readCell(const std::string& text, Structure& s)
{
    const auto lines = linesOf(text);
    bool found = false;
    auto lattice = block(lines, "lattice_cart", found);
    if (found) {
        const double scale = unitScale(lattice);
        std::array<Vec3, 3> rows{};
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) rows[r][c] = toNumber(lattice.at(r).at(c)) * scale;
        setCell(s, rows);
    } else {
        auto abc = block(lines, "lattice_abc", found);
        if (!found) throw std::runtime_error("No LATTICE_CART or LATTICE_ABC block");
        const double scale = unitScale(abc);
        const double a = toNumber(abc.at(0).at(0)) * scale, b = toNumber(abc.at(0).at(1)) * scale, c = toNumber(abc.at(0).at(2)) * scale;
        const double rad = 3.14159265358979323846 / 180;
        const double al = toNumber(abc.at(1).at(0)) * rad, be = toNumber(abc.at(1).at(1)) * rad, ga = toNumber(abc.at(1).at(2)) * rad;
        const double cx = c * std::cos(be), cy = c * (std::cos(al) - std::cos(be) * std::cos(ga)) / std::sin(ga);
        setCell(s, {Vec3{a, 0, 0}, Vec3{b * std::cos(ga), b * std::sin(ga), 0}, Vec3{cx, cy, std::sqrt(std::max(0.0, c * c - cx * cx - cy * cy))}});
    }
    auto positions = block(lines, "positions_abs", found);
    if (found) {
        const double scale = unitScale(positions);
        for (const auto& row : positions) {
            const Vec3 x = vector3(row, 1);
            s.atoms.push_back(atom(row.at(0), {x[0] * scale, x[1] * scale, x[2] * scale}));
        }
    } else {
        positions = block(lines, "positions_frac", found);
        if (!found) throw std::runtime_error("No POSITIONS_ABS or POSITIONS_FRAC block");
        for (const auto& row : positions) s.atoms.push_back(atom(row.at(0), fromFractional(s, vector3(row, 1))));
    }
    // Ionic constraints fixing single Cartesian directions.
    const auto constraints = block(lines, "ionic_constraints", found);
    std::map<std::string, std::vector<std::size_t>> bySpecies;
    for (std::size_t i = 0; i < s.atoms.size(); ++i) bySpecies[s.atoms[i].symbol].push_back(i);
    for (const auto& row : constraints) {
        if (row.size() < 6) continue;
        const std::string symbol = elementFrom(row[1]);
        const long index = std::stol(row[2]);
        const auto list = bySpecies.find(symbol);
        if (list == bySpecies.end() || index < 1 || index > static_cast<long>(list->second.size())) continue;
        const Vec3 d = vector3(row, 3);
        for (int k = 0; k < 3; ++k)
            if (std::abs(d[k]) > 0.5 && std::abs(d[(k + 1) % 3]) < 1e-9 && std::abs(d[(k + 2) % 3]) < 1e-9)
                s.atoms[list->second[static_cast<std::size_t>(index - 1)]].fixed |= static_cast<std::uint8_t>(1 << k);
    }
}

// pw.x input with ibrav = 0: CELL_PARAMETERS and ATOMIC_POSITIONS cards.
void readPwi(const std::string& text, Structure& s)
{
    const auto lines = linesOf(text);
    double alat = 0;
    long nat = -1;
    for (const auto& line : lines) {
        std::string l = lower(line);
        l.erase(std::remove(l.begin(), l.end(), ' '), l.end());
        l = l.substr(0, l.find('!'));
        const auto value = [&](const std::string& key) -> const char* {
            const auto at = l.find(key + "=");
            return at == std::string::npos ? nullptr : l.c_str() + at + key.size() + 1;
        };
        if (const char* v = value("ibrav")) {
            if (std::strtol(v, nullptr, 10) != 0) throw std::runtime_error("Only ibrav = 0 (CELL_PARAMETERS) is read; this file uses another ibrav");
        }
        if (const char* v = value("celldm(1)")) alat = toNumber(std::string(v).substr(0, std::string(v).find_first_of(",/"))) * kBohr;
        if (const char* v = value("nat")) nat = std::strtol(v, nullptr, 10);
        if (l.rfind("a=", 0) == 0 || l.find(",a=") != std::string::npos) {
            const auto at = l.rfind("a=", 0) == 0 ? 0 : l.find(",a=") + 1;
            alat = toNumber(l.substr(at + 2, l.find_first_of(",/", at + 2) - at - 2));
        }
    }
    std::array<Vec3, 3> rows{};
    bool haveCell = false;
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const auto t = tokens(lines[i]);
        if (t.empty()) continue;
        const std::string card = lower(t[0]);
        std::string unit = t.size() > 1 ? lower(t[1]) : "";
        unit.erase(std::remove_if(unit.begin(), unit.end(), [](char c) { return c == '(' || c == ')' || c == '{' || c == '}'; }), unit.end());
        if (card == "cell_parameters") {
            const double scale = unit == "bohr" ? kBohr : unit == "angstrom" ? 1.0 : (alat > 0 ? alat : 1.0);
            for (int r = 0; r < 3; ++r) {
                rows[r] = vector3(tokens(lines.at(i + 1 + r)), 0);
                for (double& x : rows[r]) x *= scale;
            }
            haveCell = true;
        }
    }
    if (!haveCell) throw std::runtime_error("No CELL_PARAMETERS card");
    setCell(s, rows);
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const auto t = tokens(lines[i]);
        if (t.empty() || lower(t[0]) != "atomic_positions") continue;
        std::string unit = t.size() > 1 ? lower(t[1]) : "alat";
        unit.erase(std::remove_if(unit.begin(), unit.end(), [](char c) { return c == '(' || c == ')' || c == '{' || c == '}'; }), unit.end());
        for (std::size_t k = i + 1; k < lines.size() && (nat < 0 || static_cast<long>(s.atoms.size()) < nat); ++k) {
            const auto a = tokens(lines[k]);
            if (a.size() < 4 || !elementFrom(a[0]).size()) break;
            Vec3 x = vector3(a, 1);
            if (unit == "crystal") x = fromFractional(s, x);
            else {
                const double scale = unit == "bohr" ? kBohr : unit == "angstrom" ? 1.0 : (alat > 0 ? alat : 1.0);
                for (double& v : x) v *= scale;
            }
            AtomSite site = atom(a[0], x);
            if (a.size() >= 7)
                for (int c = 0; c < 3; ++c)
                    if (a[4 + c] == "0") site.fixed |= static_cast<std::uint8_t>(1 << c);
            s.atoms.push_back(site);
        }
        break;
    }
    if (s.atoms.empty()) throw std::runtime_error("No ATOMIC_POSITIONS card");
}
}  // namespace

bool isSimulationFormat(const std::string& format)
{
    const std::string f = canonical(format);
    return f == "pwi" || f == "cell" || f == "gen" || f == "xsf" || f == "cp2k";
}

bool canReadSimulationFormat(const std::string& format)
{
    const std::string f = canonical(format);
    return f == "pwi" || f == "cell" || f == "gen" || f == "xsf";
}

bool writeSimulationFormat(const Structure& structure, const std::string& format, std::ostream& out, std::string& error)
{
    const std::string f = canonical(format);
    if (!structure.hasUnitCell && (f == "pwi" || f == "cell" || f == "cp2k")) {
        error = "This format needs a unit cell; give the structure one with Edit > Structure Operations > Add Vacuum or Shrink-Wrap Cell.";
        return false;
    }
    if (f == "pwi") writePwi(structure, out);
    else if (f == "cell") writeCell(structure, out);
    else if (f == "gen") writeGen(structure, out);
    else if (f == "xsf") writeXsf(structure, out);
    else if (f == "cp2k") writeCp2k(structure, out);
    else {
        error = "Unknown simulation format " + format;
        return false;
    }
    if (!out) {
        error = "Could not write the file";
        return false;
    }
    return true;
}

bool readSimulationFormat(const std::string& text, const std::string& format, Structure& structure, std::string& error)
{
    structure = Structure();
    try {
        const std::string f = canonical(format);
        if (f == "xsf") readXsf(text, structure);
        else if (f == "gen") readGen(text, structure);
        else if (f == "cell") readCell(text, structure);
        else if (f == "pwi") readPwi(text, structure);
        else throw std::runtime_error("Unknown simulation format " + format);
        if (structure.atoms.empty()) throw std::runtime_error("No atoms found");
    } catch (const std::exception& e) {
        error = std::string("Could not read the ") + format + " file: " + e.what();
        structure = Structure();
        return false;
    }
    return true;
}
}  // namespace atomforge::io
