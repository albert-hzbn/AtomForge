#include "electronic/Volume.h"
#include "util/ElementData.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace atomforge::electronic
{
namespace
{
int atomicNumberFromSymbol(const std::string& symbol)
{
    for (int z = 1; z <= 118; ++z) if (symbol == elementSymbol(z)) return z;
    return 0;
}
std::string line(std::istream& in)
{
    std::string s;
    if (!std::getline(in, s)) throw std::runtime_error("Unexpected end of volumetric file");
    return s;
}
std::vector<std::string> words(const std::string& s)
{
    std::istringstream in(s);
    std::vector<std::string> out;
    for (std::string token; in >> token;) out.push_back(token);
    return out;
}
double number(std::string s)
{
    std::replace(s.begin(), s.end(), 'D', 'E');
    std::replace(s.begin(), s.end(), 'd', 'e');
    std::size_t end = 0;
    const double value = std::stod(s, &end);
    if (end != s.size() || !std::isfinite(value)) throw std::runtime_error("Invalid finite number: " + s);
    return value;
}
int integer(const std::string& s)
{
    std::size_t end = 0;
    const int value = std::stoi(s, &end);
    if (end != s.size()) throw std::runtime_error("Invalid integer: " + s);
    return value;
}
glm::dvec3 vector(const std::vector<std::string>& w, std::size_t offset = 0)
{
    if (w.size() < offset + 3) throw std::runtime_error("Missing vector components");
    return {number(w[offset]), number(w[offset + 1]), number(w[offset + 2])};
}
void allocate(Grid& g)
{
    std::size_t n = 1;
    for (int d : g.shape)
    {
        if (d < 2 || d > 100000 || n > 100000000 / static_cast<std::size_t>(d))
            throw std::runtime_error("Invalid or excessive grid dimensions");
        n *= d;
    }
    g.values.resize(n);
}
std::string token(std::istream& in)
{
    std::string s;
    if (!(in >> s)) throw std::runtime_error("Truncated volumetric data");
    return s;
}
void setQuantity(Grid& g, const std::string& quantity, double densityScale, double potentialScale)
{
    double factor = 1;
    if (quantity == "density") { g.unit = "e/A^3"; factor = densityScale; }
    else if (quantity == "potential") { g.unit = "eV"; factor = potentialScale; }
    else if (quantity == "elf") g.unit = "1";
    else if (quantity != "raw" && quantity != "auto") throw std::invalid_argument("Quantity must be auto, raw, density, potential or elf");
    for (double& value : g.values) value *= factor;
    g.validate();
}
Volume readCube(std::istream& in, const std::string& quantity, const std::string& coordinateUnit)
{
    Volume out;
    Grid g;
    g.name = line(in);
    line(in);
    const auto header = words(line(in));
    if (header.size() < 4) throw std::runtime_error("Invalid Cube header");
    const int atoms = integer(header[0]);
    if (atoms < -10000000 || atoms > 10000000) throw std::runtime_error("Excessive Cube atom count");
    int channels = header.size() > 4 ? integer(header[4]) : 1;
    if (atoms < 0 && channels != 1) throw std::runtime_error("Orbital Cube NVAL must be 1");
    if (coordinateUnit != "bohr" && coordinateUnit != "angstrom") throw std::invalid_argument("Cube coordinate unit must be bohr or angstrom");
    const double length = coordinateUnit == "bohr" ? bohr : 1.0;
    g.origin = vector(header, 1) * length;
    for (int i = 0; i < 3; ++i)
    {
        const auto w = words(line(in));
        if (w.size() < 4) throw std::runtime_error("Invalid Cube axis");
        // Negative dimension conventions vary among writers; units are explicit.
        const int n = integer(w[0]);
        if (n == 0 || n < -100000 || n > 100000) throw std::runtime_error("Invalid Cube dimension");
        g.shape[i] = std::abs(n);
        g.cell[i] = vector(w, 1) * (length * (g.shape[i] - 1));
    }
    allocate(g);
    for (int i = 0; i < std::abs(atoms); ++i)
    {
        const auto w = words(line(in));
        if (w.size() < 5) throw std::runtime_error("Invalid Cube atom");
        out.sites.push_back({integer(w[0]), vector(w, 2) * length});
    }
    std::vector<std::string> names;
    if (atoms < 0)
    {
        channels = integer(token(in));
        if (channels < 1 || channels > 1024) throw std::runtime_error("Invalid Cube channel count");
        for (int i = 0; i < channels; ++i) names.push_back("orbital " + std::to_string(integer(token(in))));
    }
    if (channels < 1 || channels > 1024 || g.values.size() > 100000000 / static_cast<std::size_t>(channels))
        throw std::runtime_error("Excessive Cube channels or samples");
    out.fields.assign(channels, g);
    for (int c = 0; c < channels; ++c) out.fields[c].name = names.empty() ? "field " + std::to_string(c) : names[c];
    for (int x = 0; x < g.shape[0]; ++x)
        for (int y = 0; y < g.shape[1]; ++y)
            for (int z = 0; z < g.shape[2]; ++z)
                for (auto& field : out.fields) field.values[g.index(x, y, z)] = number(token(in));
    std::string extra;
    if (in >> extra) throw std::runtime_error("Extra Cube values after expected data");
    for (auto& field : out.fields) setQuantity(field, quantity, 1.0 / (bohr * bohr * bohr), hartree);
    return out;
}

Volume readVasp(std::istream& in, const std::string& quantity)
{
    Volume out;
    Grid g;
    g.periodic = true;
    g.name = line(in);
    const auto scales = words(line(in));
    if (scales.size() != 1 && scales.size() != 3) throw std::runtime_error("VASP requires one or three scale factors");
    for (int i = 0; i < 3; ++i) g.cell[i] = vector(words(line(in)));
    glm::dvec3 scaling(1.0);
    if (scales.size() == 1)
    {
        double s = number(scales[0]);
        if (s == 0 || std::abs(glm::determinant(g.cell)) < 1e-12) throw std::runtime_error("Invalid VASP lattice or scale");
        if (s < 0) s = std::cbrt(-s / std::abs(glm::determinant(g.cell)));
        scaling = glm::dvec3(s);
    }
    else
    {
        scaling = vector(scales);
        if (scaling.x <= 0 || scaling.y <= 0 || scaling.z <= 0) throw std::runtime_error("VASP component scales must be positive");
    }
    for (int i = 0; i < 3; ++i) g.cell[i] *= scaling;
    auto symbols = words(line(in));
    if (symbols.empty()) throw std::runtime_error("Missing VASP species");
    std::vector<std::string> counts;
    if (std::isdigit(static_cast<unsigned char>(symbols[0][0])))
    {
        counts = symbols;
        symbols.assign(counts.size(), "X");
    }
    else counts = words(line(in));
    if (counts.size() != symbols.size()) throw std::runtime_error("VASP species and counts differ");
    auto mode = words(line(in));
    if (!mode.empty() && (mode[0][0] == 'S' || mode[0][0] == 's')) mode = words(line(in));
    if (mode.empty()) throw std::runtime_error("Missing VASP coordinate mode");
    const char m = static_cast<char>(std::tolower(static_cast<unsigned char>(mode[0][0])));
    if (m != 'c' && m != 'k' && m != 'd') throw std::runtime_error("Unknown VASP coordinate mode");
    for (std::size_t species = 0; species < counts.size(); ++species)
    {
        const int count = integer(counts[species]);
        if (count < 0 || count > 10000000 || out.sites.size() + count > 10000000) throw std::runtime_error("Invalid VASP atom count");
        for (int j = 0; j < count; ++j)
        {
            auto p = vector(words(line(in)));
            p = m == 'd' ? g.cell * p : scaling * p;
            out.sites.push_back({atomicNumberFromSymbol(symbols[species]), p});
        }
    }
    std::string s;
    while (std::getline(in, s))
    {
        auto w = words(s);
        if (w.empty()) continue;
        if (w[0] == "augmentation")
        {
            if (w.size() < 4 || w[1] != "occupancies") throw std::runtime_error("Invalid VASP augmentation header");
            const int count = integer(w[3]);
            if (count < 0 || count > 100000000) throw std::runtime_error("Invalid augmentation count");
            for (int i = 0; i < count; ++i) number(token(in));
            if (count) std::getline(in, s);
            continue;
        }
        if (w.size() != 3) throw std::runtime_error("Unexpected VASP content after grid/atoms; trajectories are not supported");
        for (int i = 0; i < 3; ++i) g.shape[i] = integer(w[i]);
        allocate(g);
        if (!out.fields.empty() && g.shape != out.fields[0].shape) throw std::runtime_error("VASP channel dimensions differ");
        for (double& value : g.values) value = number(token(in));
        std::getline(in, s);
        if (!words(s).empty()) throw std::runtime_error("Extra values on VASP final grid line");
        setQuantity(g, quantity, 1.0 / std::abs(glm::determinant(g.cell)), 1.0);
        out.fields.push_back(g);
        if (out.fields.size() > 4) throw std::runtime_error("VASP supports at most four scalar channels");
    }
    const auto n = out.fields.size();
    if (n != 1 && n != 2 && n != 4) throw std::runtime_error("Expected one, two or four VASP grids");
    const std::vector<std::string> names = n == 4 ? std::vector<std::string>{"total", "mx", "my", "mz"} : std::vector<std::string>{"total", "magnetization"};
    for (std::size_t i = 0; i < n; ++i) out.fields[i].name = names[i];
    return out;
}

Volume readXsf(std::istream& in, const std::string& quantity)
{
    Volume out;
    std::string s;
    bool atomsSection = false;
    std::size_t totalSamples = 0;
    while (std::getline(in, s))
    {
        const auto w = words(s);
        if (w.empty() || w[0][0] == '#') continue;
        if (atomsSection)
        {
            const int z = std::isdigit(static_cast<unsigned char>(w[0][0])) ? integer(w[0]) : atomicNumberFromSymbol(w[0]);
            if (w.size() >= 4 && (z > 0 || w[0] == "0"))
            {
                if (out.sites.size() >= 10000000) throw std::runtime_error("Excessive XSF atom count");
                out.sites.push_back({z,vector(w,1)});
                continue;
            }
            atomsSection = false;
        }
        if (w[0] == "ATOMS") { atomsSection = true; continue; }
        if (w[0] == "ANIMSTEPS" || w[0].find("BANDGRID") != std::string::npos)
            throw std::runtime_error("XSF animations and reciprocal-space BANDGRID are not real-space scalar grids");
        if (w[0] == "PRIMCOORD")
        {
            const auto counts = words(line(in));
            if (counts.empty()) throw std::runtime_error("Missing XSF atom count");
            const int count = integer(counts[0]);
            if (count < 0 || count > 10000000) throw std::runtime_error("Invalid XSF atom count");
            for (int j = 0; j < count; ++j)
            {
                const auto atom = words(line(in));
                if (atom.size() < 4) throw std::runtime_error("Invalid XSF atom");
                int z = std::isdigit(static_cast<unsigned char>(atom[0][0])) ? integer(atom[0]) : atomicNumberFromSymbol(atom[0]);
                out.sites.push_back({z, vector(atom, 1)});
            }
        }
        else if (w[0].rfind("BEGIN_DATAGRID_3D", 0) == 0)
        {
            Grid g;
            g.name = w[0].size() > 17 ? w[0].substr(17) : "field";
            auto dimensions = words(line(in));
            if (dimensions.size() != 3) throw std::runtime_error("Invalid XSF dimensions");
            for (int i = 0; i < 3; ++i) g.shape[i] = integer(dimensions[i]);
            g.origin = vector(words(line(in)));
            for (int i = 0; i < 3; ++i) g.cell[i] = vector(words(line(in)));
            allocate(g);
            totalSamples += g.values.size();
            if (totalSamples > 100000000) throw std::runtime_error("Excessive total XSF samples");
            for (double& value : g.values) value = number(token(in));
            if (token(in) != "END_DATAGRID_3D") throw std::runtime_error("Missing XSF grid terminator or extra values");
            setQuantity(g, quantity, 1.0, 1.0);
            out.fields.push_back(std::move(g));
        }
    }
    if (out.fields.empty()) throw std::runtime_error("No XSF DATAGRID_3D found");
    return out;
}
void writeVector(std::ostream& out, glm::dvec3 p) { out << p.x << ' ' << p.y << ' ' << p.z << '\n'; }
void writeValues(std::ostream& out, const std::vector<double>& values, double factor)
{
    std::size_t count = 0;
    for (double v : values) out << v * factor << (++count % 5 ? ' ' : '\n');
    out << '\n';
}
}

Volume loadVolume(const std::string& path, const std::string& quantity, const std::string& coordinateUnit)
{
    std::ifstream in{std::filesystem::u8path(path)};
    if (!in) throw std::runtime_error("Cannot open volumetric file: " + path);
    std::string ext = std::filesystem::u8path(path).extension().u8string();
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (ext == ".cube" || ext == ".cub") return readCube(in, quantity, coordinateUnit);
    if (ext == ".xsf") return readXsf(in, quantity);
    std::string q = quantity;
    if (q == "auto")
    {
        const auto name = std::filesystem::u8path(path).filename().u8string();
        if (name.rfind("CHGCAR", 0) == 0 || name == "CHG") q = "density";
        else if (name.rfind("LOCPOT", 0) == 0) q = "potential";
        else if (name.rfind("ELFCAR", 0) == 0) q = "elf";
        else q = "raw";
    }
    return readVasp(in, q);
}

void saveVolume(const Volume& volume, const std::string& path, const std::string& format)
{
    if (volume.fields.empty()) throw std::invalid_argument("No fields to save");
    for (const auto& field : volume.fields) requireAligned(volume.fields[0], field);
    if (format != "cube" && format != "xsf" && format != "vasp") throw std::invalid_argument("Output format must be cube, xsf or vasp");
    const Grid& first = volume.fields[0];
    if (format == "vasp" && (!first.periodic || glm::length(first.origin) > 1e-10))
        throw std::invalid_argument("VASP output requires a periodic grid with zero origin");
    if (format == "vasp" && volume.fields.size() != 1 && volume.fields.size() != 2 && volume.fields.size() != 4)
        throw std::invalid_argument("VASP requires one, two or four fields");
    if (format == "vasp" && volume.sites.empty()) throw std::invalid_argument("VASP output requires atomic sites");
    std::ofstream out{std::filesystem::u8path(path)};
    if (!out) throw std::runtime_error("Cannot write volume: " + path);
    out << std::setprecision(17);
    if (format == "cube")
    {
        // Add periodic endpoint planes so finite-grid interchange preserves extent.
        auto shape = first.shape;
        if (first.periodic) for (int& n : shape) ++n;
        out << "AtomForge volume\nUnits: atomic units for density/potential; raw otherwise\n" << volume.sites.size() << ' ';
        const auto p = first.origin / bohr;
        out << p.x << ' ' << p.y << ' ' << p.z << ' ' << volume.fields.size() << '\n';
        for (int i = 0; i < 3; ++i) { out << shape[i] << ' '; writeVector(out, first.steps()[i] / bohr); }
        for (auto site : volume.sites) { out << site.number << " 0 "; writeVector(out, site.position / bohr); }
        for (int x = 0; x < shape[0]; ++x)
            for (int y = 0; y < shape[1]; ++y)
                for (int z = 0; z < shape[2]; ++z)
                {
                    for (const auto& f : volume.fields)
                    {
                        double factor = f.unit == "e/A^3" ? bohr * bohr * bohr : (f.unit == "eV" ? 1.0 / hartree : 1.0);
                        out << f.values[f.index(x % f.shape[0], y % f.shape[1], z % f.shape[2])] * factor << ' ';
                    }
                    out << '\n';
                }
    }
    else if (format == "xsf")
    {
        out << "CRYSTAL\nPRIMVEC\n";
        for (int i = 0; i < 3; ++i) writeVector(out, first.cell[i]);
        out << "PRIMCOORD\n" << volume.sites.size() << " 1\n";
        for (auto site : volume.sites) { out << site.number << ' '; writeVector(out, site.position); }
        out << "BEGIN_BLOCK_DATAGRID_3D\nAtomForge\n";
        for (std::size_t c = 0; c < volume.fields.size(); ++c)
        {
            const auto& f = volume.fields[c];
            auto shape = f.shape;
            if (f.periodic) for (int& n : shape) ++n;
            out << "BEGIN_DATAGRID_3D_field" << c << '\n' << shape[0] << ' ' << shape[1] << ' ' << shape[2] << '\n';
            writeVector(out, f.origin);
            for (int i = 0; i < 3; ++i) writeVector(out, f.cell[i]);
            for (int z = 0; z < shape[2]; ++z)
                for (int y = 0; y < shape[1]; ++y)
                    for (int x = 0; x < shape[0]; ++x)
                        out << f.values[f.index(x % f.shape[0], y % f.shape[1], z % f.shape[2])] << '\n';
            out << "END_DATAGRID_3D\n";
        }
        out << "END_BLOCK_DATAGRID_3D\n";
    }
    else
    {
        out << "AtomForge post-processing grid (no PAW augmentation)\n1.0\n";
        for (int i = 0; i < 3; ++i) writeVector(out, first.cell[i]);
        // Each site is its own species group, preserving site order without loss.
        for (auto site : volume.sites) out << elementSymbol(site.number) << ' ';
        out << '\n';
        for (std::size_t i = 0; i < volume.sites.size(); ++i) out << "1 ";
        out << "\nCartesian\n";
        for (auto site : volume.sites) writeVector(out, site.position);
        for (const auto& f : volume.fields)
        {
            out << '\n' << f.shape[0] << ' ' << f.shape[1] << ' ' << f.shape[2] << '\n';
            writeValues(out, f.values, f.unit == "e/A^3" ? std::abs(glm::determinant(f.cell)) : 1.0);
        }
    }
    if (!out) throw std::runtime_error("Failed while writing volume");
}
}
