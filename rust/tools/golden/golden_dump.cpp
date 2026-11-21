// Reference outputs of the C++ implementation for the Rust port's tests.
// Built against the C++ libraries (see tools/golden/build.sh); writes JSON
// files into the directory given on the command line. Each section mirrors a
// Rust test module; values are printed with full precision.

#include "math/NeighbourSearch.h"
#include "model/DisplaySettings.h"
#include "model/Structure.h"
#include "util/ElementColourSchemes.h"
#include "util/ElementData.h"

#include <cstdio>
#include <algorithm>
#include <fstream>
#include <functional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace
{
std::string num(double v)
{
    char buffer[40];
    std::snprintf(buffer, sizeof buffer, "%.17g", v);
    return buffer;
}
std::string numf(float v)
{
    char buffer[40];
    std::snprintf(buffer, sizeof buffer, "%.9g", static_cast<double>(v));
    return buffer;
}
std::string quoted(const std::string& s)
{
    std::string out = "\"";
    for (char c : s) {
        if (c == '\n') { out += "\\n"; continue; }
        if (c == '\t') { out += "\\t"; continue; }
        if (c == '"' || c == '\\') out += '\\';
        out += c;
    }
    return out + "\"";
}

void write(const std::string& dir, const std::string& name, const std::string& text)
{
    std::ofstream out(dir + "/" + name);
    out << text;
    std::printf("wrote %s/%s\n", dir.c_str(), name.c_str());
}

std::string structureJson(const Structure& s)
{
    std::ostringstream o;
    o << "{\"has_cell\":" << (s.hasUnitCell ? "true" : "false") << ",\"cell\":[";
    for (int r = 0; r < 3; ++r)
        o << (r ? "," : "") << "[" << num(s.cellVectors[r][0]) << "," << num(s.cellVectors[r][1]) << "," << num(s.cellVectors[r][2]) << "]";
    o << "],\"pbc\":[" << s.pbc[0] << "," << s.pbc[1] << "," << s.pbc[2] << "],\"offset\":[" << num(s.cellOffset[0]) << "," << num(s.cellOffset[1])
      << "," << num(s.cellOffset[2]) << "],\"atoms\":[";
    for (std::size_t i = 0; i < s.atoms.size(); ++i) {
        const AtomSite& a = s.atoms[i];
        o << (i ? "," : "") << "{\"symbol\":" << quoted(a.symbol) << ",\"z\":" << a.atomicNumber << ",\"position\":[" << num(a.x) << "," << num(a.y)
          << "," << num(a.z) << "]}";
    }
    o << "]}";
    return o.str();
}

// ---------------------------------------------------------------- elements
void dumpElements(const std::string& dir)
{
    const auto radii = makeLiteratureCovalentRadii();
    std::ostringstream o;
    o << "{\"elements\":[";
    for (int z = 0; z <= 118; ++z) {
        float r, g, b;
        getDefaultElementColor(z, r, g, b);
        o << (z ? "," : "") << "{\"z\":" << z << ",\"symbol\":" << quoted(elementSymbol(z)) << ",\"name\":" << quoted(elementName(z))
          << ",\"covalent_radius\":" << numf(radii[z]) << ",\"mass\":" << num(elementAtomicMass(z)) << ",\"colour\":[" << numf(r) << "," << numf(g)
          << "," << numf(b) << "]}";
    }
    o << "],\"schemes\":[";
    const auto& schemes = atomforge::display::colourSchemes();
    for (std::size_t k = 0; k < schemes.size(); ++k) {
        o << (k ? "," : "") << "{\"name\":" << quoted(schemes[k].name) << ",\"description\":" << quoted(schemes[k].description) << ",\"colours\":[";
        const auto colours = atomforge::display::schemeColours(static_cast<int>(k));
        for (std::size_t z = 0; z < colours.size(); ++z)
            o << (z ? "," : "") << "[" << numf(colours[z].r) << "," << numf(colours[z].g) << "," << numf(colours[z].b) << "]";
        o << "]}";
    }
    o << "]}";
    write(dir, "elements.json", o.str());
}

// ---------------------------------------------------------------- fixtures
Structure fcc(double a, int n)
{
    Structure s;
    s.hasUnitCell = true;
    s.cellVectors = {{{a * n, 0, 0}, {0, a * n, 0}, {0, 0, a * n}}};
    const double basis[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                for (const auto& b : basis) {
                    AtomSite site;
                    site.symbol = "Cu";
                    site.atomicNumber = 29;
                    site.x = (i + b[0]) * a;
                    site.y = (j + b[1]) * a;
                    site.z = (k + b[2]) * a;
                    s.atoms.push_back(site);
                }
    return s;
}

Structure randomCluster(int n, unsigned seed, double box)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0, box);
    Structure s;
    s.hasUnitCell = false;
    s.pbc = {false, false, false};
    for (int i = 0; i < n; ++i) {
        AtomSite site;
        site.symbol = i % 3 ? "C" : "O";
        site.atomicNumber = i % 3 ? 6 : 8;
        site.x = u(rng);
        site.y = u(rng);
        site.z = u(rng);
        s.atoms.push_back(site);
    }
    return s;
}

Structure triclinic(unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(0, 1);
    Structure s;
    s.hasUnitCell = true;
    s.cellVectors = {{{4.1, 0.0, 0.0}, {1.3, 3.7, 0.0}, {0.9, -0.6, 3.2}}};
    for (int i = 0; i < 9; ++i) {
        const double f[3] = {u(rng) * 1.4 - 0.2, u(rng), u(rng)};  // some outside the cell
        AtomSite site;
        site.symbol = "Si";
        site.atomicNumber = 14;
        site.x = f[0] * s.cellVectors[0][0] + f[1] * s.cellVectors[1][0] + f[2] * s.cellVectors[2][0];
        site.y = f[0] * s.cellVectors[0][1] + f[1] * s.cellVectors[1][1] + f[2] * s.cellVectors[2][1];
        site.z = f[0] * s.cellVectors[0][2] + f[1] * s.cellVectors[1][2] + f[2] * s.cellVectors[2][2];
        s.atoms.push_back(site);
    }
    return s;
}

// ---------------------------------------------------------------- neighbours
void dumpNeighbours(const std::string& dir)
{
    struct Case { std::string name; Structure s; double cutoff; bool periodic; };
    std::vector<Case> cases = {
        {"fcc_cu_2x2x2", fcc(3.615, 2), 3.0, true},
        {"fcc_cu_small_cell_long_cutoff", fcc(3.615, 1), 6.0, true},
        {"random_cluster", randomCluster(80, 7, 12.0), 3.5, false},
        {"triclinic", triclinic(11), 5.0, true},
    };
    std::ostringstream o;
    o << "{\"cases\":[";
    for (std::size_t c = 0; c < cases.size(); ++c) {
        const auto pairs = neighbours::within(cases[c].s, cases[c].cutoff, cases[c].periodic);
        o << (c ? "," : "") << "{\"name\":" << quoted(cases[c].name) << ",\"cutoff\":" << num(cases[c].cutoff)
          << ",\"periodic\":" << (cases[c].periodic ? "true" : "false") << ",\"structure\":" << structureJson(cases[c].s) << ",\"pairs\":[";
        for (std::size_t k = 0; k < pairs.size(); ++k) {
            const auto& p = pairs[k];
            o << (k ? "," : "") << "[" << p.i << "," << p.j << "," << num(p.d.x) << "," << num(p.d.y) << "," << num(p.d.z) << "," << num(p.r) << "]";
        }
        o << "]}";
    }
    o << "]}";
    write(dir, "neighbours.json", o.str());
}
// ---------------------------------------------------------------- display settings
void dumpDisplay(const std::string& dir)
{
    using namespace atomforge::display;
    std::ostringstream o;
    o << "{\"default_text\":" << quoted(toText(DisplaySettings())) << ",\"presets\":[";
    for (std::size_t i = 0; i < presets().size(); ++i) {
        DisplaySettings s;
        applyPreset(s, presets()[i]);
        o << (i ? "," : "") << "{\"name\":" << quoted(presets()[i].name) << ",\"text\":" << quoted(toText(s)) << ",\"matching\":" << matchingPreset(s) << "}";
    }
    o << "],\"parsed\":[";
    const char* inputs[] = {
        "ambient=0.3\ncontrast=9\nlight_azimuth=-30\nlight_elevation=100\nbond_colour=0.1 0.2 0.3\nbackground=theme\ncell_colour=2 -1 0.5\n"
        "metal_bonds=yes\nbond_atom_limit=1234.6\nshininess_floor=abc\n  depth_cueing = 0.25  \nlight_follows_camera=1\n",
        "atom_scale=0.01\nbond_scale=10\noutline_width=0.3\noutline_colour=1 0 0\nselection_colour=0 1 0\nbond_tolerance=1.3\n",
        "light_azimuth=725.5\nlight_elevation=-12.5\nbond_colour=atoms\nbackground=0.25 0.5 0.75\n",
    };
    for (std::size_t i = 0; i < sizeof inputs / sizeof inputs[0]; ++i) {
        const DisplaySettings s = fromText(inputs[i]);
        const glm::vec3 light = lightDirection(s);
        o << (i ? "," : "") << "{\"input\":" << quoted(inputs[i]) << ",\"text\":" << quoted(toText(s)) << ",\"matching\":" << matchingPreset(s)
          << ",\"light\":[" << numf(light.x) << "," << numf(light.y) << "," << numf(light.z) << "]}";
    }
    o << "],\"look\":";
    Look look = defaultLook();
    look.colours[29] = glm::vec3(0.1f, 0.2f, 0.3f);
    look.radii[8] = 0.7f;
    look.display.atomScale = 1.5f;
    const std::string text = lookText(look);
    const Look back = readLook(text + "radius.Xx=2\ncolour.Fe=0.5 0.5\nradius.Au=9\n");
    o << "{\"text\":" << quoted(text) << ",\"reread_text\":" << quoted(lookText(back)) << "}}";
    write(dir, "display.json", o.str());
}

// ---------------------------------------------------------------- number formatting
void dumpFormatting(const std::string& dir)
{
    const double values[] = {0, 1, -1, 0.1, 0.4f, 35.26439f, 1e-5, 1.23456789e-7, 123456789.0, 1e21, 2.5e-300, 0.000123456, 99999.95, 3.14159265358979,
                             -0.0001, 6.02214076e23, 1.0 / 3.0, 1e15, 123456.5, 5e-324};
    std::ostringstream o;
    o << "{\"cases\":[";
    bool first = true;
    for (double v : values)
        for (int precision : {1, 3, 6, 9, 17}) {
            std::ostringstream s;
            s.precision(precision);
            s << v;
            o << (first ? "" : ",") << "{\"value\":" << num(v) << ",\"precision\":" << precision << ",\"text\":" << quoted(s.str()) << "}";
            first = false;
        }
    o << "]}";
    write(dir, "formatting.json", o.str());
}
// ---------------------------------------------------------------- std::sort order with ties
void dumpSort(const std::string& dir)
{
    std::ostringstream o;
    o << "{\"cases\":[";
    std::mt19937 rng(5);
    bool first = true;
    for (int n : {0, 1, 2, 3, 5, 16, 17, 18, 31, 64, 100, 257, 1000, 5000})
        for (int distinct : {1, 2, 3, 7, 1000000})
        {
            std::vector<std::pair<float, int>> v;
            std::uniform_int_distribution<int> key(0, distinct - 1);
            for (int i = 0; i < n; ++i) v.push_back({static_cast<float>(key(rng)), i});
            const auto input = v;
            std::sort(v.begin(), v.end(), [](const std::pair<float, int>& a, const std::pair<float, int>& b) { return a.first < b.first; });
            o << (first ? "" : ",") << "{\"keys\":[";
            for (std::size_t i = 0; i < input.size(); ++i) o << (i ? "," : "") << input[i].first;
            o << "],\"order\":[";
            for (std::size_t i = 0; i < v.size(); ++i) o << (i ? "," : "") << v[i].second;
            o << "]}";
            first = false;
        }
    // An adversarial (organ-pipe) input that drives introsort to its depth limit (heap sort).
    {
        std::vector<std::pair<float, int>> v;
        const int n = 4096;
        for (int i = 0; i < n; ++i) v.push_back({static_cast<float>((i % 2) ? n - i : i) , i});
        // Median-of-3 killer pattern.
        std::vector<int> killer(n);
        for (int i = 0; i < n; ++i) killer[i] = i;
        for (int i = 0; i < n; ++i) v[i].first = static_cast<float>((i * 7919) % 97);
        const auto input = v;
        std::sort(v.begin(), v.end(), [](const std::pair<float, int>& a, const std::pair<float, int>& b) { return a.first < b.first; });
        o << ",{\"keys\":[";
        for (std::size_t i = 0; i < input.size(); ++i) o << (i ? "," : "") << input[i].first;
        o << "],\"order\":[";
        for (std::size_t i = 0; i < v.size(); ++i) o << (i ? "," : "") << v[i].second;
        o << "]}";
    }
    o << "]}";
    write(dir, "sort.json", o.str());
}
// ---------------------------------------------------------------- <random> sequences
template <typename Engine>
void randomCases(std::ostringstream& o, const char* name, std::uint64_t seed)
{
    Engine e(static_cast<typename Engine::result_type>(seed));
    o << "{\"engine\":" << quoted(name) << ",\"seed\":" << seed << ",\"raw\":[";
    for (int i = 0; i < 20; ++i) o << (i ? "," : "") << e();
    o << "],\"real\":[";
    std::uniform_real_distribution<double> real(-2.5, 7.0);
    for (int i = 0; i < 50; ++i) o << (i ? "," : "") << num(real(e));
    o << "],\"int\":[";
    bool first = true;
    for (long long b : {1LL, 2LL, 5LL, 6LL, 100LL, 1000003LL, 4294967295LL, 4294967296LL, 9007199254740993LL})
        for (int i = 0; i < 6; ++i) {
            std::uniform_int_distribution<long long> d(-3, b);
            o << (first ? "" : ",") << d(e);
            first = false;
        }
    o << "],\"uint\":[";
    first = true;
    for (unsigned long long b : {1ULL, 7ULL, 18446744073709551615ULL})
        for (int i = 0; i < 4; ++i) {
            std::uniform_int_distribution<unsigned long long> d(0, b);
            o << (first ? "" : ",") << d(e);
            first = false;
        }
    o << "],\"normal\":[";
    std::normal_distribution<double> normal(0.3, 1.7);
    for (int i = 0; i < 41; ++i) o << (i ? "," : "") << num(normal(e));
    o << "],\"shuffles\":[";
    for (int n : {0, 1, 2, 3, 4, 5, 10, 11, 100, 1001, 70000}) {
        std::vector<int> v(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) v[static_cast<std::size_t>(i)] = i;
        std::shuffle(v.begin(), v.end(), e);
        o << (n ? "," : "") << "[";
        for (int i = 0; i < n; ++i) o << (i ? "," : "") << v[static_cast<std::size_t>(i)];
        o << "]";
    }
    o << "],\"after\":" << e() << "}";
}

void dumpRandom(const std::string& dir)
{
    std::ostringstream o;
    o << "{\"cases\":[";
    randomCases<std::mt19937>(o, "mt19937", 1);
    o << ",";
    randomCases<std::mt19937>(o, "mt19937", 12345);
    o << ",";
    randomCases<std::mt19937_64>(o, "mt19937_64", 1);
    o << ",";
    randomCases<std::mt19937_64>(o, "mt19937_64", 987654321987ULL);
    o << "]}";
    write(dir, "random.json", o.str());
}
// ---------------------------------------------------------------- nth_element and partial_sort
void dumpSelect(const std::string& dir)
{
    std::ostringstream o;
    o << "{\"cases\":[";
    std::mt19937 rng(11);
    bool first = true;
    const auto byKey = [](const std::pair<float, int>& a, const std::pair<float, int>& b) { return a.first < b.first; };
    for (int n : {1, 2, 3, 4, 5, 7, 16, 17, 33, 100, 257, 1000, 4096})
        for (int distinct : {1, 2, 5, 1000000})
            for (int which = 0; which < 4; ++which)
            {
                std::vector<std::pair<float, int>> v;
                std::uniform_int_distribution<int> key(0, distinct - 1);
                for (int i = 0; i < n; ++i) v.push_back({static_cast<float>(key(rng)), i});
                if (n == 4096 && distinct == 1000000)  // the median-of-3 killer keys of the sort test
                    for (int i = 0; i < n; ++i) v[i].first = static_cast<float>((i * 7919) % 97);
                const int nth = which == 0 ? 0 : which == 1 ? n / 3 : which == 2 ? n / 2 : n - 1;
                const auto input = v;
                auto selected = v;
                std::nth_element(selected.begin(), selected.begin() + nth, selected.end(), byKey);
                auto partial = v;
                std::partial_sort(partial.begin(), partial.begin() + nth, partial.end(), byKey);
                // A tolerance comparator (not a strict weak order), as sort-atoms uses.
                auto stable = v;
                std::stable_sort(stable.begin(), stable.end(), [](const std::pair<float, int>& a, const std::pair<float, int>& b) {
                    if (std::abs(a.first - b.first) > 1.5f) return a.first < b.first;
                    return (a.second * 31) % 7 < (b.second * 31) % 7;
                });
                o << (first ? "" : ",") << "{\"nth\":" << nth << ",\"keys\":[";
                for (std::size_t i = 0; i < input.size(); ++i) o << (i ? "," : "") << input[i].first;
                o << "],\"nth_element\":[";
                for (std::size_t i = 0; i < selected.size(); ++i) o << (i ? "," : "") << selected[i].second;
                o << "],\"partial_sort\":[";
                for (std::size_t i = 0; i < partial.size(); ++i) o << (i ? "," : "") << partial[i].second;
                o << "],\"stable_sort\":[";
                for (std::size_t i = 0; i < stable.size(); ++i) o << (i ? "," : "") << stable[i].second;
                o << "]}";
                first = false;
            }
    o << "]}";
    write(dir, "select.json", o.str());
}
}  // namespace

int main(int argc, char** argv)
{
    const std::string dir = argc > 1 ? argv[1] : ".";
    dumpElements(dir);
    dumpNeighbours(dir);
    dumpDisplay(dir);
    dumpFormatting(dir);
    dumpSort(dir);
    dumpRandom(dir);
    dumpSelect(dir);
    return 0;
}
