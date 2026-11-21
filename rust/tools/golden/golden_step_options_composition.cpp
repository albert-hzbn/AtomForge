// Golden step options of Build > Substitutional Solid Solution and Build >
// Custom Structure: the dialogs' editStep() reading of a step's options and
// their stepOptions() writing them back (copied from
// ui/SubstitutionalSolidSolutionDialog.cpp and ui/CustomStructureDialog.cpp,
// whose GL/ImGui state cannot be built here), through pipeline/Options.h.
// Writes rust/tests/golden/step_options_composition.json.
#include "pipeline/Options.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace
{
const char* kSymbols[] = {"", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
                          "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"};
int atomicNumber(const std::string& s)
{
    for (int z = 1; z <= 30; ++z)
        if (s == kSymbols[z]) return z;
    return 0;
}

std::string json(const std::string& s)
{
    std::string out = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') out += '\\';
        out += c;
    }
    return out + "\"";
}

struct Entry { int atomicNumber = 26; float percent = 0; int count = 0; };

void recompute(std::vector<Entry>& entries, int total)
{
    const int n = (int)entries.size();
    if (n == 0 || total == 0) return;
    float sum = 0.0f;
    for (const auto& e : entries) sum += std::max(0.0f, e.percent);
    if (sum <= 0.0f) {
        const int base = total / n, rem = total % n;
        for (int i = 0; i < n; ++i) entries[i].count = base + (i < rem ? 1 : 0);
        return;
    }
    std::vector<float> real(n);
    for (int i = 0; i < n; ++i) real[i] = std::max(0.0f, entries[i].percent) / sum * (float)total;
    int assigned = 0;
    std::vector<int> floored(n);
    std::vector<float> remainders(n);
    for (int i = 0; i < n; ++i) {
        floored[i] = (int)std::floor(real[i]);
        remainders[i] = real[i] - (float)floored[i];
        assigned += floored[i];
    }
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return remainders[a] > remainders[b]; });
    const int leftover = total - assigned;
    for (int k = 0; k < leftover && k < n; ++k) floored[order[k]]++;
    for (int i = 0; i < n; ++i) entries[i].count = floored[i];
}

// editStep + stepOptions of the solid solution dialog; the step input is
// `total` atoms of copper (0: none).
std::string sss(const std::string& options, int total)
{
    using namespace atomforge::pipeline::options;
    Reader reader(options);
    std::vector<Entry> entries;
    std::string sourcePath;
    if (reader.has("--input")) sourcePath = reader.text("--input");  // loads (the cases' files exist)
    if (total > 0) entries.push_back({29, 100.0f, total});
    const std::string frac = reader.text("--frac");
    std::vector<Entry> parsed;
    std::istringstream stream(frac);
    for (std::string token; std::getline(stream, token, ',');) {
        const auto eq = token.find('=');
        if (eq == std::string::npos) continue;
        std::string symbol = token.substr(0, eq);
        symbol.erase(std::remove_if(symbol.begin(), symbol.end(), [](unsigned char c) { return std::isspace(c); }), symbol.end());
        const int z = atomicNumber(symbol);
        if (z <= 0) continue;
        Entry entry;
        entry.atomicNumber = z;
        try { entry.percent = std::stof(token.substr(eq + 1)) * 100.0f; } catch (const std::exception&) { continue; }
        parsed.push_back(entry);
    }
    if (!parsed.empty()) {
        entries = parsed;
        recompute(entries, total);
    }
    const int seed = std::max(0, reader.integer("--seed", 12345));
    // stepOptions
    float percentSum = 0.0f;
    for (const auto& e : entries) percentSum += std::max(0.0f, e.percent);
    std::string out;
    for (const auto& e : entries) {
        const double fraction = total > 0 ? (double)e.count / (double)total
                              : percentSum > 0.0f ? (double)std::max(0.0f, e.percent) / percentSum : 0.0;
        out += (out.empty() ? "" : ",") + std::string(kSymbols[e.atomicNumber]) + "=" + Writer::number(fraction);
    }
    Writer w;
    if (!sourcePath.empty()) w.add("--input", sourcePath);
    if (!out.empty()) w.add("--frac", out);
    w.add("--seed", seed);
    return w.str();
}

// editStep + stepOptions of the custom structure dialog (the reference is the
// step's input; --mesh files load).
std::string custom(const std::string& options)
{
    using namespace atomforge::pipeline::options;
    Reader reader(options);
    std::string modelPath = reader.has("--mesh") ? reader.text("--mesh") : "";
    float modelScale = (float)reader.number("--scale", 1.0);
    float vacuumPadding = (float)reader.number("--vacuum", 5.0);
    bool setOutputCell = !reader.has("--no-cell");
    const auto center = reader.numbers("--center");
    bool autoCenter = center.size() < 3;
    float cx = 0, cy = 0, cz = 0;
    if (center.size() >= 3) { cx = (float)center[0]; cy = (float)center[1]; cz = (float)center[2]; }
    const int repA = reader.integer("--repa", 0), repB = reader.integer("--repb", 0), repC = reader.integer("--repc", 0);
    bool autoReplicate = !(repA > 0 || repB > 0 || repC > 0);
    int ra = repA > 0 ? repA : 5, rb = repB > 0 ? repB : 5, rc = repC > 0 ? repC : 5;
    float rot[3] = {(float)reader.number("--rotx", 0.0), (float)reader.number("--roty", 0.0), (float)reader.number("--rotz", 0.0)};
    int type = 1, miller[3] = {0, 0, 1};
    bool orient;
    if (reader.has("--miller")) {
        const auto m = reader.numbers("--miller");
        orient = true;
        type = 2;
        miller[0] = m.size() > 0 ? (int)m[0] : 1;
        miller[1] = m.size() > 1 ? (int)m[1] : 0;
        miller[2] = m.size() > 2 ? (int)m[2] : 0;
    } else
        orient = rot[0] != 0.0f || rot[1] != 0.0f || rot[2] != 0.0f;
    // stepOptions
    Writer out;
    if (!modelPath.empty()) out.add("--mesh", modelPath);
    out.add("--scale", modelScale);
    if (setOutputCell) out.add("--vacuum", vacuumPadding);
    else out.flag("--no-cell");
    if (!autoCenter) out.add("--center", std::vector<double>{cx, cy, cz});
    if (!autoReplicate) out.add("--repa", std::max(1, ra)).add("--repb", std::max(1, rb)).add("--repc", std::max(1, rc));
    if (orient) {
        if (type == 2) out.add("--miller", std::to_string(miller[0]) + " " + std::to_string(miller[1]) + " " + std::to_string(miller[2]));
        else out.add("--rotx", rot[0]).add("--roty", rot[1]).add("--rotz", rot[2]);
    }
    return out.str();
}
}  // namespace

int main(int argc, char** argv)
{
    const std::string dir = argc > 1 ? argv[1] : "rust/tests/golden";
    std::ofstream f(dir + "/step_options_composition.json");
    f << "{\"sss\": [\n";
    const std::pair<std::string, int> sssCases[] = {
        {"", 0}, {"", 500}, {"--frac Cu=0.7,Ni=0.3 --seed 42", 500}, {"--frac Cu=0.7,Ni=0.3 --seed 42", 0},
        {"--frac \"Cu = 0.333, Ni=0.333,Fe=0.334\" --seed -5", 7}, {"--frac Cu=0.1x,Zn=abc,Xx=0.5,Ni=0.25", 13},
        {"--frac Cu=1e-1,Ni=2.5e-1", 3}, {"--seed 7", 4},
    };
    bool first = true;
    for (const auto& [options, total] : sssCases) {
        f << (first ? "" : ",\n") << "  {\"options\": " << json(options) << ", \"total\": " << total << ", \"out\": " << json(sss(options, total)) << "}";
        first = false;
    }
    f << "\n], \"custom\": [\n";
    const std::string customCases[] = {
        "", "--scale 12 --vacuum 5", "--scale 0.1 --no-cell", "--center 1.5,-2 3.25 --repa 3", "--repb 2 --repc 4",
        "--rotx 10.5 --rotz -30", "--miller \"1 1 0\" --rotx 45", "--miller 2", "--scale 1.2345678901 --vacuum 0.3",
        "--mesh ../../../docs/manual/examples/octahedron.obj --scale 12",
    };
    first = true;
    for (const auto& options : customCases) {
        f << (first ? "" : ",\n") << "  {\"options\": " << json(options) << ", \"out\": " << json(custom(options)) << "}";
        first = false;
    }
    f << "\n]}\n";
    return 0;
}
