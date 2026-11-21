// Reference outputs of the C++ structure loader and writers (io/StructureLoader,
// io/SimulationFormats, with Open Babel and spglib as in the desktop build) for
// the af-io tests. Usage: golden_io <fixture dir> <extra fixture dir> <output dir>.

#include "io/StructureLoader.h"
#include "io/Trajectory.h"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

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
    for (unsigned char c : s) {
        if (c == '\n') { out += "\\n"; continue; }
        if (c == '\r') { out += "\\r"; continue; }
        if (c == '\t') { out += "\\t"; continue; }
        if (c == '"' || c == '\\') out += '\\';
        if (c < 0x20) { char b[8]; std::snprintf(b, sizeof b, "\\u%04x", c); out += b; continue; }
        out += static_cast<char>(c);
    }
    return out + "\"";
}

std::string structureJson(const Structure& s)
{
    std::ostringstream o;
    o << "{\"has_cell\":" << (s.hasUnitCell ? "true" : "false") << ",\"cell\":[";
    for (int r = 0; r < 3; ++r)
        o << (r ? "," : "") << "[" << num(s.cellVectors[r][0]) << "," << num(s.cellVectors[r][1]) << "," << num(s.cellVectors[r][2]) << "]";
    o << "],\"pbc\":[" << s.pbc[0] << "," << s.pbc[1] << "," << s.pbc[2] << "],\"offset\":[" << num(s.cellOffset[0]) << "," << num(s.cellOffset[1])
      << "," << num(s.cellOffset[2]) << "],\"load_note\":" << quoted(s.loadNote) << ",\"ipf_load_status\":" << quoted(s.ipfLoadStatus) << ",\"atoms\":[";
    for (std::size_t i = 0; i < s.atoms.size(); ++i) {
        const AtomSite& a = s.atoms[i];
        o << (i ? "," : "") << "{\"symbol\":" << quoted(a.symbol) << ",\"z\":" << a.atomicNumber << ",\"position\":[" << num(a.x) << "," << num(a.y)
          << "," << num(a.z) << "],\"colour\":[" << numf(a.r) << "," << numf(a.g) << "," << numf(a.b) << "],\"fixed\":" << int(a.fixed)
          << ",\"magmom\":" << num(a.magmom) << ",\"charge\":" << num(a.charge) << ",\"velocity\":[" << num(a.velocity[0]) << "," << num(a.velocity[1])
          << "," << num(a.velocity[2]) << "]}";
    }
    o << "],\"grain_colors\":[";
    for (std::size_t i = 0; i < s.grainColors.size(); ++i)
        o << (i ? "," : "") << "[" << numf(s.grainColors[i][0]) << "," << numf(s.grainColors[i][1]) << "," << numf(s.grainColors[i][2]) << "]";
    o << "],\"grain_region_ids\":[";
    for (std::size_t i = 0; i < s.grainRegionIds.size(); ++i) o << (i ? "," : "") << s.grainRegionIds[i];
    o << "]}";
    return o.str();
}

std::string readText(const fs::path& p)
{
    std::ifstream in(p, std::ios::binary);
    std::ostringstream o;
    o << in.rdbuf();
    return o.str();
}
}  // namespace

int main(int argc, char** argv)
{
    if (argc < 4) {
        std::fprintf(stderr, "usage: golden_io <fixtures> <extra fixtures> <output>\n");
        return 2;
    }
    const fs::path fixtures = argv[1], extra = argv[2], output = argv[3];
    std::vector<fs::path> files;
    for (const auto& dir : {fixtures, extra})
        for (const auto& entry : fs::directory_iterator(dir))
            if (entry.is_regular_file() && entry.path().extension() != ".atomforge-ipf") files.push_back(entry.path());
    std::sort(files.begin(), files.end());

    std::ostringstream o;
    o << "{\"loads\":[";
    for (std::size_t k = 0; k < files.size(); ++k) {
        Structure s;
        std::string error;
        const bool ok = loadStructureFromFile(files[k].string(), s, error);
        const std::string relative = files[k].parent_path() == fixtures ? "fixtures/" + files[k].filename().string() : "extra/" + files[k].filename().string();
        o << (k ? "," : "") << "\n{\"file\":" << quoted(relative) << ",\"ok\":" << (ok ? "true" : "false") << ",\"error\":" << quoted(error)
          << ",\"structure\":" << (ok ? structureJson(s) : std::string("null")) << "}";
    }
    o << "],\n\"saves\":[";
    // Writers: each source structure in each format; the text of native writers.
    const char* sources[] = {"cell.extxyz", "nacl.vasp", "selective.vasp", "abs.cell", "charge_triclinic.data", "labels.xyz", "p1.cif"};
    const char* formats[] = {"cif", "xyz", "extxyz", "vasp", "data", "pwi", "cell", "gen", "xsf", "cp2k", "pdb", "sdf", "mol2", "gjf", "mol"};
    const fs::path temp = output / "save_tmp";
    fs::create_directories(temp);
    bool first = true;
    for (const char* source : sources) {
        Structure s;
        std::string error;
        if (!loadStructureFromFile((fixtures / source).string(), s, error)) {
            std::fprintf(stderr, "cannot load %s: %s\n", source, error.c_str());
            return 1;
        }
        for (const char* format : formats) {
            const fs::path file = temp / (std::string("out.") + format);
            fs::remove(file);
            const bool ok = saveStructure(s, file.string(), format);
            o << (first ? "" : ",") << "\n{\"source\":" << quoted(source) << ",\"format\":" << quoted(format) << ",\"ok\":" << (ok ? "true" : "false")
              << ",\"error\":" << quoted(lastSaveError()) << ",\"text\":" << quoted(ok ? readText(file) : std::string()) << "}";
            first = false;
        }
    }
    o << "],\n\"sidecar\":";
    {
        // A structure with orientation colours and regions survives a save and load.
        Structure s;
        std::string error;
        loadStructureFromFile((fixtures / "nacl.vasp").string(), s, error);
        s.grainColors.clear();
        s.grainRegionIds.clear();
        for (std::size_t i = 0; i < s.atoms.size(); ++i) {
            s.grainColors.push_back({{0.1f * static_cast<float>(i), 0.5f, 1.0f - 0.1f * static_cast<float>(i)}});
            s.grainRegionIds.push_back(static_cast<int>(i % 3));
        }
        s.ipfLoadStatus.clear();
        const fs::path file = temp / "with_ipf.xyz";
        saveStructure(s, file.string(), "xyz");
        const std::string sidecar = readText(temp / "with_ipf.atomforge-ipf");
        Structure back;
        loadStructureFromFile(file.string(), back, error);
        o << "{\"sidecar_text\":" << quoted(sidecar) << ",\"reloaded\":" << structureJson(back) << "}";
    }
    o << ",\n\"trajectories\":[";
    {
        // Trajectories: every frame, or the error; the type-only dump with and without type elements.
        const fs::path traj = fixtures.parent_path() / "traj";
        std::vector<std::pair<std::string, std::string>> cases;  // file, types
        for (const auto& entry : fs::directory_iterator(traj)) cases.push_back({entry.path().filename().string(), ""});
        std::sort(cases.begin(), cases.end());
        cases.push_back({"type_only.dump", "Cu Ni"});
        cases.push_back({"scaled_triclinic.lammpstrj", "Al,Mg"});
        bool firstCase = true;
        for (const auto& [file, types] : cases) {
            setLammpsTypeElements(types);
            std::string error;
            std::vector<Structure> frames;
            try {
                frames = loadTrajectory((traj / file).string());
            } catch (const std::exception& e) {
                error = e.what();
            }
            o << (firstCase ? "" : ",") << "\n{\"file\":" << quoted(file) << ",\"types\":" << quoted(types) << ",\"error\":" << quoted(error) << ",\"frames\":[";
            for (std::size_t k = 0; k < frames.size(); ++k) o << (k ? "," : "") << structureJson(frames[k]);
            o << "]}";
            firstCase = false;
        }
        setLammpsTypeElements("");
    }
    o << "]}\n";
    std::ofstream(output / "io.json") << o.str();
    fs::remove_all(temp);
    std::printf("wrote %s\n", (output / "io.json").string().c_str());
    return 0;
}
