// Reference grain-boundary colouring of the C++ desktop (app/EditorOps
// applyGrainBoundaryColors, in an anonymous namespace, so this file includes
// EditorOps.cpp itself) for af-render::colouring: colouring.json maps each
// fixture to the per-atom flag "boundary colour".
// Usage (from the repository root): golden_colouring <fixture dir> <output dir>.

#include "app/EditorOps.cpp"
#include "science/Json.h"
#include "util/ElementData.h"

#include <cstdio>
#include <fstream>
#include <sstream>

using atomforge::science::Json;

namespace
{
Structure readXyz(const std::string& path)
{
    std::ifstream in(path);
    std::string line;
    std::getline(in, line);
    const int n = std::stoi(line);
    std::getline(in, line);
    Structure s;
    if (line.rfind("cell ", 0) == 0) {
        std::istringstream c(line.substr(5));
        for (auto& row : s.cellVectors) for (double& v : row) c >> v;
        s.hasUnitCell = true;
        s.pbc = {true, true, true};
    }
    for (int i = 0; i < n; ++i) {
        std::getline(in, line);
        std::istringstream w(line);
        AtomSite a;
        w >> a.symbol >> a.x >> a.y >> a.z;
        int region;
        if (w >> region) s.grainRegionIds.push_back(region);
        for (int z = 1; z <= 118; ++z) if (a.symbol == elementSymbol(z)) a.atomicNumber = z;
        s.atoms.push_back(a);
    }
    return s;
}
}

int main(int argc, char** argv)
{
    if (argc < 3) { std::fprintf(stderr, "usage: golden_colouring <fixture dir> <output dir>\n"); return 1; }
    const std::string in = std::string(argv[1]) + "/";
    Json cases = Json::object();
    for (const char* name : {"cu_vacancy", "cu_cluster", "cu_perfect", "fe_bcc", "cu_jitter", "bicrystal", "amorphous", "three"}) {
        Structure s = readXyz(in + name + ".xyz");
        applyGrainBoundaryColors(s);
        Json flags = Json::array();
        for (const auto& a : s.atoms) flags.push(a.r == 1.0f);
        cases[name] = flags;
    }
    std::ofstream(std::string(argv[2]) + "/colouring.json") << cases.dump();
    std::printf("wrote %s/colouring.json (%zu cases)\n", argv[2], cases.size());
    return 0;
}

// App functions the linked app sources reference but this dumper never calls.
#include "graphics/Shader.h"
#include "ui/MeasurementOverlay.h"
#include "ui/TransformAtomsDialog.h"
GLuint createComputeProgram(const char*) { return 0; }
void MeasurementOverlayState::clearVisuals() {}
void TransformAtomsDialog::clearTransform() {}
