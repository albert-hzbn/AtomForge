// Reference measurements of the C++ desktop (ui/MeasurementOverlay
// processMeasurementRequests) for af-render::measure: distance, angle,
// dihedral and Atom Info over the instances of the render fixtures, for
// selections of instances (periodic images included) in pick order.
// Usage (from the repository root): golden_measure <fixture dir> <output dir>.

#include "app/EditorOps.h"
#include "camera/Camera.h"
#include "graphics/GraphicsBackend.h"
#include "graphics/SceneBuffers.h"
#include "graphics/StructureInstanceBuilder.h"
#include "model/DisplaySettings.h"
#include "science/Json.h"
#include "ui/MeasurementOverlay.h"
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
        std::string word;
        if (c >> word && word == "origin") c >> s.cellOffset[0] >> s.cellOffset[1] >> s.cellOffset[2];
        s.hasUnitCell = true;
        s.pbc = {true, true, true};
    }
    for (int i = 0; i < n; ++i) {
        std::getline(in, line);
        std::istringstream w(line);
        AtomSite a;
        w >> a.symbol >> a.x >> a.y >> a.z;
        for (int z = 1; z <= 118; ++z) if (a.symbol == elementSymbol(z)) a.atomicNumber = z;
        s.atoms.push_back(a);
    }
    return s;
}

// Every instance of the first `atoms` atoms, atom by atom.
std::vector<int> imagesOf(const SceneBuffers& b, int atoms)
{
    std::vector<int> out;
    for (int atom = 0; atom < atoms; ++atom)
        for (int i = 0; i < (int)b.atomIndices.size(); ++i)
            if (b.atomIndices[i] == atom) out.push_back(i);
    return out;
}
}

int main(int argc, char** argv)
{
    if (argc < 3) { std::fprintf(stderr, "usage: golden_measure <fixture dir> <output dir>\n"); return 1; }
    graphics::setBackend(graphics::Backend::Vulkan);  // CPU only: no OpenGL buffers
    const std::string in = std::string(argv[1]) + "/";
    const auto look = atomforge::display::defaultLook();
    const std::vector<float> shininess(119, 32.0f);
    Json cases = Json::array();
    for (const char* name : {"cu_fcc", "si_diamond", "nacl", "mg_hcp", "molecule", "si_super"}) {
        Structure s = readXyz(in + name + ".xyz");
        if (std::string(name) == "molecule") {
            s.atoms[0].fixed = 5;
            s.atoms[0].magmom = 1.23456;
            s.atoms[0].charge = -0.5;
            s.atoms[0].velocity = {0.001, -2.5, 1e-7};
        }
        const int identity[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        StructureInstanceData data = buildStructureInstanceData(s, false, identity, look.radii, shininess);
        SceneBuffers buffers;
        buffers.upload(data, false, {});
        const int n = (int)buffers.atomPositions.size();
        std::vector<std::vector<int>> selections = {{}, {0}, {1}, {0, 1}, {1, 0}, {0, 1, 2}, {2, 0, 1}, {0, 1, 2, 3}, {3, 2, 1, 0},
                                                    {n - 1, 0}, {0, n / 2, n - 1}, {0, n / 3, 2 * n / 3, n - 1}, {0, 1, 2, 3, 4}};
        for (int k = 1; k <= 4; ++k) selections.push_back(imagesOf(buffers, k));
        Json entry = Json::object();
        entry["name"] = name;
        Json results = Json::array();
        for (const auto& selection : selections) {
            MeasurementOverlayState state;
            processMeasurementRequests(state, true, true, true, true, selection, buffers, s, look.radii);
            Json r = Json::object();
            Json sel = Json::array();
            for (int i : selection) sel.push(i);
            r["selection"] = sel;
            r["distance"] = Json::array({std::string(state.distanceValue), std::string(state.distanceMessage)});
            r["angle"] = Json::array({std::string(state.angleValue), std::string(state.angleMessage)});
            r["dihedral"] = Json::array({std::string(state.dihedralValue), std::string(state.dihedralMessage)});
            r["atom_info"] = std::string(state.atomInfoMessage);
            Json lines = Json::array();
            if (state.showDistanceLine) { lines.push(state.distanceLineIdx0); lines.push(state.distanceLineIdx1); }
            r["distance_lines"] = lines;
            Json angleLines = Json::array();
            if (state.showAngleLines) { angleLines.push(state.angleLineIdx0); angleLines.push(state.angleLineIdx1); angleLines.push(state.angleLineIdx2); }
            r["angle_lines"] = angleLines;
            results.push(r);
        }
        entry["results"] = results;
        cases.push(entry);
    }
    std::ofstream(std::string(argv[2]) + "/measure.json") << cases.dump();
    // The toolbar's Rotate (rotateCrystalAroundAxis): cameras after rotations.
    Json rotations = Json::array();
    const float starts[][3] = {{45.0f, 35.2643897f, 0.0f}, {0, 0, 0}, {120, -30, 15}, {-170, 80, -45}, {90, 0, 0}, {0, 90, 0}};
    for (const auto& start : starts)
        for (int axis = 0; axis < 3; ++axis)
            for (double angle : {5.0, -30.0, 90.0, 180.0, 0.1, 37.5}) {
                Camera camera;
                camera.yaw = start[0]; camera.pitch = start[1]; camera.roll = start[2];
                Json r = Json::object();
                r["start"] = Json::array({start[0], start[1], start[2]});
                r["axis"] = axis;
                r["angle"] = angle;
                Json chain = Json::array();
                for (int k = 0; k < 4; ++k) {
                    rotateCrystalAroundAxis(camera, axis, angle);
                    chain.push(Json::array({camera.yaw, camera.pitch, camera.roll}));
                }
                r["cameras"] = chain;
                rotations.push(r);
            }
    std::ofstream(std::string(argv[2]) + "/camera_rotate.json") << rotations.dump();
    std::printf("wrote %s/measure.json (%zu cases)\n", argv[2], cases.size());
    return 0;
}

// App functions the linked app sources reference but this dumper never calls.
#include "graphics/Shader.h"
#include "ui/TransformAtomsDialog.h"
GLuint createComputeProgram(const char*) { return 0; }
void TransformAtomsDialog::clearTransform() {}
