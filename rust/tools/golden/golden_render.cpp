// Reference scene data of the C++ renderer (render.json) for the Rust
// af-render tests: atom instances with periodic boundary images, ambient
// occlusion, bonds, the default view and the frame matrices.
// Usage (from the repository root): golden_render <fixture dir> <output dir>.

#include "app/EditorOps.h"
#include "app/SceneView.h"
#include "camera/Camera.h"
#include "graphics/AmbientOcclusion.h"
#include "graphics/GraphicsBackend.h"
#include "graphics/SceneBuffers.h"
#include "graphics/StructureInstanceBuilder.h"
#include "model/DisplaySettings.h"
#include "science/Json.h"
#include "util/ElementData.h"


#include <cstdio>
#include <fstream>
#include <sstream>

using atomforge::science::Json;

namespace
{
Json f3(const glm::vec3& v) { return Json::array({v.x, v.y, v.z}); }
Json vecs(const std::vector<glm::vec3>& list)
{
    Json out = Json::array();
    for (const auto& v : list) out.push(f3(v));
    return out;
}
Json floats(const std::vector<float>& list)
{
    Json out = Json::array();
    for (float v : list) out.push(v);
    return out;
}
Json ints(const std::vector<int>& list)
{
    Json out = Json::array();
    for (int v : list) out.push(v);
    return out;
}
Json mat(const glm::mat4& m)
{
    Json out = Json::array();
    for (int c = 0; c < 4; ++c)
        for (int r = 0; r < 4; ++r) out.push(m[c][r]);
    return out;
}
Json instancesJson(const StructureInstanceData& d)
{
    Json o = Json::object();
    o["positions"] = vecs(d.positions);
    o["colors"] = vecs(d.colors);
    o["scales"] = floats(d.scales);
    o["bond_radii"] = floats(d.bondRadii);
    o["shininess"] = floats(d.shininess);
    o["atomic_numbers"] = ints(d.atomicNumbers);
    o["box_lines"] = vecs(d.boxLines);
    o["atom_indices"] = ints(d.atomIndices);
    o["orbit_center"] = f3(d.orbitCenter);
    return o;
}
Json bondsJson(const SceneBuffers& b)
{
    Json o = Json::object();
    o["starts"] = vecs(b.bondStarts);
    o["ends"] = vecs(b.bondEnds);
    o["colors_a"] = vecs(b.bondColorsA);
    o["colors_b"] = vecs(b.bondColorsB);
    o["radii"] = floats(b.bondRadiiCpu);
    return o;
}
Json frameJson(const FrameView& f)
{
    Json o = Json::object();
    o["projection"] = mat(f.projection);
    o["view"] = mat(f.view);
    o["light_mvp"] = mat(f.lightMVP);
    o["light_position"] = f3(f.lightPosition);
    o["camera_position"] = f3(f.cameraPosition);
    o["depth"] = Json::array({f.depthNear, f.depthFar});
    return o;
}
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
}

int main(int argc, char** argv)
{
    if (argc < 3) { std::fprintf(stderr, "usage: golden_render <fixture dir> <output dir>\n"); return 1; }
    graphics::setBackend(graphics::Backend::Vulkan);  // CPU only: no OpenGL buffers
    const std::string in = std::string(argv[1]) + "/";
    const auto look = atomforge::display::defaultLook();
    const std::vector<float> shininess(119, 32.0f);
    Json cases = Json::array();
    for (const char* name : {"cu_fcc", "si_diamond", "nacl", "mg_hcp", "molecule", "si_super"}) {
        Structure s = readXyz(in + name + ".xyz");
        for (auto& a : s.atoms) {
            const glm::vec3 c = look.colours[static_cast<std::size_t>(a.atomicNumber)];
            a.r = c.r; a.g = c.g; a.b = c.b;
        }
        Json entry = Json::object();
        entry["name"] = name;
        const int identity[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        const int doubled[3][3] = {{2, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        const int skew[3][3] = {{1, 1, 0}, {-1, 1, 0}, {0, 0, 1}};
        StructureInstanceData data = buildStructureInstanceData(s, false, identity, look.radii, shininess);
        entry["instances"] = instancesJson(data);
        entry["doubled"] = instancesJson(buildStructureInstanceData(s, true, doubled, look.radii, shininess));
        entry["skew"] = instancesJson(buildStructureInstanceData(s, true, skew, look.radii, shininess));
        for (float& scale : data.scales) scale *= 1.3f;
        atomforge::display::darkenOccluded(data.positions, data.scales, data.colors, 0.6f);
        entry["occluded"] = vecs(data.colors);
        Json bonds = Json::array();
        for (int style = 0; style < 3; ++style) {
            SceneBuffers buffers;
            buffers.bondStyle.metalBonds = style == 1;
            buffers.bondStyle.uniformColour = style == 2;
            buffers.bondStyle.thickness = style == 2 ? 1.7f : 1.0f;
            buffers.bondStyle.tolerance = style == 2 ? 1.3f : 1.18f;
            buffers.upload(data, false, {});
            bonds.push(bondsJson(buffers));
        }
        entry["bonds"] = bonds;
        SceneBuffers buffers;
        buffers.upload(data, false, {});
        Json views = Json::array();
        for (const auto& [w, h] : std::vector<std::pair<int, int>>{{1600, 1200}, {64, 48}, {300, 900}}) {
            Camera camera;
            applyDefaultView(camera, buffers, w, h, true);
            Json view = Json::object();
            view["distance"] = camera.distance;
            for (int variant = 0; variant < 3; ++variant) {
                Camera c = camera;
                atomforge::display::DisplaySettings display;
                if (variant == 1) { c.yaw = 20; c.pitch = -15; c.roll = 30; display.lightFollowsCamera = !display.lightFollowsCamera; }
                if (variant == 2) { c.distance = 7.5f; display.lightAzimuth = 200; display.lightElevation = 80; }
                FrameView frame;
                frame.framebufferWidth = w;
                frame.framebufferHeight = h;
                buildFrameView(c, buffers, variant == 2, frame, display);
                view[std::string("frame") + std::to_string(variant)] = frameJson(frame);
            }
            views.push(view);
        }
        entry["views"] = views;
        cases.push(entry);
    }
    std::ofstream(std::string(argv[2]) + "/render.json") << cases.dump();
    std::printf("wrote %s/render.json (%zu cases)\n", argv[2], cases.size());
    return 0;
}

// App functions the linked app sources reference but this dumper never calls.
#include "graphics/Shader.h"
#include "ui/MeasurementOverlay.h"
#include "ui/TransformAtomsDialog.h"
GLuint createComputeProgram(const char*) { return 0; }
void MeasurementOverlayState::clearVisuals() {}
void TransformAtomsDialog::clearTransform() {}
