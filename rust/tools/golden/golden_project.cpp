// Reference AtomForge projects (model/Workspace) for af-electronic::project:
// two tabs with every archived field set, saved in formats 1-7 as
// project_v<N>.afproject. Usage: golden_project <output dir>.

#include "model/Workspace.h"

#include <cmath>
#include <cstdio>
#include <limits>

using namespace atomforge;

namespace
{
Structure structure(int n, double shift)
{
    Structure s;
    s.hasUnitCell = true;
    s.cellVectors = {{{4.05, 0.0, 0.0}, {0.1, 4.05, 0.0}, {0.0, 0.2, 4.05}}};
    s.cellOffset = {0.5, -0.25, shift};
    s.pbc = {true, false, true};
    s.pbcBoundaryTol = 0.125f;
    for (int i = 0; i < n; ++i) {
        AtomSite a;
        a.symbol = i % 2 ? "Al" : "Cu";
        a.atomicNumber = i % 2 ? 13 : 29;
        a.x = 0.1 * i + shift; a.y = 1.0 / 3.0 * i; a.z = -0.7 * i;
        a.r = 0.25f * (i % 4); a.g = 0.5f; a.b = 1.0f / 3.0f;
        a.fixed = static_cast<std::uint8_t>(i % 8);
        a.magmom = 0.1 + i; a.charge = -0.3 * i; a.velocity = {1e-3 * i, -2.5, 1.0 / 7.0};
        s.atoms.push_back(a);
        s.grainColors.push_back({0.1f * i, 0.2f, 0.3f});
        s.grainRegionIds.push_back(i / 2);
        s.atomProperty.push_back(i == 1 ? std::numeric_limits<double>::quiet_NaN() : 1.5 * i);
        s.frozenProperty.push_back(-1.0 * i);
    }
    s.atomPropertyName = "coordination";
    s.frozenPropertyName = "frozen energy";
    s.ipfLoadStatus = "IPF from sidecar";
    s.dislocationLoopPoints = {glm::vec3(1, 2, 3), glm::vec3(-1.5f, 0.25f, 9)};
    s.dislocationDetectionDone = true;
    return s;
}
}

int main(int argc, char** argv)
{
    if (argc < 2) { std::fprintf(stderr, "usage: golden_project <output dir>\n"); return 1; }
    std::vector<Workspace> tabs(2);
    Workspace& w = tabs[0];
    w.structure = structure(5, 0.0);
    electronic::Grid g;
    g.shape = {2, 3, 2};
    g.cell = glm::dmat3(glm::dvec3(4, 0, 0), glm::dvec3(0, 5, 0), glm::dvec3(0.5, 0, 6));
    g.origin = glm::dvec3(0.1, 0.2, 0.3);
    g.periodic = true;
    g.name = "density";
    g.unit = "e/A^3";
    for (int i = 0; i < 12; ++i) g.values.push_back(0.5 * i - 1.0 / 3.0);
    w.volume.sites = {{29, glm::dvec3(0, 0, 0)}, {13, glm::dvec3(1, 2, 3)}};
    w.volume.fields = {g};
    w.reference.fields = {g};
    w.surface.vertices = {glm::dvec3(0, 0, 0), glm::dvec3(1, 0, 0), glm::dvec3(0, 1, 0)};
    w.surface.colors = {0.1, 0.2, 0.3};
    w.table = {1, 2, 3, 4, 5, 6};
    w.columns = 3;
    w.heading = "x,y,z";
    w.settings = {{"isovalue", 0.05}, {"alpha", 1}};
    w.sliceSettings = {{"plane", 2}};
    w.camera = {30, -20, 5, 12.5, 0.1, 0.2, 0.3};
    w.sourcePath = "C:/data/CHGCAR";
    w.referencePath = "C:/data/AECCAR0";
    w.title = "Cu Al tab";
    w.history = {"Opened", "Replicate 2 1 1"};
    w.science = "{\"tool\":\"rdf\"}";
    w.pipeline = "{\"steps\":[]}";
    w.hasPipelineInput = true;
    w.pipelineInput = structure(3, 1.0);
    w.display = "ambient 0.3\n";
    tabs[1].structure = structure(2, 2.0);
    tabs[1].structure.atomProperty.clear();
    tabs[1].structure.frozenProperty.clear();
    tabs[1].title = "second";
    for (int format = 1; format <= 7; ++format) {
        const std::string path = std::string(argv[1]) + "/project_v" + std::to_string(format) + ".afproject";
        saveWorkspace(tabs, path, format);
    }
    std::printf("wrote %s/project_v1-7.afproject\n", argv[1]);
    return 0;
}
