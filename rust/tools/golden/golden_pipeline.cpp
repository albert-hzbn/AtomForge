// Reference outputs of the C++ pipeline (pipeline.json) for the Rust
// af-pipeline tests: text and JSON forms, expression values and the outputs
// and stage results of pipelines on several input structures.
// Usage (from the repository root): golden_pipeline <output dir>.

#include "pipeline/Expression.h"
#include "pipeline/ModifierSupport.h"
#include "pipeline/Pipeline.h"
#include "science/ScienceData.h"
#include "util/ElementData.h"

#include <cmath>

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using namespace atomforge::pipeline;
using atomforge::science::Json;

namespace
{
template <typename F>
Json attempt(F&& f)
{
    try {
        return f();
    } catch (const std::exception& e) {
        Json error = Json::object();
        error["error"] = std::string(e.what());
        return error;
    }
}

Json vec3(const std::array<double, 3>& v) { return Json::array({v[0], v[1], v[2]}); }

Json structureDump(const Structure& s)
{
    Json o = Json::object();
    Json atoms = Json::array();
    for (const auto& a : s.atoms) {
        Json atom = Json::object();
        atom["symbol"] = a.symbol;
        atom["z"] = a.atomicNumber;
        atom["position"] = Json::array({a.x, a.y, a.z});
        atom["colour"] = Json::array({double(a.r), double(a.g), double(a.b)});
        atom["fixed"] = int(a.fixed);
        atom["magmom"] = a.magmom;
        atom["charge"] = a.charge;
        atom["velocity"] = Json::array({a.velocity[0], a.velocity[1], a.velocity[2]});
        atoms.push(atom);
    }
    o["atoms"] = atoms;
    o["has_cell"] = s.hasUnitCell;
    o["cell"] = Json::array({vec3(s.cellVectors[0]), vec3(s.cellVectors[1]), vec3(s.cellVectors[2])});
    o["pbc"] = Json::array({s.pbc[0], s.pbc[1], s.pbc[2]});
    o["offset"] = vec3(s.cellOffset);
    Json colours = Json::array();
    for (const auto& c : s.grainColors) colours.push(Json::array({double(c[0]), double(c[1]), double(c[2])}));
    o["grain_colors"] = colours;
    Json regions = Json::array();
    for (int r : s.grainRegionIds) regions.push(r);
    o["grain_region_ids"] = regions;
    Json property = Json::array(), frozen = Json::array();
    for (double v : s.atomProperty) property.push(std::isnan(v) ? Json("nan") : Json(v));
    for (double v : s.frozenProperty) frozen.push(std::isnan(v) ? Json("nan") : Json(v));
    o["atom_property"] = property;
    o["frozen_property"] = frozen;
    o["atom_property_name"] = s.atomPropertyName;
    o["frozen_property_name"] = s.frozenPropertyName;
    o["loop_points"] = int(s.dislocationLoopPoints.size());
    o["detection_done"] = s.dislocationDetectionDone;
    return o;
}

AtomSite atom(const char* symbol, double x, double y, double z)
{
    AtomSite a;
    a.symbol = symbol;
    a.atomicNumber = atomforge::science::atomicNumber(symbol);
    a.x = x; a.y = y; a.z = z;
    getDefaultElementColor(a.atomicNumber, a.r, a.g, a.b);
    return a;
}

std::vector<Structure> inputs()
{
    std::vector<Structure> out;
    // fcc Cu with an O interstitial, per-atom metadata and a property.
    Structure fcc;
    fcc.hasUnitCell = true;
    fcc.cellVectors = {{{3.61, 0, 0}, {0, 3.61, 0}, {0, 0, 3.61}}};
    for (const auto& p : std::vector<std::array<double, 3>>{{0, 0, 0}, {0, 1.805, 1.805}, {1.805, 0, 1.805}, {1.805, 1.805, 0}})
        fcc.atoms.push_back(atom("Cu", p[0], p[1], p[2]));
    fcc.atoms.push_back(atom("O", 1.805, 1.805, 1.805));
    fcc.atoms[1].magmom = 0.5;
    fcc.atoms[2].charge = -0.25;
    fcc.atoms[3].velocity = {0.01, 0.02, -0.03};
    fcc.atoms[4].fixed = 5;
    fcc.grainColors = {{0.1f, 0.2f, 0.3f}, {0.4f, 0.5f, 0.6f}, {0.7f, 0.8f, 0.9f}, {1, 0, 0}, {0, 1, 0}};
    fcc.grainRegionIds = {0, 1, 1, 2, 3};
    fcc.atomProperty = {0.1, 0.5, std::nan(""), 0.9, 1.5};
    fcc.atomPropertyName = "test";
    out.push_back(fcc);
    // A triclinic cell with an offset, mixed elements and atoms outside the cell.
    Structure tri;
    tri.hasUnitCell = true;
    tri.cellVectors = {{{4.2, 0, 0}, {1.3, 3.9, 0}, {0.7, -0.4, 5.1}}};
    tri.cellOffset = {0.5, -0.25, 0.1};
    tri.pbc = {true, true, false};
    const char* symbols[] = {"Na", "Cl", "Fe", "O", "H", "Cu", "Zn", "O"};
    for (int i = 0; i < 8; ++i)
        tri.atoms.push_back(atom(symbols[i], -1.0 + 0.9 * i, 0.3 + 0.55 * ((i * 3) % 7), -0.8 + 0.85 * ((i * 5) % 8)));
    tri.frozenProperty = {1, 2, 3, 4, 5, 6, 7, 8};
    out.push_back(tri);
    // A molecule without a cell.
    Structure water;
    water.atoms = {atom("O", 0, 0, 0.1173), atom("H", 0, 0.7572, -0.4692), atom("H", 0, -0.7572, -0.4692)};
    out.push_back(water);
    // An empty structure.
    out.push_back(Structure{});
    return out;
}

const char* const kTexts[] = {
    "",
    "   ",
    "replicate 2 2 1",
    "replicate 2,2,1 | select-expression fz > 0.5 && element == Cu | delete-selected",
    "select-expression \"x > 1 || y < 2\" mode=add",
    "select-element Cu, O mode=subtract",
    "select-element elements=Cu O mode=intersect",
    "select-slab c 0.2 0.7",
    "select-sphere center=1,2,3 radius=2.5",
    "select-sphere 1 2 3 2.5 add",
    "displace vector=0.1 0 0",
    "displace 0.1 0 0 selected true false",
    "set-cell a=5,0,0 b=0,5,0 c=0,0,5 scale_atoms=no",
    "supercell \"1 1 0 -1 1 0 0 0 2\"",
    "supercell matrix='2 0 0 0 2 0 0 0 1'",
    "merge file=\"my file.xyz\" offset=1,2,3",
    "replicate 2 2",
    "replicate 2 2 2 2",
    "replicate count=2,2,2",
    "replicat 2 2 2",
    "select-random fraction=x",
    "select-random fraction=0.5 seed=1.5",
    "fix-atoms axes=w",
    "set-pbc a=maybe",
    "replicate 2 2 2 |",
    "| replicate 2 2 2",
    "replicate 2 2 2 || wrap",
    "select-expression \"x > 1",
    "wrap | center | shrink-wrap 1.5 | add-vacuum b 12 false",
    "strain normal=0.01,-0.02,0 shear=0,0,0.05 | transform 1 0 0 0 1 0 0 0.5 1 transform_cell=false",
    "select-expression x<=1 | select-expression index==2 mode=add | select-expression x>=2 mode=subtract",
    "select-expression a=1",
    "random-displacement 0.1 7 selected",
    "slice normal=1,1,1 distance=2.5 width=1 invert=true action=select",
    "add-atom element=fe position=0.5,0.5,0.5 fractional=true select=false",
    "assign-element Zn target=all | free-atoms all",
    "expand-selection 2.5 2 nearest=6",
    "select-property minimum=-1e300 maximum=1e300",
    "replicate counts=1e400,1,1",
    "replicate counts=0x2,1,1",
};

struct Run { const char* text; std::vector<int> disabled; };
const Run kRuns[] = {
    {"select-element Cu O | delete-selected", {}},
    {"select-expression \"fz > 0.5 && element == Cu\" | invert-selection", {}},
    {"select-slab axis=c 0.2 0.7 | fix-atoms axes=xy | select-expression fixed >= 2", {}},
    {"select-slab axis=z -0.5 1.9 mode=add | select-slab x 0 1 mode=intersect", {}},
    {"select-sphere 1 1 1 2.5 | expand-selection cutoff=3 steps=2 nearest=4", {}},
    {"select-sphere 0 0 0 1.5 | expand-selection cutoff=2.6 steps=1", {}},
    {"select-random fraction=0.4 seed=7 | random-displacement amplitude=0.1 seed=3 target=selected", {}},
    {"select-random fraction=1 seed=-3 | random-displacement 0.05 -1", {}},
    {"replicate 2 1 3 | select-property 0 1", {}},
    {"supercell '1 1 0 -1 1 0 0 0 2' | wrap", {}},
    {"supercell '0 1 0 1 0 0 0 0 -1'", {}},
    {"transform 1 0.1 0 0 1 0 0 0 1 | strain normal=0.01,0,-0.02 shear=0,0.03,0", {}},
    {"select-element O | transform 0 1 0 1 0 0 0 0 1 target=selected", {}},
    {"set-cell 5 0 0 0 5 0 0 0 5 | center", {}},
    {"set-cell 5 0 0 0 5 0 0 0 5 scale_atoms=false | center", {}},
    {"add-vacuum c 10 | shrink-wrap padding=2", {}},
    {"add-vacuum a 3.5 false", {}},
    {"add-atom O 0.5 0.5 0.5 fractional=true | assign-element Fe target=selected", {}},
    {"add-atom N 1 2 3 select=false | select-expression Z == 7", {}},
    {"slice 0 0 1 distance=2 width=1 invert=true | slice normal=1,1,0 distance=1 action=select", {}},
    {"slice 0 0 1 distance=0.5", {}},
    {"set-pbc true false true | displace 0.25 0 0 fractional=true wrap=true", {}},
    {"set-pbc false false false | wrap", {}},
    {"displace 0.3 -0.2 0.1 wrap=true", {}},
    {"merge rust/tests/fixtures/science/structure.json offset=1,2,3 | select-element Fe", {}},
    {"merge rust/tests/fixtures/science/molecule.json select=false use_cell=false", {}},
    {"merge missing_file.xyz", {}},
    {"clear-selection | free-atoms target=all | select-expression \"index % 2 == 0 or Z == 8\" mode=add", {}},
    {"replicate 2 2 2 | select-property 0 1 | delete-selected", {}},
    {"replicate 2 1 1 | select-element Cu | delete-selected | replicate 1 2 1", {1}},
    {"select-element Cu | delete-selected | add-atom H 0 0 0", {0}},
    {"fix-atoms | free-atoms", {}},
    {"select-expression \"element == 'Cu' && frozen > 2\" | invert-selection | fix-atoms x target=selected", {}},
    {"select-expression property", {}},
    {"select-expression \"element\"", {}},
    {"select-expression \"element > Cu\"", {}},
    {"select-expression \"sqrt(x*x + y*y) < 2.6 and not (z > 1)\"", {}},
    {"set-pbc false true true | select-sphere 0.1 0.1 0.1 2", {}},
    {"replicate 101 1 1", {}},
    {"insert-interstitials H kind=octahedral count=0", {}},
    {"replicate 2 2 2 | insert-interstitials C count=3 order=random seed=4", {}},
    {"insert-interstitials He kind=tetrahedral count=2 resolution=10 clearance=0.4 separation=0.7", {}},
    {"insert-interstitials H kind=irregular count=0 | select-element H | invert-selection", {}},
    {"insert-interstitials Xx", {}},
    {"insert-interstitials H count=-1", {}},
    // EditSteps
    {"mirror 1 1 0", {}},
    {"select-element O | mirror normal=0,0,1 center=point point=1,1,1 target=selected", {}},
    {"select-element O | mirror 0 1 0 center=cell target=selected wrap=true", {}},
    {"mirror 0 0 0", {}},
    {"invert center=origin wrap=true", {}},
    {"invert center=cell", {}},
    {"rotate axis=1,1,1 angle=30", {}},
    {"rotate mode=euler euler=10,20,30 convention=zxz rotate_cell=false center=origin", {}},
    {"rotate mode=euler euler=5,-40,70 convention=xyz", {}},
    {"rotate mode=euler euler=5,-40,70 convention=zyz target=all center=centroid rotate_cell=false", {}},
    {"align-direction 1 1 1 lattice=true axis=z second_direction=1,-1,0 second_axis=x", {}},
    {"align-direction direction=0,1,1 lattice=false axis=x", {}},
    {"align-direction 1 0 0 axis=z second_axis=z", {}},
    {"align-direction direction=1,0,0 lattice=false axis=y second_direction=2,0,0 second_axis=x", {}},
    {"sort-atoms key=z then=element tolerance=0.5", {}},
    {"sort-atoms key=element element_order='O Cu Cl' descending=true", {}},
    {"sort-atoms key=distance then=x descending=true", {}},
    {"sort-atoms key=property", {}},
    {"sort-atoms key=mass then=selected", {}},
    {"select-element O | sort-atoms key=selected then=y", {}},
    {"sort-atoms key=x tolerance=1.0 then=z", {}},
    {"assign-element C target=all | saturate-hydrogen", {}},
    {"assign-element Si target=all | saturate-hydrogen hybridization=sp3 terminator=F max_per_atom=2", {}},
    {"assign-element N target=all | saturate-hydrogen bond_scale=1.0", {}},
    {"select-element O | saturate-hydrogen elements=O target=selected bond_length=1.0 select=false", {}},
    {"clear-selection | saturate-hydrogen target=selected", {}},
    {"saturate-hydrogen elements='Cu O' terminator=Cl", {}},
    {"swap-elements Cu O", {}},
    {"swap-elements first=H second=O", {}},
    {"measure-geometry atoms='0 1'", {}},
    {"measure-geometry atoms='0 1 2'", {}},
    {"measure-geometry atoms='0 1 2 3'", {}},
    {"measure-geometry atoms='0 x'", {}},
    {"set-geometry atoms='0 1' value=1.2", {}},
    {"set-geometry atoms='1 0 2' value=120 move=both", {}},
    {"set-geometry atoms='0 1 2' value=95 move=atom", {}},
    {"set-geometry atoms='0 1 2 3' value=60 move=atom", {}},
    {"set-geometry atoms='0 1 2 3' value=-120 move=both", {}},
    {"set-geometry atoms='0 1 2' value=200", {}},
    {"set-geometry atoms='0 0'", {}},
    {"set-geometry atoms='1 2' value=2 move=both", {}},
    {"select-element O | duplicate-selected offset=0,0,2 copies=2", {}},
    {"select-element H | duplicate-selected select=false", {}},
    {"clear-selection | duplicate-selected", {}},
    {"niggli-reduce", {}},
    {"supercell '1 1 0 -1 1 0 0 0 1' | niggli-reduce", {}},
    {"set-cell 3 0 0 2.9 1 0 0.5 0.3 4 | niggli-reduce", {}},
    {"scale-cell 1.05", {}},
    {"scale-cell mode=volume value=50 scale_atoms=false", {}},
    {"scale-cell mode=density value=8", {}},
    {"scale-cell mode=axes factors=1.1,0.9,1", {}},
    {"scale-cell mode=foo", {}},
    {"unwrap-molecules", {}},
    {"displace 3 0 0 wrap=true | unwrap-molecules bond_scale=1.5", {}},
    {"select-fragment rule=size min_atoms=2", {}},
    {"select-element H | select-fragment mode=add", {}},
    {"set-magnetic-moments 2.5", {}},
    {"set-magnetic-moments mode=elements moments='Cu 1.5 O -0.5'", {}},
    {"set-magnetic-moments mode=layered axis=x period=1 layer_tolerance=0.3", {}},
    {"replicate 2 2 1 | set-magnetic-moments mode=layered axis=z period=2", {}},
    {"set-magnetic-moments mode=checkerboard", {}},
    {"replicate 2 1 1 | set-magnetic-moments mode=random seed=5", {}},
    {"set-magnetic-moments mode=elements moments='Cu'", {}},
    {"set-magnetic-moments mode=elements moments='Cu 1x'", {}},
    {"set-charges 0.5 neutralize=shift", {}},
    {"set-charges mode=elements charges='Cu 2 O -2' neutralize=scale", {}},
    {"set-charges mode=elements charges='Na 1 Cl -1 Fe 2 H 0.4' neutralize=scale", {}},
    {"set-charges -1 neutralize=scale", {}},
    {"set-charges 1 neutralize=bad", {}},
    {"set-velocities 300 seed=3", {}},
    {"set-velocities temperature=500 zero_rotation=true seed=11", {}},
    {"set-velocities 0", {}},
    {"set-velocities -1", {}},
    {"stack-layers file=rust/tests/fixtures/science/structure.json gap=2", {}},
    {"stack-layers file=rust/tests/fixtures/science/molecule.json copies=2 shift=0.5,0.5,0 twist=15", {}},
    {"stack-layers file=rust/tests/fixtures/science/structure.json match=average", {}},
    {"stack-layers file=rust/tests/fixtures/science/structure.json match=substrate copies=2 vacuum=5", {}},
    {"stack-layers", {}},
    {"scale-cell value=1.05", {}},
    // ShapeSteps
    {"build-shape script='sphere(r = 5);'", {}},
    {"build-shape script='R = 6; sphere(r = R); region(\"Au\") sphere(r = 0.5 * R);' vacuum=4", {}},
    {"build-shape script='difference() { box(size = [10, 10, 10]); sphere(r = 3); }'", {}},
    {"build-shape script='shell(t = 2) sphere(r = 6);' anchor=cell-centre", {}},
    {"build-shape script='facets(hkl = [[1, 1, 1], [1, 0, 0]], d = [6, 7], family = true);'", {}},
    {"build-shape script='sphere(r = 2.5);' mode=keep", {}},
    {"build-shape script='sphere(r = 2);' mode=delete anchor=centroid", {}},
    {"build-shape script='box(size = [3, 3, 3]);' mode=select combine=add position=1,1,1", {}},
    {"build-shape script='sphere(r = 3); region(\"none\") halfspace(normal = [0, 0, -1], d = 0);' mode=keep anchor=origin", {}},
    {"build-shape script='sphere(r = 3); region(\"Ag\") halfspace(normal = [0, 0, -1], d = 0);' mode=keep anchor=origin", {}},
    {"build-shape script='halfspace(normal = [0, 0, 1], d = 0);'", {}},
    {"build-shape script='sphere(r = ;'", {}},
    {"build-shape script='sphere(r = 50);' max_atoms=100", {}},
    {"build-shape script='sphere(r = 4);' anchor=origin mode=select", {}},
    {"replicate 2 2 2 | build-shape script='torus(R = 4, r = 1.5);' mode=keep", {}},
    {"build-shape script='sphere(r = 0.1);' position=0.3,0.3,0.3", {}},
    // GeneratorSteps
    {"build-monolayer", {}},
    {"build-monolayer kind=hbn na=2 nb=2 layers=2 stacking=AB", {}},
    {"build-monolayer kind=mx2-2h layers=2 stacking=2H rectangular=true strain=1.5", {}},
    {"build-monolayer kind=mx2-1t elements='W Se' layers=3 stacking=ABC periodic_c=false", {}},
    {"build-monolayer kind=silicene na=2", {}},
    {"build-monolayer kind=phosphorene layers=2 stacking=AB", {}},
    {"build-monolayer kind=phosphorene layers=2 stacking=AA nb=2", {}},
    {"build-monolayer kind=stanene thickness=1 a=4.8", {}},
    {"build-monolayer kind=germanene elements='Ge Ge Ge'", {}},
    {"build-monolayer vacuum=-1", {}},
    {"build-nanotube", {}},
    {"build-nanotube n=8 m=0 length=2", {}},
    {"build-nanotube n=6 m=3 elements='B N'", {}},
    {"build-nanotube n=5 m=5 walls=2 axis=x", {}},
    {"build-nanotube n=8 m=0 finite=true length_A=12 end_cap=H axis=y", {}},
    {"build-nanotube n=4 m=4 finite=true end_cap=none", {}},
    {"build-nanotube n=4 m=2 elements='Si C' bond=0", {}},
    {"build-nanotube n=3 m=5", {}},
    {"build-nanotube n=6 m=3 walls=2", {}},
    {"build-nanoribbon", {}},
    {"build-nanoribbon edge=zigzag width=4 length=2", {}},
    {"build-nanoribbon width=5 finite=true length=3 termination=F", {}},
    {"build-nanoribbon edge=zigzag width=3 elements='B N' passivate=false", {}},
    {"build-nanoribbon edge=zigzag width=2 finite=true length=2 termination=Cl termination_bond=1.7", {}},
    {"build-nanoribbon width=8 termination=none", {}},
    {"build-twisted-bilayer index=2", {}},
    {"build-twisted-bilayer angle=13 stacking=AB", {}},
    {"build-twisted-bilayer kind=hbn top_kind=graphene index=1", {}},
    {"build-twisted-bilayer kind=mx2-2h index=1 top_kind=mx2-1t top_elements='W Se'", {}},
    {"build-twisted-bilayer kind=phosphorene", {}},
    {"build-twisted-bilayer index=0", {}},
    {"build-cluster", {}},
    {"build-cluster kind=decahedron shells=4", {}},
    {"build-cluster kind=cuboctahedron shells=3 second_element=Ag fraction=0.3 distribution=core", {}},
    {"build-cluster kind=octahedron shells=4 second_element=Pd distribution=shell fraction=0.5", {}},
    {"build-cluster kind=icosahedron shells=3 second_element=Au distribution=random seed=9 fraction=0.25", {}},
    {"build-cluster kind=decahedron shells=3 second_element=Ni distribution=janus", {}},
    {"build-cluster kind=decahedron shells=1", {}},
    {"build-cluster shells=1 element=pt", {}},
    {"build-cluster second_element=Ag fraction=1.5", {}},
    {"build-molecule", {}},
    {"build-molecule name=benzene orientation=principal", {}},
    {"build-molecule name=C60 vacuum=0", {}},
    {"build-molecule name=CH3OH orientation=random seed=4 periodic=true", {}},
    {"build-molecule name=Fe", {}},
    {"build-molecule name=rust/tests/fixtures/science/molecule.json orientation=principal", {}},
    {"build-molecule name=xyz123", {}},
    {"build-molecule name='carbon dioxide'", {}},
    {"build-molecule name=NH3 orientation=principal", {}},
    {"build-molecule name=C2H6 orientation=principal", {}},
    {"build-molecule name=H2CO orientation=random seed=2", {}},
    {"build-packing molecules='H2O 8' density=1", {}},
    {"build-packing molecules='H2O 4 CO2 2' box=10,10,10 seed=3", {}},
    {"build-packing molecules='CH4 3' box=12,12,12 region=sphere region_center=6,6,6 region_radius=4", {}},
    {"build-packing molecules='N2 0' density=0.8 box=8,8,8", {}},
    {"build-packing molecules='Ar 0' density=1.4 box=10,10,10 region=slab region_low=0,0,2 region_high=0,0,8 min_distance=3", {}},
    {"build-packing molecules='H2O 2' mode=solvate min_distance=1.5 solute_distance=1.2", {}},
    {"build-packing molecules='H2O x'", {}},
    {"build-packing molecules='H2O 200' box=5,5,5 attempts=30", {}},
    {"build-packing molecules='H2O 3' fixed_orientation=true box=9,9,9 region=cylinder region_center=4.5,4.5,0 region_radius=3 region_axis=z", {}},
    {"build-packing molecules='He 4' box=8,8,8 region=outside-sphere region_center=4,4,4 region_radius=2", {}},
    {"build-packing molecules='He 4' box=8,8,8 region=box region_low=1,1,1 region_high=5,5,5", {}},
    {"add-vacuum c 10 | add-adsorbate", {}},
    {"add-vacuum c 10 | replicate 2 2 1 | add-adsorbate species=O site=fcc coverage=all", {}},
    {"add-vacuum c 10 | replicate 2 2 1 | add-adsorbate species=O site=hcp coverage=all surface=bottom", {}},
    {"add-vacuum c 10 | replicate 2 2 1 | add-adsorbate species=H site=bridge coverage=fraction fraction=0.5 seed=3 surface=both min_separation=1", {}},
    {"add-vacuum c 10 | add-adsorbate species=H2O anchor=0 tilt=30 azimuth=45 principal=true height=2.1 position=0.25,0.25,0", {}},
    {"add-vacuum c 10 | add-adsorbate species=CO anchor=5", {}},
    {"add-vacuum c 10 | add-adsorbate species=NH3 site=bridge select=false", {}},
    {"near-cubic-supercell target_atoms=40 range=1", {}},
    {"near-cubic-supercell target_atoms=20 shape=fcc range=1 tolerance=20", {}},
    {"near-cubic-supercell target_atoms=30 shape=image range=1", {}},
    {"near-cubic-supercell shape=bcc range=1 target_atoms=10", {}},
    {"near-cubic-supercell target_atoms=32 range=2", {}},
    {"near-cubic-supercell range=4", {}},
    {"interpolate file=rust/tests/fixtures/pipeline/fcc_end.xyz fraction=0.3", {}},
    {"interpolate file=rust/tests/fixtures/pipeline/fcc_end.xyz method=idpp fraction=0.5 idpp_steps=50", {}},
    {"interpolate file=rust/tests/fixtures/pipeline/water_end.xyz method=idpp idpp_cutoff=1.2 fraction=0.5", {}},
    {"interpolate file=rust/tests/fixtures/pipeline/water_end.xyz method=idpp fraction=0.25 idpp_steps=200", {}},
    {"interpolate", {}},
    {"interpolate file=rust/tests/fixtures/science/molecule.json", {}},
    {"create-antisite elements='Cu O'", {}},
    {"replicate 2 2 2 | create-antisite elements='Cu O' count=3 min_separation=2 partner=random seed=5", {}},
    {"create-antisite elements='Na Cl' count=2", {}},
    {"create-antisite elements='Cu Cu'", {}},
    {"replicate 2 2 2 | create-frenkel-pair min_distance=3", {}},
    {"replicate 2 2 2 | create-frenkel-pair interstitial=dumbbell direction=1,0,0 min_distance=2 seed=4", {}},
    {"replicate 2 2 2 | create-frenkel-pair interstitial=dumbbell lattice_direction=false direction=1,1,0 dumbbell_length=2 max_distance=4", {}},
    {"create-frenkel-pair element=Zn", {}},
    // AnalysisSteps
    {"select-overlapping distance=2 | invert-selection", {}},
    {"select-overlapping 2.6 both mode=add", {}},
    {"select-overlapping 0", {}},
    {"select-surface probe_radius=2", {}},
    {"replicate 2 2 2 | select-surface 3 mode=add", {}},
    {"select-region regions='region r1 sphere 1.8 1.8 1.8 2 units box'", {}},
    {"select-region regions='region b block INF INF INF INF 1 INF units lattice' lattice_scale=2 mode=add", {}},
    {"assign-color 0.2 0.4 0.6 all", {}},
    {"select-element O | assign-color color=1,0,0", {}},
    {"assign-color 1.5 0 0", {}},
    {"ambient-occlusion intensity=0.5 directions=12", {}},
    {"ambient-occlusion intensity=2", {}},
    {"clear-property | reduce-property", {}},
    {"freeze-property | compute-expression 'property*2' name=double", {}},
    {"clear-property | freeze-property", {}},
    {"select-element Cu | compute-expression 'z*2' target=selected | reduce-property target=selected bins=4", {}},
    {"compute-expression 'x+y' | reduce-property bins=3", {}},
    {"compute-property structure-type", {}},
    {"replicate 2 2 2 | compute-property cna", {}},
    {"replicate 2 2 2 | compute-property ptm", {}},
    {"replicate 2 2 2 | compute-property orientation", {}},
    {"replicate 2 2 2 | compute-property ptm-ordering", {}},
    {"compute-property vorotop", {}},
    {"compute-property diamond", {}},
    {"replicate 2 2 2 | compute-property centrosymmetry", {}},
    {"compute-property coordination cutoff=3", {}},
    {"compute-property q6 cutoff=2.7", {}},
    {"compute-property bond-valence", {}},
    {"compute-property voronoi-volume", {}},
    {"compute-property voronoi-coordination", {}},
    {"compute-property surface", {}},
    {"replicate 2 2 2 | compute-property planar-faults", {}},
    {"compute-property chill-plus", {}},
    {"compute-property local-entropy cutoff=4", {}},
    {"compute-property rings cutoff=2.6", {}},
    {"compute-property bogus", {}},
    {"replicate 2 2 2 | compute-property cna | freeze-property | random-displacement 0.4 5 | compute-property cna | select-expression 'property != frozen'", {}},
    {"replicate 2 2 2 | elastic-strain", {}},
    {"replicate 2 2 2 | elastic-strain component=xy a=3.6", {}},
    {"replicate 2 2 2 | elastic-strain component=shear lattice=bcc a=2.87", {}},
    {"displacements", {}},
    {"displacements reference=rust/tests/fixtures/pipeline/fcc_end.xyz component=x", {}},
    {"displacements reference=rust/tests/fixtures/pipeline/fcc_end.xyz affine_mapping=to_reference", {}},
    {"bin-property bins=1,1,4", {}},
    {"bin-property 2 1 1 count", {}},
    {"compute-expression z | bin-property 1 1 3 operation=max", {}},
    {"replicate 3 3 3 | segment-grains min_size=5", {}},
    {"replicate 2 2 2 | segment-grains threshold=10 min_size=1 adopt=false", {}},
    // scienceStepTypes
    {"relax calculator='{\"potential\": \"EMT\"}' steps=5", {}},
    {"relax calculator='{\"potential\": \"LJ\", \"sigma\": 2.3}' steps=3 relax_cell=true fmax=0.05", {}},
    {"nvt-dynamics calculator='{\"potential\": \"EMT\"}' steps=4 sample_interval=2 seed=1", {}},
    {"npt-dynamics calculator='{\"potential\": \"EMT\"}' steps=3 sample_interval=1", {}},
    {"standardize-cell", {}},
    {"standardize-cell cell=conventional symprec_A=0.1", {}},
    {"standardize-cell cell=symmetrized", {}},
    {"relax", {}},
    {"relax calculator='{bad'", {}},
    // build-shape surfaces (OBJ files compared byte for byte by the Rust test)
    {"build-shape script='sphere(r = 5);' mesh_file=rust/tests/golden/shape_sphere.obj mesh_spacing=0.8", {}},
    {"build-shape script='difference() { box(size = [8, 8, 8]); sphere(r = 3); }' mode=keep mesh_file=rust/tests/golden/shape_void.obj mesh_spacing=1", {}},
    {"build-shape script='torus(R = 5, r = 2);' mode=select mesh_file=rust/tests/golden/shape_torus.obj mesh_spacing=0.9", {}},
    {"build-shape script='halfspace(normal = [0, 0, 1], d = 0);' mode=select mesh_file=rust/tests/golden/shape_unbounded.obj", {}},
    {"build-shape script='sphere(r = 5);' mesh_file=rust/tests/golden/shape_bad.obj mesh_spacing=0", {}},
};

const char* const kExpressions[] = {
    "x", "x + y * z - fx / 2", "-x ^ 2", "2 ^ -1", "(x + 1) % 3", "-index % 2", "x % 0",
    "fz > 0.5 && element == Cu", "element != \"O\" or Z >= 26", "not selected", "!(x < 1)", "x != 1",
    "abs(y - 2) + sqrt(abs(z)) + exp(-x) + log(1 + index)", "sin(x) + cos(y) + tan(z)", "floor(x) + ceil(-y)",
    "min(x, y) - max(z, index)", "pi * 2", "count - index", "property", "frozen > 2", "fixed", "magmom + charge",
    "vx + vy + vz", "atomic_number == Z", "type == cu", "element == FE", "1e-3 * x", ".5 + 1.", "0x1p2 + 0X10",
    "x > 1 > 2", "x >", "(x", "foo", "elemnt == Cu", "abs(x, y)", "min(x)", "bar(x)", "1 +", "\"Cu\"",
    "element + 1", "1 ? 2", "'unterminated", "x andy", "x and y", "x or y", "x &&& y", "1e999", "1e-400",
    "x == y == z", "--x", "+-x", "not not x", "!=x", "",
};
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::fprintf(stderr, "usage: golden_pipeline <output dir>\n");
        return 1;
    }
    Json root = Json::object();

    Json types = Json::array();
    for (const auto& t : modifierTypes()) types.push(std::string(t.id));
    root["types"] = types;

    Json texts = Json::array();
    for (const char* text : kTexts) {
        Json item = Json::object();
        item["text"] = std::string(text);
        item["parse"] = attempt([&] {
            const Pipeline p = Pipeline::parse(text);
            Json r = Json::object();
            r["json"] = p.toJson();
            r["text"] = p.toText();
            r["round_trip"] = Pipeline::fromJson(Json::parse(p.toJson().dump())).toText();
            return r;
        });
        texts.push(item);
    }
    root["texts"] = texts;

    // JSON forms: values of other types, unknown parameters, labels.
    Json jsons = Json::array();
    for (const char* text : {
             R"({"modifiers": [{"type": "replicate", "parameters": {"counts": [3, 1, 2]}, "label": "big", "enabled": false}]})",
             R"({"modifiers": [{"type": "replicate", "parameters": {"counts": "1 2 3"}}]})",
             R"({"modifiers": [{"type": "select-slab", "parameters": {"minimum": 0.25, "axis": "b", "mode": "add"}}]})",
             R"({"modifiers": [{"type": "select-slab", "parameters": {"axis": "w"}}]})",
             R"({"modifiers": [{"type": "select-slab", "parameters": {"minimun": 1}}]})",
             R"({"modifiers": [{"type": "set-pbc", "parameters": {"a": false, "b": "no"}}]})",
             R"({"modifiers": [{"type": "set-pbc", "parameters": {"a": 1}}]})",
             R"({"modifiers": [{"type": "merge", "parameters": {"file": 5}}]})",
             R"({"modifiers": [{"type": "random-displacement", "parameters": {"seed": 2.5}}]})",
             R"({"modifiers": [{"type": "nonsense"}]})",
             R"({"modifiers": [{"parameters": {}}]})",
             R"({"modifier": []})",
             R"([1, 2])",
             R"({"modifiers": [{"type": "wrap", "enabled": "yes"}]})",
             R"({"modifiers": [{"type": "displace", "parameters": {"vector": [1, 2]}}]})",
         }) {
        Json item = Json::object();
        item["text"] = std::string(text);
        item["result"] = attempt([&] {
            const Pipeline p = Pipeline::fromJson(Json::parse(text));
            Json r = Json::object();
            r["json"] = p.toJson();
            r["text"] = p.toText();
            return r;
        });
        jsons.push(item);
    }
    root["jsons"] = jsons;

    const auto structures = inputs();
    Json inputList = Json::array();
    for (const auto& s : structures) inputList.push(structureDump(s));
    root["inputs"] = inputList;

    Json expressions = Json::array();
    for (const char* text : kExpressions) {
        Json item = Json::object();
        item["text"] = std::string(text);
        Json values = Json::array();
        for (const auto& s : structures) {
            values.push(attempt([&] {
                const Expression e(text);
                PipelineData data;
                data.structure = s;
                data.clearSelection();
                if (!data.selected.empty()) data.selected[0] = 1;
                Json list = Json::array();
                for (double v : support::evaluateForAtoms(e, data)) list.push(std::isnan(v) ? Json("nan") : Json(v));
                return list;
            }));
        }
        item["values"] = values;
        expressions.push(item);
    }
    root["expressions"] = expressions;

    Json runs = Json::array();
    for (const Run& run : kRuns) {
        Json item = Json::object();
        item["text"] = std::string(run.text);
        Json disabled = Json::array();
        for (int d : run.disabled) disabled.push(d);
        item["disabled"] = disabled;
        Json outputs = Json::array();
        for (const auto& s : structures) {
            outputs.push(attempt([&] {
                Pipeline p = Pipeline::parse(run.text);
                for (int d : run.disabled) p.modifiers[static_cast<std::size_t>(d)].enabled = false;
                std::vector<StageResult> stages;
                const PipelineData out = p.evaluate(s, &stages);
                Json r = Json::object();
                r["structure"] = structureDump(out.structure);
                Json selected = Json::array();
                for (char c : out.selected) selected.push(int(c));
                r["selected"] = selected;
                Json notes = Json::array();
                for (const auto& n : out.notes) notes.push(n);
                r["notes"] = notes;
                Json stageList = Json::array();
                for (const auto& st : stages) {
                    Json o = Json::object();
                    o["atoms"] = st.atoms;
                    o["selected"] = st.selected;
                    o["error"] = st.error;
                    Json stageNotes = Json::array();
                    for (const auto& n : st.notes) stageNotes.push(n);
                    o["notes"] = stageNotes;
                    o["skipped"] = st.skipped;
                    stageList.push(o);
                }
                r["stages"] = stageList;
                return r;
            }));
        }
        item["outputs"] = outputs;
        runs.push(item);
    }
    root["runs"] = runs;

    std::ofstream(std::string(argv[1]) + "/pipeline.json") << root.dump(1);
    std::printf("wrote %s/pipeline.json\n", argv[1]);
    return 0;
}
