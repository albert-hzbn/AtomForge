// Reference results of the C++ chemistry library (chemistry_ops.json) for the
// Rust af-chemistry tests: the UFF force field, point groups, editing,
// builders, Z-matrices, quantum-chemistry inputs and distance-geometry
// embedding (graph, stereo and energy only: the coordinates come from a
// different eigensolver in Rust).
// Usage (from the repository root): golden_chemistry_ops <fixture dir> <output dir>.

#include "chemistry/Builders.h"
#include "chemistry/Editing.h"
#include "chemistry/Embedding.h"
#include "chemistry/ForceField.h"
#include "chemistry/Molecule.h"
#include "chemistry/PointGroup.h"
#include "chemistry/Smiles.h"
#include "science/Json.h"

#include <cstdio>
#include <fstream>
#include <functional>
#include <sstream>

using namespace atomforge::chemistry;
using atomforge::science::Json;

namespace
{
Json vec3(const Vec3& v) { return Json::array({v[0], v[1], v[2]}); }
Json positions(const std::vector<Vec3>& x)
{
    Json out = Json::array();
    for (const auto& v : x) out.push(vec3(v));
    return out;
}
Json graphJson(const Molecule& m)
{
    Json o = Json::object();
    Json atoms = Json::array();
    for (const auto& a : m.atoms) atoms.push(Json::array({a.z, a.charge, a.hydrogens, a.aromatic, a.dummy}));
    o["atoms"] = atoms;
    Json bonds = Json::array();
    for (const auto& b : m.bonds) bonds.push(Json::array({b.a, b.b, b.order, b.aromatic, b.dative}));
    o["bonds"] = bonds;
    Json tetra = Json::array();
    for (const auto& t : m.tetrahedral) tetra.push(Json::array({t.centre, t.n[0], t.n[1], t.n[2], t.n[3], t.sign}));
    o["tetrahedral"] = tetra;
    Json doubles = Json::array();
    for (const auto& d : m.doubleBonds) doubles.push(Json::array({d.a, d.x, d.y, d.b, d.cis}));
    o["double_bonds"] = doubles;
    o["formula"] = formula(m);
    return o;
}
Json moleculeJson(const Molecule& m)
{
    Json o = graphJson(m);
    std::vector<Vec3> x;
    for (const auto& a : m.atoms) x.push_back(a.position);
    o["positions"] = positions(x);
    Json images = Json::array();
    for (const auto& b : m.bonds) images.push(vec3(b.image));
    o["images"] = images;
    return o;
}
Json ints(const std::vector<int>& v)
{
    Json out = Json::array();
    for (int x : v) out.push(x);
    return out;
}
Json doubles(const std::vector<double>& v)
{
    Json out = Json::array();
    for (double x : v) out.push(x);
    return out;
}
Json termsJson(const EnergyTerms& t)
{
    return Json::array({t.bond, t.angle, t.torsion, t.inversion, t.vdw, t.electrostatic, t.restraint, t.total()});
}
Json minimizeJson(const MinimizeResult& r) { return Json::array({r.initialEnergy, r.energy, r.steps, r.converged, r.maxForce}); }
Json groupJson(const PointGroupResult& g)
{
    Json o = Json::object();
    o["name"] = g.name;
    o["order"] = g.order;
    o["rotational"] = g.rotationalSymmetryNumber;
    Json ops = Json::array();
    for (const auto& m : g.operations)
        ops.push(Json::array({m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2]}));
    o["operations"] = ops;
    o["centre"] = vec3(g.centre);
    o["principal"] = vec3(g.principalAxis);
    o["secondary"] = vec3(g.secondaryAxis);
    o["deviation"] = g.maxDeviation;
    return o;
}
Json cases = Json::array();
void attempt(const std::string& name, const std::function<Json()>& work)
{
    Json entry = Json::object();
    entry["name"] = name;
    try { entry["result"] = work(); }
    catch (const std::exception& e) { entry["error"] = std::string(e.what()); }
    cases.push(entry);
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
        s.hasUnitCell = true;
        s.pbc = {true, true, true};
    }
    for (int i = 0; i < n; ++i) {
        std::getline(in, line);
        std::istringstream w(line);
        AtomSite a;
        w >> a.symbol >> a.x >> a.y >> a.z;
        a.atomicNumber = atomicNumber(a.symbol);
        s.atoms.push_back(a);
    }
    return s;
}
std::vector<Vec3> coordinates(const Molecule& m)
{
    std::vector<Vec3> x;
    for (const auto& a : m.atoms) x.push_back(a.position);
    return x;
}
}

int main(int argc, char** argv)
{
    if (argc < 3) { std::fprintf(stderr, "usage: golden_chemistry_ops <fixture dir> <output dir>\n"); return 1; }
    const std::string in = std::string(argv[1]) + "/";
    auto load = [&](const char* name) { return readXyz(in + name + ".xyz"); };

    // ---- force field
    for (const char* name : {"aspirin", "caffeine", "benzene", "ethanol", "chfclbr", "z_butene", "ammonium_acetate", "polyethylene", "butane", "water", "sf6"}) {
        const Structure s = load(name);
        attempt(std::string("uff ") + name, [&] {
            const Molecule m = perceive(s);
            const Structure* cell = s.hasUnitCell ? &s : nullptr;
            Json o = Json::object();
            UffForceField field(m, cell);
            Json types = Json::array();
            for (const auto& t : field.types()) types.push(t);
            o["types"] = types;
            std::vector<Vec3> x = coordinates(m), g;
            EnergyTerms terms;
            o["energy"] = field.energy(x, &g, &terms);
            o["terms"] = termsJson(terms);
            o["gradient"] = positions(g);
            ForceFieldOptions eo;
            eo.electrostatics = true;
            eo.dielectric = 4.0;
            eo.cutoff = 9.0;
            UffForceField charged(m, cell, eo);
            EnergyTerms ct;
            std::vector<Vec3> cg;
            o["charged_energy"] = charged.energy(x, &cg, &ct);
            o["charged_terms"] = termsJson(ct);
            o["charged_gradient"] = positions(cg);
            UffForceField restrained(m, cell);
            const int n = static_cast<int>(m.atoms.size());
            Restraint r;
            r.kind = 2; r.atoms = {0, 1, 0, 0}; r.target = 1.7; r.k = 300;
            restrained.restraints.push_back(r);
            if (n >= 3) { r.kind = 3; r.atoms = {0, 1, 2, 0}; r.target = 100; r.k = 50; restrained.restraints.push_back(r); }
            if (n >= 4) { r.kind = 4; r.atoms = {0, 1, 2, 3}; r.target = 30; r.k = 20; restrained.restraints.push_back(r); }
            EnergyTerms rt;
            std::vector<Vec3> rg;
            o["restrained_energy"] = restrained.energy(x, &rg, &rt);
            o["restrained_terms"] = termsJson(rt);
            o["restrained_gradient"] = positions(rg);
            std::vector<Vec3> y = x;
            o["minimize"] = minimizeJson(minimize(field, y, {}, 300, 0.01));
            o["minimized"] = positions(y);
            std::vector<Vec3> z = x;
            std::vector<std::uint8_t> fixed(n, 0);
            fixed[0] = 7;
            fixed[1] = 1;
            o["minimize_fixed"] = minimizeJson(minimize(restrained, z, fixed, 100, 0.05));
            o["minimized_fixed"] = positions(z);
            return o;
        });
    }
    attempt("uff tables", [&] {
        Json o = Json::object();
        Json has = Json::array();
        for (int z = 0; z <= 120; ++z) has.push(uffHasElement(z));
        o["has"] = has;
        Json rows = Json::array();
        for (const char* t : {"C_3", "C_R", "C_2", "C_1", "N_3", "N_R", "O_3", "O_2", "H_", "S_3+2", "Fe3+2", "Cl", "P_3+5", "Si3"})
            rows.push(Json::array({t, uffAngle(t), uffVdwRadius(t), uffBondLength(t, "C_3", 1.0), uffBondLength(t, "O_2", 2.0), uffBondLength(t, "C_R", 1.5)}));
        o["rows"] = rows;
        return o;
    });
    attempt("uff unknown type", [&] { return Json(uffBondLength("Qq", "C_3", 1.0)); });
    attempt("uff implicit hydrogens", [&] { UffForceField f(parseSmiles("CCO")); return Json(0); });

    // ---- point groups
    for (const char* name : {"methane", "sf6", "water", "ammonia", "pcl5", "ethane", "trans_dichloroethylene", "allene", "co2", "hcn", "benzene",
                             "benzene_distorted", "ethanol", "aspirin", "caffeine", "chfclbr"}) {
        const Structure s = load(name);
        attempt(std::string("point group ") + name, [&] {
            std::vector<int> z;
            std::vector<std::array<double, 3>> x;
            for (const auto& a : s.atoms) { z.push_back(a.atomicNumber); x.push_back({a.x, a.y, a.z}); }
            Json o = Json::object();
            const PointGroupResult g = pointGroup(z, x, 0.1);
            o["group"] = groupJson(g);
            o["tight"] = groupJson(pointGroup(z, x, 0.01));
            std::vector<std::array<double, 3>> sym = x;
            o["symmetrized_group"] = groupJson(symmetrize(z, sym, 0.1));
            o["symmetrized"] = positions(sym);
            std::vector<std::array<double, 3>> oriented = x;
            standardOrientation(z, oriented, g);
            o["oriented"] = positions(oriented);
            return o;
        });
    }
    attempt("point group single", [&] { return groupJson(pointGroup({6}, {{1, 2, 3}}, 0.1)); });
    attempt("point group empty", [&] { return groupJson(pointGroup({}, {}, 0.1)); });

    // ---- editing
    for (const char* name : {"aspirin", "caffeine", "ethanol", "ammonium_acetate", "z_butene", "benzene"}) {
        const Structure s = load(name);
        attempt(std::string("gasteiger ") + name, [&] {
            const Molecule m = perceive(s);
            Json o = Json::object();
            o["six"] = doubles(gasteigerCharges(m));
            o["two"] = doubles(gasteigerCharges(m, 2));
            return o;
        });
    }
    attempt("rmsd", [&] {
        const Molecule a = perceive(load("chfclbr")), b = perceive(load("chfclbr_mirror"));
        Json o = Json::object();
        std::vector<Vec3> aligned;
        o["mirror"] = rmsd(coordinates(a), coordinates(b), {}, &aligned);
        o["aligned"] = positions(aligned);
        o["subset"] = rmsd(coordinates(a), coordinates(b), {0, 1, 2}, nullptr);
        const Molecule c = perceive(load("butane"));
        std::vector<Vec3> turned = coordinates(c);
        for (auto& p : turned) p = {0.3 * p[0] - 0.954 * p[1] + 1, 0.954 * p[0] + 0.3 * p[1] - 2, p[2] + 0.5};
        o["turned"] = rmsd(coordinates(c), turned, {}, &aligned);
        o["turned_aligned"] = positions(aligned);
        return o;
    });
    attempt("butane torsions", [&] {
        Molecule m = perceive(load("butane"));
        Json o = Json::object();
        o["side12"] = ints(sideOf(m, 1, 2));
        o["side21"] = ints(sideOf(m, 2, 1));
        o["dihedral"] = dihedralAngle(m.atoms[0].position, m.atoms[1].position, m.atoms[2].position, m.atoms[3].position);
        setDihedral(m, 0, 1, 2, 3, 60);
        o["set60"] = positions(coordinates(m));
        setDihedral(m, 0, 1, 2, 3, -75.5);
        o["set-75"] = positions(coordinates(m));
        return o;
    });
    attempt("dihedral in ring", [&] {
        Molecule m = perceive(load("benzene"));
        setDihedral(m, 0, 1, 2, 3, 10);
        return Json(0);
    });
    auto editCase = [&](const std::string& name, const char* file, const std::function<Json(Molecule&)>& work) {
        const Structure s = load(file);
        attempt(name, [&] {
            Molecule m = perceive(s);
            Json o = Json::object();
            o["value"] = work(m);
            o["molecule"] = moleculeJson(m);
            return o;
        });
    };
    editCase("bond order 2", "ethanol", [](Molecule& m) { return Json(setBondOrder(m, 0, 1, 2)); });
    editCase("bond order 3", "ethanol", [](Molecule& m) { return Json(setBondOrder(m, 0, 1, 3)); });
    editCase("bond order 0", "ethanol", [](Molecule& m) { return Json(setBondOrder(m, 0, 1, 0)); });
    editCase("bond order new", "ethanol", [](Molecule& m) { return Json(setBondOrder(m, 0, 2, 1)); });
    editCase("bond order hydrogen", "ethanol", [](Molecule& m) { int h = -1; for (int x : m.neighbours(0)) if (m.atoms[x].z == 1) h = x; return Json(setBondOrder(m, 0, h, 2)); });
    editCase("bond order same", "ethanol", [](Molecule& m) { return Json(setBondOrder(m, 1, 1, 1)); });
    editCase("bond order butene", "z_butene", [](Molecule& m) { return Json(setBondOrder(m, 0, 1, 1)); });
    editCase("add hydrogens", "aspirin", [](Molecule& m) {
        std::vector<char> keep;
        for (const auto& a : m.atoms) keep.push_back(a.z != 1);
        m.removeAtoms(keep);
        return Json(addHydrogens(m, {}));
    });
    editCase("add hydrogens masked", "caffeine", [](Molecule& m) {
        std::vector<char> keep;
        for (const auto& a : m.atoms) keep.push_back(a.z != 1);
        m.removeAtoms(keep);
        std::vector<char> mask(m.atoms.size(), 0);
        mask[0] = mask[2] = 1;
        return Json(addHydrogens(m, mask));
    });
    editCase("remove hydrogens", "aspirin", [](Molecule& m) { return Json(removeHydrogens(m, {})); });
    editCase("remove hydrogens masked", "ethanol", [](Molecule& m) { std::vector<char> mask(m.atoms.size(), 0); mask[0] = 1; return Json(removeHydrogens(m, mask)); });
    editCase("protonate ammonia", "ammonia", [](Molecule& m) { protonate(m, 0); return Json(0); });
    editCase("protonate water", "water", [](Molecule& m) { protonate(m, 0); return Json(0); });
    editCase("protonate hydrogen", "water", [](Molecule& m) { protonate(m, 1); return Json(0); });
    editCase("deprotonate ethanol", "ethanol", [](Molecule& m) { int o = -1; for (int i = 0; i < (int)m.atoms.size(); ++i) if (m.atoms[i].z == 8) o = i; deprotonate(m, o); return Json(o); });
    editCase("deprotonate none", "sf6", [](Molecule& m) { deprotonate(m, 0); return Json(0); });
    editCase("ionize aspirin", "aspirin", [](Molecule& m) { return Json(ionizeAtNeutralPH(m)); });
    editCase("ionize caffeine", "caffeine", [](Molecule& m) { return Json(ionizeAtNeutralPH(m)); });
    editCase("invert chfclbr", "chfclbr", [](Molecule& m) { invertStereocentre(m, 0); return Json(stereoMatches(m)); });
    editCase("invert three bonds", "ammonia", [](Molecule& m) { invertStereocentre(m, 0); return Json(0); });
    editCase("flip butene", "z_butene", [](Molecule& m) { flipDoubleBond(m, 0, 1); return Json(0); });
    editCase("flip unbonded", "z_butene", [](Molecule& m) { flipDoubleBond(m, 2, 3); return Json(0); });
    editCase("flip ring", "benzene", [](Molecule& m) { flipDoubleBond(m, 0, 1); return Json(0); });
    for (const char* s : {"NCC(=O)O", "CCN", "NC(N)=N", "OS(=O)(=O)c1ccccc1", "OP(=O)(O)O", "CN(C)C", "Nc1ccccc1", "CC(N)=O", "C[NH3+]", "OCCN"})
        attempt(std::string("ionize ") + s, [&] {
            Molecule m = parseSmiles(s);
            makeHydrogensExplicit(m, true);
            Json o = Json::object();
            o["text"] = ionizeAtNeutralPH(m);
            o["molecule"] = moleculeJson(m);
            return o;
        });
    attempt("group smiles", [&] {
        Json o = Json::array();
        for (const char* g : {"Me", "methyl", "tBu", "TBU", "ph", "COOH", "coo-", "*CC", "nitro", "Tert_Butyl", "hydrogen", "b(oh)2", "OAc"}) o.push(groupSmiles(g));
        Json groups = Json::array();
        for (const auto& [n, s] : functionalGroups()) groups.push(Json::array({n, s}));
        o.push(groups);
        return o;
    });
    attempt("group unknown", [&] { return Json(groupSmiles("unobtainium")); });
    for (const auto& [file, atom, group] : std::vector<std::tuple<const char*, int, const char*>>{
             {"ethanol", 3, "methyl"}, {"ethanol", 3, "phenyl"}, {"ethanol", 4, "*C(F)(F)F"}, {"benzene", 6, "nitro"}, {"benzene", 7, "COOH"},
             {"ethanol", 0, "methyl"}, {"ethanol", 3, "unknown-group"}, {"ethanol", 3, "*C*"}, {"ethanol", 3, "*=C"}, {"ethanol", 99, "methyl"}}) {
        const Structure s = load(file);
        attempt(std::string("substitute ") + file + " " + std::to_string(atom) + " " + group, [&] {
            Molecule m = perceive(s);
            const SubstituteResult r = substitute(m, atom, group);
            Json o = Json::object();
            o["new_atoms"] = ints(r.newAtoms);
            o["attachment"] = r.attachment;
            o["graph"] = graphJson(m);
            return o;
        });
    }
    for (const auto& [file, method, trials] : std::vector<std::tuple<const char*, const char*, int>>{
             {"butane", "random", 12}, {"butane", "systematic", 0}, {"ethanol", "random", 8}, {"aspirin", "random", 4}, {"benzene", "random", 5}}) {
        const Structure s = load(file);
        attempt(std::string("conformers ") + file + " " + method, [&] {
            const Molecule m = perceive(s);
            ConformerOptions co;
            co.method = method;
            co.trials = trials;
            co.seed = 3;
            co.optimizeSteps = 150;
            co.maxConformers = 4;
            co.rmsdThreshold = 0.3;
            int rotors = -1;
            const auto found = conformerSearch(m, co, &rotors);
            Json o = Json::object();
            o["rotors"] = rotors;
            Json list = Json::array();
            for (const auto& c : found) list.push(Json::array({c.energy, positions(c.positions)}));
            o["conformers"] = list;
            return o;
        });
    }
    for (const auto& [atoms, start, stop, step, relaxed] : std::vector<std::tuple<std::vector<int>, double, double, double, bool>>{
             {{0, 1, 2, 3}, 0, 180, 60, true}, {{0, 1, 2, 3}, 180, -180, 90, false}, {{0, 1}, 1.4, 1.7, 0.1, false}, {{0, 1}, 1.4, 1.7, 0.15, true},
             {{0, 1, 2}, 100, 120, 10, false}, {{0, 1, 2}, 100, 120, 10, true}, {{0}, 0, 1, 1, false}, {{0, 1}, 0, 1, 0, false}}) {
        const Structure s = load("butane");
        std::string name = "scan";
        for (int a : atoms) name += " " + std::to_string(a);
        name += std::string(relaxed ? " relaxed " : " rigid ") + std::to_string(start) + " " + std::to_string(step);
        attempt(name, [&] {
            const Molecule m = perceive(s);
            Json list = Json::array();
            for (const auto& p : geometryScan(m, atoms, start, stop, step, relaxed, 80)) list.push(Json::array({p.value, p.energy, positions(p.positions)}));
            return list;
        });
    }

    // ---- builders
    attempt("compound library", [&] {
        Json o = Json::array();
        for (const auto& c : compoundLibrary()) o.push(Json::array({c.name, c.category, c.smiles}));
        return o;
    });
    attempt("compound smiles", [&] {
        Json o = Json::array();
        for (const char* n : {"water", "Ethanol", "diethyl-ether", "DIETHYL_ETHER", "meoh", "Tylenol", "bpy", "lys", "nothing", "", "1,4-dioxane", "14dioxane"})
            o.push(compoundSmiles(n));
        return o;
    });
    attempt("polymer library", [&] {
        Json o = Json::array();
        for (const auto& [n, s] : polymerLibrary()) o.push(Json::array({n, s}));
        return o;
    });
    attempt("ligand library", [&] {
        Json o = Json::array();
        for (const auto& l : ligandLibrary()) o.push(Json::array({l.name, l.smiles, ints(l.donors), l.description}));
        return o;
    });
    struct PeptideCase { const char* sequence; const char* conformation; const char* termini; bool ionize; bool d; double phi, psi, omega; };
    for (const PeptideCase& pc : std::vector<PeptideCase>{
             {"ACDEFGHIKLMNPQRSTVWY", "alpha", "charged", true, false, 0, 0, 0}, {"Ala-Gly-Ser", "beta", "neutral", true, false, 0, 0, 0},
             {"GPG", "3-10", "capped", true, false, 0, 0, 0}, {"KRDE", "extended", "charged", false, false, 0, 0, 0},
             {"AW", "custom", "neutral", true, false, -60, -40, 175}, {"AA", "ppii", "charged", true, true, 0, 0, 0},
             {"HIE-CYX,ash glh LYN", "pi", "capped", true, false, 0, 0, 0}, {"y q", "parallel-beta", "Charged", true, false, 0, 0, 0},
             {"", "alpha", "charged", true, false, 0, 0, 0}, {"AXZ", "alpha", "charged", true, false, 0, 0, 0},
             {"AG", "spiral", "charged", true, false, 0, 0, 0}, {"AG", "alpha", "zwitter", true, false, 0, 0, 0}}) {
        attempt(std::string("peptide ") + pc.sequence + " " + pc.conformation + " " + pc.termini + (pc.d ? " D" : ""), [&] {
            PeptideOptions o;
            o.sequence = pc.sequence;
            o.conformation = pc.conformation;
            o.termini = pc.termini;
            o.ionizeSideChains = pc.ionize;
            o.dAmino = pc.d;
            if (std::string(pc.conformation) == "custom") { o.phi = pc.phi; o.psi = pc.psi; o.omega = pc.omega; }
            std::vector<std::string> labels;
            const Molecule m = buildPeptide(o, &labels);
            Json r = Json::object();
            r["molecule"] = moleculeJson(m);
            Json l = Json::array();
            for (const auto& s : labels) l.push(s);
            r["labels"] = l;
            return r;
        });
    }
    struct PolymerCase { const char* monomer; int units; const char* tacticity; const char* conformation; const char* head; const char* tail; bool periodic; unsigned seed; };
    for (const PolymerCase& pc : std::vector<PolymerCase>{
             {"PE", 5, "isotactic", "trans", "H", "H", false, 1}, {"PS", 4, "syndiotactic", "trans", "CH3", "H", false, 1},
             {"pp", 6, "atactic", "trans", "H", "phenyl", false, 7}, {"PE", 6, "isotactic", "random", "H", "H", false, 5},
             {"PE", 4, "isotactic", "trans", "H", "H", true, 1}, {"*CC(Cl)*", 3, "isotactic", "trans", "methyl", "methyl", false, 2},
             {"PS", 3, "syndiotactic", "trans", "H", "H", true, 1}, {"*CC", 3, "isotactic", "trans", "H", "H", false, 1},
             {"PE", 0, "isotactic", "trans", "H", "H", false, 1}, {"PE", 3, "random", "trans", "H", "H", false, 1},
             {"PE", 3, "isotactic", "helix", "H", "H", false, 1}, {"PET", 2, "isotactic", "trans", "H", "H", false, 1}}) {
        attempt(std::string("polymer ") + pc.monomer + " " + std::to_string(pc.units) + " " + pc.tacticity + " " + pc.conformation + " " + pc.head + " " +
                    pc.tail + (pc.periodic ? " periodic" : ""),
                [&] {
                    PolymerOptions o;
                    o.monomer = pc.monomer;
                    o.units = pc.units;
                    o.tacticity = pc.tacticity;
                    o.conformation = pc.conformation;
                    o.headGroup = pc.head;
                    o.tailGroup = pc.tail;
                    o.periodic = pc.periodic;
                    o.seed = pc.seed;
                    const PolymerResult r = buildPolymer(o);
                    Json j = Json::object();
                    j["graph"] = graphJson(r.molecule);
                    j["monomer"] = r.monomerSmiles;
                    return j;
                });
    }
    struct ComplexCase { const char* metal; int ox; const char* geometry; const char* ligands; const char* spin; double distance; };
    for (const ComplexCase& cc : std::vector<ComplexCase>{
             {"Fe", 2, "octahedral", "NH3 6", "high", 0}, {"Co", 3, "octahedral", "en 3", "low", 0}, {"Ni", 2, "square-planar", "Cl 4", "high", 0},
             {"Cu", 2, "tetrahedral", "Cl 4", "low", 0}, {"Pt", 2, "square planar", "NH3 2 Cl 2", "high", 2.05}, {"Ru", 2, "octahedral", "bipy 3", "low", 0},
             {"Ni", 2, "octahedral", "smiles:NCCN:0,3 2 H2O 2", "high", 0}, {"Zn", 2, "tetrahedral", "py 2 Cl 2", "high", 0},
             {"Fe", 3, "trigonal-bipyramidal", "CO 5", "high", 0}, {"V", 4, "square-pyramidal", "O 1 F 4", "high", 0},
             {"Na", 1, "linear", "H2O 2", "high", 0}, {"Cr", 3, "octahedral", "acac 3", "high", 0},
             {"Xx", 2, "octahedral", "NH3 6", "high", 0}, {"Fe", 2, "octahedral", "NH3 5", "high", 0}, {"Fe", 2, "hexagonal", "NH3 6", "high", 0},
             {"Fe", 2, "octahedral", "foo 6", "high", 0}, {"Fe", 2, "octahedral", "smiles:NCCN", "high", 0}, {"Fe", 2, "octahedral", "smiles:NCC(N)CN:0,4,5 2", "high", 0},
             {"Fe", 2, "octahedral", "smiles:NCCN:0,9 3", "high", 0}}) {
        attempt(std::string("complex ") + cc.metal + " " + cc.geometry + " " + cc.ligands + " " + cc.spin, [&] {
            ComplexOptions o;
            o.metal = cc.metal;
            o.oxidationState = cc.ox;
            o.geometry = cc.geometry;
            o.ligands = cc.ligands;
            o.spin = cc.spin;
            o.distance = cc.distance;
            const ComplexResult r = buildComplex(o);
            Json j = Json::object();
            j["graph"] = graphJson(r.molecule);
            j["charge"] = r.charge;
            j["d"] = r.dElectrons;
            j["multiplicity"] = r.multiplicity;
            return j;
        });
    }
    const char* zmatrices[] = {
        "O\nH 1 0.96\nH 1 0.96 2 104.5\n",
        "O\nH 1 r\nH 1 r 2 a\n\nr 0.97\na=105.0\n",
        "C\nC 1 1.54\nH 1 1.09 2 109.5\nH 1 1.09 2 109.5 3 120\nH 1 1.09 2 109.5 3 -120\nH 2 1.09 1 109.5 3 60\nH 2 1.09 1 109.5 6 120\nH 2 1.09 1 109.5 6 -120\n",
        "%chk=x.chk\n# HF/STO-3G\n\ntitle line\n\n0 1\nC1\nO2 C1 rco\nH3 C1 rch O2 ang\nH4 C1 rch O2 ang H3 180.0\nVariables:\nrco 1.2\nrch 1.1\nang 121.0\n",
        "X\nN 1 1.0\nH 2 1.01 1 90.0\nH 2 1.01 1 90.0 3 120.0\nH 2 1.01 1 90.0 3 -120.0 ! comment\n",
        "C\nO 1 1.2\nH 0.5 0.9 0.1\nCl 1 1.8 2 -r 3 90\n\nr 120\n",
        "O\nH 1 r\n",
        "O\nH 2 0.96\n",
        "Q\nH 1 1\n",
        "O\nH 1\n",
        "",
        "  \n\n",
        "Fe\nCL 1 2.2\nCa 1 2.4 2 90\n",
        "O\nH 1 r\n\nr abc\n",
        "C\r\nH 1 1.09\r\nH 1 1.09 2 109.47\r\nH 1 1.09 2 109.47 3 120\r\n",
        "C\nC 1 1.5\nC 2 1.5 1 180\nC 3 1.5 2 120 1 0\n",
    };
    for (std::size_t k = 0; k < sizeof(zmatrices) / sizeof(zmatrices[0]); ++k)
        attempt("zmatrix " + std::to_string(k), [&] {
            const Molecule m = parseZMatrix(zmatrices[k]);
            Json o = Json::object();
            o["molecule"] = moleculeJson(m);
            o["written"] = writeZMatrix(m, false);
            return o;
        });
    for (const char* name : {"aspirin", "ethanol", "water", "co2", "hcn", "sf6"})
        attempt(std::string("write zmatrix ") + name, [&] {
            const Molecule m = perceive(load(name));
            return Json::array({writeZMatrix(m, false), writeZMatrix(m, true)});
        });
    attempt("quantum formats", [&] {
        Json o = Json::array();
        for (const char* f : {"orca", "ORCA", "nwchem", "psi4", "qchem", "mopac", "turbomole", "gzmat", "gamess", "zmat", "xyz", ""}) o.push(isQuantumFormat(f));
        return o;
    });
    {
        Structure ions = load("ammonium_acetate");
        ions.atoms[0].charge = 1;
        ions.atoms[8].charge = -1;
        Structure partial = load("ethanol");
        partial.atoms[0].charge = 0.25;
        Structure radical = load("water");
        radical.atoms.pop_back();
        Structure empty;
        Structure unknown = load("water");
        unknown.atoms[0].atomicNumber = 0;
        unknown.atoms[0].symbol = "O";
        const std::vector<std::pair<std::string, Structure>> structures = {
            {"ethanol", load("ethanol")}, {"polyethylene", load("polyethylene")}, {"ammonium_acetate", load("ammonium_acetate")},
            {"ions", ions}, {"partial", partial}, {"radical", radical}, {"empty", empty}, {"unknown", unknown}};
        for (const auto& [name, s] : structures) {
            attempt("charge " + name, [&] {
                int q = 0, mult = 0;
                chargeAndMultiplicity(s, q, mult);
                return Json::array({q, mult});
            });
            for (const char* f : {"orca", "nwchem", "psi4", "qchem", "mopac", "turbomole", "gzmat", "gamess", "zmat", "cp2k"})
                attempt("quantum " + name + " " + f, [&] {
                    std::ostringstream out;
                    std::string error;
                    const bool ok = writeQuantumInput(s, f, out, error);
                    return Json::array({ok, out.str(), error});
                });
        }
    }

    // ---- embedding (graph, stereo and energy; coordinates differ by design)
    for (const char* s : {"CCO", "c1ccccc1", "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O", "F/C=C/F", "F/C=C\\F", "CC(C)(C)C", "OC[C@@H](O)[C@H](O)CO",
                          "C1CCCCC1", "CC(=O)Oc1ccccc1C(=O)O", "C", "[Na+].[Cl-]", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "C1CC1", "N#N", "C=C=C",
                          "CC(N)C(C)O", "C/C=C/C=C\\C", "[C@H](F)(Cl)Br", "*CC*", "[Fe+2]", "C1CC2CCC1C2", "[Xx]", "C("}) {
        for (bool optimize : {true, false})
            attempt(std::string("embed ") + s + (optimize ? "" : " plain"), [&] {
                EmbedOptions eo;
                eo.attempts = 3;
                eo.seed = 11;
                eo.optimize = optimize;
                EmbedResult r;
                const Molecule m = moleculeFromSmiles(s, eo, &r);
                Json o = Json::object();
                o["graph"] = graphJson(m);
                o["energy"] = r.energy;
                o["succeeded"] = r.succeeded;
                o["unspecified"] = r.unspecifiedStereo;
                o["stereo"] = stereoMatches(m);
                return o;
            });
    }
    attempt("embed implicit hydrogens", [&] { Molecule m = parseSmiles("CCO"); embed(m); return Json(0); });
    attempt("embed empty", [&] { Molecule m; embed(m); return Json(0); });

    std::ofstream(std::string(argv[2]) + "/chemistry_ops.json") << cases.dump();
    std::printf("wrote %s/chemistry_ops.json (%zu cases)\n", argv[2], cases.size());
    return 0;
}
