// Reference results of the C++ chemistry library (chemistry.json) for the Rust
// af-chemistry tests: SMILES reading and writing, molecular properties,
// perception from coordinates, CIP labels and hydrogen placement.
// Usage (from the repository root): golden_chemistry <fixture dir> <output dir>.

#include "chemistry/Molecule.h"
#include "chemistry/Smiles.h"
#include "science/Json.h"
#include "util/ElementData.h"

#include <cstdio>
#include <fstream>
#include <functional>
#include <sstream>

using namespace atomforge::chemistry;
using atomforge::science::Json;

namespace
{
Json moleculeJson(const Molecule& m)
{
    Json o = Json::object();
    Json atoms = Json::array();
    for (const auto& a : m.atoms)
        atoms.push(Json::array({a.z, a.charge, a.hydrogens, a.aromatic, a.isotope, a.mapNumber, a.dummy, a.position[0], a.position[1], a.position[2]}));
    o["atoms"] = atoms;
    Json bonds = Json::array();
    for (const auto& b : m.bonds) bonds.push(Json::array({b.a, b.b, b.order, b.aromatic, b.dative, b.image[0], b.image[1], b.image[2]}));
    o["bonds"] = bonds;
    Json tetra = Json::array();
    for (const auto& t : m.tetrahedral) tetra.push(Json::array({t.centre, t.n[0], t.n[1], t.n[2], t.n[3], t.sign}));
    o["tetrahedral"] = tetra;
    Json doubles = Json::array();
    for (const auto& d : m.doubleBonds) doubles.push(Json::array({d.a, d.x, d.y, d.b, d.cis}));
    o["double_bonds"] = doubles;
    return o;
}
Json ints(const std::vector<int>& v)
{
    Json out = Json::array();
    for (int x : v) out.push(x);
    return out;
}
Json rings(const std::vector<std::vector<int>>& list)
{
    Json out = Json::array();
    for (const auto& r : list) out.push(ints(r));
    return out;
}
Json properties(const Molecule& m)
{
    Json o = Json::object();
    o["formula"] = formula(m);
    o["molar_mass"] = molarMass(m);
    o["monoisotopic_mass"] = monoisotopicMass(m);
    o["charge"] = totalCharge(m);
    o["electrons"] = electronCount(m);
    o["multiplicity"] = multiplicity(m);
    o["fragments"] = rings(fragments(m));
    o["rotatable"] = ints(rotatableBonds(m));
    o["ranks"] = ints(canonicalRanks(m));
    const auto r = smallestRings(m);
    o["rings"] = rings(r);
    o["ring_sizes"] = ints(smallestRingSize(m, r));
    std::vector<int> hyb, valences, hcounts;
    for (int i = 0; i < static_cast<int>(m.atoms.size()); ++i) { hyb.push_back(hybridisation(m, i)); valences.push_back(m.valence(i)); hcounts.push_back(m.hydrogenCount(i)); }
    o["hybridisation"] = ints(hyb);
    o["valences"] = ints(valences);
    o["hydrogen_counts"] = ints(hcounts);
    Json labels = Json::array();
    for (const auto& l : cipLabels(m)) labels.push(Json::array({l.atom, l.other, l.label}));
    o["cip"] = labels;
    o["smiles"] = writeSmiles(m, false);
    o["smiles_stereo"] = writeSmiles(m, true);
    return o;
}
Json cases = Json::array();
void attempt(const std::string& name, const std::function<Json()>& work)
{
    Json entry = Json::object();
    entry["name"] = name;
    try { entry["result"] = work(); }
    catch (const SmilesError& e) { entry["error"] = std::string(e.what()); entry["position"] = static_cast<double>(e.position); }
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
}

int main(int argc, char** argv)
{
    if (argc < 3) { std::fprintf(stderr, "usage: golden_chemistry <fixture dir> <output dir>\n"); return 1; }
    const std::string in = std::string(argv[1]) + "/";
    const char* smiles[] = {
        "CCO", "c1ccccc1", "C1=CC=CC=C1", "c1ccncc1", "c1cc[nH]c1", "c1ccnc1", "o1cccc1", "s1cccc1", "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O",
        "N[C@@H](C)C(=O)O", "F/C=C/F", "F/C=C\\F", "C/C=C(/C)C", "[NH4+]", "CC(=O)[O-]", "[Na+].[Cl-]", "[13CH4]", "[CH3:1]C", "C%10CC%10",
        "c1ccc2ccccc2c1", "O=C1NC(=O)C=C1", "[Fe+2]", "*CC*", "CC(C)(C)C", "C#N", "O=S(=O)(O)O", "P(=O)(O)(O)O", "C[N+](=O)[O-]",
        "[O-][N+](=O)c1ccccc1", "OC[C@@H](O)[C@H](O)CO", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "C1CCCCCCC1", "C/1=C/CCCCCC1", "CC.CC",
        "C-1CC1", "[2H]O", "C$C", "c1ccccc1-c2ccccc2", "[C@H](F)(Cl)Br", "CCO ethanol", "[CH]", "[C]", "[O]", "O=O", "[Cu+2].[O-2]",
        "[se]1cccc1", "c1ccc[se]1", "[as]1cccc1", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "[nH]1cccc1", "C[S@](=O)c1ccccc1", "[C@@](F)(Cl)(Br)I",
        "N#[N+][O-]", "[B-](F)(F)(F)F", "[Si](C)(C)(C)C", "c1cccc1", "C1CC", "C(", "C)", "(C)", "[Xx]", "[c]1ccccc1", "Q", "%1", "C%1",
        "C%1a", "1C", "C11", "[CH3", "[C@SP1H4]", "CC==C", "", " ", "C.", "C1.C1", "C=1CC-1", "C=1CCC#1", "C[NH3+]", "[O-]C(=O)C",
    };
    for (const char* s : smiles)
        attempt(std::string("smiles ") + s, [&] {
            Molecule m = parseSmiles(s);
            Json o = Json::object();
            o["molecule"] = moleculeJson(m);
            o["properties"] = properties(m);
            Molecule h = m;
            o["explicit_added"] = makeHydrogensExplicit(h, true);
            o["explicit"] = moleculeJson(h);
            o["explicit_smiles"] = writeSmiles(h, true);
            return o;
        });
    for (const char* name : {"aspirin", "caffeine", "benzene", "ethanol", "chfclbr", "chfclbr_mirror", "z_butene", "e_butene", "ammonium_acetate", "polyethylene"}) {
        const Structure s = readXyz(in + name + ".xyz");
        attempt(std::string("perceive ") + name, [&] {
            Json o = Json::object();
            const Molecule m = perceive(s);
            o["molecule"] = moleculeJson(m);
            o["properties"] = properties(m);
            PerceiveOptions plain;
            plain.bondOrders = false;
            plain.bondScale = 1.3;
            o["plain"] = moleculeJson(perceive(s, plain));
            // Hydrogens removed, counted back and placed by hybridisation.
            Molecule heavy = m;
            std::vector<char> keep;
            for (const auto& a : heavy.atoms) keep.push_back(a.z != 1);
            o["map"] = ints(heavy.removeAtoms(keep));
            assignImplicitHydrogens(heavy, {});
            o["placed_added"] = makeHydrogensExplicit(heavy, true);
            o["placed"] = moleculeJson(heavy);
            Json structure = Json::array();
            for (const auto& a : toStructure(m, true).atoms) structure.push(Json::array({a.symbol, a.atomicNumber, a.x, a.y, a.z, a.charge, a.r, a.g, a.b}));
            o["structure"] = structure;
            Json vectors = Json::array();
            for (int j = 1; j < static_cast<int>(s.atoms.size()) && j < 4; ++j) {
                const auto v = bondVector(s, 0, j);
                vectors.push(Json::array({v[0], v[1], v[2]}));
            }
            o["bond_vectors"] = vectors;
            return o;
        });
    }
    attempt("elements", [&] {
        Json o = Json::array();
        for (int z = 0; z <= 118; ++z) {
            int valences[8];
            for (int s = 0; s < 8; ++s) valences[s] = defaultValence(z, s % 3 - 1, s);
            Json row = Json::array({symbolOf(z), isMetal(z), covalentRadius(z), hydrogenBondLength(z), atomicNumber(symbolOf(z))});
            for (int v : valences) row.push(v);
            o.push(row);
        }
        for (const char* s : {"cl", "CL", "D", "T", "Xx", ""}) o.push(Json::array({s, atomicNumber(s)}));
        return o;
    });
    std::ofstream(std::string(argv[2]) + "/chemistry.json") << cases.dump();
    std::printf("wrote %s/chemistry.json (%zu cases)\n", argv[2], cases.size());
    return 0;
}
