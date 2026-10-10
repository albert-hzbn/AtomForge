// Non-destructive modification pipeline: expressions, text syntax, every
// modifier, caching, error handling and that the input is never changed.
#include "pipeline/Expression.h"
#include "pipeline/Options.h"
#include "pipeline/Pipeline.h"
#include "pipeline/PipelineEditor.h"
#include "science/ScienceCore.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace atomforge::pipeline;

namespace
{
int failures = 0;

void check(bool condition, const std::string& what)
{
    if (!condition) throw std::runtime_error(what);
}

void close(double a, double b, double tolerance, const std::string& what)
{
    if (!(std::abs(a - b) <= tolerance)) throw std::runtime_error(what + ": " + std::to_string(a) + " vs " + std::to_string(b));
}

void expectError(const std::function<void()>& work, const std::string& what)
{
    try { work(); } catch (const std::exception&) { return; }
    throw std::runtime_error("expected an error: " + what);
}

void test(const std::string& name, const std::function<void()>& body)
{
    try { body(); std::cout << "  ok   " << name << "\n"; }
    catch (const std::exception& error) { ++failures; std::cout << "  FAIL " << name << ": " << error.what() << "\n"; }
}

AtomSite atom(const std::string& symbol, int z, double x, double y, double w)
{
    AtomSite a;
    a.symbol = symbol; a.atomicNumber = z; a.x = x; a.y = y; a.z = w;
    return a;
}

// Conventional fcc cell(s) of copper, a = 3.615.
Structure copper(int repeat = 1)
{
    Structure s;
    const double a = 3.615;
    s.hasUnitCell = true;
    s.cellVectors = {{{a * repeat, 0, 0}, {0, a * repeat, 0}, {0, 0, a * repeat}}};
    for (int i = 0; i < repeat; ++i) for (int j = 0; j < repeat; ++j) for (int k = 0; k < repeat; ++k)
        for (const auto& f : {std::array<double, 3>{0, 0, 0}, {0, .5, .5}, {.5, 0, .5}, {.5, .5, 0}})
            s.atoms.push_back(atom("Cu", 29, (i + f[0]) * a, (j + f[1]) * a, (k + f[2]) * a));
    return s;
}

PipelineData run(const std::string& text, const Structure& input)
{
    std::vector<StageResult> stages;
    const auto out = Pipeline::parse(text).evaluate(input, &stages);
    for (const auto& stage : stages) if (!stage.error.empty()) throw std::runtime_error(stage.error);
    return out;
}

// Shortest interatomic distance, over the neighbouring periodic images.
double nearestDistance(const Structure& s)
{
    double best = std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < s.atoms.size(); ++i)
        for (std::size_t j = 0; j < s.atoms.size(); ++j)
            for (int a = -1; a <= 1; ++a)
                for (int b = -1; b <= 1; ++b)
                    for (int c = -1; c <= 1; ++c) {
                        if (i == j && !a && !b && !c) continue;
                        double d2 = 0;
                        for (int k = 0; k < 3; ++k) {
                            const double pi[3] = {s.atoms[i].x, s.atoms[i].y, s.atoms[i].z}, pj[3] = {s.atoms[j].x, s.atoms[j].y, s.atoms[j].z};
                            const double d = pj[k] + a * s.cellVectors[0][k] + b * s.cellVectors[1][k] + c * s.cellVectors[2][k] - pi[k];
                            d2 += d * d;
                        }
                        best = std::min(best, std::sqrt(d2));
                    }
    return best;
}

double evaluate(const std::string& expression, const AtomVariables& v = {}) { return Expression(expression).evaluate(v); }
}

int main()
{
    std::cout << "Structure pipeline regressions\n";

    test("selection expressions", [] {
        AtomVariables v;
        v.x = 2; v.y = -1; v.z = 4; v.fz = 0.75; v.element = "O"; v.atomicNumber = 8; v.index = 3;
        close(evaluate("1 + 2 * 3 ^ 2", v), 19, 0, "precedence");
        close(evaluate("-2 ^ 2"), -4, 0, "unary minus binds weaker than power");
        close(evaluate("(1 + 2) * 3"), 9, 0, "parentheses");
        close(evaluate("x * y + z", v), 2, 0, "variables");
        check(evaluate("fz > 0.5 && element == O", v) == 1, "element without quotes");
        check(evaluate("element == \"O\" and Z == 8", v) == 1, "quoted element and keywords");
        check(evaluate("element != 'O' or x > 100", v) == 0, "single quotes and or");
        check(evaluate("!(x > 1)", v) == 0 && evaluate("not x > 5", v) == 1, "negation");
        check(evaluate("x != 2", v) == 0, "!= is not negation");
        close(evaluate("abs(y) + sqrt(16) + max(1, 3) + min(1, 3) + floor(2.7)", v), 11, 0, "functions");
        check(evaluate("index >= 3 && index <= 3", v) == 1, "comparisons");
        check(std::isnan(evaluate("property", v)) == false, "property defaults to zero in an empty variable set");
        expectError([] { Expression("unknown(1)"); }, "unknown function");
        expectError([] { Expression("(x > 1"); }, "missing parenthesis");
        expectError([] { Expression("x > "); }, "incomplete");
        expectError([] { Expression(""); }, "empty");
        expectError([&] { Expression("element").evaluate(v); }, "text result is not a condition");
        expectError([&] { Expression("element < 3").evaluate(v); }, "ordering element names");
    });

    test("text syntax, JSON and round trips", [] {
        const Pipeline p = Pipeline::parse("replicate 2 2 1 | select-expression fz > 0.5 || element == O | delete-selected");
        check(p.modifiers.size() == 3, "three steps (|| is not a pipe)");
        check(p.modifiers[0].parameters.at("counts").items()[2].number() == 1, "positional vector");
        check(p.modifiers[1].parameters.at("expression").string() == "fz > 0.5 || element == O", "expression takes the rest of its step");
        const Pipeline keyed = Pipeline::parse("select-slab axis=c minimum=0.25 maximum=0.75 mode=add | replicate counts=3,1,1");
        check(keyed.modifiers[0].parameters.at("axis").string() == "c" && keyed.modifiers[0].parameters.at("mode").string() == "add", "key=value");
        check(keyed.modifiers[1].parameters.at("counts").items()[0].number() == 3, "comma vector");
        const Pipeline quoted = Pipeline::parse("select-expression 'element == \"O\" && x > 1' mode=intersect");
        check(quoted.modifiers[0].parameters.at("expression").string() == "element == \"O\" && x > 1", "quoted expression");
        // Text and JSON round trips.
        for (const Pipeline& original : {p, keyed, quoted}) {
            const Pipeline again = Pipeline::parse(original.toText());
            check(again.toJson().dump() == original.toJson().dump(), "text round trip of " + original.toText());
            check(Pipeline::fromJson(original.toJson()).toJson().dump() == original.toJson().dump(), "JSON round trip");
        }
        check(Pipeline::parse("").modifiers.empty(), "empty pipeline");
        expectError([] { Pipeline::parse("no-such-modifier"); }, "unknown modifier");
        expectError([] { Pipeline::parse("replicate 2 2"); }, "incomplete vector");
        expectError([] { Pipeline::parse("replicate 2.5 1 1"); }, "integer counts");
        expectError([] { Pipeline::parse("select-slab axis=w"); }, "choice values");
        expectError([] { Pipeline::parse("wrap | | wrap"); }, "empty step");
        expectError([] { Pipeline::parse("clear-selection extra"); }, "too many values");
        expectError([] { Pipeline::fromJson(Json::parse(R"({"modifiers": [{"type": "replicate", "parameters": {"bogus": 1}}]})")); }, "unknown parameter in JSON");
        check(findModifierType("replicate") && !findModifierType("nope"), "type lookup");
        for (const auto& type : modifierTypes()) {
            const Modifier m = makeModifier(type.id);  // every default parses
            check(m.parameters.size() == type.parameters.size(), std::string("defaults of ") + type.id);
        }
    });

    test("the input is never modified; disabled steps and errors", [] {
        const Structure input = copper(2);
        Pipeline p = Pipeline::parse("select-element Cu | delete-selected");
        std::vector<StageResult> stages;
        const auto out = p.evaluate(input, &stages);
        check(out.structure.atoms.empty() && input.atoms.size() == 32, "input untouched");
        check(stages[0].selected == 32 && stages[1].atoms == 0, "stage results");
        p.modifiers[1].enabled = false;
        const auto kept = p.evaluate(input, &stages);
        check(kept.structure.atoms.size() == 32 && kept.selectedCount() == 32 && stages[1].skipped, "disabled step skipped");
        // A failing step stops the pipeline; the output is the data before it.
        const Pipeline broken = Pipeline::parse("replicate 2 1 1 | select-property 0 1 | delete-selected");
        const auto partial = broken.evaluate(input, &stages);
        check(!stages[1].error.empty() && stages[2].skipped, "error stops later steps");
        check(partial.structure.atoms.size() == 64, "output before the failing step");
    });

    test("per-stage cache matches full evaluation after edits and moves", [] {
        const Structure input = copper(2);
        Pipeline p = Pipeline::parse("replicate 2 1 1 | select-random 0.25 seed=3 | random-displacement 0.1 | delete-selected");
        PipelineCache cache;
        p.evaluate(input, cache);
        // Edit the last-but-one step; only it and later steps re-run.
        p.modifiers[2].parameters["amplitude"] = 0.2;
        cache.invalidateFrom(2);
        const auto cached = p.evaluate(input, cache);
        const auto fresh = p.evaluate(input);
        check(cached.structure.atoms.size() == fresh.structure.atoms.size(), "same atoms");
        for (std::size_t i = 0; i < fresh.structure.atoms.size(); ++i) close(cached.structure.atoms[i].x, fresh.structure.atoms[i].x, 0, "same positions");
        // Moving a step invalidates from the earlier of the two positions.
        std::swap(p.modifiers[0], p.modifiers[1]);
        cache.invalidateFrom(0);
        check(p.evaluate(input, cache).structure.atoms.size() == p.evaluate(input).structure.atoms.size(), "after a move");
        check(cache.stages.size() == 4, "one stage result per step");
    });

    test("selection modifiers", [] {
        const Structure fcc = copper(2);  // 32 atoms, 7.23 A cell
        check(run("select-expression fz < 0.5", fcc).selectedCount() == 16, "half the cell by expression");
        check(run("select-slab z 0 3", fcc).selectedCount() == 16, "Cartesian slab");
        check(run("select-slab c 0.49 1", fcc).selectedCount() == 16, "fractional slab");
        // The sphere uses the minimum image: around a corner atom it reaches the other corners.
        check(run("select-sphere 0 0 0 2.6", fcc).selectedCount() == 13, "corner atom and its 12 neighbours");
        check(run("select-sphere center=0,0,0 radius=0.1 | expand-selection 2.6", fcc).selectedCount() == 13, "expand to neighbours");
        check(run("select-sphere 0 0 0 0.1 | expand-selection 2.6 steps=2", fcc).selectedCount() == 1 + 12 + 6 + 24 - 12, "second shell growth");
        const auto random = run("select-random 0.25 seed=7", fcc);
        check(random.selectedCount() == 8 && run("select-random 0.25 seed=7", fcc).selected == random.selected, "exact, reproducible random fraction");
        check(run("select-random 0.25 seed=7", fcc).selected != run("select-random 0.25 seed=8", fcc).selected, "seed changes the selection");
        check(run("select-expression fz < 0.5 | select-expression fx < 0.5 mode=add", fcc).selectedCount() == 24, "add");
        check(run("select-expression fz < 0.5 | select-expression fx < 0.5 mode=subtract", fcc).selectedCount() == 8, "subtract");
        check(run("select-expression fz < 0.5 | select-expression fx < 0.5 mode=intersect", fcc).selectedCount() == 8, "intersect");
        check(run("select-expression fz < 0.5 | invert-selection", fcc).selectedCount() == 16, "invert");
        check(run("select-element Cu | clear-selection", fcc).selectedCount() == 0, "clear");
        Structure oxide = fcc;
        oxide.atoms[0].symbol = "O"; oxide.atoms[0].atomicNumber = 8;
        oxide.atoms[5].symbol = "Ni"; oxide.atoms[5].atomicNumber = 28;
        check(run("select-element O Ni", oxide).selectedCount() == 2 && run("select-element elements=O,Ni", oxide).selectedCount() == 2, "element lists");
        expectError([&] { run("select-element Xx", oxide); }, "unknown element");
        expectError([&] { run("select-slab c 0.6 0.2", fcc); }, "inverted slab");
    });

    test("modification modifiers", [] {
        const Structure fcc = copper(2);
        const auto deleted = run("select-expression index < 4 | delete-selected", fcc);
        check(deleted.structure.atoms.size() == 28 && deleted.selectedCount() == 0, "delete");
        close(deleted.structure.atoms[0].x, fcc.atoms[4].x, 0, "remaining atoms keep their order");
        const auto alloy = run("select-random 0.5 seed=1 | assign-element Ni", fcc);
        check(std::count_if(alloy.structure.atoms.begin(), alloy.structure.atoms.end(), [](const AtomSite& a) { return a.symbol == "Ni" && a.atomicNumber == 28; }) == 16, "assign element to the selection");
        expectError([&] { run("assign-element Ni", fcc); }, "assign without a selection");
        check(run("assign-element Ni target=all", fcc).structure.atoms[31].symbol == "Ni", "assign to all");
        const auto moved = run("select-expression index == 0 | displace 1 2 3 target=selected", fcc);
        close(moved.structure.atoms[0].z, 3, 1e-12, "displace selected");
        close(moved.structure.atoms[1].z, fcc.atoms[1].z, 0, "unselected atoms stay");
        // Random displacement: Gaussian with the given standard deviation, reproducible.
        const auto noisy = run("replicate 4 4 4 | random-displacement 0.1 seed=5", fcc);
        double sum2 = 0;
        const auto big = run("replicate 4 4 4", fcc);
        for (std::size_t i = 0; i < big.structure.atoms.size(); ++i) {
            const double d = noisy.structure.atoms[i].x - big.structure.atoms[i].x;
            sum2 += d * d;
        }
        close(std::sqrt(sum2 / static_cast<double>(big.structure.atoms.size())), 0.1, 0.005, "displacement standard deviation");
        check(run("replicate 4 4 4 | random-displacement 0.1 seed=5", fcc).structure.atoms[7].y == noisy.structure.atoms[7].y, "reproducible");
        // Slice: half-space, slab, inverted, select-only.
        check(run("slice 0 0 1 3.0", fcc).structure.atoms.size() == 16, "half-space above z = 3 removed");
        check(run("slice normal=0,0,1 distance=3 invert=true", fcc).structure.atoms.size() == 16, "inverted half-space");
        check(run("slice 0 0 1 3.615 width=1", fcc).structure.atoms.size() == 24, "slab of one (001) plane removed");
        const auto sliced = run("slice 0 0 1 3.0 action=select", fcc);
        check(sliced.structure.atoms.size() == 32 && sliced.selectedCount() == 16, "slice can select instead");
    });

    test("cell modifiers", [] {
        Structure fcc = copper(1);
        fcc.atomProperty = {1, 2, 3, 4};
        fcc.atomPropertyName = "test";
        const auto big = run("replicate 2 3 1", fcc);
        check(big.structure.atoms.size() == 24, "replicated atoms");
        close(big.structure.cellVectors[1][1], 3 * 3.615, 1e-12, "replicated cell");
        check(big.structure.atomProperty.size() == 24 && big.structure.atomProperty[5] == 2, "per-atom property follows the copies");
        expectError([&] { Structure open = fcc; open.hasUnitCell = false; run("replicate 2 2 2", open); }, "replicate needs a cell");
        const auto doubled = run("transform 2 0 0 0 2 0 0 0 2", fcc);
        close(atomforge::science::cellVolume({doubled.structure.cellVectors[0], doubled.structure.cellVectors[1], doubled.structure.cellVectors[2]}), 8 * std::pow(3.615, 3), 1e-9, "transform scales the cell");
        close(doubled.structure.atoms[1].y, 3.615, 1e-12, "and the atoms");
        expectError([&] { run("transform 1 0 0 1 0 0 0 0 1", fcc); }, "singular matrix");
        const auto strained = run("strain 0.01 0 0", fcc);
        close(strained.structure.cellVectors[0][0], 3.615 * 1.01, 1e-12, "normal strain");
        const auto sheared = run("strain shear=0,0,0.02", fcc);
        close(sheared.structure.cellVectors[0][1], 3.615 * 0.01, 1e-12, "engineering shear is twice the tensor component");
        Structure outside = fcc;
        outside.atoms[0].x = -0.5; outside.atoms[1].z = 3.615 + 0.25;
        const auto wrapped = run("wrap", outside);
        close(wrapped.structure.atoms[0].x, 3.115, 1e-12, "wrapped below");
        close(wrapped.structure.atoms[1].z, 0.25, 1e-12, "wrapped above");
        const auto slab = run("add-vacuum c 10", fcc);
        close(slab.structure.cellVectors[2][2], 13.615, 1e-12, "vacuum lengthens c");
        close(slab.structure.atoms[0].z, 5, 1e-12, "atoms moved to the middle");
        const auto centred = run("add-vacuum c 10 center=false | center", fcc);
        const double top = 0.5 * 3.615, middle = 0.5 * 13.615;
        close(centred.structure.atoms[0].z, middle - top / 2, 1e-9, "center puts the atoms' extent in the middle");
        close(centred.structure.atoms[1].z, middle + top / 2, 1e-9, "both ends of the extent");
    });

    test("edit modifiers: add atom, set cell, merge", [] {
        const Structure fcc = copper(1);
        const auto added = run("add-atom O 0.5 0.5 0.5 fractional=true", fcc);
        check(added.structure.atoms.size() == 5 && added.structure.atoms[4].symbol == "O" && added.selectedCount() == 1, "atom added and selected");
        close(added.structure.atoms[4].x, 1.8075, 1e-12, "fractional position");
        const auto bigger = run("set-cell 4 0 0 0 4 0 0 0 4", fcc);
        close(bigger.structure.atoms[1].y, 2.0, 1e-12, "atoms scale with the cell");
        const auto fixed = run("set-cell a=4,0,0 b=0,4,0 c=0,0,4 scale_atoms=false", fcc);
        close(fixed.structure.atoms[1].y, 1.8075, 1e-12, "atoms stay without scaling");
        expectError([&] { run("set-cell 1 0 0 2 0 0 0 0 1", fcc); }, "coplanar cell");
        const auto file = std::filesystem::temp_directory_path() / "atomforge_merge.extxyz";
        { std::ofstream out(file); out << "2\nmerge test\nNi 0 0 0\nNi 1 0 0\n"; }
        const auto merged = run("merge file=" + file.u8string() + " offset=10,0,0", fcc);
        check(merged.structure.atoms.size() == 6 && merged.structure.atoms[5].symbol == "Ni" && merged.selectedCount() == 2, "merged and selected");
        close(merged.structure.atoms[5].x, 11, 1e-12, "offset");
        std::filesystem::remove(file);
        expectError([&] { run("merge", fcc); }, "merge needs a file");
    });

    test("insert-interstitials fills fcc voids", [] {
        const Structure fcc = copper(2);
        const auto octahedral = run("insert-interstitials H kind=octahedral count=0", fcc);
        const auto tetrahedral = run("insert-interstitials H kind=tetrahedral count=0", fcc);
        // An fcc cell has as many octahedral sites as atoms and twice as many tetrahedral sites.
        check(octahedral.structure.atoms.size() == 32 + 32, "octahedral sites: " + std::to_string(octahedral.structure.atoms.size() - 32));
        check(tetrahedral.structure.atoms.size() == 32 + 64, "tetrahedral sites: " + std::to_string(tetrahedral.structure.atoms.size() - 32));
        const auto three = run("insert-interstitials C count=3 order=random seed=4", fcc);
        check(three.structure.atoms.size() == 35 && three.selectedCount() == 3 && three.structure.atoms.back().symbol == "C", "three new, selected atoms");
        check(run("insert-interstitials C count=3 order=random seed=4", fcc).structure.atoms.back().x == three.structure.atoms.back().x, "reproducible");
    });

    test("supercell (Edit > Transform Structure)", [] {
        const Structure fcc = copper();
        const auto doubled = run("supercell \"2 0 0 0 2 0 0 0 2\"", fcc);
        check(doubled.structure.atoms.size() == 32, "2x2x2 supercell: " + std::to_string(doubled.structure.atoms.size()));
        close(doubled.structure.cellVectors[0][0], 2 * 3.615, 1e-9, "doubled a");
        // A rotated cell of twice the volume: (a+b, -a+b, c).
        const auto rotated = run("supercell matrix=1,1,0,-1,1,0,0,0,1", fcc);
        check(rotated.structure.atoms.size() == 8, "rotated supercell: " + std::to_string(rotated.structure.atoms.size()));
        close(rotated.structure.cellVectors[1][0], -3.615, 1e-9, "second vector is -a+b");
        const auto shrunk = run("supercell \"1 0 0 0 1 0 0 0 1\" | select-random 0.5 seed=2", fcc);
        check(shrunk.structure.atoms.size() == 4 && shrunk.selectedCount() == 2, "identity keeps the cell");
        expectError([&] { run("supercell \"1 0 0 1 0 0 0 0 1\"", fcc); }, "singular matrix");
        expectError([&] { run("supercell \"1.5 0 0 0 1 0 0 0 1\"", fcc); }, "non-integer matrix");
    });

    test("scientific tools as steps", [] {
        Structure fcc = copper();
        fcc.atoms[0].x += 0.1;
        const auto relaxed = run("relax steps=200 fmax=0.005", fcc);
        check(relaxed.structure.atoms.size() == 4, "relaxation keeps the atoms");
        // Relaxing restores the nearest-neighbour distance shortened by the displacement.
        const double ideal = 3.615 / std::sqrt(2.0), before = nearestDistance(fcc), after = nearestDistance(relaxed.structure);
        check(before < ideal - 0.05 && std::abs(after - ideal) < 0.5 * std::abs(before - ideal), "nearest distance " + std::to_string(before) + " -> " + std::to_string(after));
        check(findModifierType("nvt-dynamics") && findModifierType("npt-dynamics"), "dynamics steps");
#ifdef ATOMS_ENABLE_SPGLIB
        const auto primitive = run("standardize-cell cell=primitive", copper(2));
        check(primitive.structure.atoms.size() == 1, "primitive fcc cell: " + std::to_string(primitive.structure.atoms.size()));
        check(run("standardize-cell cell=conventional", copper(2)).structure.atoms.size() == 4, "conventional fcc cell");
#endif
    });

    test("command-line options of Build steps", [] {
        namespace opt = atomforge::pipeline::options;
        const opt::Reader reader("--count 4 --frac \"Cu=0.7,Zn=0.3\" --shift -0.5 --conventional --atom \"Cu 0 0 0\" --atom 'O 0.5 0.5 0.5' --axis 1,1,0");
        check(reader.integer("--count", 0) == 4 && reader.text("--frac") == "Cu=0.7,Zn=0.3", "values");
        close(reader.number("--shift", 0), -0.5, 1e-12, "negative numbers are values");
        check(reader.has("--conventional") && reader.values("--conventional").empty(), "flags without a value");
        check(reader.values("--atom").size() == 2 && reader.values("--atom")[1] == "O 0.5 0.5 0.5", "repeated flags");
        check(reader.numbers("--axis") == std::vector<double>({1, 1, 0}), "number lists");
        check(reader.number("--missing", 2.5) == 2.5, "defaults");
        opt::Writer writer;
        writer.add("--count", 4).add("--frac", "Cu=0.7,Zn=0.3").add("--atom", "Cu 0 0 0").add("--a", 3.615).flag("--conventional").add("--axis", std::vector<double>{0, 0, 1});
        check(writer.str() == "--count 4 --frac Cu=0.7,Zn=0.3 --atom \"Cu 0 0 0\" --a 3.615 --conventional --axis \"0 0 1\"", "written: " + writer.str());
        const opt::Reader back(writer.str());
        check(back.text("--atom") == "Cu 0 0 0" && back.has("--conventional"), "written options read back");
    });

    test("compute-property feeds select-property", [] {
        Structure fcc = copper(3);
        fcc.atoms.erase(fcc.atoms.begin());  // a vacancy
        const auto types = run("compute-property structure-type", fcc);
        check(types.structure.atomProperty.size() == fcc.atoms.size() && !types.structure.atomPropertyName.empty(), "structure type per atom");
        // Ackland-Jones tolerates one missing neighbour, so every atom stays fcc here.
        check(*std::min_element(types.structure.atomProperty.begin(), types.structure.atomProperty.end()) == 1, "all fcc despite the vacancy");
        const auto defects = run("compute-property coordination cutoff=3 | select-property 0 11.5", fcc);
        check(defects.selectedCount() == 12, "the vacancy's 12 neighbours have coordination 11");
        const auto core = run("compute-property centrosymmetry cutoff=4 | select-property 0.5 1000 | delete-selected", fcc);
        check(core.structure.atoms.size() == fcc.atoms.size() - 12, "delete atoms next to the vacancy by centrosymmetry");
        const auto coordination = run("compute-property coordination cutoff=3", fcc);
        check(*std::min_element(coordination.structure.atomProperty.begin(), coordination.structure.atomProperty.end()) == 11, "coordination 11 next to the vacancy");
    });

    test("pipeline editor: add, move, duplicate, disable, delete, bake", [] {
        PipelineEditor editor;
        check(!editor.active(), "inactive until it has an input");
        expectError([&] { editor.evaluate(); }, "evaluation needs an input");
        const Structure input = copper(2);
        editor.setInput(input);
        const auto same = [&](const std::string& what) {
            const auto& cached = editor.evaluate();
            const auto fresh = editor.pipeline().evaluate(input);
            check(cached.structure.atoms.size() == fresh.structure.atoms.size(), what + ": atom count");
            for (std::size_t i = 0; i < fresh.structure.atoms.size(); ++i) {
                check(cached.structure.atoms[i].symbol == fresh.structure.atoms[i].symbol, what + ": species");
                close(cached.structure.atoms[i].x, fresh.structure.atoms[i].x, 0, what + ": positions");
            }
            check(cached.selected == fresh.selected, what + ": selection");
        };
        check(editor.evaluate().structure.atoms.size() == 32, "empty pipeline passes the input through");
        const std::size_t replicate = editor.add("replicate");
        editor.setParameter(replicate, "counts", Json::array({2.0, 1.0, 1.0}));
        same("add replicate");
        check(editor.evaluate().structure.atoms.size() == 64, "replicated");
        const std::size_t sel = editor.add("select-expression");
        editor.setParameter(sel, "expression", Json(std::string("fx < 0.25")));
        editor.add("delete-selected");
        same("select and delete");
        const std::size_t afterDelete = editor.evaluate().structure.atoms.size();
        // Insert at the front, move, duplicate, disable and delete, checking against fresh runs.
        editor.add("assign-element", 0);
        editor.setParameter(0, "element", Json(std::string("Ni")));
        editor.setParameter(0, "target", Json(std::string("all")));
        same("insert at the front");
        check(editor.evaluate().structure.atoms[0].symbol == "Ni", "assigned before replication");
        editor.move(0, 3);
        same("move to the end");
        check(editor.modifier(3).type == "assign-element" && editor.modifier(0).type == "replicate", "order after the move");
        editor.duplicate(0);
        same("duplicate replicate");
        check(editor.evaluate().structure.atoms.size() == 2 * afterDelete, "duplicated replicate doubles the cell again");
        editor.setEnabled(1, false);
        same("disable the copy");
        check(editor.evaluate().structure.atoms.size() == afterDelete, "bypassed step has no effect");
        editor.remove(1);
        same("remove the copy");
        check(editor.size() == 4, "four steps");
        expectError([&] { editor.setParameter(0, "bogus", Json(1.0)); }, "unknown parameter");
        expectError([&] { editor.move(0, 9); }, "move out of range");
        expectError([&] { editor.remove(9); }, "remove out of range");
        // Change detection and the input stay separate from the output.
        const auto output = editor.evaluate().structure;
        check(editor.isCurrentOutput(output), "output recognised");
        Structure edited = output;
        edited.atoms[0].x += 0.5;
        check(!editor.isCurrentOutput(edited), "outside edits detected");
        check(editor.input().atoms.size() == 32, "input kept");
        // JSON round trip of the editor's pipeline.
        PipelineEditor copy;
        copy.fromJson(editor.toJson());
        copy.setInput(input);
        check(copy.evaluate().structure.atoms.size() == output.atoms.size(), "restored pipeline gives the same output");
        const Structure baked = editor.bake();
        check(baked.atoms.size() == output.atoms.size() && !editor.active() && editor.size() == 0, "bake ends the pipeline with its output");
    });

    if (failures) {
        std::cout << failures << " pipeline regression(s) failed\n";
        return 1;
    }
    std::cout << "All pipeline regressions passed\n";
    return 0;
}
