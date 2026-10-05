// Validates the dislocation-lines tool (src/science/DislocationLines.h) on
// dislocations inserted by the dislocation builder, whose line direction,
// Burgers vector and position are known exactly.
#include "algorithms/DislocationBuilder.h"
#include "science/ScienceTools.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace atomforge;
using science::Json;

namespace
{
int failures = 0;

void check(bool condition, const std::string& what)
{
    if (!condition) throw std::runtime_error(what);
}

void test(const std::string& name, const std::function<void()>& body)
{
    try { body(); std::cout << "  ok   " << name << '\n'; }
    catch (const std::exception& error) { ++failures; std::cout << "  FAIL " << name << ": " << error.what() << '\n'; }
}

Structure fcc(double a, int n)
{
    const double frac[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
    Structure s;
    s.hasUnitCell = true;
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) s.cellVectors[i][j] = i == j ? a * n : 0.0;
    for (int x = 0; x < n; ++x) for (int y = 0; y < n; ++y) for (int z = 0; z < n; ++z)
        for (const auto& f : frac) {
            AtomSite atom;
            atom.symbol = "Cu"; atom.atomicNumber = 29;
            atom.x = (x + f[0]) * a; atom.y = (y + f[1]) * a; atom.z = (z + f[2]) * a;
            s.atoms.push_back(atom);
        }
    return s;
}

std::string describe(const Json& result)
{
    std::string text;
    for (const auto& line : result.at("lines").items())
        text += " [" + std::to_string(static_cast<int>(line.at("atoms").number())) + " atoms at " + line.at("centroid_A").dump() +
                " b=" + line.at("burgers_vector_A").dump() + " t=" + line.at("direction").dump() + "]";
    return text;
}

Json lines(const Structure& reference, const Structure& deformed)
{
    // Cylinder cuts are finite clusters: analyse without periodicity.
    Structure ref = reference, def = deformed;
    ref.hasUnitCell = def.hasUnitCell = false;
    Json request = Json::object();
    request["reference"] = science::structureJson(ref);
    request["structure"] = science::structureJson(def);
    request["reference"]["cell"] = Json();
    request["structure"]["cell"] = Json();
    return science::runTool("dislocation-lines", request).result;
}

glm::dvec3 vec(const Json& value)
{
    return {value.items()[0].number(), value.items()[1].number(), value.items()[2].number()};
}

DislocationResult build(DislocationCharacter character, bool dipole)
{
    DislocationParams params;
    params.character = character;
    params.shape = DislocationShape::Cylinder;
    // Apply the elastic field to the whole crystal: a finite cylinder cut
    // would leave a displacement discontinuity (real Nye content) at its surface.
    params.cylinderRadius = 1.0e6f;
    params.coreRadius = 1.2f;
    params.cutoffRadius = 0.0f;
    params.dipole = dipole;
    params.dipoleOffset = glm::vec2(16.0f, 0.0f);
    const auto result = buildDislocation(fcc(3.61, dipole ? 14 : 11), params);
    check(result.success, "builder: " + result.message);
    return result;
}
}

int main()
{
    std::cout << "Dislocation line extraction regressions\n";
    const Structure reference13 = fcc(3.61, 11), reference16 = fcc(3.61, 14);

    test("edge dislocation: line, Burgers vector and character", [&] {
        const auto edge = build(DislocationCharacter::Edge, false);
        const Json result = lines(reference13, edge.output);
        check(result.at("line_count").number() == 1, "one line, got" + describe(result));
        const Json& line = result.at("lines").items()[0];
        const glm::dvec3 t = vec(line.at("direction")), b = vec(line.at("burgers_vector_A"));
        check(std::abs(glm::dot(t, glm::normalize(glm::dvec3(edge.lineDirection)))) > 0.97, "line direction");
        check(std::abs(glm::dot(glm::normalize(b), glm::normalize(glm::dvec3(edge.burgersDirection)))) > 0.9, "Burgers direction");
        // The circuit closure is an exact lattice vector.
        check(std::abs(glm::length(b) - edge.burgersMagnitude) < 1e-3, "Burgers magnitude equals a/sqrt(2) (got " + std::to_string(glm::length(b)) + ")");
        check(line.at("character_angle_deg").number() > 70, "edge character");
        const glm::dvec3 offset = vec(line.at("centroid_A")) - glm::dvec3(edge.linePoint);
        const double distance = glm::length(offset - glm::dot(offset, glm::dvec3(edge.lineDirection)) * glm::dvec3(edge.lineDirection));
        check(distance < 3.0, "line position within 3 A of the inserted core (" + std::to_string(distance) + ")");
    });

    test("screw dislocation character", [&] {
        const auto screw = build(DislocationCharacter::Screw, false);
        const Json result = lines(reference13, screw.output);
        check(result.at("line_count").number() == 1, "one screw line");
        const Json& line = result.at("lines").items()[0];
        check(line.at("character_angle_deg").number() < 20, "screw character (angle " + std::to_string(line.at("character_angle_deg").number()) + ")");
        check(std::abs(line.at("burgers_magnitude_A").number() - screw.burgersMagnitude) < 1e-3, "screw Burgers magnitude");
    });

    test("dipole: two lines with opposite Burgers vectors", [&] {
        const auto dipole = build(DislocationCharacter::Edge, true);
        const Json result = lines(reference16, dipole.output);
        check(result.at("line_count").number() == 2, "two lines, got" + describe(result));
        const Json& a = result.at("lines").items()[0];
        const Json& b = result.at("lines").items()[1];
        const glm::dvec3 ba = vec(a.at("burgers_vector_A")), bb = vec(b.at("burgers_vector_A"));
        check(glm::length(ba + bb) < 1e-3, "exactly opposite Burgers vectors");
        const glm::dvec3 d = vec(a.at("centroid_A")) - vec(b.at("centroid_A"));
        const glm::dvec3 t(dipole.lineDirection);
        const double separation = glm::length(d - glm::dot(d, t) * t);
        check(std::abs(separation - 16.0) < 3.0, "separation matches the dipole offset (" + std::to_string(separation) + " A)");
    });

    test("perfect crystal has no dislocation cores", [&] {
        bool rejected = false;
        try { lines(reference13, reference13); } catch (const std::exception&) { rejected = true; }
        check(rejected, "zero Nye tensor rejected");
    });

    if (failures) { std::cout << failures << " dislocation line regression(s) failed\n"; return 1; }
    std::cout << "All dislocation line regressions passed\n";
    return 0;
}
