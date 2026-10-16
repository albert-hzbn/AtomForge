#include "pipeline/ModifierSupport.h"
#include "util/ElementData.h"

#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace atomforge::pipeline::support
{
namespace
{
const ModifierType& self(const char* id) { return *findModifierType(id); }
}

double number(const Json& p, const char* id, const char* name) { return parameter(p, self(id), name).number(); }
bool flag(const Json& p, const char* id, const char* name) { return parameter(p, self(id), name).boolean(); }
std::string text(const Json& p, const char* id, const char* name) { return parameter(p, self(id), name).string(); }
Vec3 vector(const Json& p, const char* id, const char* name)
{
    const Json v = parameter(p, self(id), name);
    return {v.items()[0].number(), v.items()[1].number(), v.items()[2].number()};
}

Mat3 cellOf(const Structure& s) { Mat3 m{}; for (int r = 0; r < 3; ++r) m[r] = s.cellVectors[r]; return m; }

void requireCell(const PipelineData& data, const char* what)
{
    if (!data.structure.hasUnitCell) throw std::runtime_error(std::string(what) + " needs a structure with a unit cell");
}

void rebuild(PipelineData& data, const std::vector<std::size_t>& keep)
{
    Structure& s = data.structure;
    const std::size_t n = s.atoms.size();
    Structure out = s;
    out.atoms.clear();
    out.grainColors.clear();
    out.grainRegionIds.clear();
    out.atomProperty.clear();
    out.frozenProperty.clear();
    std::vector<char> selected;
    for (std::size_t i : keep) {
        out.atoms.push_back(s.atoms[i]);
        if (s.grainColors.size() == n) out.grainColors.push_back(s.grainColors[i]);
        if (s.grainRegionIds.size() == n) out.grainRegionIds.push_back(s.grainRegionIds[i]);
        if (s.atomProperty.size() == n) out.atomProperty.push_back(s.atomProperty[i]);
        if (s.frozenProperty.size() == n) out.frozenProperty.push_back(s.frozenProperty[i]);
        selected.push_back(i < data.selected.size() ? data.selected[i] : 0);
    }
    out.dislocationLoopPoints.clear();
    out.dislocationDetectionDone = false;
    data.structure = std::move(out);
    data.selected = std::move(selected);
}

std::vector<std::size_t> targets(const PipelineData& data, bool selectedOnly, const char* what)
{
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < data.structure.atoms.size(); ++i)
        if (!selectedOnly || data.selected[i]) out.push_back(i);
    if (selectedOnly && out.empty()) throw std::runtime_error(std::string(what) + ": no atoms are selected; add a selection modifier above");
    return out;
}

void select(PipelineData& data, const std::vector<char>& mask, const std::string& mode)
{
    for (std::size_t i = 0; i < mask.size(); ++i) {
        char& s = data.selected[i];
        if (mode == "replace") s = mask[i];
        else if (mode == "add") s = s || mask[i];
        else if (mode == "subtract") s = s && !mask[i];
        else s = s && mask[i];  // intersect
    }
    data.notes.push_back(std::to_string(data.selectedCount()) + " of " + std::to_string(mask.size()) + " atoms selected");
}

Vec3 minimumImage(const Structure& s, Vec3 d)
{
    if (!s.hasUnitCell) return d;
    const Mat3 cell = cellOf(s);
    Vec3 f = science::fractional(d, cell);
    for (std::size_t k = 0; k < 3; ++k)
        if (s.pbc[k]) f[k] -= std::round(f[k]);
    return science::rowTimes(f, cell);
}

std::string standardElement(const std::string& name)
{
    std::string lower = name;
    for (char& c : lower) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    for (int z = 1; z <= 118; ++z) {
        std::string symbol = elementSymbol(z);
        for (char& c : symbol) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (symbol == lower) return elementSymbol(z);
    }
    throw std::runtime_error("Unknown element " + name);
}

std::vector<double> evaluateForAtoms(const Expression& expression, const PipelineData& data)
{
    const Structure& s = data.structure;
    const Mat3 inv = s.hasUnitCell ? science::inverse(cellOf(s)) : Mat3{};
    std::vector<double> values(s.atoms.size());
    AtomVariables v;
    v.count = static_cast<double>(s.atoms.size());
    for (std::size_t i = 0; i < s.atoms.size(); ++i) {
        const AtomSite& a = s.atoms[i];
        v.x = a.x; v.y = a.y; v.z = a.z;
        if (s.hasUnitCell) {
            const Vec3 f = science::rowTimes(positionOf(a), inv);
            v.fx = f[0]; v.fy = f[1]; v.fz = f[2];
        }
        v.index = static_cast<double>(i);
        v.property = s.atomProperty.size() == s.atoms.size() ? s.atomProperty[i] : std::numeric_limits<double>::quiet_NaN();
        v.frozen = s.frozenProperty.size() == s.atoms.size() ? s.frozenProperty[i] : std::numeric_limits<double>::quiet_NaN();
        v.selected = data.selected[i] ? 1 : 0;
        v.atomicNumber = a.atomicNumber;
        v.fixed = static_cast<double>((a.fixed & 1) + ((a.fixed >> 1) & 1) + ((a.fixed >> 2) & 1));
        v.element = a.symbol;
        values[i] = expression.evaluate(v);
    }
    return values;
}
}
