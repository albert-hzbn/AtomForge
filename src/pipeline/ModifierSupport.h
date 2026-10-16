#pragma once
// Helpers shared by the modifier implementations (Modifiers.cpp,
// AnalysisSteps.cpp, PythonScriptStep.cpp): parameter access, atom
// bookkeeping and selections. Not part of the public pipeline interface.

#include "pipeline/Expression.h"
#include "pipeline/Pipeline.h"
#include "science/ScienceCore.h"

#include <string>
#include <vector>

namespace atomforge::pipeline::support
{
using science::Mat3;
using science::Vec3;

// ---------------------------------------------------------------- parameters
// A parameter of the step `id` (its value, or its default).
double number(const Json& p, const char* id, const char* name);
bool flag(const Json& p, const char* id, const char* name);
std::string text(const Json& p, const char* id, const char* name);
Vec3 vector(const Json& p, const char* id, const char* name);

// How a selection step combines with the current selection.
inline const char* const kModes = "replace|add|subtract|intersect";
inline const ModifierParameter kModeParameter = {"mode", "Combine with the current selection", "choice", "replace", kModes};

// ---------------------------------------------------------------- atoms
Mat3 cellOf(const Structure& s);
inline Vec3 positionOf(const AtomSite& a) { return {a.x, a.y, a.z}; }
inline void place(AtomSite& a, const Vec3& x) { a.x = x[0]; a.y = x[1]; a.z = x[2]; }
void requireCell(const PipelineData& data, const char* what);

// Keeps atoms by index (in the given order), carrying aligned per-atom metadata.
void rebuild(PipelineData& data, const std::vector<std::size_t>& keep);

// Indices the operation applies to: the selection, or every atom.
std::vector<std::size_t> targets(const PipelineData& data, bool selectedOnly, const char* what);

// Combines a new selection with the current one by `mode`.
void select(PipelineData& data, const std::vector<char>& mask, const std::string& mode);

// The shortest image of a vector along the periodic cell vectors.
Vec3 minimumImage(const Structure& s, Vec3 d);

// An element symbol in any case ("cu", "CU" -> "Cu"); throws for an unknown one.
std::string standardElement(const std::string& name);

// The value of an expression for every atom (variables as in select-expression).
std::vector<double> evaluateForAtoms(const Expression& expression, const PipelineData& data);
}
