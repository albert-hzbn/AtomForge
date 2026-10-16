#pragma once

#include "science/AtomProperties.h"

#include "science/ScienceData.h"

#include <filesystem>
#include <string>

namespace atomforge::science
{
// Runs one catalog tool natively. Request keys are validated against the
// catalog; data references are resolved relative to the base directory.
ToolOutput runTool(const std::string& tool, const Json& request, const std::filesystem::path& base = ".",
                   const StructureReader& reader = {});
// Periodic axes of a structure: none without a cell.
Pbc periodicAxes(const Structure& structure);

// Positions, cell and periodic axes of a structure as tool parameters.
Json atomsRequest(const Structure& structure, const Pbc& pbc);

// The request of a per-atom analysis of a structure: its atoms, the species
// for tools that take them, and the neighbour cutoff only for tools that need
// one; tools with their own default or automatic range choose it themselves.
Json perAtomRequest(const std::string& tool, const Structure& structure, const Pbc& pbc, double cutoffA);

// Per-atom analysis of one structure through a catalog tool that reports
// per-atom properties (structure-type, centrosymmetry, bond-order, ...); returns
// the property named `property` (the first one when empty).
AtomProperty analysePerAtom(const Structure& structure, const std::string& tool, double cutoffA, const std::string& property = "");

// Tool ids with a runner (for registry consistency checks).
std::vector<std::string> runnableTools();

// Compact human-readable summary; full arrays remain in the result JSON.
std::string resultReport(const std::string& tool, const Json& result);

// {"tool", "parameters", "result"} document, as written by the CLI and desktop.
Json resultDocument(const std::string& tool, const Json& request, const ToolOutput& output);
}
