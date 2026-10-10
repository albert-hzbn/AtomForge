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
