#pragma once

#include "science/ScienceData.h"

#include <filesystem>
#include <string>

namespace atomforge::science
{
// Runs one catalog tool natively. Request keys are validated against the
// catalog; data references are resolved relative to the base directory.
ToolOutput runTool(const std::string& tool, const Json& request, const std::filesystem::path& base = ".",
                   const StructureReader& reader = {});
// Tool ids with a runner (for registry consistency checks).
std::vector<std::string> runnableTools();

// Compact human-readable summary; full arrays remain in the result JSON.
std::string resultReport(const std::string& tool, const Json& result);

// {"tool", "parameters", "result"} document, as written by the CLI and desktop.
Json resultDocument(const std::string& tool, const Json& request, const ToolOutput& output);
}
