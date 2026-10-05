#pragma once

#include "science/ScienceData.h"

#include <filesystem>
#include <string>

namespace atomforge::science
{
// Runs one tool many times. The batch object holds
//   "tool":    catalog tool id
//   "base":    parameters shared by every run (optional)
//   "sweep":   {"param": [values...], ...}  Cartesian product of values
//   "files":   {"param": "pattern"}         one run per matching file; the
//              pattern may contain * and ? in its file name (sorted matches)
//   "collect": ["key", "key.sub", ...]      result fields tabulated per run
//   "continue_on_error": true (default)     record failures and go on
// Relative paths resolve against base. Returns {"tool", "runs": [...],
// "table": {"columns": [...], "rows": [...]}, "failures"}.
Json runBatch(const Json& batch, const std::filesystem::path& base, const StructureReader& reader = {});

// CSV text of a batch table (parameter columns then collected values).
std::string batchCsv(const Json& batchResult);

// Value at a dotted path ("a.b.0.c") in a result, or null.
Json valueAt(const Json& value, const std::string& path);
}
