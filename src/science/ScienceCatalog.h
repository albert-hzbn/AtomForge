// Desktop and CLI metadata for the native scientific tools (src/science).
// Parameter names match the JSON request keys of the atomforge.science Python API.
//
// Adding a tool: describe it in ScienceCatalog.cpp, add its runner to the
// table in ScienceTools.cpp and, optionally, plots (ResultPlots.cpp) and
// per-atom properties (AtomProperties.cpp). Registry tests check that every
// table names a catalog tool and that every catalog tool can run.
#pragma once
#include <string>
#include <vector>

struct ScienceParameterDef { const char* name; const char* label; const char* value; const char* kind; bool required; };

// How the desktop presents a tool: per-atom analyses beside the viewport,
// trajectory analyses, property calculators, plot workspaces, simulations
// and input-file generators each have their own window layout.
enum class ScienceToolView { Atoms, Trajectory, Properties, Plot, Simulation, Generator };

struct ScienceToolDef
{
    const char* id;
    const char* title;
    const char* help;      // first line: one-line description; the rest: method notes (CLI catalog)
    const char* category;  // Analysis menu section
    std::vector<ScienceParameterDef> parameters;
    ScienceToolView view;
};

const std::vector<ScienceToolDef>& scienceToolCatalog();
// The catalog entry with this id, or nullptr.
const ScienceToolDef* findScienceTool(const std::string& id);
