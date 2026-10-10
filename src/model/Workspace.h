#pragma once

#include "model/Structure.h"
#include "electronic/Volume.h"
#include <string>
#include <vector>
#include <map>

namespace atomforge
{
// Portable project data. GPU resources and running computations are never archived.
struct Workspace
{
    Structure structure;
    electronic::Volume volume, reference;
    electronic::Mesh surface;
    std::vector<double> table, camera;
    std::map<std::string,double> settings, sliceSettings;
    int columns = 0;
    std::string heading, sourcePath, referencePath, title;
    std::vector<std::string> history;
    // Scientific tools dialog state (selected tool, inputs and last result) as JSON.
    std::string science;
    // Structure pipeline (JSON) and, while it is active, its unmodified input.
    std::string pipeline;
    bool hasPipelineInput = false;
    Structure pipelineInput;
};

// Projects are written in format 3 (format 2 added per-atom properties and
// scientific tool state, format 3 the structure pipeline); formats 1 and 2 still
// load. formatVersion exists for compatibility tests.
void saveWorkspace(const std::vector<Workspace>& tabs, const std::string& path, int formatVersion = 3);
std::vector<Workspace> loadWorkspace(const std::string& path);
}
