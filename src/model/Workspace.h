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
};

void saveWorkspace(const std::vector<Workspace>& tabs, const std::string& path);
std::vector<Workspace> loadWorkspace(const std::string& path);
}
