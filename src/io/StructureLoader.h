#pragma once
#include <vector>
#include <string>

struct AtomSite
{
    std::string symbol;
    double x,y,z;
};

struct Structure
{
    std::vector<AtomSite> atoms;
};

Structure loadStructure(const std::string& filename);