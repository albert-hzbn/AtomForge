#pragma once
#include <vector>
#include <string>

struct AtomSite
{
    std::string symbol;
    int atomicNumber;

    double x;
    double y;
    double z;

    float r;
    float g;
    float b;
};

struct Structure
{
    std::vector<AtomSite> atoms;
};

Structure loadStructure(const std::string& filename);