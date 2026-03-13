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

void getDefaultElementColor(int atomicNumber, float& r, float& g, float& b);

Structure loadStructure(const std::string& filename);