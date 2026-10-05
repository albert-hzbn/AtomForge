#pragma once

#include "model/Structure.h"
#include "science/ScienceCore.h"

#include <filesystem>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace atomforge::science
{
// Optional reader for structure formats beyond the built-in XYZ/extXYZ,
// VASP POSCAR/XDATCAR, LAMMPS dump and JSON readers (the desktop passes Open Babel).
using StructureReader = std::function<Structure(const std::filesystem::path&)>;

// One frame of a structure or trajectory file, with optional dynamics data.
struct FrameData
{
    Structure structure;
    Pbc pbc{};
    std::vector<Vec3> velocities;  // Angstrom/fs; empty when absent
    std::vector<double> masses;    // amu; empty when absent
    std::optional<double> timeFs;
};

std::vector<FrameData> readFrames(const std::filesystem::path& path, const StructureReader& reader = {});

// Standard atomic masses (amu) for chemical symbols.
double atomicMass(const std::string& symbol);
int atomicNumber(const std::string& symbol);

// A structure input together with its periodicity.
struct StructureInput
{
    Structure structure;
    Pbc pbc{};
};

// Resolved tool parameters. Data references ({"file": ...}) are loaded when
// constructed, relative to the base directory; files are never evaluated.
class Parameters
{
public:
    Parameters() = default;
    Parameters(const std::string& tool, const Json& request, const std::filesystem::path& base,
               const StructureReader& reader = {});

    bool has(const std::string& name) const;
    const Json& json(const std::string& name) const;
    NdArray array(const std::string& name) const;
    double number(const std::string& name) const;
    double number(const std::string& name, double fallback) const;
    bool boolean(const std::string& name, bool fallback) const;
    Pbc pbc(const std::string& name, Pbc fallback) const;
    const StructureInput& structure(const std::string& name) const;

private:
    std::map<std::string, Json> m_json;
    std::map<std::string, NdArray> m_arrays;
    std::map<std::string, StructureInput> m_structures;
};

// A tool's numerical result plus any structural frames it produced.
struct ToolOutput
{
    Json result;
    std::vector<Structure> frames;               // images, trajectory frames or a primitive cell
    std::vector<std::vector<Vec3>> velocities;   // Angstrom/fs per frame (dynamics only)
    std::vector<double> times;                   // fs per frame (dynamics only)
};

// Writes frames as extXYZ (ASE-compatible momenta and time_fs metadata).
void writeExtxyz(const std::filesystem::path& path, const std::vector<Structure>& frames,
                 const std::vector<std::vector<Vec3>>& velocities, const std::vector<double>& times);

Json structureJson(const Structure& structure);
}
