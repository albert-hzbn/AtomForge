#include "cli/CliArgs.h"

#include <cctype>
#include <cmath>
#include <cstring>
#include <sstream>
#include <stdexcept>

namespace cli
{
const char* findArg(int argc, char* argv[], const char* flag)
{
    for (int i = 1; i < argc - 1; ++i)
        if (std::strcmp(argv[i], flag) == 0)
            return argv[i + 1];
    return nullptr;
}

bool hasFlag(int argc, char* argv[], const char* flag)
{
    for (int i = 1; i < argc; ++i)
        if (std::strcmp(argv[i], flag) == 0)
            return true;
    return false;
}

// Collect all values after repeated occurrences of `flag`.
std::vector<std::string> findAllArgs(int argc, char* argv[], const char* flag)
{
    std::vector<std::string> out;
    for (int i = 1; i < argc - 1; ++i)
        if (std::strcmp(argv[i], flag) == 0)
            out.emplace_back(argv[i + 1]);
    return out;
}

double argDouble(int argc, char* argv[], const char* flag, double def)
{
    const char* v = findArg(argc, argv, flag);
    if (!v)
        return def;

    try
    {
        std::size_t parsed = 0;
        const std::string text(v);
        const double value = std::stod(text, &parsed);
        if (parsed != text.size() || !std::isfinite(value))
            throw std::invalid_argument("not a finite number");
        return value;
    }
    catch (const std::exception&)
    {
        throw std::invalid_argument(std::string("invalid numeric value '") + v + "' for " + flag);
    }
}

int argInt(int argc, char* argv[], const char* flag, int def)
{
    const char* v = findArg(argc, argv, flag);
    if (!v)
        return def;

    try
    {
        std::size_t parsed = 0;
        const std::string text(v);
        const int value = std::stoi(text, &parsed);
        if (parsed != text.size())
            throw std::invalid_argument("trailing characters");
        return value;
    }
    catch (const std::exception&)
    {
        throw std::invalid_argument(std::string("invalid integer value '") + v + "' for " + flag);
    }
}

// Detect Open Babel format string from output filename extension.
std::string detectFormat(const std::string& path)
{
    const auto dot = path.rfind('.');
    if (dot == std::string::npos)
        return "cif";
    std::string ext = path.substr(dot + 1);
    // lower-case
    for (char& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));

    if (ext == "cif")              return "cif";
    if (ext == "xyz")              return "xyz";
    if (ext == "vasp" || ext == "poscar" || ext == "contcar") return "vasp";
    if (ext == "lmp"  || ext == "lammps") return "lammps";
    if (ext == "pdb")              return "pdb";
    if (ext == "mol2")             return "mol2";
    if (ext == "extxyz")           return "extxyz";
    return ext; // pass-through for anything else
}

bool parseIvec3(const char* text, glm::ivec3& out)
{
    if (!text)
        return false;
    std::istringstream iss(text);
    int x = 0, y = 0, z = 0;
    if (!(iss >> x >> y >> z))
        return false;
    out = glm::ivec3(x, y, z);
    return true;
}

bool parseVec3(const char* text, glm::vec3& out)
{
    if (!text)
        return false;
    std::istringstream iss(text);
    float x = 0.0f, y = 0.0f, z = 0.0f;
    if (!(iss >> x >> y >> z))
        return false;
    out = glm::vec3(x, y, z);
    return true;
}
}
