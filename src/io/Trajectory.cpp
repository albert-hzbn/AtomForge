#include "io/Trajectory.h"
#include "util/ElementData.h"
#include "util/TaskControl.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

std::vector<Structure> loadXyzTrajectory(const std::string& path)
{
    std::ifstream input(std::filesystem::u8path(path));
    if (!input) throw std::runtime_error("Cannot open trajectory");
    std::vector<Structure> frames;
    std::string line;
    std::size_t total=0;
    while (std::getline(input,line)) {
        atomforge::taskCheckpoint();
        if (line.find_first_not_of(" \t\r")==std::string::npos) continue;
        std::istringstream countLine(line);
        int count=0; std::string trailing;
        if (!(countLine>>count) || count<1 || count>1000000 || countLine>>trailing)
            throw std::runtime_error("Invalid XYZ frame atom count");
        total+=count;
        if (total>10000000 || frames.size()>=100000) throw std::runtime_error("Trajectory exceeds memory limit");
        if (!std::getline(input,line)) throw std::runtime_error("Missing XYZ frame comment");
        Structure frame;
        const auto lattice=line.find("Lattice=\"");
        if (lattice!=std::string::npos) {
            std::istringstream cell(line.substr(lattice+9));
            for (auto& vector:frame.cellVectors) for (auto& value:vector)
                if (!(cell>>value) || !std::isfinite(value)) throw std::runtime_error("Invalid XYZ lattice");
            frame.hasUnitCell=true;
        }
        const auto properties=line.find("Properties=");
        const std::string expected="Properties=species:S:1:pos:R:3";
        if (properties!=std::string::npos && line.compare(properties,expected.size(),expected)!=0)
            throw std::runtime_error("Desktop trajectory expects species:S:1:pos:R:3 first; convert other layouts using Python");
        for (int i=0;i<count;++i) {
            atomforge::taskCheckpoint();
            if (!std::getline(input,line)) throw std::runtime_error("Truncated XYZ trajectory");
            std::istringstream row(line); AtomSite atom;
            if (!(row>>atom.symbol>>atom.x>>atom.y>>atom.z) || !std::isfinite(atom.x) ||
                !std::isfinite(atom.y) || !std::isfinite(atom.z)) throw std::runtime_error("Invalid XYZ atom row");
            for (int z=1;z<=118;++z)
                if (atom.symbol==elementSymbol(z)) { atom.atomicNumber=z; break; }
            if (!atom.atomicNumber) throw std::runtime_error("Unknown XYZ element: "+atom.symbol);
            getDefaultElementColor(atom.atomicNumber,atom.r,atom.g,atom.b);
            frame.atoms.push_back(atom);
        }
        frames.push_back(std::move(frame));
    }
    if (frames.empty()) throw std::runtime_error("Trajectory contains no frames");
    return frames;
}
