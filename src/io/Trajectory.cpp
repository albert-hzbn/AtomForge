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

std::vector<Structure> loadLammpsTrajectory(const std::string& path)
{
    std::ifstream input(std::filesystem::u8path(path));
    if (!input) throw std::runtime_error("Cannot open trajectory");
    std::vector<Structure> frames;
    std::string line;
    std::size_t totalAtoms=0;
    while (std::getline(input,line)) {
        atomforge::taskCheckpoint();
        if (line.find_first_not_of(" \t\r")==std::string::npos) continue;
        if (line.rfind("ITEM: TIMESTEP",0)!=0) throw std::runtime_error("Expected 'ITEM: TIMESTEP', got: "+line);
        if (!std::getline(input,line)) throw std::runtime_error("Truncated LAMMPS dump: missing timestep value");
        if (!std::getline(input,line) || line.rfind("ITEM: NUMBER OF ATOMS",0)!=0)
            throw std::runtime_error("Expected 'ITEM: NUMBER OF ATOMS'");
        if (!std::getline(input,line)) throw std::runtime_error("Truncated LAMMPS dump: missing atom count");
        std::istringstream countLine(line);
        long long count=0;
        if (!(countLine>>count) || count<1 || count>1000000) throw std::runtime_error("Invalid LAMMPS atom count");
        totalAtoms+=static_cast<std::size_t>(count);
        if (totalAtoms>10000000 || frames.size()>=100000) throw std::runtime_error("Trajectory exceeds memory limit");
        if (!std::getline(input,line) || line.rfind("ITEM: BOX BOUNDS",0)!=0)
            throw std::runtime_error("Expected 'ITEM: BOX BOUNDS'");
        const bool triclinic = line.find("xy xz yz")!=std::string::npos;
        double lo[3]{}, hi[3]{}, tilt[3]{};
        for (int axis=0;axis<3;++axis) {
            if (!std::getline(input,line)) throw std::runtime_error("Truncated LAMMPS box bounds");
            std::istringstream row(line);
            const bool ok = triclinic ? static_cast<bool>(row>>lo[axis]>>hi[axis]>>tilt[axis])
                                       : static_cast<bool>(row>>lo[axis]>>hi[axis]);
            if (!ok || !std::isfinite(lo[axis]) || !std::isfinite(hi[axis]) || !std::isfinite(tilt[axis]))
                throw std::runtime_error("Invalid LAMMPS box bounds");
        }
        if (!std::getline(input,line) || line.rfind("ITEM: ATOMS ",0)!=0)
            throw std::runtime_error("Expected 'ITEM: ATOMS <columns>'");
        std::vector<std::string> columns;
        { std::istringstream header(line.substr(std::string("ITEM: ATOMS ").size())); for (std::string token; header>>token;) columns.push_back(token); }
        int elementCol=-1, xCol=-1, yCol=-1, zCol=-1;
        for (std::size_t i=0;i<columns.size();++i) {
            if (columns[i]=="element") elementCol=static_cast<int>(i);
            else if (columns[i]=="x" || columns[i]=="xu") xCol=static_cast<int>(i);
            else if (columns[i]=="y" || columns[i]=="yu") yCol=static_cast<int>(i);
            else if (columns[i]=="z" || columns[i]=="zu") zCol=static_cast<int>(i);
        }
        if (elementCol<0)
            throw std::runtime_error("LAMMPS dump has no 'element' column; add 'dump_modify <id> element <symbols...>' "
                                      "to your LAMMPS input so atom types map to chemical symbols");
        if (xCol<0 || yCol<0 || zCol<0) throw std::runtime_error("LAMMPS dump is missing x/y/z (or xu/yu/zu) position columns");
        Structure frame;
        // Convert LAMMPS's bounding-box convention to true (tilted) cell
        // edges: LAMMPS manual, "Triclinic simulation boxes".
        const double xyT=triclinic?tilt[0]:0, xzT=triclinic?tilt[1]:0, yzT=triclinic?tilt[2]:0;
        const double xloB=lo[0]-std::min({0.0,xyT,xzT,xyT+xzT}), xhiB=hi[0]-std::max({0.0,xyT,xzT,xyT+xzT});
        const double yloB=lo[1]-std::min(0.0,yzT), yhiB=hi[1]-std::max(0.0,yzT);
        frame.cellOffset={xloB,yloB,lo[2]};
        frame.cellVectors={{ {xhiB-xloB,0,0}, {xyT,yhiB-yloB,0}, {xzT,yzT,hi[2]-lo[2]} }};
        frame.hasUnitCell=true;
        const int maxColumn=std::max({elementCol,xCol,yCol,zCol});
        for (long long i=0;i<count;++i) {
            atomforge::taskCheckpoint();
            if (!std::getline(input,line)) throw std::runtime_error("Truncated LAMMPS trajectory");
            std::istringstream row(line);
            std::vector<std::string> values;
            for (std::string token; row>>token;) values.push_back(token);
            if (static_cast<int>(values.size())<=maxColumn) throw std::runtime_error("LAMMPS atom row has too few columns");
            AtomSite atom;
            atom.symbol=values[static_cast<std::size_t>(elementCol)];
            if (!(std::istringstream(values[static_cast<std::size_t>(xCol)])>>atom.x) || !std::isfinite(atom.x))
                throw std::runtime_error("Invalid LAMMPS x coordinate");
            if (!(std::istringstream(values[static_cast<std::size_t>(yCol)])>>atom.y) || !std::isfinite(atom.y))
                throw std::runtime_error("Invalid LAMMPS y coordinate");
            if (!(std::istringstream(values[static_cast<std::size_t>(zCol)])>>atom.z) || !std::isfinite(atom.z))
                throw std::runtime_error("Invalid LAMMPS z coordinate");
            for (int z=1;z<=118;++z) if (atom.symbol==elementSymbol(z)) { atom.atomicNumber=z; break; }
            if (!atom.atomicNumber) throw std::runtime_error("Unknown LAMMPS element: "+atom.symbol);
            getDefaultElementColor(atom.atomicNumber,atom.r,atom.g,atom.b);
            frame.atoms.push_back(atom);
        }
        frames.push_back(std::move(frame));
    }
    if (frames.empty()) throw std::runtime_error("Trajectory contains no frames");
    return frames;
}

std::vector<Structure> loadTrajectory(const std::string& path)
{
    std::ifstream probe(std::filesystem::u8path(path));
    if (!probe) throw std::runtime_error("Cannot open trajectory");
    std::string line;
    while (std::getline(probe,line))
        if (line.find_first_not_of(" \t\r")!=std::string::npos) break;
    if (line.rfind("ITEM:",0)==0) return loadLammpsTrajectory(path);
    return loadXyzTrajectory(path);
}
