#include "io/Trajectory.h"
#include "util/ElementData.h"
#include "util/TaskControl.h"
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <cctype>
#include <mutex>
#include <stdexcept>

namespace
{
// An element symbol in standard case ("CU", "cu" -> "Cu"); labels that are
// not purely letters stay as written.
std::string standardSymbol(std::string symbol)
{
    bool letters = !symbol.empty();
    for (char c : symbol) letters = letters && std::isalpha(static_cast<unsigned char>(c));
    if (!letters) return symbol;
    symbol[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(symbol[0])));
    for (std::size_t k = 1; k < symbol.size(); ++k) symbol[k] = static_cast<char>(std::tolower(static_cast<unsigned char>(symbol[k])));
    return symbol;
}
}

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
        // extXYZ Properties=name:type:count... ; locate the species and position columns.
        std::size_t speciesColumn=0, positionColumn=1, columnCount=4;
        const auto properties=line.find("Properties=");
        if (properties!=std::string::npos) {
            std::string spec=line.substr(properties+11);
            spec=spec.substr(0,spec.find_first_of(" \t\r"));
            std::istringstream fields(spec);
            std::vector<std::string> parts;
            for (std::string part; std::getline(fields,part,':');) parts.push_back(part);
            if (parts.size()%3) throw std::runtime_error("Invalid extXYZ Properties");
            bool hasSpecies=false, hasPositions=false;
            columnCount=0;
            for (std::size_t p=0;p<parts.size();p+=3) {
                const int width=std::atoi(parts[p+2].c_str());
                if (width<1) throw std::runtime_error("Invalid extXYZ Properties");
                if (parts[p]=="species" && width==1) { speciesColumn=columnCount; hasSpecies=true; }
                if ((parts[p]=="pos" || parts[p]=="positions") && width==3) { positionColumn=columnCount; hasPositions=true; }
                columnCount+=static_cast<std::size_t>(width);
            }
            if (!hasSpecies || !hasPositions) throw std::runtime_error("extXYZ frames need species and pos columns");
        }
        for (int i=0;i<count;++i) {
            atomforge::taskCheckpoint();
            if (!std::getline(input,line)) throw std::runtime_error("Truncated XYZ trajectory");
            std::istringstream row(line); AtomSite atom;
            std::vector<std::string> tokens;
            for (std::string token; row>>token;) tokens.push_back(token);
            if (tokens.size()<columnCount) throw std::runtime_error("Invalid XYZ atom row");
            atom.symbol=standardSymbol(tokens[speciesColumn]);
            try {
                atom.x=std::stod(tokens[positionColumn]);
                atom.y=std::stod(tokens[positionColumn+1]);
                atom.z=std::stod(tokens[positionColumn+2]);
            } catch (const std::exception&) { throw std::runtime_error("Invalid XYZ atom row"); }
            if (!std::isfinite(atom.x) || !std::isfinite(atom.y) || !std::isfinite(atom.z)) throw std::runtime_error("Invalid XYZ atom row");
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

namespace
{
std::mutex g_typeMutex;
std::vector<std::string> g_typeElements;
std::string g_typeText;
}

void setLammpsTypeElements(const std::string& symbols)
{
    std::string text = symbols;
    for (char& c : text) if (c == ',' || c == ';') c = ' ';
    std::istringstream words(text);
    std::vector<std::string> elements;
    for (std::string word; words >> word;) {
        const std::string symbol = standardSymbol(word);
        int z = 0;
        for (int k = 1; k <= 118 && !z; ++k) if (symbol == elementSymbol(k)) z = k;
        if (!z) throw std::runtime_error("Unknown element in the LAMMPS types: " + word);
        elements.push_back(symbol);
    }
    std::lock_guard<std::mutex> lock(g_typeMutex);
    g_typeElements = elements;
    g_typeText = symbols;
}

std::string lammpsTypeElements()
{
    std::lock_guard<std::mutex> lock(g_typeMutex);
    return g_typeText;
}

std::vector<Structure> loadLammpsTrajectory(const std::string& path)
{
    std::ifstream input(std::filesystem::u8path(path));
    if (!input) throw std::runtime_error("Cannot open trajectory");
    return loadLammpsFrames(input, static_cast<std::size_t>(-1));
}

std::vector<Structure> loadLammpsFrames(std::istream& input, std::size_t maxFrames)
{
    std::vector<Structure> frames;
    std::string line;
    std::size_t totalAtoms=0;
    while (frames.size()<maxFrames && std::getline(input,line)) {
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
        // Cartesian (x, xu) or scaled (xs, xsu: fractions of the box, as
        // LAMMPS's default "dump atom" writes) positions; atoms ordered by id.
        int elementCol=-1, xCol=-1, yCol=-1, zCol=-1, idCol=-1, typeCol=-1;
        bool scaled=false;
        for (std::size_t i=0;i<columns.size();++i) {
            const std::string& c=columns[i];
            if (c=="element") elementCol=static_cast<int>(i);
            else if (c=="id") idCol=static_cast<int>(i);
            else if (c=="type") typeCol=static_cast<int>(i);
            else if (c=="x" || c=="xu") xCol=static_cast<int>(i);
            else if (c=="y" || c=="yu") yCol=static_cast<int>(i);
            else if (c=="z" || c=="zu") zCol=static_cast<int>(i);
            else if (c=="xs" || c=="xsu") { xCol=static_cast<int>(i); scaled=true; }
            else if (c=="ys" || c=="ysu") yCol=static_cast<int>(i);
            else if (c=="zs" || c=="zsu") zCol=static_cast<int>(i);
        }
        std::vector<std::string> typeElements;
        { std::lock_guard<std::mutex> lock(g_typeMutex); typeElements=g_typeElements; }
        if (elementCol<0 && (typeCol<0 || typeElements.empty()))
            throw std::runtime_error("LAMMPS dump has no 'element' column: name the elements of types 1, 2, ... "
                                      "(LAMMPS types in the Trajectory window, --types on the command line), "
                                      "or add 'dump_modify <id> element <symbols...>' to the LAMMPS input");
        if (xCol<0 || yCol<0 || zCol<0) throw std::runtime_error("LAMMPS dump is missing x/y/z (or xu, xs, xsu) position columns");
        Structure frame;
        // Convert LAMMPS's bounding-box convention to true (tilted) cell
        // edges: LAMMPS manual, "Triclinic simulation boxes".
        const double xyT=triclinic?tilt[0]:0, xzT=triclinic?tilt[1]:0, yzT=triclinic?tilt[2]:0;
        const double xloB=lo[0]-std::min({0.0,xyT,xzT,xyT+xzT}), xhiB=hi[0]-std::max({0.0,xyT,xzT,xyT+xzT});
        const double yloB=lo[1]-std::min(0.0,yzT), yhiB=hi[1]-std::max(0.0,yzT);
        frame.cellOffset={xloB,yloB,lo[2]};
        frame.cellVectors={{ {xhiB-xloB,0,0}, {xyT,yhiB-yloB,0}, {xzT,yzT,hi[2]-lo[2]} }};
        frame.hasUnitCell=true;
        const int maxColumn=std::max({elementCol,xCol,yCol,zCol,idCol,typeCol});
        std::vector<std::pair<long long, AtomSite>> rows;
        rows.reserve(static_cast<std::size_t>(count));
        for (long long i=0;i<count;++i) {
            atomforge::taskCheckpoint();
            if (!std::getline(input,line)) throw std::runtime_error("Truncated LAMMPS trajectory");
            std::istringstream row(line);
            std::vector<std::string> values;
            for (std::string token; row>>token;) values.push_back(token);
            if (static_cast<int>(values.size())<=maxColumn) throw std::runtime_error("LAMMPS atom row has too few columns");
            AtomSite atom;
            if (elementCol>=0) atom.symbol=standardSymbol(values[static_cast<std::size_t>(elementCol)]);
            else {
                long type=0;
                if (!(std::istringstream(values[static_cast<std::size_t>(typeCol)])>>type) || type<1 ||
                    type>static_cast<long>(typeElements.size()))
                    throw std::runtime_error("LAMMPS atom type "+values[static_cast<std::size_t>(typeCol)]+" has no element in the LAMMPS types ("+
                                             std::to_string(typeElements.size())+" given)");
                atom.symbol=typeElements[static_cast<std::size_t>(type-1)];
            }
            if (!(std::istringstream(values[static_cast<std::size_t>(xCol)])>>atom.x) || !std::isfinite(atom.x))
                throw std::runtime_error("Invalid LAMMPS x coordinate");
            if (!(std::istringstream(values[static_cast<std::size_t>(yCol)])>>atom.y) || !std::isfinite(atom.y))
                throw std::runtime_error("Invalid LAMMPS y coordinate");
            if (!(std::istringstream(values[static_cast<std::size_t>(zCol)])>>atom.z) || !std::isfinite(atom.z))
                throw std::runtime_error("Invalid LAMMPS z coordinate");
            for (int z=1;z<=118;++z) if (atom.symbol==elementSymbol(z)) { atom.atomicNumber=z; break; }
            if (!atom.atomicNumber) throw std::runtime_error("Unknown LAMMPS element: "+atom.symbol);
            getDefaultElementColor(atom.atomicNumber,atom.r,atom.g,atom.b);
            if (scaled) {
                const auto& cv=frame.cellVectors;
                const double fa=atom.x, fb=atom.y, fc=atom.z;
                atom.x=frame.cellOffset[0]+fa*cv[0][0]+fb*cv[1][0]+fc*cv[2][0];
                atom.y=frame.cellOffset[1]+fa*cv[0][1]+fb*cv[1][1]+fc*cv[2][1];
                atom.z=frame.cellOffset[2]+fa*cv[0][2]+fb*cv[1][2]+fc*cv[2][2];
            }
            long long id=i;
            if (idCol>=0 && !(std::istringstream(values[static_cast<std::size_t>(idCol)])>>id))
                throw std::runtime_error("Invalid LAMMPS atom id");
            rows.push_back({id,atom});
        }
        // Parallel runs write atoms in a different order in every frame.
        std::stable_sort(rows.begin(),rows.end(),[](const auto& a,const auto& b){ return a.first<b.first; });
        frame.atoms.reserve(rows.size());
        for (auto& row:rows) frame.atoms.push_back(std::move(row.second));
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
