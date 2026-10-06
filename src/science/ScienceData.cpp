#include "science/ScienceData.h"
#include "science/ScienceCatalog.h"
#include "io/Trajectory.h"
#include "util/ElementData.h"
#include "util/TaskControl.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#ifdef ATOMFORGE_SCIENCE_ZLIB
#include <zlib.h>
#endif
#ifdef ATOMFORGE_SCIENCE_BZIP2
#include <bzlib.h>
#endif
#ifdef ATOMFORGE_SCIENCE_LZMA
#include <lzma.h>
#endif

namespace atomforge::science
{
namespace
{
// ASE time unit conversion (CODATA 2014): 1 fs in sqrt(amu*A^2/eV).
const double kAseFemtosecond = 1e-5 * std::sqrt(1.6021766208e-19 / 1.660539040e-27);

std::string upper(std::string text)
{
    for (char& c : text) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return text;
}

std::string lower(std::string text)
{
    for (char& c : text) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return text;
}

std::string readFile(const std::filesystem::path& path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open " + path.u8string());
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

std::vector<std::string> split(const std::string& line)
{
    std::istringstream stream(line);
    std::vector<std::string> tokens;
    for (std::string token; stream >> token;) tokens.push_back(token);
    return tokens;
}

double parseDouble(const std::string& token, const char* what)
{
    try {
        std::size_t used = 0;
        const double value = std::stod(token, &used);
        if (used != token.size() || !std::isfinite(value)) throw std::invalid_argument(token);
        return value;
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("Invalid number in ") + what + ": " + token);
    }
}

AtomSite makeAtom(const std::string& symbol, const Vec3& position)
{
    AtomSite atom;
    atom.symbol = symbol;
    atom.atomicNumber = atomicNumber(symbol);
    if (!atom.atomicNumber) throw std::runtime_error("Unknown element: " + symbol);
    atom.x = position[0]; atom.y = position[1]; atom.z = position[2];
    getDefaultElementColor(atom.atomicNumber, atom.r, atom.g, atom.b);
    return atom;
}

// extXYZ comment-line key=value pairs, honouring quoted values.
std::map<std::string, std::string> parseComment(const std::string& line)
{
    std::map<std::string, std::string> result;
    std::size_t i = 0;
    while (i < line.size()) {
        while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
        const std::size_t keyStart = i;
        while (i < line.size() && line[i] != '=' && !std::isspace(static_cast<unsigned char>(line[i]))) ++i;
        const std::string key = line.substr(keyStart, i - keyStart);
        if (key.empty()) { ++i; continue; }
        if (i >= line.size() || line[i] != '=') { result[lower(key)] = "T"; continue; }
        ++i;
        std::string value;
        if (i < line.size() && (line[i] == '"' || line[i] == '{')) {
            const char close = line[i] == '"' ? '"' : '}';
            const std::size_t end = line.find(close, i + 1);
            value = line.substr(i + 1, (end == std::string::npos ? line.size() : end) - i - 1);
            i = end == std::string::npos ? line.size() : end + 1;
        } else {
            const std::size_t start = i;
            while (i < line.size() && !std::isspace(static_cast<unsigned char>(line[i]))) ++i;
            value = line.substr(start, i - start);
        }
        result[lower(key)] = value;
    }
    return result;
}

}

std::vector<FrameData> readExtxyzFrames(std::istream& input, std::size_t maxFrames)
{
    std::vector<FrameData> frames;
    std::string line;
    std::size_t total = 0;
    while (frames.size() < maxFrames && std::getline(input, line)) {
        taskCheckpoint();
        if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
        const auto countTokens = split(line);
        if (countTokens.size() != 1) throw std::runtime_error("Invalid XYZ frame atom count");
        const long long count = static_cast<long long>(parseDouble(countTokens[0], "XYZ atom count"));
        if (count < 0 || count > 10000000) throw std::runtime_error("Invalid XYZ frame atom count");
        total += static_cast<std::size_t>(count);
        if (total > 50000000) throw std::runtime_error("Trajectory exceeds memory limit");
        if (!std::getline(input, line)) throw std::runtime_error("Missing XYZ frame comment");
        const auto info = parseComment(line);
        FrameData frame;
        if (const auto lattice = info.find("lattice"); lattice != info.end()) {
            const auto tokens = split(lattice->second);
            if (tokens.size() != 9) throw std::runtime_error("Invalid XYZ lattice");
            for (int k = 0; k < 9; ++k) frame.structure.cellVectors[k / 3][k % 3] = parseDouble(tokens[k], "XYZ lattice");
            frame.structure.hasUnitCell = true;
            frame.pbc = {true, true, true};
        }
        if (const auto pbc = info.find("pbc"); pbc != info.end()) {
            const auto tokens = split(pbc->second);
            if (tokens.size() == 3)
                for (int k = 0; k < 3; ++k) frame.pbc[k] = upper(tokens[k])[0] == 'T';
        }
        if (const auto time = info.find("time_fs"); time != info.end()) frame.timeFs = parseDouble(time->second, "time_fs");
        // Columns: default plain XYZ layout species x y z.
        struct Column { std::string name; char type; int count; };
        std::vector<Column> columns = {{"species", 'S', 1}, {"pos", 'R', 3}};
        if (const auto properties = info.find("properties"); properties != info.end()) {
            columns.clear();
            std::vector<std::string> parts;
            std::stringstream stream(properties->second);
            for (std::string part; std::getline(stream, part, ':');) parts.push_back(part);
            if (parts.size() % 3) throw std::runtime_error("Invalid extXYZ Properties");
            for (std::size_t p = 0; p < parts.size(); p += 3)
                columns.push_back({lower(parts[p]), static_cast<char>(std::toupper(static_cast<unsigned char>(parts[p + 1][0]))),
                                   static_cast<int>(parseDouble(parts[p + 2], "Properties"))});
        }
        std::vector<Vec3> momenta;
        bool hasMomenta = false, hasVelocities = false, hasMasses = false;
        for (long long a = 0; a < count; ++a) {
            if (!std::getline(input, line)) throw std::runtime_error("Truncated XYZ trajectory");
            const auto tokens = split(line);
            std::size_t cursor = 0;
            std::string symbol;
            Vec3 position{}, velocity{}, momentum{};
            double mass = 0;
            for (const auto& column : columns) {
                if (cursor + static_cast<std::size_t>(column.count) > tokens.size()) throw std::runtime_error("Invalid XYZ atom row");
                if (column.name == "species") symbol = tokens[cursor];
                else if (column.name == "z" && column.type == 'I' && symbol.empty())
                    symbol = elementSymbol(static_cast<int>(parseDouble(tokens[cursor], "atomic number")));
                else if ((column.name == "pos" || column.name == "positions") && column.count == 3)
                    for (int k = 0; k < 3; ++k) position[k] = parseDouble(tokens[cursor + k], "XYZ positions");
                else if (column.name == "momenta" && column.count == 3) {
                    hasMomenta = true;
                    for (int k = 0; k < 3; ++k) momentum[k] = parseDouble(tokens[cursor + k], "XYZ momenta");
                } else if ((column.name == "velocities" || column.name == "vel") && column.count == 3) {
                    hasVelocities = true;
                    for (int k = 0; k < 3; ++k) velocity[k] = parseDouble(tokens[cursor + k], "XYZ velocities");
                } else if (column.name == "masses" && column.count == 1) {
                    hasMasses = true;
                    mass = parseDouble(tokens[cursor], "XYZ masses");
                }
                cursor += static_cast<std::size_t>(column.count);
            }
            if (symbol.empty()) throw std::runtime_error("extXYZ frame has no species column");
            frame.structure.atoms.push_back(makeAtom(symbol, position));
            frame.velocities.push_back(velocity);
            momenta.push_back(momentum);
            frame.masses.push_back(hasMasses ? mass : atomicMass(symbol));
        }
        if (hasMomenta && !hasVelocities)
            for (std::size_t a = 0; a < momenta.size(); ++a)
                frame.velocities[a] = scale(momenta[a], kAseFemtosecond / frame.masses[a]);
        if (!hasMomenta && !hasVelocities) frame.velocities.clear();
        if (!hasMasses) frame.masses.clear();
        frames.push_back(std::move(frame));
    }
    if (frames.empty()) throw std::runtime_error("Structure/trajectory file contains no frames");
    return frames;
}

namespace
{
std::vector<FrameData> readExtxyz(const std::filesystem::path& path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open " + path.u8string());
    return readExtxyzFrames(input, static_cast<std::size_t>(-1));
}

// VASP POSCAR/CONTCAR and XDATCAR (constant or repeated variable-cell headers).
}

bool VaspReader::readHeader(bool skipComment)
{
    std::string line;
    {
        if (!skipComment && !std::getline(in, line)) return false;
        if (!std::getline(in, line)) throw std::runtime_error("Truncated VASP header");
        const auto scaleTokens = split(line);
        if (scaleTokens.empty()) throw std::runtime_error("Missing VASP scale factor");
        double factor = parseDouble(scaleTokens[0], "VASP scale");
        for (auto& row : cell) {
            if (!std::getline(in, line)) throw std::runtime_error("Truncated VASP lattice");
            const auto tokens = split(line);
            if (tokens.size() < 3) throw std::runtime_error("Invalid VASP lattice");
            for (int k = 0; k < 3; ++k) row[k] = parseDouble(tokens[k], "VASP lattice");
        }
        if (factor < 0) factor = std::cbrt(-factor / cellVolume(cell));
        for (auto& row : cell) row = scale(row, factor);
        if (!std::getline(in, line)) throw std::runtime_error("Truncated VASP species");
        auto names = split(line);
        if (names.empty() || std::isdigit(static_cast<unsigned char>(names[0][0])))
            throw std::runtime_error("VASP files must list element symbols (VASP 5 format)");
        if (!std::getline(in, line)) throw std::runtime_error("Truncated VASP counts");
        const auto counts = split(line);
        if (counts.size() != names.size()) throw std::runtime_error("VASP species and counts differ");
        symbols.clear();
        for (std::size_t s = 0; s < names.size(); ++s) {
            std::string name = names[s].substr(0, names[s].find_first_of("/_"));
            const int count = static_cast<int>(parseDouble(counts[s], "VASP counts"));
            for (int c = 0; c < count; ++c) symbols.push_back(name);
        }
        return true;
    }
}

FrameData VaspReader::readPositions(bool direct)
{
    std::string line;
    {
        FrameData frame;
        frame.pbc = {true, true, true};
        frame.structure.hasUnitCell = true;
        for (int r = 0; r < 3; ++r) frame.structure.cellVectors[r] = cell[r];
        for (const auto& symbol : symbols) {
            if (!std::getline(in, line)) throw std::runtime_error("Truncated VASP coordinates");
            const auto tokens = split(line);
            if (tokens.size() < 3) throw std::runtime_error("Invalid VASP coordinates");
            Vec3 value{parseDouble(tokens[0], "VASP coordinates"), parseDouble(tokens[1], "VASP coordinates"), parseDouble(tokens[2], "VASP coordinates")};
            frame.structure.atoms.push_back(makeAtom(symbol, direct ? rowTimes(value, cell) : value));
        }
        return frame;
    }
}

namespace
{
std::vector<FrameData> readVasp(const std::filesystem::path& path, bool trajectory)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open " + path.u8string());
    std::vector<FrameData> frames;
    std::string line;
    VaspReader vasp(input);
    auto readHeader = [&](bool skipComment) { return vasp.readHeader(skipComment); };
    auto readPositions = [&](bool direct) { return vasp.readPositions(direct); };
    if (!readHeader(false)) throw std::runtime_error("Empty VASP file");
    if (!trajectory) {
        if (!std::getline(input, line)) throw std::runtime_error("Truncated POSCAR");
        if (!line.empty() && (line[0] == 'S' || line[0] == 's') && !std::getline(input, line)) throw std::runtime_error("Truncated POSCAR");
        const char mode = line.empty() ? 'D' : static_cast<char>(std::toupper(static_cast<unsigned char>(line.find_first_not_of(" \t") == std::string::npos ? 'D' : line[line.find_first_not_of(" \t")])));
        frames.push_back(readPositions(mode != 'C' && mode != 'K'));
        return frames;
    }
    while (std::getline(input, line)) {
        taskCheckpoint();
        if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
        const std::string head = lower(line.substr(line.find_first_not_of(" \t")));
        if (head.rfind("direct", 0) == 0 || head.rfind("cartesian", 0) == 0) {
            frames.push_back(readPositions(head[0] == 'd'));
        } else {
            readHeader(true);  // Variable-cell XDATCAR repeats the header.
        }
    }
    if (frames.empty()) throw std::runtime_error("XDATCAR contains no configurations");
    return frames;
}

Structure structureFromJson(const Json& value)
{
    if (!value.isObject()) throw std::runtime_error("Structure requires a file reference or symbols, positions and optional cell");
    for (const auto& member : value.members())
        if (member.first != "symbols" && member.first != "positions" && member.first != "cell")
            throw std::runtime_error("Structure requires a file reference or symbols, positions and optional cell");
    Structure structure;
    const auto& symbols = value.at("symbols").items();
    const auto positions = points(toArray(value.at("positions"), "positions"), "positions");
    if (symbols.size() != positions.size()) throw std::runtime_error("Structure symbols and positions differ in length");
    for (std::size_t i = 0; i < symbols.size(); ++i) structure.atoms.push_back(makeAtom(symbols[i].string(), positions[i]));
    if (const Json* cell = value.find("cell"); cell && !cell->isNull()) {
        const Mat3 matrix = cellMatrix(toArray(*cell, "cell"));
        for (int r = 0; r < 3; ++r) structure.cellVectors[r] = matrix[r];
        structure.hasUnitCell = true;
    }
    return structure;
}

NdArray loadNpy(const std::filesystem::path& path)
{
    const std::string data = readFile(path);
    if (data.size() < 10 || data.compare(0, 6, "\x93NUMPY") != 0) throw std::runtime_error("Invalid NPY file");
    const int major = static_cast<unsigned char>(data[6]);
    std::size_t headerLength = 0, offset = 0;
    if (major == 1) { headerLength = static_cast<unsigned char>(data[8]) | (static_cast<unsigned char>(data[9]) << 8); offset = 10; }
    else {
        if (data.size() < 12) throw std::runtime_error("Invalid NPY file");
        for (int b = 0; b < 4; ++b) headerLength |= static_cast<std::size_t>(static_cast<unsigned char>(data[8 + b])) << (8 * b);
        offset = 12;
    }
    if (offset + headerLength > data.size()) throw std::runtime_error("Invalid NPY header");
    const std::string header = data.substr(offset, headerLength);
    auto field = [&](const std::string& key) {
        const auto position = header.find("'" + key + "'");
        if (position == std::string::npos) throw std::runtime_error("NPY header lacks " + key);
        return header.substr(header.find(':', position) + 1);
    };
    std::string descr = field("descr");
    descr = descr.substr(descr.find('\'') + 1);
    descr = descr.substr(0, descr.find('\''));
    const bool fortran = field("fortran_order").find("True") < field("fortran_order").find(',');
    std::string shapeText = field("shape");
    shapeText = shapeText.substr(shapeText.find('(') + 1);
    shapeText = shapeText.substr(0, shapeText.find(')'));
    NdArray result;
    std::stringstream shapeStream(shapeText);
    for (std::string token; std::getline(shapeStream, token, ',');)
        if (token.find_first_not_of(" ") != std::string::npos) result.shape.push_back(static_cast<std::size_t>(std::stoull(token)));
    std::size_t count = 1;
    for (std::size_t extent : result.shape) count *= extent;
    if (descr.size() < 3 || descr[0] == '>') throw std::runtime_error("Only little-endian numeric NPY arrays are supported");
    const char kind = descr[1];
    const std::size_t width = static_cast<std::size_t>(std::stoul(descr.substr(2)));
    const std::size_t start = offset + headerLength;
    if (start + count * width > data.size()) throw std::runtime_error("Truncated NPY data");
    std::vector<double> values(count);
    for (std::size_t i = 0; i < count; ++i) {
        const char* bytes = data.data() + start + i * width;
        if (kind == 'f' && width == 8) { double v; std::memcpy(&v, bytes, 8); values[i] = v; }
        else if (kind == 'f' && width == 4) { float v; std::memcpy(&v, bytes, 4); values[i] = v; }
        else if (kind == 'i' && width == 8) { std::int64_t v; std::memcpy(&v, bytes, 8); values[i] = static_cast<double>(v); }
        else if (kind == 'i' && width == 4) { std::int32_t v; std::memcpy(&v, bytes, 4); values[i] = v; }
        else if (kind == 'i' && width == 2) { std::int16_t v; std::memcpy(&v, bytes, 2); values[i] = v; }
        else if (kind == 'u' && width == 1) values[i] = static_cast<unsigned char>(bytes[0]);
        else if (kind == 'b' && width == 1) values[i] = bytes[0] ? 1.0 : 0.0;
        else throw std::runtime_error("Unsupported NPY data type " + descr);
    }
    if (fortran && result.shape.size() > 1) {
        // Reorder column-major storage into row-major.
        result.values.assign(count, 0.0);
        std::vector<std::size_t> index(result.shape.size(), 0);
        for (std::size_t linear = 0; linear < count; ++linear) {
            std::size_t c = 0, stride = 1;
            for (std::size_t axis = result.shape.size(); axis-- > 0;) { c += index[axis] * stride; stride *= result.shape[axis]; }
            result.values[c] = values[linear];
            for (std::size_t axis = 0; axis < index.size(); ++axis) {
                if (++index[axis] < result.shape[axis]) break;
                index[axis] = 0;
            }
        }
    } else result.values = std::move(values);
    return result;
}

NdArray loadTable(const std::filesystem::path& path, bool comma)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open " + path.u8string());
    std::vector<std::vector<double>> rows;
    for (std::string line; std::getline(input, line);) {
        line = line.substr(0, line.find('#'));
        if (comma) std::replace(line.begin(), line.end(), ',', ' ');
        const auto tokens = split(line);
        if (tokens.empty()) continue;
        std::vector<double> row;
        for (const auto& token : tokens) row.push_back(parseDouble(token, "numeric table"));
        if (!rows.empty() && row.size() != rows[0].size()) throw std::runtime_error("Numeric table rows have different lengths");
        rows.push_back(std::move(row));
    }
    if (rows.empty()) throw std::runtime_error("Numeric table is empty");
    NdArray result({rows.size(), rows[0].size()});
    for (std::size_t i = 0; i < rows.size(); ++i)
        for (std::size_t j = 0; j < rows[i].size(); ++j) result(i, j) = rows[i][j];
    return result;
}

// Reads a text file, decompressing .gz, .bz2 and .xz files when the build has the codec.
std::string readMaybeCompressedImpl(const std::filesystem::path& path)
{
    const std::string extension = lower(path.extension().u8string());
    if (extension != ".gz" && extension != ".bz2" && extension != ".xz") return readFile(path);
    const std::string packed = readFile(path);
    std::string text;
    if (extension == ".gz") {
#ifdef ATOMFORGE_SCIENCE_ZLIB
        z_stream stream{};
        if (inflateInit2(&stream, 15 + 32) != Z_OK) throw std::runtime_error("Cannot initialise gzip decompression");
        stream.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(packed.data()));
        stream.avail_in = static_cast<uInt>(packed.size());
        char buffer[65536];
        int status = Z_OK;
        do {
            stream.next_out = reinterpret_cast<Bytef*>(buffer);
            stream.avail_out = sizeof(buffer);
            status = inflate(&stream, Z_NO_FLUSH);
            if (status != Z_OK && status != Z_STREAM_END) { inflateEnd(&stream); throw std::runtime_error("Invalid gzip file " + path.u8string()); }
            text.append(buffer, sizeof(buffer) - stream.avail_out);
        } while (status != Z_STREAM_END);
        inflateEnd(&stream);
        return text;
#endif
    } else if (extension == ".bz2") {
#ifdef ATOMFORGE_SCIENCE_BZIP2
        bz_stream stream{};
        if (BZ2_bzDecompressInit(&stream, 0, 0) != BZ_OK) throw std::runtime_error("Cannot initialise bzip2 decompression");
        stream.next_in = const_cast<char*>(packed.data());
        stream.avail_in = static_cast<unsigned>(packed.size());
        char buffer[65536];
        int status = BZ_OK;
        do {
            stream.next_out = buffer;
            stream.avail_out = sizeof(buffer);
            status = BZ2_bzDecompress(&stream);
            if (status != BZ_OK && status != BZ_STREAM_END) { BZ2_bzDecompressEnd(&stream); throw std::runtime_error("Invalid bzip2 file " + path.u8string()); }
            text.append(buffer, sizeof(buffer) - stream.avail_out);
        } while (status != BZ_STREAM_END);
        BZ2_bzDecompressEnd(&stream);
        return text;
#endif
    } else {
#ifdef ATOMFORGE_SCIENCE_LZMA
        lzma_stream stream = LZMA_STREAM_INIT;
        if (lzma_stream_decoder(&stream, UINT64_MAX, 0) != LZMA_OK) throw std::runtime_error("Cannot initialise xz decompression");
        stream.next_in = reinterpret_cast<const uint8_t*>(packed.data());
        stream.avail_in = packed.size();
        uint8_t buffer[65536];
        lzma_ret status = LZMA_OK;
        do {
            stream.next_out = buffer;
            stream.avail_out = sizeof(buffer);
            status = lzma_code(&stream, LZMA_FINISH);
            if (status != LZMA_OK && status != LZMA_STREAM_END) { lzma_end(&stream); throw std::runtime_error("Invalid xz file " + path.u8string()); }
            text.append(reinterpret_cast<const char*>(buffer), sizeof(buffer) - stream.avail_out);
        } while (status != LZMA_STREAM_END);
        lzma_end(&stream);
        return text;
#endif
    }
    throw std::runtime_error("This build cannot decompress " + extension + " files; decompress " + path.u8string() + " first");
}

// VASP EIGENVAL -> (spin, k, band) energies.
NdArray loadEigenval(const std::filesystem::path& path)
{
    std::istringstream input(readMaybeCompressedImpl(path));
    std::string line;
    std::getline(input, line);
    const auto first = split(line);
    if (first.size() < 4) throw std::runtime_error("Invalid EIGENVAL header");
    const std::size_t spins = static_cast<std::size_t>(parseDouble(first[3], "EIGENVAL"));
    for (int skip = 0; skip < 4; ++skip) std::getline(input, line);
    std::getline(input, line);
    const auto sizes = split(line);
    if (sizes.size() < 3 || (spins != 1 && spins != 2)) throw std::runtime_error("Invalid EIGENVAL header");
    const std::size_t kpoints = static_cast<std::size_t>(parseDouble(sizes[1], "EIGENVAL"));
    const std::size_t bands = static_cast<std::size_t>(parseDouble(sizes[2], "EIGENVAL"));
    NdArray result({spins, kpoints, bands});
    for (std::size_t k = 0; k < kpoints; ++k) {
        do { if (!std::getline(input, line)) throw std::runtime_error("Truncated EIGENVAL"); } while (split(line).empty());
        for (std::size_t b = 0; b < bands; ++b) {
            if (!std::getline(input, line)) throw std::runtime_error("Truncated EIGENVAL");
            const auto tokens = split(line);
            if (tokens.size() < 1 + spins) throw std::runtime_error("Invalid EIGENVAL band row");
            for (std::size_t s = 0; s < spins; ++s) result(s, k, b) = parseDouble(tokens[1 + s], "EIGENVAL");
        }
    }
    return result;
}

Json selectField(Json value, const std::string& field)
{
    std::stringstream stream(field);
    for (std::string component; std::getline(stream, component, '.');) {
        if (value.isArray()) {
            const long long index = std::stoll(component);
            if (index < 0 || static_cast<std::size_t>(index) >= value.size()) throw std::runtime_error("JSON field index out of range: " + component);
            value = Json(value.items()[static_cast<std::size_t>(index)]);
        } else {
            value = Json(value.at(component));
        }
    }
    return value;
}

NdArray selectColumn(const NdArray& table, std::optional<std::size_t> column)
{
    if (!column) return table;
    if (table.ndim() != 2 || *column >= table.shape[1]) throw std::runtime_error("Column is outside the two-dimensional data table");
    NdArray result({table.shape[0]});
    for (std::size_t i = 0; i < table.shape[0]; ++i) result(i) = table(i, *column);
    return result;
}

bool isVaspStructure(const std::filesystem::path& path)
{
    const std::string name = upper(path.filename().u8string());
    const std::string extension = lower(path.extension().u8string());
    // Common names: POSCAR, CONTCAR, POSCAR.relax, Si_POSCAR, *.vasp.
    return name.find("POSCAR") != std::string::npos || name.find("CONTCAR") != std::string::npos ||
           extension == ".vasp" || extension == ".poscar";
}

bool isLammpsDump(const std::filesystem::path& path)
{
    const std::string extension = lower(path.extension().u8string());
    if (extension == ".lammpstrj" || extension == ".dump") return true;
    std::ifstream input(path);
    std::string line;
    while (std::getline(input, line))
        if (line.find_first_not_of(" \t\r") != std::string::npos) return line.rfind("ITEM:", 0) == 0;
    return false;
}
}

std::string readTextFile(const std::filesystem::path& path)
{
    return readMaybeCompressedImpl(path);
}

int atomicNumber(const std::string& symbol)
{
    for (int z = 1; z <= 118; ++z)
        if (symbol == elementSymbol(z)) return z;
    return 0;
}

double atomicMass(const std::string& symbol)
{
    const int z = atomicNumber(symbol);
    if (!z) throw std::runtime_error("Unknown element: " + symbol);
    return elementAtomicMass(z);
}

std::vector<FrameData> readFrames(const std::filesystem::path& path, const StructureReader& reader)
{
    if (!std::filesystem::is_regular_file(path)) throw std::runtime_error("File not found: " + path.u8string());
    const std::string name = upper(path.filename().u8string());
    const std::string extension = lower(path.extension().u8string());
    if (name.rfind("XDATCAR", 0) == 0) return readVasp(path, true);
    if (isVaspStructure(path)) return readVasp(path, false);
    if (extension == ".xyz" || extension == ".extxyz") return readExtxyz(path);
    if (extension == ".json") {
        FrameData frame;
        frame.structure = structureFromJson(Json::parse(readFile(path)));
        if (frame.structure.hasUnitCell) frame.pbc = {true, true, true};
        return {frame};
    }
    if (isLammpsDump(path)) {
        std::vector<FrameData> frames;
        for (auto& structure : loadLammpsTrajectory(path.u8string())) {
            FrameData frame;
            frame.pbc = {structure.hasUnitCell, structure.hasUnitCell, structure.hasUnitCell};
            frame.structure = std::move(structure);
            frames.push_back(std::move(frame));
        }
        return frames;
    }
    if (!reader) throw std::runtime_error("Unsupported structure file: " + path.u8string());
    FrameData frame;
    frame.structure = reader(path);
    if (frame.structure.atoms.empty()) throw std::runtime_error("Structure file contains no atoms");
    if (frame.structure.hasUnitCell) frame.pbc = {true, true, true};
    return {frame};
}

TrajectoryStream::TrajectoryStream(const std::filesystem::path& path, const StructureReader& reader) : m_path(path)
{
    if (!std::filesystem::is_regular_file(path)) throw std::runtime_error("File not found: " + path.u8string());
    const std::string name = upper(path.filename().u8string());
    const std::string extension = lower(path.extension().u8string());
    // Binary mode keeps stream offsets exact on every platform; parsers accept CRLF.
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open " + path.u8string());
    std::string line;
    const auto skip = [&](long long lines) {
        for (long long k = 0; k < lines; ++k)
            if (!std::getline(input, line)) throw std::runtime_error("Truncated trajectory frame in " + path.u8string());
    };
    if (name.rfind("XDATCAR", 0) == 0) {
        m_format = Format::Xdatcar;
        VaspReader vasp(input);
        std::streamoff header = 0;
        if (!vasp.readHeader(false)) throw std::runtime_error("Empty VASP file");
        for (std::streamoff offset = input.tellg(); std::getline(input, line); offset = input.tellg()) {
            taskCheckpoint();
            if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
            const std::string head = lower(line.substr(line.find_first_not_of(" \t")));
            if (head.rfind("direct", 0) == 0 || head.rfind("cartesian", 0) == 0) {
                m_entries.push_back({input.tellg(), header, head[0] == 'd'});
                skip(static_cast<long long>(vasp.symbols.size()));
            } else {
                header = offset;  // a variable-cell XDATCAR repeats the header
                vasp.readHeader(true);
            }
        }
    } else if (extension == ".xyz" || extension == ".extxyz") {
        m_format = Format::Extxyz;
        for (std::streamoff offset = input.tellg(); std::getline(input, line); offset = input.tellg()) {
            taskCheckpoint();
            if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
            const auto tokens = split(line);
            if (tokens.size() != 1) throw std::runtime_error("Invalid XYZ frame atom count");
            const double count = parseDouble(tokens[0], "XYZ atom count");
            if (count < 0 || count > 10000000) throw std::runtime_error("Invalid XYZ frame atom count");
            m_entries.push_back({offset, 0, true});
            skip(static_cast<long long>(count) + 1);
        }
    } else if (!isVaspStructure(path) && extension != ".json" && isLammpsDump(path)) {
        m_format = Format::Lammps;
        for (std::streamoff offset = input.tellg(); std::getline(input, line); offset = input.tellg()) {
            taskCheckpoint();
            if (line.rfind("ITEM: TIMESTEP", 0) == 0) m_entries.push_back({offset, 0, true});
        }
    } else {
        m_kept = readFrames(path, reader);
        return;
    }
    if (m_entries.empty()) throw std::runtime_error("Structure/trajectory file contains no frames");
}

std::size_t TrajectoryStream::size() const
{
    return m_format == Format::Kept ? m_kept.size() : m_entries.size();
}

FrameData TrajectoryStream::frame(std::size_t index) const
{
    if (index >= size()) throw std::out_of_range("Trajectory frame " + std::to_string(index) + " of " + std::to_string(size()));
    if (m_format == Format::Kept) return m_kept[index];
    std::ifstream input(m_path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open " + m_path.u8string());
    const Entry& entry = m_entries[index];
    if (m_format == Format::Xdatcar) {
        VaspReader vasp(input);
        input.seekg(entry.header);
        vasp.readHeader(false);
        input.seekg(entry.start);
        return vasp.readPositions(entry.direct);
    }
    input.seekg(entry.start);
    if (m_format == Format::Extxyz) return readExtxyzFrames(input, 1).front();
    FrameData frame;
    frame.structure = loadLammpsFrames(input, 1).front();
    frame.pbc = {frame.structure.hasUnitCell, frame.structure.hasUnitCell, frame.structure.hasUnitCell};
    return frame;
}

namespace
{
NdArray framesField(const std::vector<FrameData>& frames, const std::string& field, std::optional<double> timestep)
{
    if (frames.empty()) throw std::runtime_error("Structure/trajectory file contains no frames");
    if (timestep && (field == "positions" || field == "velocities")) {
        bool any = false, all = true;
        for (const auto& frame : frames) { any = any || frame.timeFs.has_value(); all = all && frame.timeFs.has_value(); }
        if (any) {
            if (!all) throw std::runtime_error("Trajectory has incomplete time_fs metadata");
            for (std::size_t f = 1; f < frames.size(); ++f) {
                const double interval = *frames[f].timeFs - *frames[f - 1].timeFs;
                if (!(interval > 0) || std::abs(interval - *timestep) > 1e-10 + 1e-8 * std::abs(*timestep))
                    throw std::runtime_error("Trajectory sampling must be uniform and match timestep_fs; resample or trim the final partial interval");
            }
        }
    }
    const auto& first = frames[0].structure.atoms;
    for (const auto& frame : frames) {
        const auto& atoms = frame.structure.atoms;
        bool same = atoms.size() == first.size();
        for (std::size_t a = 0; same && a < atoms.size(); ++a) same = atoms[a].symbol == first[a].symbol;
        if (!same) throw std::runtime_error("Trajectory atom species/order must remain consistent");
    }
    const std::size_t atoms = first.size();
    if (field == "velocities") {
        NdArray result({frames.size(), atoms, 3});
        for (std::size_t f = 0; f < frames.size(); ++f) {
            if (frames[f].velocities.size() != atoms) throw std::runtime_error("Trajectory has no velocities; supply explicit A/fs data");
            for (std::size_t a = 0; a < atoms; ++a)
                for (std::size_t k = 0; k < 3; ++k) result(f, a, k) = frames[f].velocities[a][k];
        }
        return result;
    }
    if (field == "cell") {
        if (!frames[0].structure.hasUnitCell) throw std::runtime_error("Structure file has no cell");
        for (const auto& frame : frames)
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    if (!frame.structure.hasUnitCell || std::abs(frame.structure.cellVectors[r][c] - frames[0].structure.cellVectors[r][c]) > 1e-10)
                        throw std::runtime_error("This analysis requires a fixed cell");
        NdArray result({3, 3});
        for (std::size_t r = 0; r < 3; ++r)
            for (std::size_t c = 0; c < 3; ++c) result(r, c) = frames[0].structure.cellVectors[r][c];
        return result;
    }
    if (field == "masses") {
        NdArray result({atoms});
        for (std::size_t a = 0; a < atoms; ++a)
            result(a) = frames[0].masses.size() == atoms ? frames[0].masses[a] : atomicMass(first[a].symbol);
        return result;
    }
    if (field != "positions") throw std::runtime_error("Structure field must be positions, velocities, cell or masses");
    NdArray result({frames.size(), atoms, 3});
    for (std::size_t f = 0; f < frames.size(); ++f)
        for (std::size_t a = 0; a < atoms; ++a) {
            const auto& atom = frames[f].structure.atoms[a];
            result(f, a, 0) = atom.x; result(f, a, 1) = atom.y; result(f, a, 2) = atom.z;
        }
    if (frames.size() == 1) result.shape.erase(result.shape.begin());
    return result;
}

StructureInput resolveStructure(const Json& value, const std::filesystem::path& base, const StructureReader& reader)
{
    if (value.isObject() && value.size() == 1 && value.contains("file")) {
        const auto path = (base / std::filesystem::u8path(value.at("file").string())).lexically_normal();
        if (lower(path.extension().u8string()) == ".json") {
            const Json parsed = Json::parse(readFile(path));
            return resolveStructure(parsed, path.parent_path(), reader);
        }
        auto frames = readFrames(path, reader);
        // A trajectory supplies its final frame, as when reading one structure.
        return {frames.back().structure, frames.back().pbc};
    }
    StructureInput input;
    input.structure = structureFromJson(value);
    const bool periodic = input.structure.hasUnitCell;
    input.pbc = {periodic, periodic, periodic};
    return input;
}
}

namespace
{
// Structure inputs are the catalog parameters of kind "structure" (for internal
// requests without a catalog entry: structure, initial and final).
bool isStructureParameter(const std::string& tool, const std::string& name)
{
    if (const ScienceToolDef* definition = findScienceTool(tool)) {
        for (const auto& parameter : definition->parameters)
            if (name == parameter.name) return std::strcmp(parameter.kind, "structure") == 0;
        return false;
    }
    return name == "structure" || name == "initial" || name == "final";
}
}

Parameters::Parameters(const std::string& tool, const Json& request, const std::filesystem::path& base,
                       const StructureReader& reader)
{
    if (!request.isObject()) throw std::runtime_error("parameters must be an object");
    std::optional<double> timestep;
    if (tool == "msd" || tool == "vacf" || tool == "vibrational-spectrum")
        if (const Json* value = request.find("timestep_fs"); value && value->isNumber()) timestep = value->number();
    for (const auto& [name, value] : request.members()) {
        taskCheckpoint();
        if (isStructureParameter(tool, name)) {
            m_structures[name] = resolveStructure(value, base, reader);
            continue;
        }
        if (name == "restart_images") {
            if (!value.isObject() || !value.contains("file")) throw std::runtime_error("restart_images must reference a structure/trajectory file");
            m_frames[name] = readFrames((base / std::filesystem::u8path(value.at("file").string())).lexically_normal(), reader);
            continue;
        }
        if (name.size() > 5 && name.compare(name.size() - 5, 5, "_file") == 0) {
            const std::string path = value.isString() ? value.string() : value.at("file").string();
            m_json[name] = (base / std::filesystem::u8path(path)).lexically_normal().u8string();
            continue;
        }
        if (!value.isObject() || !value.contains("file") || name == "calculator" || name == "calculator_factory") {
            m_json[name] = value;
            continue;
        }
        for (const auto& member : value.members())
            if (member.first != "file" && member.first != "field" && member.first != "column")
                throw std::runtime_error("Unknown data-reference options for " + name);
        const auto path = (base / std::filesystem::u8path(value.at("file").string())).lexically_normal();
        std::optional<std::string> field;
        if (const Json* f = value.find("field"); f && !f->isNull()) field = f->string();
        std::optional<std::size_t> column;
        if (const Json* c = value.find("column"); c && !c->isNull()) {
            if (!c->isNumber() || c->number() < 0 || c->number() != std::floor(c->number()))
                throw std::runtime_error("Column must be a nonnegative integer");
            column = static_cast<std::size_t>(c->number());
        }
        const std::string extension = lower(path.extension().u8string());
        const std::string filename = upper(path.filename().u8string());
        if (!std::filesystem::is_regular_file(path)) throw std::runtime_error("File not found: " + path.u8string());
        if (extension == ".json") {
            Json data = Json::parse(readFile(path));
            if (field) data = selectField(data, *field);
            // A JSON structure file without a field selects its positions.
            if (!field && data.isObject() && data.contains("positions") && data.contains("symbols")) data = Json(data.at("positions"));
            if (column) m_arrays[name] = selectColumn(toArray(data, name), column);
            else m_json[name] = data;
        } else if (filename.rfind("EIGENVAL", 0) == 0 && name == "energies_eV") {
            if (field || column) throw std::runtime_error("EIGENVAL supplies all spin/k-point/band energies; field/column selection is unsupported");
            m_arrays[name] = loadEigenval(path);
        } else if (extension == ".npy") {
            if (field) throw std::runtime_error("NPY arrays do not have named fields");
            m_arrays[name] = selectColumn(loadNpy(path), column);
        } else if (extension == ".csv" || extension == ".txt" || extension == ".dat") {
            if (field) throw std::runtime_error("Numeric text tables do not have named fields; select a column instead");
            NdArray table = loadTable(path, extension == ".csv");
            if (!column && table.shape[1] == 1) table.shape = {table.shape[0]};
            m_arrays[name] = selectColumn(table, column);
        } else {
            if (column) throw std::runtime_error("Column selection is supported only for JSON/NPY/numeric tables");
            std::string selected = field ? *field : "positions";
            if (!field) {
                if (name == "velocities") selected = "velocities";
                else if (name == "masses") selected = "masses";
                else if (name == "cell" || name == "reference_cell" || name == "current_cell") selected = "cell";
            }
            m_arrays[name] = framesField(readFrames(path, reader), selected, timestep);
        }
    }
}

bool Parameters::has(const std::string& name) const
{
    // An explicit null selects the documented default, as None does in Python.
    if (const auto found = m_json.find(name); found != m_json.end()) return !found->second.isNull();
    return m_arrays.count(name) || m_structures.count(name) || m_frames.count(name);
}

const Json& Parameters::json(const std::string& name) const
{
    const auto found = m_json.find(name);
    if (found == m_json.end()) throw std::runtime_error("Provide " + name);
    return found->second;
}

NdArray Parameters::array(const std::string& name) const
{
    if (const auto found = m_arrays.find(name); found != m_arrays.end()) return found->second;
    return toArray(json(name), name);
}

double Parameters::number(const std::string& name) const
{
    if (const auto found = m_json.find(name); found != m_json.end() && found->second.isNumber()) {
        if (!std::isfinite(found->second.number())) throw std::runtime_error(name + " must be finite");
        return found->second.number();
    }
    const NdArray value = array(name);
    if (value.size() != 1 || !std::isfinite(value.values[0])) throw std::runtime_error(name + " must be one finite number");
    return value.values[0];
}

double Parameters::number(const std::string& name, double fallback) const
{
    if (!has(name)) return fallback;
    return number(name);
}

bool Parameters::boolean(const std::string& name, bool fallback) const
{
    if (!has(name)) return fallback;
    const Json& value = json(name);
    if (!value.isBool()) throw std::runtime_error(name + " must be a boolean, not a string or number");
    return value.boolean();
}

Pbc Parameters::pbc(const std::string& name, Pbc fallback) const
{
    if (!has(name)) return fallback;
    return parsePbc(json(name));
}

const std::vector<FrameData>& Parameters::frames(const std::string& name) const
{
    const auto found = m_frames.find(name);
    if (found == m_frames.end()) throw std::runtime_error("Provide " + name);
    return found->second;
}

const StructureInput& Parameters::structure(const std::string& name) const
{
    const auto found = m_structures.find(name);
    if (found == m_structures.end()) throw std::runtime_error("Provide " + name);
    return found->second;
}

void writeExtxyz(const std::filesystem::path& path, const std::vector<Structure>& frames,
                 const std::vector<std::vector<Vec3>>& velocities, const std::vector<double>& times)
{
    if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Cannot write " + path.u8string());
    output << std::setprecision(12);
    for (std::size_t f = 0; f < frames.size(); ++f) {
        const auto& frame = frames[f];
        const bool withMomenta = f < velocities.size() && velocities[f].size() == frame.atoms.size();
        output << frame.atoms.size() << '\n';
        if (frame.hasUnitCell) {
            output << "Lattice=\"";
            for (int k = 0; k < 9; ++k) output << (k ? " " : "") << frame.cellVectors[k / 3][k % 3];
            output << "\" ";
        }
        output << "Properties=species:S:1:pos:R:3" << (withMomenta ? ":momenta:R:3" : "");
        if (f < times.size()) output << " time_fs=" << times[f];
        output << " pbc=\"" << (frame.hasUnitCell ? "T T T" : "F F F") << "\"\n";
        for (std::size_t a = 0; a < frame.atoms.size(); ++a) {
            const auto& atom = frame.atoms[a];
            output << atom.symbol << ' ' << atom.x << ' ' << atom.y << ' ' << atom.z;
            if (withMomenta) {
                const double mass = atomicMass(atom.symbol);
                for (int k = 0; k < 3; ++k) output << ' ' << velocities[f][a][k] * mass / kAseFemtosecond;
            }
            output << '\n';
        }
    }
    if (!output) throw std::runtime_error("Cannot write " + path.u8string());
}

Json structureJson(const Structure& structure)
{
    Json result = Json::object();
    Json symbols = Json::array(), positions = Json::array();
    for (const auto& atom : structure.atoms) {
        symbols.push(atom.symbol);
        positions.push(Json::array({atom.x, atom.y, atom.z}));
    }
    result["symbols"] = symbols;
    result["positions"] = positions;
    if (structure.hasUnitCell) {
        Mat3 cell{};
        for (int r = 0; r < 3; ++r) cell[r] = structure.cellVectors[r];
        result["cell"] = toJson(cell);
    } else result["cell"] = Json();
    return result;
}
}
