#pragma once
// Out-of-core atom clouds (.afcloud) for structures far larger than memory:
// billions of atoms, as in large molecular-dynamics snapshots.
//
// Atoms are stored by space in chunks of up to ~0.5 M atoms. Each atom takes
// 8 bytes: x, y, z quantized to 16 bits within its chunk's box and a 16-bit
// species index. Inside a chunk the atoms are in random order, so any prefix
// of a chunk is a uniform random subsample of it: a renderer reads and draws
// only as much of each chunk as its size on screen needs (level of detail
// without storing extra levels).
//
// Layout: "AFCLOUD1", a little-endian uint64 header length, the header as
// JSON (species, chunks with offsets and boxes), then the chunk data.
#include <array>
#include <atomic>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace atomforge::cloud
{
struct Species
{
    std::string name;            // element symbol or "type N"
    std::array<float, 3> color{0.7f, 0.7f, 0.7f};
    float radius = 1.3f;         // display radius (Angstrom)
};

struct PackedAtom { std::uint16_t x, y, z, species; };
static_assert(sizeof(PackedAtom) == 8, "packed atoms are 8 bytes");

struct Chunk
{
    std::uint64_t offset = 0;    // file offset of its first atom
    std::uint64_t count = 0;
    std::array<double, 3> lower{}, upper{};  // box of its atoms (quantization frame)
};

struct CloudInfo
{
    std::uint64_t atoms = 0;
    std::array<double, 3> lower{}, upper{};
    bool hasCell = false;
    std::array<std::array<double, 3>, 3> cell{};
    std::array<double, 3> cellOrigin{};
    std::vector<Species> species;
    std::vector<Chunk> chunks;
    std::string source;          // what it was built from
};

// A read-only .afcloud file. Reads are thread-safe.
class CloudFile
{
public:
    explicit CloudFile(const std::filesystem::path& path);
    const CloudInfo& info() const { return m_info; }
    const std::filesystem::path& path() const { return m_path; }
    // Reads atoms [first, first + count) of a chunk.
    void read(std::size_t chunk, std::uint64_t first, std::uint64_t count, PackedAtom* out) const;
    // Decodes a packed atom of a chunk to a position.
    std::array<double, 3> position(std::size_t chunk, const PackedAtom& atom) const;

private:
    std::filesystem::path m_path;
    CloudInfo m_info;
    mutable std::mutex m_mutex;
    mutable std::ifstream m_file;
};

// One atom handed to the builder.
struct SourceAtom { double x, y, z; std::uint16_t species; };

// Progress of long operations: fraction (0..1, or <0 when unknown) and a
// message; return false to cancel.
using Progress = std::function<bool(double fraction, const std::string& message)>;

struct BuildOptions
{
    std::uint64_t chunkAtoms = 1u << 19;     // target atoms per chunk
    std::filesystem::path temporary;         // folder for the bucket files (default: next to the output)
    unsigned threads = 0;                    // 0 = hardware concurrency
    std::uint64_t seed = 12345;              // shuffle seed (reproducible files)
};

// Streams atoms into a .afcloud without holding them in memory: atoms are
// first spread over spatial buckets on disk, then each bucket is split into
// chunks, shuffled, quantized and written.
class CloudBuilder
{
public:
    // `lower`/`upper` bound the atoms (atoms outside are clamped into the grid);
    // `expectedAtoms` sizes the bucket grid (an estimate is fine).
    CloudBuilder(const std::filesystem::path& output, const std::array<double, 3>& lower, const std::array<double, 3>& upper,
                 std::uint64_t expectedAtoms, BuildOptions options = {});
    ~CloudBuilder();
    // Species index for a name (added on first use; colours and radii from the element table).
    std::uint16_t species(const std::string& name);
    void setSpecies(std::vector<Species> species) { m_species = std::move(species); }
    void setCell(const std::array<std::array<double, 3>, 3>& cell, const std::array<double, 3>& origin);
    void setSource(const std::string& source) { m_source = source; }
    void add(const SourceAtom& atom);
    void add(const SourceAtom* atoms, std::size_t count);
    std::uint64_t added() const { return m_added; }
    // Writes the file; returns its description.
    CloudInfo finish(const Progress& progress = {});

private:
    struct Bucket;
    void flush(Bucket& bucket);
    std::filesystem::path m_output, m_folder;
    BuildOptions m_options;
    std::array<double, 3> m_lower, m_upper;
    std::array<int, 3> m_grid{1, 1, 1};
    std::vector<std::unique_ptr<Bucket>> m_buckets;
    std::vector<Species> m_species;
    bool m_hasCell = false;
    std::array<std::array<double, 3>, 3> m_cell{};
    std::array<double, 3> m_cellOrigin{};
    std::string m_source;
    std::uint64_t m_added = 0;
    bool m_finished = false;
};

// Converters. Each streams its input (memory use does not grow with the file).
// LAMMPS text dumps (first frame unless `frame` is given): columns x y z, xs ys zs,
// xu yu zu or xsu ysu zsu, with element or type; `typeNames` names LAMMPS types
// ("Cu,Ni" for types 1 and 2).
CloudInfo buildFromLammpsDump(const std::filesystem::path& input, const std::filesystem::path& output,
                              const std::vector<std::string>& typeNames = {}, std::size_t frame = 0,
                              const Progress& progress = {}, BuildOptions options = {});
// XYZ and extended XYZ (first frame unless `frame` is given).
CloudInfo buildFromXyz(const std::filesystem::path& input, const std::filesystem::path& output, std::size_t frame = 0,
                       const Progress& progress = {}, BuildOptions options = {});
// A crystal of nx x ny x nz cells: lattice fcc, bcc, sc, diamond or hcp with
// lattice constant a (c for hcp, ideal when 0) and one element (or two, alternating
// on the basis, for "B2"/"L12"-like test data with `second`).
CloudInfo generateCrystal(const std::filesystem::path& output, const std::string& lattice, double a, double c,
                          std::array<std::uint64_t, 3> cells, const std::string& element, const std::string& second = "",
                          const Progress& progress = {}, BuildOptions options = {});

// Whether a file is (by its first bytes) an .afcloud.
bool isCloudFile(const std::filesystem::path& path);
// Atom count of a LAMMPS dump or XYZ file from its first header (0 if unknown):
// lets the desktop decide whether a file needs the large-data path.
std::uint64_t peekAtomCount(const std::filesystem::path& path);
}
