#include "cloud/AtomCloud.h"

#include "science/Json.h"
#include "util/ElementData.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <thread>

namespace atomforge::cloud
{
using science::Json;

namespace
{
constexpr char kMagic[8] = {'A', 'F', 'C', 'L', 'O', 'U', 'D', '1'};
constexpr std::uint64_t kInMemoryLimit = 8u << 20;   // records a bucket may hold in memory

// A bucket record on disk (16 bytes).
struct Record { float x, y, z; std::uint16_t species, pad; };
static_assert(sizeof(Record) == 16, "bucket records are 16 bytes");

int elementNumber(const std::string& symbol)
{
    for (int z = 1; z <= 118; ++z)
        if (symbol == elementSymbol(z)) return z;
    // Lower/upper-case variants ("CU", "cu").
    if (symbol.empty() || symbol.size() > 3) return 0;
    std::string fixed = symbol;
    fixed[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(fixed[0])));
    for (std::size_t i = 1; i < fixed.size(); ++i) fixed[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(fixed[i])));
    for (int z = 1; z <= 118; ++z)
        if (fixed == elementSymbol(z)) return z;
    return 0;
}

Species speciesFor(const std::string& name, std::size_t index)
{
    Species s;
    s.name = name;
    if (const int z = elementNumber(name)) {
        static const auto radii = makeLiteratureCovalentRadii();
        float r = 0.7f, g = 0.7f, b = 0.7f;
        getDefaultElementColor(z, r, g, b);
        s.color = {r, g, b};
        s.radius = std::max(0.3f, radii[static_cast<std::size_t>(z)]);
    } else {
        // Distinct colours for LAMMPS types without element names.
        static const float palette[][3] = {{0.85f, 0.53f, 0.20f}, {0.25f, 0.55f, 0.90f}, {0.35f, 0.75f, 0.35f}, {0.85f, 0.30f, 0.35f},
                                           {0.65f, 0.45f, 0.85f}, {0.90f, 0.80f, 0.30f}, {0.40f, 0.80f, 0.80f}, {0.60f, 0.60f, 0.60f}};
        const auto& p = palette[index % 8];
        s.color = {p[0], p[1], p[2]};
        s.radius = 1.3f;
    }
    return s;
}

Json vec(const std::array<double, 3>& v) { return Json::array({v[0], v[1], v[2]}); }
std::array<double, 3> vec(const Json& j) { return {j.items().at(0).number(), j.items().at(1).number(), j.items().at(2).number()}; }

void writeU64(std::ostream& out, std::uint64_t value)
{
    unsigned char bytes[8];
    for (int i = 0; i < 8; ++i) bytes[i] = static_cast<unsigned char>(value >> (8 * i));
    out.write(reinterpret_cast<const char*>(bytes), 8);
}
std::uint64_t readU64(std::istream& in)
{
    unsigned char bytes[8] = {};
    in.read(reinterpret_cast<char*>(bytes), 8);
    std::uint64_t value = 0;
    for (int i = 0; i < 8; ++i) value |= static_cast<std::uint64_t>(bytes[i]) << (8 * i);
    return value;
}

// Fast whitespace-separated field parsing for large text files.
inline const char* skipSpace(const char* p, const char* end)
{
    while (p < end && (*p == ' ' || *p == '\t' || *p == '\r')) ++p;
    return p;
}
inline const char* skipField(const char* p, const char* end)
{
    while (p < end && *p != ' ' && *p != '\t' && *p != '\r' && *p != '\n') ++p;
    return p;
}
bool parseDouble(const char* begin, const char* end, double& value)
{
    const auto result = std::from_chars(begin, end, value);
    return result.ec == std::errc() && result.ptr == end;
}

// Reads a text file line by line through a large buffer, tracking the bytes consumed.
class LineReader
{
public:
    explicit LineReader(const std::filesystem::path& path) : m_file(path, std::ios::binary)
    {
        if (!m_file) throw std::runtime_error("Cannot read " + path.u8string());
        std::error_code ignored;
        m_size = std::filesystem::file_size(path, ignored);
        m_buffer.resize(64u << 20);
    }
    // The next line without its newline; false at the end of the file.
    bool next(const char*& begin, const char*& end)
    {
        for (;;) {
            const char* data = m_buffer.data();
            const char* newline = static_cast<const char*>(std::memchr(data + m_position, '\n', m_filled - m_position));
            if (newline) {
                begin = data + m_position;
                end = newline;
                m_position = static_cast<std::size_t>(newline - data) + 1;
                m_consumed += static_cast<std::uint64_t>(end - begin) + 1;
                if (end > begin && end[-1] == '\r') --end;
                return true;
            }
            if (m_eof) {
                if (m_position >= m_filled) return false;
                begin = data + m_position;
                end = data + m_filled;
                m_consumed += m_filled - m_position;
                m_position = m_filled;
                if (end > begin && end[-1] == '\r') --end;
                return true;
            }
            // Move the partial line to the front and refill.
            const std::size_t rest = m_filled - m_position;
            if (rest == m_buffer.size()) m_buffer.resize(m_buffer.size() * 2);
            std::memmove(m_buffer.data(), m_buffer.data() + m_position, rest);
            m_filled = rest;
            m_position = 0;
            m_file.read(m_buffer.data() + m_filled, static_cast<std::streamsize>(m_buffer.size() - m_filled));
            m_filled += static_cast<std::size_t>(m_file.gcount());
            if (m_file.gcount() == 0 || !m_file) m_eof = true;
        }
    }
    std::string line()
    {
        const char *b = nullptr, *e = nullptr;
        if (!next(b, e)) throw std::runtime_error("Unexpected end of file");
        return std::string(b, e);
    }
    double fraction() const { return m_size ? static_cast<double>(m_consumed) / static_cast<double>(m_size) : -1.0; }

private:
    std::ifstream m_file;
    std::vector<char> m_buffer;
    std::size_t m_filled = 0, m_position = 0;
    std::uint64_t m_consumed = 0, m_size = 0;
    bool m_eof = false;
};

std::vector<std::string> words(const std::string& line)
{
    std::istringstream in(line);
    std::vector<std::string> out;
    for (std::string w; in >> w;) out.push_back(w);
    return out;
}

unsigned threadCount(unsigned requested)
{
    const unsigned hardware = std::max(1u, std::thread::hardware_concurrency());
    return requested ? requested : std::max(1u, std::min(hardware, 16u));
}
}

// ---------------------------------------------------------------- reading
CloudFile::CloudFile(const std::filesystem::path& path) : m_path(path), m_file(path, std::ios::binary)
{
    if (!m_file) throw std::runtime_error("Cannot open " + path.u8string());
    char magic[8] = {};
    m_file.read(magic, 8);
    if (!m_file || std::memcmp(magic, kMagic, 8) != 0) throw std::runtime_error(path.u8string() + " is not an AtomForge atom cloud");
    const std::uint64_t headerOffset = readU64(m_file), headerLength = readU64(m_file);
    if (!headerOffset || headerLength > (1ull << 32)) throw std::runtime_error(path.u8string() + " is incomplete (its build did not finish)");
    std::string text(headerLength, '\0');
    m_file.seekg(static_cast<std::streamoff>(headerOffset));
    m_file.read(text.data(), static_cast<std::streamsize>(headerLength));
    if (!m_file) throw std::runtime_error("Cannot read the header of " + path.u8string());
    const Json header = Json::parse(text);
    m_info.atoms = static_cast<std::uint64_t>(header.at("atoms").number());
    m_info.lower = vec(header.at("lower"));
    m_info.upper = vec(header.at("upper"));
    if (const Json* source = header.find("source")) m_info.source = source->string();
    if (const Json* cell = header.find("cell")) {
        m_info.hasCell = true;
        for (int r = 0; r < 3; ++r) m_info.cell[static_cast<std::size_t>(r)] = vec(cell->items().at(static_cast<std::size_t>(r)));
        m_info.cellOrigin = vec(header.at("cell_origin"));
    }
    for (const auto& s : header.at("species").items()) {
        Species species;
        species.name = s.at("name").string();
        const auto c = vec(s.at("color"));
        species.color = {static_cast<float>(c[0]), static_cast<float>(c[1]), static_cast<float>(c[2])};
        species.radius = static_cast<float>(s.at("radius").number());
        m_info.species.push_back(species);
    }
    for (const auto& c : header.at("chunks").items()) {
        Chunk chunk;
        chunk.offset = static_cast<std::uint64_t>(c.at("offset").number());
        chunk.count = static_cast<std::uint64_t>(c.at("count").number());
        chunk.lower = vec(c.at("lower"));
        chunk.upper = vec(c.at("upper"));
        m_info.chunks.push_back(chunk);
    }
}

void CloudFile::read(std::size_t chunk, std::uint64_t first, std::uint64_t count, PackedAtom* out) const
{
    const Chunk& c = m_info.chunks.at(chunk);
    if (first + count > c.count) throw std::out_of_range("Atom cloud read past the end of a chunk");
    std::lock_guard<std::mutex> lock(m_mutex);
    m_file.clear();
    m_file.seekg(static_cast<std::streamoff>(c.offset + first * sizeof(PackedAtom)));
    m_file.read(reinterpret_cast<char*>(out), static_cast<std::streamsize>(count * sizeof(PackedAtom)));
    if (!m_file) throw std::runtime_error("Cannot read atoms from " + m_path.u8string());
}

std::array<double, 3> CloudFile::position(std::size_t chunk, const PackedAtom& atom) const
{
    const Chunk& c = m_info.chunks.at(chunk);
    const std::uint16_t q[3] = {atom.x, atom.y, atom.z};
    std::array<double, 3> p{};
    for (int k = 0; k < 3; ++k)
        p[static_cast<std::size_t>(k)] = c.lower[static_cast<std::size_t>(k)] +
            (c.upper[static_cast<std::size_t>(k)] - c.lower[static_cast<std::size_t>(k)]) * q[k] / 65535.0;
    return p;
}

bool isCloudFile(const std::filesystem::path& path)
{
    std::ifstream in(path, std::ios::binary);
    char magic[8] = {};
    in.read(magic, 8);
    return in && std::memcmp(magic, kMagic, 8) == 0;
}

std::uint64_t peekAtomCount(const std::filesystem::path& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) return 0;
    std::string line;
    for (int i = 0; i < 64 && std::getline(in, line); ++i) {
        if (i == 0) {
            // XYZ: the first line is the atom count.
            std::istringstream first(line);
            unsigned long long n = 0;
            std::string rest;
            if (first >> n && !(first >> rest)) return n;
        }
        if (line.rfind("ITEM: NUMBER OF ATOMS", 0) == 0 && std::getline(in, line)) {
            try { return std::stoull(line); } catch (const std::exception&) { return 0; }
        }
    }
    return 0;
}

// ---------------------------------------------------------------- building
struct CloudBuilder::Bucket
{
    std::filesystem::path path;
    std::vector<Record> pending;
    std::uint64_t count = 0;
};

CloudBuilder::CloudBuilder(const std::filesystem::path& output, const std::array<double, 3>& lower, const std::array<double, 3>& upper,
                           std::uint64_t expectedAtoms, BuildOptions options)
    : m_output(output), m_options(std::move(options)), m_lower(lower), m_upper(upper)
{
    if (m_options.chunkAtoms < 1024) m_options.chunkAtoms = 1024;
    m_folder = m_options.temporary.empty() ? std::filesystem::path(output).replace_extension(".afcloud-build") : m_options.temporary;
    std::filesystem::create_directories(m_folder);
    // Buckets of a few chunks each, laid out in proportion to the box.
    const double buckets = std::clamp(std::ceil(static_cast<double>(expectedAtoms) / (4.0 * static_cast<double>(m_options.chunkAtoms))), 1.0, 4096.0);
    std::array<double, 3> extent{};
    double volume = 1.0;
    int used = 0;
    for (int k = 0; k < 3; ++k) {
        extent[static_cast<std::size_t>(k)] = std::max(0.0, upper[static_cast<std::size_t>(k)] - lower[static_cast<std::size_t>(k)]);
        if (extent[static_cast<std::size_t>(k)] > 1e-9) { volume *= extent[static_cast<std::size_t>(k)]; ++used; }
    }
    const double perAxis = used ? std::pow(buckets / volume, 1.0 / used) : 0.0;
    std::size_t total = 1;
    for (int k = 0; k < 3; ++k) {
        const double e = extent[static_cast<std::size_t>(k)];
        m_grid[static_cast<std::size_t>(k)] = e > 1e-9 ? std::clamp(static_cast<int>(std::lround(e * perAxis)), 1, 64) : 1;
        total *= static_cast<std::size_t>(m_grid[static_cast<std::size_t>(k)]);
    }
    m_buckets.reserve(total);
    for (std::size_t i = 0; i < total; ++i) {
        auto bucket = std::make_unique<Bucket>();
        bucket->path = m_folder / ("bucket-" + std::to_string(i) + ".bin");
        std::error_code ignored;
        std::filesystem::remove(bucket->path, ignored);
        m_buckets.push_back(std::move(bucket));
    }
}

CloudBuilder::~CloudBuilder()
{
    std::error_code ignored;
    std::filesystem::remove_all(m_folder, ignored);
}

std::uint16_t CloudBuilder::species(const std::string& name)
{
    for (std::size_t i = 0; i < m_species.size(); ++i)
        if (m_species[i].name == name) return static_cast<std::uint16_t>(i);
    if (m_species.size() >= 65535) throw std::runtime_error("Too many species in an atom cloud");
    m_species.push_back(speciesFor(name, m_species.size()));
    return static_cast<std::uint16_t>(m_species.size() - 1);
}

void CloudBuilder::setCell(const std::array<std::array<double, 3>, 3>& cell, const std::array<double, 3>& origin)
{
    m_hasCell = true;
    m_cell = cell;
    m_cellOrigin = origin;
}

void CloudBuilder::flush(Bucket& bucket)
{
    if (bucket.pending.empty()) return;
    std::FILE* file = std::fopen(bucket.path.u8string().c_str(), "ab");
    if (!file) throw std::runtime_error("Cannot write " + bucket.path.u8string());
    const std::size_t written = std::fwrite(bucket.pending.data(), sizeof(Record), bucket.pending.size(), file);
    std::fclose(file);
    if (written != bucket.pending.size()) throw std::runtime_error("Cannot write " + bucket.path.u8string() + " (disk full?)");
    bucket.pending.clear();
}

void CloudBuilder::add(const SourceAtom& atom)
{
    std::size_t index = 0, stride = 1;
    const double p[3] = {atom.x, atom.y, atom.z};
    for (int k = 0; k < 3; ++k) {
        const auto K = static_cast<std::size_t>(k);
        const double e = m_upper[K] - m_lower[K];
        int cell = e > 1e-9 ? static_cast<int>((p[k] - m_lower[K]) / e * m_grid[K]) : 0;
        cell = std::clamp(cell, 0, m_grid[K] - 1);
        index += static_cast<std::size_t>(cell) * stride;
        stride *= static_cast<std::size_t>(m_grid[K]);
    }
    Bucket& bucket = *m_buckets[index];
    bucket.pending.push_back({static_cast<float>(atom.x), static_cast<float>(atom.y), static_cast<float>(atom.z), atom.species, 0});
    ++bucket.count;
    ++m_added;
    if (bucket.pending.size() >= 16384) flush(bucket);
}

void CloudBuilder::add(const SourceAtom* atoms, std::size_t count)
{
    for (std::size_t i = 0; i < count; ++i) add(atoms[i]);
}

namespace
{
// Chunk data produced from one bucket.
struct BuiltChunk { std::vector<PackedAtom> atoms; Chunk chunk; };

void boundsOf(const Record* r, std::size_t n, std::array<double, 3>& lower, std::array<double, 3>& upper)
{
    lower = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    upper = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
    for (std::size_t i = 0; i < n; ++i) {
        const double p[3] = {r[i].x, r[i].y, r[i].z};
        for (int k = 0; k < 3; ++k) {
            lower[static_cast<std::size_t>(k)] = std::min(lower[static_cast<std::size_t>(k)], p[k]);
            upper[static_cast<std::size_t>(k)] = std::max(upper[static_cast<std::size_t>(k)], p[k]);
        }
    }
}

// Splits records into chunks of at most `limit` atoms (median cuts along the longest axis).
void split(Record* r, std::size_t n, std::uint64_t limit, std::mt19937_64& random, std::vector<BuiltChunk>& out)
{
    std::array<double, 3> lower{}, upper{};
    boundsOf(r, n, lower, upper);
    if (n > limit) {
        int axis = 0;
        for (int k = 1; k < 3; ++k)
            if (upper[static_cast<std::size_t>(k)] - lower[static_cast<std::size_t>(k)] > upper[static_cast<std::size_t>(axis)] - lower[static_cast<std::size_t>(axis)]) axis = k;
        const std::size_t half = n / 2;
        std::nth_element(r, r + half, r + n, [axis](const Record& a, const Record& b) {
            return (axis == 0 ? a.x : axis == 1 ? a.y : a.z) < (axis == 0 ? b.x : axis == 1 ? b.y : b.z);
        });
        split(r, half, limit, random, out);
        split(r + half, n - half, limit, random, out);
        return;
    }
    // Random order: any prefix is a uniform subsample.
    std::shuffle(r, r + n, random);
    BuiltChunk built;
    built.chunk.count = n;
    built.chunk.lower = lower;
    built.chunk.upper = upper;
    built.atoms.resize(n);
    double scale[3];
    for (int k = 0; k < 3; ++k) {
        const double e = upper[static_cast<std::size_t>(k)] - lower[static_cast<std::size_t>(k)];
        scale[k] = e > 0 ? 65535.0 / e : 0.0;
    }
    for (std::size_t i = 0; i < n; ++i) {
        const double p[3] = {r[i].x, r[i].y, r[i].z};
        std::uint16_t q[3];
        for (int k = 0; k < 3; ++k)
            q[k] = static_cast<std::uint16_t>(std::clamp(std::lround((p[k] - lower[static_cast<std::size_t>(k)]) * scale[k]), 0L, 65535L));
        built.atoms[i] = {q[0], q[1], q[2], r[i].species};
    }
    out.push_back(std::move(built));
}

std::vector<Record> readRecords(const std::filesystem::path& path, std::uint64_t count)
{
    std::vector<Record> records(count);
    if (!count) return records;
    std::FILE* file = std::fopen(path.u8string().c_str(), "rb");
    if (!file) throw std::runtime_error("Cannot read " + path.u8string());
    const std::size_t got = std::fread(records.data(), sizeof(Record), records.size(), file);
    std::fclose(file);
    if (got != records.size()) throw std::runtime_error("Short read from " + path.u8string());
    return records;
}

// Splits a bucket file too large for memory into eight files by octant, streaming.
std::vector<std::pair<std::filesystem::path, std::uint64_t>> splitFile(const std::filesystem::path& path, std::uint64_t count)
{
    // The bucket's bounds, then its octants.
    std::array<double, 3> lower{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    std::array<double, 3> upper{std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
    std::vector<Record> block(1u << 20);
    const auto stream = [&](const std::function<void(const Record*, std::size_t)>& use) {
        std::FILE* file = std::fopen(path.u8string().c_str(), "rb");
        if (!file) throw std::runtime_error("Cannot read " + path.u8string());
        for (std::size_t got; (got = std::fread(block.data(), sizeof(Record), block.size(), file)) > 0;) use(block.data(), got);
        std::fclose(file);
    };
    stream([&](const Record* r, std::size_t n) {
        std::array<double, 3> lo{}, hi{};
        boundsOf(r, n, lo, hi);
        for (std::size_t k = 0; k < 3; ++k) { lower[k] = std::min(lower[k], lo[k]); upper[k] = std::max(upper[k], hi[k]); }
    });
    const double mid[3] = {0.5 * (lower[0] + upper[0]), 0.5 * (lower[1] + upper[1]), 0.5 * (lower[2] + upper[2])};
    std::vector<std::pair<std::filesystem::path, std::uint64_t>> parts(8);
    std::vector<std::vector<Record>> pending(8);
    for (int i = 0; i < 8; ++i) {
        parts[static_cast<std::size_t>(i)].first = path.string() + "." + std::to_string(i);
        std::error_code ignored;
        std::filesystem::remove(parts[static_cast<std::size_t>(i)].first, ignored);
    }
    const auto flushPart = [&](int i) {
        auto& p = pending[static_cast<std::size_t>(i)];
        if (p.empty()) return;
        std::FILE* file = std::fopen(parts[static_cast<std::size_t>(i)].first.u8string().c_str(), "ab");
        if (!file) throw std::runtime_error("Cannot write a bucket part");
        std::fwrite(p.data(), sizeof(Record), p.size(), file);
        std::fclose(file);
        p.clear();
    };
    stream([&](const Record* r, std::size_t n) {
        for (std::size_t i = 0; i < n; ++i) {
            const int octant = (r[i].x > mid[0] ? 1 : 0) | (r[i].y > mid[1] ? 2 : 0) | (r[i].z > mid[2] ? 4 : 0);
            pending[static_cast<std::size_t>(octant)].push_back(r[i]);
            ++parts[static_cast<std::size_t>(octant)].second;
            if (pending[static_cast<std::size_t>(octant)].size() >= 65536) flushPart(octant);
        }
    });
    for (int i = 0; i < 8; ++i) flushPart(i);
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
    // All atoms at one point cannot be split by space: fall back to halves in file order.
    for (const auto& part : parts)
        if (part.second == count) throw std::runtime_error("Too many atoms at one position to split into chunks");
    return parts;
}
}

CloudInfo CloudBuilder::finish(const Progress& progress)
{
    if (m_finished) throw std::logic_error("CloudBuilder::finish called twice");
    m_finished = true;
    for (auto& bucket : m_buckets) flush(*bucket);

    std::ofstream out(m_output, std::ios::binary | std::ios::trunc);
    if (!out) throw std::runtime_error("Cannot write " + m_output.u8string());
    out.write(kMagic, 8);
    writeU64(out, 0);  // header offset and length, written at the end
    writeU64(out, 0);
    std::uint64_t offset = 24;

    // Buckets (and parts of oversized ones) to turn into chunks.
    struct Job { std::filesystem::path path; std::uint64_t count; std::size_t order; };
    std::vector<Job> jobs;
    for (std::size_t i = 0; i < m_buckets.size(); ++i)
        if (m_buckets[i]->count) jobs.push_back({m_buckets[i]->path, m_buckets[i]->count, i});
    std::vector<Chunk> chunks;
    std::vector<std::pair<std::size_t, std::size_t>> chunkOrder;  // (job order, index within job) for a stable layout
    std::mutex writing;
    std::uint64_t done = 0;
    std::atomic<bool> cancelled{false};
    std::exception_ptr failure;
    std::size_t next = 0;
    std::mutex jobsMutex;
    const auto worker = [&] {
        try {
            for (;;) {
                Job job;
                {
                    std::lock_guard<std::mutex> lock(jobsMutex);
                    if (next >= jobs.size() || cancelled || failure) return;
                    job = jobs[next++];
                }
                if (job.count > kInMemoryLimit) {
                    // Too large for memory: split on disk and queue the parts.
                    auto parts = splitFile(job.path, job.count);
                    std::lock_guard<std::mutex> lock(jobsMutex);
                    for (std::size_t p = 0; p < parts.size(); ++p)
                        if (parts[p].second) jobs.push_back({parts[p].first, parts[p].second, job.order * 8 + p + (1u << 20)});
                    continue;
                }
                auto records = readRecords(job.path, job.count);
                std::error_code ignored;
                std::filesystem::remove(job.path, ignored);
                std::mt19937_64 random(m_options.seed ^ (0x9E3779B97F4A7C15ull * (job.order + 1)));
                std::vector<BuiltChunk> built;
                split(records.data(), records.size(), m_options.chunkAtoms, random, built);
                records = {};
                std::lock_guard<std::mutex> lock(writing);
                for (std::size_t c = 0; c < built.size(); ++c) {
                    built[c].chunk.offset = offset;
                    out.write(reinterpret_cast<const char*>(built[c].atoms.data()), static_cast<std::streamsize>(built[c].atoms.size() * sizeof(PackedAtom)));
                    if (!out) throw std::runtime_error("Cannot write " + m_output.u8string() + " (disk full?)");
                    offset += built[c].atoms.size() * sizeof(PackedAtom);
                    chunks.push_back(built[c].chunk);
                    chunkOrder.push_back({job.order, c});
                    done += built[c].chunk.count;
                }
                if (progress && !progress(m_added ? static_cast<double>(done) / static_cast<double>(m_added) : 1.0, "Writing chunks")) cancelled = true;
            }
        } catch (...) {
            std::lock_guard<std::mutex> lock(jobsMutex);
            if (!failure) failure = std::current_exception();
        }
    };
    // Workers pull jobs until none are left (a worker that splits a bucket queues and takes its parts).
    for (;;) {
        std::vector<std::thread> threads;
        for (unsigned t = 0; t < threadCount(m_options.threads); ++t) threads.emplace_back(worker);
        for (auto& t : threads) t.join();
        if (failure) std::rethrow_exception(failure);
        if (cancelled) throw std::runtime_error("Cancelled");
        if (next >= jobs.size()) break;
    }

    // A stable chunk order (by bucket) regardless of thread timing.
    std::vector<std::size_t> order(chunks.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return chunkOrder[a] < chunkOrder[b]; });

    CloudInfo info;
    info.atoms = done;
    info.species = m_species;
    info.hasCell = m_hasCell;
    info.cell = m_cell;
    info.cellOrigin = m_cellOrigin;
    info.source = m_source;
    info.lower = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    info.upper = {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
    Json chunkList = Json::array();
    for (std::size_t i : order) {
        const Chunk& c = chunks[i];
        info.chunks.push_back(c);
        for (std::size_t k = 0; k < 3; ++k) { info.lower[k] = std::min(info.lower[k], c.lower[k]); info.upper[k] = std::max(info.upper[k], c.upper[k]); }
        Json entry = Json::object();
        entry["offset"] = static_cast<double>(c.offset);
        entry["count"] = static_cast<double>(c.count);
        entry["lower"] = vec(c.lower);
        entry["upper"] = vec(c.upper);
        chunkList.push(entry);
    }
    if (info.chunks.empty()) info.lower = info.upper = {0, 0, 0};
    Json header = Json::object();
    header["format"] = "atomforge-cloud";
    header["version"] = 1;
    header["atoms"] = static_cast<double>(info.atoms);
    header["lower"] = vec(info.lower);
    header["upper"] = vec(info.upper);
    header["source"] = info.source;
    if (m_hasCell) {
        header["cell"] = Json::array({vec(m_cell[0]), vec(m_cell[1]), vec(m_cell[2])});
        header["cell_origin"] = vec(m_cellOrigin);
    }
    Json speciesList = Json::array();
    for (const auto& s : m_species) {
        Json entry = Json::object();
        entry["name"] = s.name;
        entry["color"] = Json::array({s.color[0], s.color[1], s.color[2]});
        entry["radius"] = s.radius;
        speciesList.push(entry);
    }
    header["species"] = speciesList;
    header["chunks"] = chunkList;
    const std::string text = header.dump();
    out.write(text.data(), static_cast<std::streamsize>(text.size()));
    out.seekp(8);
    writeU64(out, offset);
    writeU64(out, text.size());
    out.close();
    if (!out) throw std::runtime_error("Cannot write " + m_output.u8string());
    if (progress) progress(1.0, "Done");
    return info;
}

// ---------------------------------------------------------------- converters
CloudInfo buildFromLammpsDump(const std::filesystem::path& input, const std::filesystem::path& output,
                              const std::vector<std::string>& typeNames, std::size_t frame, const Progress& progress, BuildOptions options)
{
    LineReader reader(input);
    for (std::size_t current = 0;; ++current) {
        // Frame header.
        std::string line;
        do { line = reader.line(); } while (line.rfind("ITEM: TIMESTEP", 0) != 0);
        reader.line();
        if (reader.line().rfind("ITEM: NUMBER OF ATOMS", 0) != 0) throw std::runtime_error("Not a LAMMPS dump: " + input.u8string());
        const std::uint64_t atoms = std::stoull(reader.line());
        const std::string boxLine = reader.line();
        if (boxLine.rfind("ITEM: BOX BOUNDS", 0) != 0) throw std::runtime_error("LAMMPS dump without BOX BOUNDS");
        const bool triclinic = boxLine.find("xy") != std::string::npos;
        double bounds[3][3] = {};
        for (auto& row : bounds) {
            const auto w = words(reader.line());
            if (w.size() < 2) throw std::runtime_error("Bad BOX BOUNDS line");
            row[0] = std::stod(w[0]); row[1] = std::stod(w[1]); row[2] = w.size() > 2 ? std::stod(w[2]) : 0.0;
        }
        const auto columns = words(reader.line());
        if (columns.size() < 3 || columns[0] != "ITEM:" || columns[1] != "ATOMS") throw std::runtime_error("LAMMPS dump without an ATOMS section");
        if (current < frame) {
            for (std::uint64_t i = 0; i < atoms; ++i) reader.line();
            continue;
        }
        // Box: orthogonal or triclinic (LAMMPS bounding-box convention).
        const double xy = triclinic ? bounds[0][2] : 0, xz = triclinic ? bounds[1][2] : 0, yz = triclinic ? bounds[2][2] : 0;
        const double xlo = bounds[0][0] - std::min({0.0, xy, xz, xy + xz}), xhi = bounds[0][1] - std::max({0.0, xy, xz, xy + xz});
        const double ylo = bounds[1][0] - std::min(0.0, yz), yhi = bounds[1][1] - std::max(0.0, yz);
        const double zlo = bounds[2][0], zhi = bounds[2][1];
        // Columns of the positions, element and type.
        int cx = -1, cy = -1, cz = -1, element = -1, type = -1;
        bool scaled = false;
        for (std::size_t i = 2; i < columns.size(); ++i) {
            const std::string& c = columns[i];
            const int index = static_cast<int>(i) - 2;
            if (c == "x" || c == "xu") cx = index;
            else if (c == "y" || c == "yu") cy = index;
            else if (c == "z" || c == "zu") cz = index;
            else if (c == "element") element = index;
            else if (c == "type") type = index;
        }
        if (cx < 0) {
            for (std::size_t i = 2; i < columns.size(); ++i) {
                const std::string& c = columns[i];
                const int index = static_cast<int>(i) - 2;
                if (c == "xs" || c == "xsu") cx = index;
                else if (c == "ys" || c == "ysu") cy = index;
                else if (c == "zs" || c == "zsu") cz = index;
            }
            scaled = true;
        }
        if (cx < 0 || cy < 0 || cz < 0) throw std::runtime_error("The dump has no x y z (or xs ys zs) columns");
        CloudBuilder builder(output, {bounds[0][0], bounds[1][0], bounds[2][0]}, {bounds[0][1], bounds[1][1], bounds[2][1]}, atoms, options);
        builder.setCell({{{xhi - xlo, 0, 0}, {xy, yhi - ylo, 0}, {xz, yz, zhi - zlo}}}, {xlo, ylo, zlo});
        builder.setSource("LAMMPS dump " + input.filename().u8string() + ", frame " + std::to_string(frame));
        // Type and element names map to species indices through small caches.
        std::vector<int> typeSpecies;
        std::vector<std::pair<std::string, std::uint16_t>> elementSpecies;
        const int needed = std::max({cx, cy, cz, element, type});
        std::vector<std::pair<const char*, const char*>> fields(static_cast<std::size_t>(needed) + 1);
        std::vector<SourceAtom> batch;
        batch.reserve(1u << 16);
        for (std::uint64_t i = 0; i < atoms; ++i) {
            const char *b = nullptr, *e = nullptr;
            if (!reader.next(b, e)) throw std::runtime_error("The dump ends after " + std::to_string(i) + " of " + std::to_string(atoms) + " atoms");
            const char* p = b;
            for (int f = 0; f <= needed; ++f) {
                p = skipSpace(p, e);
                const char* start = p;
                p = skipField(p, e);
                if (start == p) throw std::runtime_error("Short atom line in the dump: " + std::string(b, e));
                fields[static_cast<std::size_t>(f)] = {start, p};
            }
            double r[3];
            const int col[3] = {cx, cy, cz};
            for (int k = 0; k < 3; ++k) {
                const auto& f = fields[static_cast<std::size_t>(col[k])];
                if (!parseDouble(f.first, f.second, r[k])) throw std::runtime_error("Bad coordinate in the dump: " + std::string(f.first, f.second));
            }
            if (scaled) {
                const double s0 = r[0], s1 = r[1], s2 = r[2];
                r[0] = xlo + s0 * (xhi - xlo) + s1 * xy + s2 * xz;
                r[1] = ylo + s1 * (yhi - ylo) + s2 * yz;
                r[2] = zlo + s2 * (zhi - zlo);
            }
            std::uint16_t s = 0;
            if (element >= 0) {
                const auto& f = fields[static_cast<std::size_t>(element)];
                const std::size_t length = static_cast<std::size_t>(f.second - f.first);
                bool found = false;
                for (const auto& [name, index] : elementSpecies)
                    if (name.size() == length && std::memcmp(name.data(), f.first, length) == 0) { s = index; found = true; break; }
                if (!found) { std::string name(f.first, f.second); s = builder.species(name); elementSpecies.push_back({name, s}); }
            } else if (type >= 0) {
                const auto& f = fields[static_cast<std::size_t>(type)];
                int t = 0;
                std::from_chars(f.first, f.second, t);
                if (t < 0 || t > 65535) t = 0;
                if (static_cast<std::size_t>(t) >= typeSpecies.size()) typeSpecies.resize(static_cast<std::size_t>(t) + 1, -1);
                if (typeSpecies[static_cast<std::size_t>(t)] < 0)
                    typeSpecies[static_cast<std::size_t>(t)] = builder.species(t >= 1 && static_cast<std::size_t>(t) <= typeNames.size() ? typeNames[static_cast<std::size_t>(t) - 1] : "type " + std::to_string(t));
                s = static_cast<std::uint16_t>(typeSpecies[static_cast<std::size_t>(t)]);
            } else s = builder.species("X");
            batch.push_back({r[0], r[1], r[2], s});
            if (batch.size() == batch.capacity()) {
                builder.add(batch.data(), batch.size());
                batch.clear();
                if (progress && !progress(0.5 * reader.fraction(), "Reading " + input.filename().u8string())) throw std::runtime_error("Cancelled");
            }
        }
        builder.add(batch.data(), batch.size());
        return builder.finish([&](double f, const std::string& m) { return !progress || progress(0.5 + 0.5 * f, m); });
    }
}

CloudInfo buildFromXyz(const std::filesystem::path& input, const std::filesystem::path& output, std::size_t frame,
                       const Progress& progress, BuildOptions options)
{
    // Pass 1: the frame's bounds (XYZ has no box); pass 2: the atoms.
    std::array<double, 3> lower{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    std::array<double, 3> upper{std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest()};
    std::uint64_t atoms = 0;
    std::string comment;
    int species = 0, position = 1;
    const auto scan = [&](const std::function<void(const char* symbol, std::size_t length, const double* r)>& use, double share, double base) {
        LineReader reader(input);
        for (std::size_t current = 0;; ++current) {
            const auto countLine = reader.line();
            atoms = std::stoull(countLine);
            comment = reader.line();
            // Extended XYZ: species and position columns from Properties=.
            species = 0; position = 1;
            if (const auto at = comment.find("Properties="); at != std::string::npos) {
                std::string spec = comment.substr(at + 11, comment.find_first_of(" \t", at + 11) - at - 11);
                std::vector<std::string> parts;
                std::stringstream split(spec);
                for (std::string part; std::getline(split, part, ':');) parts.push_back(part);
                int column = 0;
                for (std::size_t i = 0; i + 2 < parts.size(); i += 3) {
                    const int width = std::stoi(parts[i + 2]);
                    if (parts[i] == "species") species = column;
                    if (parts[i] == "pos") position = column;
                    column += width;
                }
            }
            if (current < frame) { for (std::uint64_t i = 0; i < atoms; ++i) reader.line(); continue; }
            const int needed = std::max(species, position + 2);
            std::vector<std::pair<const char*, const char*>> fields(static_cast<std::size_t>(needed) + 1);
            for (std::uint64_t i = 0; i < atoms; ++i) {
                const char *b = nullptr, *e = nullptr;
                if (!reader.next(b, e)) throw std::runtime_error("The XYZ file ends early");
                const char* p = b;
                for (int f = 0; f <= needed; ++f) {
                    p = skipSpace(p, e);
                    const char* start = p;
                    p = skipField(p, e);
                    fields[static_cast<std::size_t>(f)] = {start, p};
                }
                double r[3];
                for (int k = 0; k < 3; ++k) {
                    const auto& f = fields[static_cast<std::size_t>(position + k)];
                    if (!parseDouble(f.first, f.second, r[k])) throw std::runtime_error("Bad coordinate in the XYZ file: " + std::string(b, e));
                }
                const auto& s = fields[static_cast<std::size_t>(species)];
                use(s.first, static_cast<std::size_t>(s.second - s.first), r);
                if ((i & 0xFFFF) == 0 && progress && !progress(base + share * reader.fraction(), "Reading " + input.filename().u8string()))
                    throw std::runtime_error("Cancelled");
            }
            return;
        }
    };
    scan([&](const char*, std::size_t, const double* r) {
        for (std::size_t k = 0; k < 3; ++k) { lower[k] = std::min(lower[k], r[k]); upper[k] = std::max(upper[k], r[k]); }
    }, 0.25, 0.0);
    CloudBuilder builder(output, lower, upper, atoms, options);
    builder.setSource("XYZ " + input.filename().u8string() + ", frame " + std::to_string(frame));
    // Extended XYZ lattice.
    if (const auto at = comment.find("Lattice=\""); at != std::string::npos) {
        std::istringstream values(comment.substr(at + 9, comment.find('"', at + 9) - at - 9));
        std::array<std::array<double, 3>, 3> cell{};
        bool ok = true;
        for (auto& row : cell) for (double& v : row) ok = ok && static_cast<bool>(values >> v);
        if (ok) builder.setCell(cell, {0, 0, 0});
    }
    std::vector<std::pair<std::string, std::uint16_t>> names;
    scan([&](const char* symbol, std::size_t length, const double* r) {
        std::uint16_t s = 0;
        bool found = false;
        for (const auto& [name, index] : names)
            if (name.size() == length && std::memcmp(name.data(), symbol, length) == 0) { s = index; found = true; break; }
        if (!found) { std::string name(symbol, length); s = builder.species(name); names.push_back({name, s}); }
        builder.add({r[0], r[1], r[2], s});
    }, 0.25, 0.25);
    return builder.finish([&](double f, const std::string& m) { return !progress || progress(0.5 + 0.5 * f, m); });
}

CloudInfo generateCrystal(const std::filesystem::path& output, const std::string& lattice, double a, double c,
                          std::array<std::uint64_t, 3> cells, const std::string& element, const std::string& second,
                          const Progress& progress, BuildOptions options)
{
    if (a <= 0) throw std::invalid_argument("The lattice constant must be positive");
    for (auto n : cells) if (n < 1) throw std::invalid_argument("Cell counts must be at least 1");
    // Cell vectors and basis (fractional).
    std::array<std::array<double, 3>, 3> cell{{{a, 0, 0}, {0, a, 0}, {0, 0, a}}};
    std::vector<std::array<double, 3>> basis;
    if (lattice == "fcc") basis = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
    else if (lattice == "bcc") basis = {{0, 0, 0}, {0.5, 0.5, 0.5}};
    else if (lattice == "sc") basis = {{0, 0, 0}};
    else if (lattice == "diamond") basis = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}, {0.25, 0.25, 0.25}, {0.75, 0.75, 0.25}, {0.75, 0.25, 0.75}, {0.25, 0.75, 0.75}};
    else if (lattice == "hcp") {
        // Orthohexagonal cell (a, sqrt(3) a, c) with four atoms.
        const double cc = c > 0 ? c : a * std::sqrt(8.0 / 3.0);
        cell = {{{a, 0, 0}, {0, std::sqrt(3.0) * a, 0}, {0, 0, cc}}};
        basis = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 1.0 / 6.0, 0.5}, {0, 2.0 / 3.0, 0.5}};
    } else throw std::invalid_argument("Unknown lattice " + lattice + " (fcc, bcc, sc, diamond or hcp)");
    const double total = static_cast<double>(cells[0]) * static_cast<double>(cells[1]) * static_cast<double>(cells[2]) * static_cast<double>(basis.size());
    if (total > 2.0e10) throw std::invalid_argument("More than 20 billion atoms");
    std::array<std::array<double, 3>, 3> super{};
    for (std::size_t r = 0; r < 3; ++r)
        for (std::size_t k = 0; k < 3; ++k) super[r][k] = cell[r][k] * static_cast<double>(cells[r]);
    const std::array<double, 3> upper{super[0][0], super[1][1], super[2][2]};
    CloudBuilder builder(output, {0, 0, 0}, upper, static_cast<std::uint64_t>(total), options);
    builder.setCell(super, {0, 0, 0});
    std::ostringstream description;
    description << lattice << " " << element << (second.empty() ? "" : "/" + second) << ", a = " << a << ", " << cells[0] << " x " << cells[1] << " x " << cells[2] << " cells";
    builder.setSource(description.str());
    const std::uint16_t first = builder.species(element);
    const std::uint16_t other = second.empty() ? first : builder.species(second);
    std::vector<SourceAtom> batch;
    batch.reserve(1u << 16);
    for (std::uint64_t k = 0; k < cells[2]; ++k) {
        for (std::uint64_t j = 0; j < cells[1]; ++j)
            for (std::uint64_t i = 0; i < cells[0]; ++i)
                for (std::size_t b = 0; b < basis.size(); ++b) {
                    const double f[3] = {static_cast<double>(i) + basis[b][0], static_cast<double>(j) + basis[b][1], static_cast<double>(k) + basis[b][2]};
                    batch.push_back({f[0] * cell[0][0] + f[1] * cell[1][0] + f[2] * cell[2][0],
                                     f[0] * cell[0][1] + f[1] * cell[1][1] + f[2] * cell[2][1],
                                     f[0] * cell[0][2] + f[1] * cell[1][2] + f[2] * cell[2][2], b % 2 ? other : first});
                    if (batch.size() == batch.capacity()) { builder.add(batch.data(), batch.size()); batch.clear(); }
                }
        if (progress && !progress(0.5 * static_cast<double>(k + 1) / static_cast<double>(cells[2]), "Generating atoms")) throw std::runtime_error("Cancelled");
    }
    builder.add(batch.data(), batch.size());
    return builder.finish([&](double f, const std::string& m) { return !progress || progress(0.5 + 0.5 * f, m); });
}
}
