#include "science/GifWriter.h"

#include <algorithm>
#include <array>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

namespace atomforge::science
{
namespace
{
void put16(std::ofstream& out, int value)
{
    out.put(static_cast<char>(value & 0xFF));
    out.put(static_cast<char>((value >> 8) & 0xFF));
}

// LSB-first bit packer emitting GIF data sub-blocks.
class BitWriter
{
public:
    explicit BitWriter(std::ofstream& out) : m_out(out) {}
    void write(int code, int bits)
    {
        m_accumulator |= static_cast<std::uint32_t>(code) << m_bits;
        m_bits += bits;
        while (m_bits >= 8) {
            m_block.push_back(static_cast<unsigned char>(m_accumulator & 0xFF));
            m_accumulator >>= 8;
            m_bits -= 8;
            if (m_block.size() == 255) flushBlock();
        }
    }
    void finish()
    {
        if (m_bits > 0) m_block.push_back(static_cast<unsigned char>(m_accumulator & 0xFF));
        m_bits = 0;
        m_accumulator = 0;
        flushBlock();
        m_out.put(0);
    }
private:
    void flushBlock()
    {
        if (m_block.empty()) return;
        m_out.put(static_cast<char>(m_block.size()));
        m_out.write(reinterpret_cast<const char*>(m_block.data()), static_cast<std::streamsize>(m_block.size()));
        m_block.clear();
    }
    std::ofstream& m_out;
    std::uint32_t m_accumulator = 0;
    int m_bits = 0;
    std::vector<unsigned char> m_block;
};

void encode(std::ofstream& out, const std::vector<unsigned char>& indices)
{
    const int clear = 256, end = 257;
    out.put(8);  // minimum code size
    BitWriter writer(out);
    std::unordered_map<std::uint32_t, int> dictionary;
    dictionary.reserve(8192);
    int codeSize = 9, next = 258;
    writer.write(clear, codeSize);
    if (indices.empty()) { writer.write(end, codeSize); writer.finish(); return; }
    int prefix = indices[0];
    for (std::size_t i = 1; i < indices.size(); ++i) {
        const int c = indices[i];
        const std::uint32_t key = (static_cast<std::uint32_t>(prefix) << 8) | static_cast<std::uint32_t>(c);
        const auto found = dictionary.find(key);
        if (found != dictionary.end()) { prefix = found->second; continue; }
        writer.write(prefix, codeSize);
        if (next < 4096) {
            dictionary[key] = next++;
            // The decoder widens codes one entry later than it is added here.
            if (next > (1 << codeSize) && codeSize < 12) ++codeSize;
        } else {
            writer.write(clear, codeSize);
            dictionary.clear();
            codeSize = 9;
            next = 258;
        }
        prefix = c;
    }
    writer.write(prefix, codeSize);
    writer.write(end, codeSize);
    writer.finish();
}
}

GifWriter::GifWriter(const std::filesystem::path& path, int width, int height, int delayCentiseconds)
    : m_out(path, std::ios::binary), m_width(width), m_height(height), m_delay(delayCentiseconds)
{
    if (!m_out) throw std::runtime_error("Cannot write " + path.u8string());
    if (width < 1 || height < 1 || width > 65535 || height > 65535) throw std::runtime_error("GIF size must be 1-65535 pixels");
    m_out.write("GIF89a", 6);
    put16(m_out, width);
    put16(m_out, height);
    m_out.put(0);  // no global colour table
    m_out.put(0);
    m_out.put(0);
    // Loop forever (NETSCAPE2.0 application extension).
    const unsigned char loop[] = {0x21, 0xFF, 0x0B, 'N', 'E', 'T', 'S', 'C', 'A', 'P', 'E', '2', '.', '0', 0x03, 0x01, 0x00, 0x00, 0x00};
    m_out.write(reinterpret_cast<const char*>(loop), sizeof(loop));
}

GifWriter::~GifWriter()
{
    try { close(); } catch (...) {}
}

void GifWriter::addFrame(const std::vector<unsigned char>& pixels, int channels)
{
    if (m_closed) throw std::runtime_error("GIF already closed");
    const std::size_t count = static_cast<std::size_t>(m_width) * static_cast<std::size_t>(m_height);
    if ((channels != 3 && channels != 4) || pixels.size() < count * static_cast<std::size_t>(channels))
        throw std::runtime_error("GIF frame has the wrong size");
    // Popularity palette on a 15-bit colour grid; each entry is its bin's mean colour.
    std::vector<std::uint32_t> histogram(32768, 0);
    std::vector<std::array<std::uint64_t, 3>> sums(32768, {0, 0, 0});
    std::vector<std::uint16_t> keys(count);
    for (std::size_t i = 0; i < count; ++i) {
        const unsigned char* p = &pixels[i * static_cast<std::size_t>(channels)];
        const std::uint16_t key = static_cast<std::uint16_t>(((p[0] >> 3) << 10) | ((p[1] >> 3) << 5) | (p[2] >> 3));
        keys[i] = key;
        ++histogram[key];
        for (int k = 0; k < 3; ++k) sums[key][static_cast<std::size_t>(k)] += p[k];
    }
    std::vector<int> bins;
    for (int key = 0; key < 32768; ++key) if (histogram[static_cast<std::size_t>(key)]) bins.push_back(key);
    std::sort(bins.begin(), bins.end(), [&](int a, int b) { return histogram[static_cast<std::size_t>(a)] > histogram[static_cast<std::size_t>(b)]; });
    if (bins.size() > 256) bins.resize(256);
    std::array<std::array<unsigned char, 3>, 256> palette{};
    for (std::size_t i = 0; i < bins.size(); ++i)
        for (int k = 0; k < 3; ++k)
            palette[i][static_cast<std::size_t>(k)] = static_cast<unsigned char>(sums[static_cast<std::size_t>(bins[i])][static_cast<std::size_t>(k)] / histogram[static_cast<std::size_t>(bins[i])]);
    std::vector<int> mapping(32768, -1);
    std::vector<unsigned char> indices(count);
    for (std::size_t i = 0; i < count; ++i) {
        int& index = mapping[keys[i]];
        if (index < 0) {
            const unsigned char* p = &pixels[i * static_cast<std::size_t>(channels)];
            int best = 0;
            long bestDistance = -1;
            for (std::size_t c = 0; c < bins.size(); ++c) {
                long distance = 0;
                for (int k = 0; k < 3; ++k) { const long d = static_cast<long>(p[k]) - palette[c][static_cast<std::size_t>(k)]; distance += d * d; }
                if (bestDistance < 0 || distance < bestDistance) { bestDistance = distance; best = static_cast<int>(c); }
            }
            index = best;
        }
        indices[i] = static_cast<unsigned char>(index);
    }
    // Graphic control extension: delay, no transparency.
    const unsigned char control[] = {0x21, 0xF9, 0x04, 0x04};
    m_out.write(reinterpret_cast<const char*>(control), sizeof(control));
    put16(m_out, m_delay);
    m_out.put(0);
    m_out.put(0);
    // Image descriptor with a 256-entry local colour table.
    m_out.put(0x2C);
    put16(m_out, 0);
    put16(m_out, 0);
    put16(m_out, m_width);
    put16(m_out, m_height);
    m_out.put(static_cast<char>(0x87));
    for (const auto& colour : palette) m_out.write(reinterpret_cast<const char*>(colour.data()), 3);
    encode(m_out, indices);
    ++m_frames;
    if (!m_out) throw std::runtime_error("Writing the GIF failed");
}

void GifWriter::close()
{
    if (m_closed) return;
    m_closed = true;
    m_out.put(0x3B);
    m_out.close();
}

GifImage readGif(const std::filesystem::path& path)
{
    std::ifstream in(path, std::ios::binary);
    std::vector<unsigned char> data((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    if (data.size() < 13 || std::memcmp(data.data(), "GIF8", 4) != 0) throw std::runtime_error("Not a GIF file");
    GifImage image;
    image.width = data[6] | (data[7] << 8);
    image.height = data[8] | (data[9] << 8);
    std::size_t pos = 13;
    std::vector<std::array<unsigned char, 3>> global;
    if (data[10] & 0x80) {
        const std::size_t size = std::size_t(2) << (data[10] & 7);
        for (std::size_t i = 0; i < size; ++i) global.push_back({data[pos + 3 * i], data[pos + 3 * i + 1], data[pos + 3 * i + 2]});
        pos += 3 * size;
    }
    int delay = 0;
    auto subBlocks = [&](std::vector<unsigned char>* out) {
        while (pos < data.size()) {
            const std::size_t length = data[pos++];
            if (length == 0) return;
            if (out) out->insert(out->end(), data.begin() + static_cast<std::ptrdiff_t>(pos), data.begin() + static_cast<std::ptrdiff_t>(pos + length));
            pos += length;
        }
    };
    while (pos < data.size()) {
        const unsigned char tag = data[pos++];
        if (tag == 0x3B) break;
        if (tag == 0x21) {
            const unsigned char label = data[pos++];
            if (label == 0xF9) { delay = data[pos + 2] | (data[pos + 3] << 8); }
            if (label == 0xFF && pos + 12 <= data.size() && std::memcmp(&data[pos + 1], "NETSCAPE2.0", 11) == 0) image.looping = true;
            subBlocks(nullptr);
            continue;
        }
        if (tag != 0x2C) throw std::runtime_error("Unexpected GIF block");
        const int w = data[pos + 4] | (data[pos + 5] << 8), h = data[pos + 6] | (data[pos + 7] << 8);
        const unsigned char flags = data[pos + 8];
        pos += 9;
        if (flags & 0x40) throw std::runtime_error("Interlaced GIF frames are not supported");
        std::vector<std::array<unsigned char, 3>> palette = global;
        if (flags & 0x80) {
            palette.clear();
            const std::size_t size = std::size_t(2) << (flags & 7);
            for (std::size_t i = 0; i < size; ++i) palette.push_back({data[pos + 3 * i], data[pos + 3 * i + 1], data[pos + 3 * i + 2]});
            pos += 3 * size;
        }
        const int minimum = data[pos++];
        std::vector<unsigned char> stream;
        subBlocks(&stream);
        // LZW decode.
        const int clear = 1 << minimum, end = clear + 1;
        std::vector<std::vector<unsigned char>> table;
        auto reset = [&] { table.assign(static_cast<std::size_t>(end + 1), {}); for (int i = 0; i < clear; ++i) table[static_cast<std::size_t>(i)] = {static_cast<unsigned char>(i)}; };
        reset();
        int codeSize = minimum + 1;
        std::vector<unsigned char> indices;
        std::size_t bit = 0;
        int previous = -1;
        while (bit + static_cast<std::size_t>(codeSize) <= stream.size() * 8) {
            int code = 0;
            for (int b = 0; b < codeSize; ++b, ++bit) code |= ((stream[bit / 8] >> (bit % 8)) & 1) << b;
            if (code == clear) { reset(); codeSize = minimum + 1; previous = -1; continue; }
            if (code == end) break;
            std::vector<unsigned char> entry;
            if (code < static_cast<int>(table.size()) && !table[static_cast<std::size_t>(code)].empty()) entry = table[static_cast<std::size_t>(code)];
            else if (previous >= 0) { entry = table[static_cast<std::size_t>(previous)]; entry.push_back(entry[0]); }
            else throw std::runtime_error("Corrupt GIF LZW stream");
            indices.insert(indices.end(), entry.begin(), entry.end());
            if (previous >= 0 && table.size() < 4096) {
                auto added = table[static_cast<std::size_t>(previous)];
                added.push_back(entry[0]);
                table.push_back(added);
                if (table.size() >= (std::size_t(1) << codeSize) && codeSize < 12) ++codeSize;
            }
            previous = code;
        }
        if (indices.size() < static_cast<std::size_t>(w) * static_cast<std::size_t>(h)) throw std::runtime_error("Truncated GIF frame");
        std::vector<unsigned char> rgb(static_cast<std::size_t>(w) * static_cast<std::size_t>(h) * 3);
        for (std::size_t i = 0; i < static_cast<std::size_t>(w) * static_cast<std::size_t>(h); ++i)
            for (int k = 0; k < 3; ++k) rgb[3 * i + static_cast<std::size_t>(k)] = palette[indices[i]][static_cast<std::size_t>(k)];
        image.frames.push_back(rgb);
        image.delays.push_back(delay);
    }
    return image;
}
}
