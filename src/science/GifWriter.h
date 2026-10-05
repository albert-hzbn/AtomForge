#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <vector>

namespace atomforge::science
{
// Animated GIF89a writer (looping). Each frame gets its own 256-colour
// palette chosen by colour popularity on a 5-bit-per-channel grid, which
// suits rendered scenes; pixels map to the nearest palette colour.
class GifWriter
{
public:
    GifWriter(const std::filesystem::path& path, int width, int height, int delayCentiseconds);
    ~GifWriter();
    // RGBA (or RGB with channels = 3), rows top to bottom.
    void addFrame(const std::vector<unsigned char>& pixels, int channels = 4);
    void close();
    int frames() const { return m_frames; }

private:
    std::ofstream m_out;
    int m_width, m_height, m_delay, m_frames = 0;
    bool m_closed = false;
};

// Decodes a GIF written by GifWriter (or any GIF without interlacing) into
// RGB frames; used to verify encoding.
struct GifImage
{
    int width = 0, height = 0;
    std::vector<std::vector<unsigned char>> frames;  // RGB
    std::vector<int> delays;                         // centiseconds
    bool looping = false;
};
GifImage readGif(const std::filesystem::path& path);
}
