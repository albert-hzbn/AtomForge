// Reference GIFs of the C++ GifWriter for the Rust af-science gif tests:
// synthetic frames with many colours (palette truncation and ties).
// Usage (from the repository root): golden_gif <output dir>.

#include "science/GifWriter.h"

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

// The frames of rust/crates/af-science/tests/golden_gif.rs.
std::vector<unsigned char> frame(int width, int height, int index, int channels)
{
    std::vector<unsigned char> pixels(static_cast<std::size_t>(width) * height * channels);
    std::uint32_t state = 12345u + static_cast<std::uint32_t>(index) * 977u;
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x) {
            state = state * 1664525u + 1013904223u;
            unsigned char* p = &pixels[(static_cast<std::size_t>(y) * width + x) * channels];
            p[0] = static_cast<unsigned char>((x * 255) / (width - 1));
            p[1] = static_cast<unsigned char>((y * 255) / (height - 1));
            p[2] = static_cast<unsigned char>((state >> 24) & (index == 2 ? 0x0F : 0xFF));
            if (channels == 4) p[3] = 255;
        }
    return pixels;
}

int main(int argc, char** argv)
{
    if (argc < 2) { std::fprintf(stderr, "usage: golden_gif <output dir>\n"); return 1; }
    const std::string out = argv[1];
    {
        atomforge::science::GifWriter gif(out + "/many_colours.gif", 64, 48, 8);
        for (int i = 0; i < 3; ++i) gif.addFrame(frame(64, 48, i, 4), 4);
        gif.close();
    }
    {
        atomforge::science::GifWriter gif(out + "/rgb.gif", 17, 5, 3);
        gif.addFrame(frame(17, 5, 2, 3), 3);
        gif.close();
    }
    std::printf("wrote many_colours.gif and rgb.gif\n");
    return 0;
}
