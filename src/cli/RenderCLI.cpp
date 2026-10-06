#include "cli/RenderCLI.h"
#include "cli/CliArgs.h"

#include "app/EditorOps.h"
#include "app/ImageExport.h"
#include "app/SceneView.h"
#include "camera/Camera.h"
#include "graphics/BillboardMesh.h"
#include "graphics/CylinderMesh.h"
#include "graphics/LowPolyMesh.h"
#include "graphics/Renderer.h"
#include "graphics/SceneBuffers.h"
#include "graphics/ShadowMap.h"
#include "graphics/SphereMesh.h"
#include "graphics/StructureInstanceBuilder.h"
#include "io/StructureLoader.h"
#include "model/Structure.h"
#include "science/GifWriter.h"
#include "science/ScienceData.h"
#include "ui/FileBrowser.h"
#include "util/ElementData.h"
#include "util/PathUtils.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cctype>
#include <cstring>
#include <filesystem>
#include <memory>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

using namespace cli;

namespace
{

// Apply "SYMBOL r g b" overrides (0..1 components) onto a size-119 element
// color table. Pulled out as a free function so flag semantics are testable
// without an OpenGL context.
void applyColorOverrides(const std::vector<std::string>& specs, std::vector<glm::vec3>& colors)
{
    for (const auto& spec : specs)
    {
        std::istringstream iss(spec);
        std::string symbol;
        float r = 0.0f, g = 0.0f, b = 0.0f;
        if (!(iss >> symbol >> r >> g >> b))
            throw std::invalid_argument("Cannot parse --color value '" + spec + "'. Expected: \"SYMBOL r g b\"");
        const int z = atomicNumberFromSymbol(symbol);
        if (z <= 0)
            throw std::invalid_argument("Unknown element symbol in --color: '" + symbol + "'");
        colors[static_cast<std::size_t>(z)] = glm::vec3(r, g, b);
    }
}

// Apply "SYMBOL radiusAngstrom" overrides onto a size-119 covalent radius table.
void applyRadiusOverrides(const std::vector<std::string>& specs, std::vector<float>& radii)
{
    for (const auto& spec : specs)
    {
        std::istringstream iss(spec);
        std::string symbol;
        float radius = 0.0f;
        if (!(iss >> symbol >> radius) || !std::isfinite(radius) || radius <= 0.0f)
            throw std::invalid_argument("Cannot parse --radius value '" + spec + "'. Expected: \"SYMBOL radiusAngstrom\"");
        const int z = atomicNumberFromSymbol(symbol);
        if (z <= 0)
            throw std::invalid_argument("Unknown element symbol in --radius: '" + symbol + "'");
        radii[static_cast<std::size_t>(z)] = radius;
    }
}

// Standard PNG chunk CRC-32 (polynomial 0xEDB88320, reflected).
uint32_t pngCrc32(const unsigned char* data, std::size_t length)
{
    uint32_t crc = 0xFFFFFFFFu;
    for (std::size_t i = 0; i < length; ++i)
    {
        crc ^= data[i];
        for (int bit = 0; bit < 8; ++bit)
            crc = (crc & 1u) ? (crc >> 1) ^ 0xEDB88320u : (crc >> 1);
    }
    return crc ^ 0xFFFFFFFFu;
}

void appendBigEndian32(std::vector<unsigned char>& buffer, uint32_t value)
{
    buffer.push_back(static_cast<unsigned char>((value >> 24) & 0xFF));
    buffer.push_back(static_cast<unsigned char>((value >> 16) & 0xFF));
    buffer.push_back(static_cast<unsigned char>((value >> 8) & 0xFF));
    buffer.push_back(static_cast<unsigned char>(value & 0xFF));
}

// stb_image_write has no facility for embedding physical resolution, so a
// pHYs chunk (pixels-per-metre, PNG's DPI equivalent) is spliced into the
// already-written file immediately after the mandatory IHDR chunk.
bool embedPngDpi(const std::string& path, double dpi, std::string& errorMessage)
{
    std::ifstream in(path, std::ios::binary);
    std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    in.close();

    static const unsigned char signature[8] = {0x89, 'P', 'N', 'G', '\r', '\n', 0x1A, '\n'};
    const std::size_t afterIhdr = 8 + 25; // signature + (length+type+13-byte data+crc)
    if (bytes.size() < afterIhdr || !std::equal(signature, signature + 8, bytes.begin())
        || bytes[12] != 'I' || bytes[13] != 'H' || bytes[14] != 'D' || bytes[15] != 'R')
    {
        errorMessage = "Unexpected PNG layout; cannot embed DPI metadata.";
        return false;
    }

    const uint32_t pixelsPerMetre = static_cast<uint32_t>(std::lround(dpi / 0.0254));

    std::vector<unsigned char> typeAndData;
    typeAndData.push_back('p'); typeAndData.push_back('H'); typeAndData.push_back('Y'); typeAndData.push_back('s');
    appendBigEndian32(typeAndData, pixelsPerMetre);
    appendBigEndian32(typeAndData, pixelsPerMetre);
    typeAndData.push_back(1); // unit specifier: metre

    std::vector<unsigned char> chunk;
    appendBigEndian32(chunk, 9); // pHYs data length
    chunk.insert(chunk.end(), typeAndData.begin(), typeAndData.end());
    appendBigEndian32(chunk, pngCrc32(typeAndData.data(), typeAndData.size()));

    bytes.insert(bytes.begin() + static_cast<std::ptrdiff_t>(afterIhdr), chunk.begin(), chunk.end());

    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out)
    {
        errorMessage = "Failed to rewrite '" + path + "' with DPI metadata.";
        return false;
    }
    out.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    if (!out.good())
    {
        errorMessage = "Failed to rewrite '" + path + "' with DPI metadata.";
        return false;
    }
    return true;
}

void printHelp()
{
    std::cout <<
"AtomForge --render - headless structure snapshot (no GUI required)\n"
"\n"
"Usage:\n"
"  AtomForge --render --input FILE --output FILE.png [options]\n"
"\n"
"  --width N              Image width in pixels          (default: 1600)\n"
"  --height N             Image height in pixels          (default: 1200)\n"
"  --yaw D                Camera yaw, degrees              (default: 0)\n"
"  --pitch D              Camera pitch, degrees            (default: 0)\n"
"  --roll D               Camera roll, degrees             (default: 0)\n"
"  --distance A           Camera distance, Angstrom       (default: auto-fit)\n"
"  --orthographic         Use an orthographic projection  (default: perspective)\n"
"  --background R G B     Background color, 0..1          (default: 1 1 1)\n"
"  --no-bonds             Hide bonds\n"
"  --no-box               Hide the unit-cell/bounding box\n"
"  --color SYMBOL R G B   Override one element's color, 0..1 (repeatable)\n"
"  --radius SYMBOL VALUE  Override one element's covalent radius, Angstrom (repeatable)\n"
"  --radius-scale FACTOR  Global atom-size multiplier      (default: 1.0)\n"
"  --dpi N                Embed a physical resolution (dots per inch) in the\n"
"                         saved PNG's metadata; does not change pixel size\n"
"                         (default: none written, matching plain stb PNGs)\n"
"  --frames N             Turntable: render N frames (writes output-000.png, ...)\n"
"  --yaw-step D           Turntable: yaw increment per frame, degrees\n"
"  --fps N                Animated GIF frame rate          (default: 12)\n"
"\n"
"Animated GIF: with --output FILE.gif a trajectory input (XYZ/extXYZ, XDATCAR,\n"
"LAMMPS dump) becomes one frame per trajectory frame, using every Nth frame with\n"
"--every N. A single structure becomes a looping turntable of --frames frames\n"
"(default 36) turning by --yaw-step degrees (default 360/frames).\n"
"\n"
"Example:\n"
"  AtomForge --render --input cu_fcc.cif --output cu.png ^\n"
"            --yaw 30 --pitch 20 --color Cu 0.9 0.5 0.2 --output cu_view.png\n"
<< std::endl;
}

GLFWwindow* createHeadlessContext()
{
    if (!glfwInit())
        return nullptr;

    const int versions[][2] = {{4, 6}, {4, 3}, {3, 3}};
    for (const auto& version : versions)
    {
        glfwDefaultWindowHints();
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, version[0]);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, version[1]);
        glfwWindowHint(GLFW_OPENGL_PROFILE, version[0] >= 3 ? GLFW_OPENGL_CORE_PROFILE : GLFW_OPENGL_ANY_PROFILE);
        glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
        if (GLFWwindow* window = glfwCreateWindow(64, 64, "AtomForge (headless)", nullptr, nullptr))
            return window;
    }
    return nullptr;
}

} // namespace

int runRenderCLI(int argc, char* argv[])
{
    if (hasFlag(argc, argv, "--help") || hasFlag(argc, argv, "-h"))
    {
        printHelp();
        return 0;
    }

    try
    {
        const char* inputPath = findArg(argc, argv, "--input");
        const char* outputPath = findArg(argc, argv, "--output");
        if (!inputPath || !outputPath)
            throw std::invalid_argument("--input FILE and --output FILE.png are required");

        const int width = argInt(argc, argv, "--width", 1600);
        const int height = argInt(argc, argv, "--height", 1200);
        if (width <= 0 || height <= 0)
            throw std::invalid_argument("--width/--height must be positive");

        const float yaw = static_cast<float>(argDouble(argc, argv, "--yaw", 0.0));
        const float pitch = static_cast<float>(argDouble(argc, argv, "--pitch", 0.0));
        const float roll = static_cast<float>(argDouble(argc, argv, "--roll", 0.0));
        const bool hasExplicitDistance = findArg(argc, argv, "--distance") != nullptr;
        const float distance = static_cast<float>(argDouble(argc, argv, "--distance", 0.0));
        const bool orthographic = hasFlag(argc, argv, "--orthographic");
        const bool showBonds = !hasFlag(argc, argv, "--no-bonds");
        const bool showBox = !hasFlag(argc, argv, "--no-box");
        const float radiusScale = static_cast<float>(argDouble(argc, argv, "--radius-scale", 1.0));
        if (!(radiusScale > 0.0f))
            throw std::invalid_argument("--radius-scale must be positive");

        const bool hasDpi = findArg(argc, argv, "--dpi") != nullptr;
        const double dpi = argDouble(argc, argv, "--dpi", 0.0);
        if (hasDpi && !(dpi > 0.0))
            throw std::invalid_argument("--dpi must be positive");

        glm::vec4 background(1.0f, 1.0f, 1.0f, 1.0f);
        if (const char* bg = findArg(argc, argv, "--background"))
        {
            std::istringstream iss(bg);
            float r = 1.0f, g = 1.0f, b = 1.0f;
            if (!(iss >> r >> g >> b))
                throw std::invalid_argument("Cannot parse --background \"r g b\"");
            background = glm::vec4(r, g, b, 1.0f);
        }

        std::string outputLower(outputPath);
        std::transform(outputLower.begin(), outputLower.end(), outputLower.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        const bool gifOutput = outputLower.size() > 4 && outputLower.compare(outputLower.size() - 4, 4, ".gif") == 0;
        const double fps = argDouble(argc, argv, "--fps", 12.0);
        if (!(fps > 0.0) || fps > 100.0)
            throw std::invalid_argument("--fps must be between 0 and 100");
        const int every = argInt(argc, argv, "--every", 1);
        if (every <= 0)
            throw std::invalid_argument("--every must be positive");

        // ---- Load structure (or trajectory frames for a GIF) ----
        std::vector<Structure> structures;
        if (gifOutput)
        {
            // Trajectory formats first; anything else is read as one structure.
            try
            {
                const auto frames = atomforge::science::readFrames(std::filesystem::u8path(inputPath));
                for (std::size_t i = 0; i < frames.size(); i += static_cast<std::size_t>(every))
                    structures.push_back(frames[i].structure);
            }
            catch (const std::exception&)
            {
                structures.clear();
            }
        }
        if (structures.empty())
        {
            Structure loaded;
            std::string loadError;
            if (!loadStructureFromFile(inputPath, loaded, loadError))
                throw std::runtime_error("Error loading input structure: " + loadError);
            structures.push_back(std::move(loaded));
        }
        if (structures.empty() || structures.front().atoms.empty())
            throw std::runtime_error("Error loading input structure: no atoms found");
        const bool trajectoryGif = gifOutput && structures.size() > 1;

        const bool hasFrames = findArg(argc, argv, "--frames") != nullptr;
        const int frameCount = trajectoryGif ? static_cast<int>(structures.size())
                                             : argInt(argc, argv, "--frames", gifOutput ? 36 : 1);
        if (frameCount <= 0 || (trajectoryGif && hasFrames))
            throw std::invalid_argument(trajectoryGif ? "--frames does not apply to trajectory GIFs; use --every"
                                                      : "--frames must be positive");
        const float yawStep = static_cast<float>(argDouble(argc, argv, "--yaw-step",
            gifOutput && !trajectoryGif ? 360.0 / frameCount : 0.0));

        std::vector<glm::vec3> elementColors = makeDefaultElementColors();
        std::vector<float> elementRadii = makeLiteratureCovalentRadii();
        applyColorOverrides(findAllArgs(argc, argv, "--color"), elementColors);
        applyRadiusOverrides(findAllArgs(argc, argv, "--radius"), elementRadii);
        for (float& radius : elementRadii)
            radius *= radiusScale;

        for (auto& structure : structures)
            for (auto& atom : structure.atoms)
            {
                if (atom.atomicNumber >= 0 && atom.atomicNumber < (int)elementColors.size())
                {
                    const glm::vec3& color = elementColors[static_cast<std::size_t>(atom.atomicNumber)];
                    atom.r = color.r; atom.g = color.g; atom.b = color.b;
                }
            }
        const Structure& structure = structures.front();

        // ---- Headless OpenGL bootstrap ----
        GLFWwindow* window = createHeadlessContext();
        if (!window)
        {
            std::cerr << "Error: cannot create an OpenGL context; --render needs a display/GPU driver.\n";
            glfwTerminate();
            return 1;
        }
        glfwMakeContextCurrent(window);
        glewExperimental = GL_TRUE;
        if (glewInit() != GLEW_OK)
        {
            std::cerr << "Error: failed to initialize OpenGL function loading (GLEW).\n";
            glfwDestroyWindow(window);
            glfwTerminate();
            return 1;
        }
        while (glGetError() != GL_NO_ERROR) {}

        {
            SphereMesh sphere(40, 40);
            LowPolyMesh lowPolyMesh;
            BillboardMesh billboardMesh;
            CylinderMesh cylinder(32);

            Renderer renderer;
            renderer.init();
            ShadowMap shadow = createShadowMap(2048, 2048);

            SceneBuffers sceneBuffers;
            sceneBuffers.init(
                sphere.vbo, sphere.ebo, sphere.indexCount,
                lowPolyMesh.vbo, lowPolyMesh.ebo, lowPolyMesh.indexCount,
                billboardMesh.vbo, billboardMesh.ebo, billboardMesh.indexCount,
                cylinder.vbo, cylinder.vertexCount);

            const int identity[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
            const std::vector<float> elementShininess(119, 32.0f);
            const StructureInstanceData instanceData = buildStructureInstanceData(
                structure, false, identity, elementRadii, elementShininess);
            sceneBuffers.upload(instanceData, false, {});

            Camera camera;
            camera.yaw = yaw;
            camera.pitch = pitch;
            camera.roll = roll;
            if (hasExplicitDistance)
                camera.distance = distance;
            else
                applyDefaultView(camera, sceneBuffers, width, height, true);

            ImageExportRequest request;
            request.format = ImageExportFormat::Png;
            request.includeBackground = true;
            request.resolutionScale = 1;
            request.includeGizmo = false;

            std::unique_ptr<atomforge::science::GifWriter> gif;
            if (gifOutput)
                gif = std::make_unique<atomforge::science::GifWriter>(
                    std::filesystem::u8path(outputPath), width, height,
                    std::max(2, static_cast<int>(std::lround(100.0 / fps))));

            for (int frame = 0; frame < frameCount; ++frame)
            {
                if (trajectoryGif && frame > 0)
                {
                    const StructureInstanceData frameData = buildStructureInstanceData(
                        structures[static_cast<std::size_t>(frame)], false, identity, elementRadii, elementShininess);
                    sceneBuffers.upload(frameData, false, {});
                }
                FrameView frameView;
                frameView.framebufferWidth = width;
                frameView.framebufferHeight = height;
                buildFrameView(camera, sceneBuffers, orthographic, frameView);

                ImageExportView view;
                view.width = width;
                view.height = height;
                view.projection = frameView.projection;
                view.view = frameView.view;
                view.lightMVP = frameView.lightMVP;
                view.lightPosition = frameView.lightPosition;
                view.cameraPosition = frameView.cameraPosition;

                if (gif)
                {
                    std::vector<unsigned char> pixels;
                    std::string renderError;
                    if (!renderSceneToRgba(view, background, showBonds, true, showBox,
                                           sceneBuffers, renderer, shadow, pixels, renderError))
                        throw std::runtime_error("Error rendering frame: " + renderError);
                    gif->addFrame(pixels);
                    camera.yaw += yawStep;
                    continue;
                }

                if (frameCount > 1)
                {
                    char suffix[32];
                    std::snprintf(suffix, sizeof(suffix), "-%03d", frame);
                    const std::string path(outputPath);
                    const auto dot = path.rfind('.');
                    request.outputPath = (dot == std::string::npos)
                        ? path + suffix
                        : path.substr(0, dot) + suffix + path.substr(dot);
                }
                else
                {
                    request.outputPath = outputPath;
                }

                std::string exportError;
                if (!exportStructureImage(request, view, background, showBonds, true, showBox,
                                          sceneBuffers, renderer, shadow, structure, exportError))
                {
                    throw std::runtime_error("Error writing image: " + exportError);
                }
                if (hasDpi)
                {
                    std::string dpiError;
                    if (!embedPngDpi(request.outputPath, dpi, dpiError))
                        throw std::runtime_error(dpiError);
                }
                std::cout << "Saved to: " << request.outputPath << "\n";

                camera.yaw += yawStep;
            }

            if (gif)
            {
                gif->close();
                std::cout << "Saved " << gif->frames() << " frames to: " << outputPath << "\n";
            }
            sceneBuffers.destroy();
        }

        glfwDestroyWindow(window);
        glfwTerminate();
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
