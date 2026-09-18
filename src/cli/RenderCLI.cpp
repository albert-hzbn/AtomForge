#include "cli/RenderCLI.h"

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
#include "ui/FileBrowser.h"
#include "util/ElementData.h"
#include "util/PathUtils.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace
{

const char* findArg(int argc, char* argv[], const char* flag)
{
    for (int i = 1; i < argc - 1; ++i)
        if (std::strcmp(argv[i], flag) == 0)
            return argv[i + 1];
    return nullptr;
}

bool hasFlag(int argc, char* argv[], const char* flag)
{
    for (int i = 1; i < argc; ++i)
        if (std::strcmp(argv[i], flag) == 0)
            return true;
    return false;
}

std::vector<std::string> findAllArgs(int argc, char* argv[], const char* flag)
{
    std::vector<std::string> out;
    for (int i = 1; i < argc - 1; ++i)
        if (std::strcmp(argv[i], flag) == 0)
            out.emplace_back(argv[i + 1]);
    return out;
}

double argDouble(int argc, char* argv[], const char* flag, double def)
{
    const char* v = findArg(argc, argv, flag);
    if (!v) return def;
    std::size_t parsed = 0;
    const std::string text(v);
    const double value = std::stod(text, &parsed);
    if (parsed != text.size() || !std::isfinite(value))
        throw std::invalid_argument(std::string("invalid numeric value '") + v + "' for " + flag);
    return value;
}

int argInt(int argc, char* argv[], const char* flag, int def)
{
    const char* v = findArg(argc, argv, flag);
    if (!v) return def;
    std::size_t parsed = 0;
    const std::string text(v);
    const int value = std::stoi(text, &parsed);
    if (parsed != text.size())
        throw std::invalid_argument(std::string("invalid integer value '") + v + "' for " + flag);
    return value;
}

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
"  --frames N             Turntable: render N frames (writes output-000.png, ...)\n"
"  --yaw-step D           Turntable: yaw increment per frame, degrees\n"
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

        glm::vec4 background(1.0f, 1.0f, 1.0f, 1.0f);
        if (const char* bg = findArg(argc, argv, "--background"))
        {
            std::istringstream iss(bg);
            float r = 1.0f, g = 1.0f, b = 1.0f;
            if (!(iss >> r >> g >> b))
                throw std::invalid_argument("Cannot parse --background \"r g b\"");
            background = glm::vec4(r, g, b, 1.0f);
        }

        const int frameCount = argInt(argc, argv, "--frames", 1);
        const float yawStep = static_cast<float>(argDouble(argc, argv, "--yaw-step", 0.0));
        if (frameCount <= 0)
            throw std::invalid_argument("--frames must be positive");

        // ---- Load structure and resolve per-element appearance ----
        Structure structure;
        std::string loadError;
        if (!loadStructureFromFile(inputPath, structure, loadError))
            throw std::runtime_error("Error loading input structure: " + loadError);

        std::vector<glm::vec3> elementColors = makeDefaultElementColors();
        std::vector<float> elementRadii = makeLiteratureCovalentRadii();
        applyColorOverrides(findAllArgs(argc, argv, "--color"), elementColors);
        applyRadiusOverrides(findAllArgs(argc, argv, "--radius"), elementRadii);
        for (float& radius : elementRadii)
            radius *= radiusScale;

        for (auto& atom : structure.atoms)
        {
            if (atom.atomicNumber >= 0 && atom.atomicNumber < (int)elementColors.size())
            {
                const glm::vec3& color = elementColors[static_cast<std::size_t>(atom.atomicNumber)];
                atom.r = color.r; atom.g = color.g; atom.b = color.b;
            }
        }

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

            for (int frame = 0; frame < frameCount; ++frame)
            {
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
                std::cout << "Saved to: " << request.outputPath << "\n";

                camera.yaw += yawStep;
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
