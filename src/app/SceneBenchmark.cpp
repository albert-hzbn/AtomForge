#include "app/SceneBenchmark.h"

#include "camera/Camera.h"
#include "graphics/FrameCapture.h"
#include "graphics/GraphicsBackend.h"
#include "graphics/SceneBuffers.h"
#ifdef ATOMFORGE_VULKAN
#include "graphics/vulkan/VulkanSceneRenderer.h"
#endif

#include "third_party/stb_image_write.h"

#include <GLFW/glfw3.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

namespace
{
struct View { const char* name; float distance; bool orbit; };
// Fractions of the fitted distance; the last fixed view is close to a corner.
const View kViews[] = {{"overview", 1.0f, false}, {"half distance", 0.6f, false}, {"close-up", 0.4f, false},
                       {"orbiting overview", 1.0f, true}};
constexpr int kWarmup = 40;
constexpr int kMeasured = 150;
}

std::unique_ptr<SceneBenchmark> SceneBenchmark::fromEnvironment()
{
    const char* report = std::getenv("ATOMFORGE_SCENE_BENCHMARK");
    if (!report || !*report) return nullptr;
    auto benchmark = std::make_unique<SceneBenchmark>();
    benchmark->m_report = report;
    return benchmark;
}

void SceneBenchmark::beforeFrame(Camera& camera, const SceneBuffers& buffers)
{
    if (m_done) return;
    if (m_distance == 0.0f)
    {
        // The structure and its fitted view first.
        if (buffers.atomCount == 0 || ++m_wait < 20) return;
        graphics::setVsync(false);
        m_distance = camera.distance;
        m_yaw = camera.yaw;
        m_pitch = camera.pitch;
        m_frame = 0;
    }
    const View& view = kViews[m_view];
    camera.distance = m_distance * view.distance;
    camera.pitch = m_pitch;
    camera.yaw = view.orbit ? m_yaw + 2.0f * static_cast<float>(m_frame) : m_yaw;
}

void SceneBenchmark::afterFrame(const SceneBuffers& buffers)
{
    if (m_done || m_distance == 0.0f) return;
    ++m_frame;
    double gpu = -1.0;
    std::string culling;
#ifdef ATOMFORGE_VULKAN
    if (graphics::usingVulkan())
    {
        const vkr::SceneStats& stats = vkr::sceneRenderer().stats();
        if (stats.valid)
        {
            gpu = stats.gpuMilliseconds;
            uint64_t drawn = 0;
            char levels[160];
            std::snprintf(levels, sizeof(levels), "[%u, %u, %u, %u, %u, %u]", stats.atomLevels[0], stats.atomLevels[1], stats.atomLevels[2],
                          stats.atomLevels[3], stats.atomLevels[4], stats.atomLevels[5]);
            for (uint32_t n : stats.atomLevels) drawn += n;
            char text[600];
            std::snprintf(text, sizeof(text), ", \"drawn_atoms\": %llu, \"atom_levels\": %s, \"drawn_bonds\": %u, \"triangles\": %llu, \"shadow_atoms\": %u, \"first_phase\": %u, \"occlusion\": %s, \"stages_ms\": [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f]",
                          static_cast<unsigned long long>(drawn), levels, stats.bondLevels[0] + stats.bondLevels[1] + stats.bondLevels[2],
                          static_cast<unsigned long long>(stats.triangles), stats.shadowAtoms, stats.firstPhaseAtoms,
                          stats.occlusion ? "true" : "false", stats.stageMilliseconds[0], stats.stageMilliseconds[1],
                          stats.stageMilliseconds[2], stats.stageMilliseconds[3], stats.stageMilliseconds[4], stats.stageMilliseconds[5],
                          stats.stageMilliseconds[6]);
            culling = text;
        }
    }
#endif
    if (m_frame == kWarmup)
    {
        graphics::waitIdle();
        m_start = glfwGetTime();
        m_gpuSum = 0;
        return;
    }
    if (m_frame <= kWarmup) return;
    m_gpuSum += gpu;
    if (m_frame < kWarmup + kMeasured) return;
    graphics::waitIdle();
    const double seconds = glfwGetTime() - m_start;
    if (const char* shots = std::getenv("ATOMFORGE_SCENE_SHOTS"); shots && *shots)
    {
        // The view's last frame (ATOMFORGE_SCENE_SHOTS=FOLDER; its copies make frames slower).
        std::vector<unsigned char> rgb;
        int w = 0, h = 0;
        if (graphics::lastFrame(rgb, w, h))
            stbi_write_png((std::string(shots) + "/" + std::to_string(m_view + 1) + ".png").c_str(), w, h, 3, rgb.data(), w * 3);
    }
    const double frameMs = 1000.0 * seconds / kMeasured;
    char entry[1024];
    std::snprintf(entry, sizeof(entry), "  {\"view\": \"%s\", \"renderer\": \"%s\", \"atoms\": %zu, \"frames\": %d, \"frame_ms\": %.3f, \"fps\": %.1f, \"gpu_ms\": %.3f%s}",
                  kViews[m_view].name, graphics::backendName(graphics::backend()), buffers.atomCount, kMeasured, frameMs, 1000.0 / frameMs,
                  gpu < 0 ? -1.0 : m_gpuSum / kMeasured, culling.c_str());
    m_entries.push_back(entry);
    std::cout << "[SceneBenchmark] " << entry << std::endl;
    m_frame = 0;
    if (++m_view == static_cast<int>(std::size(kViews))) finish();
}

void SceneBenchmark::finish()
{
    std::ofstream out(m_report);
    out << "[\n";
    for (std::size_t i = 0; i < m_entries.size(); ++i) out << m_entries[i] << (i + 1 < m_entries.size() ? ",\n" : "\n");
    out << "]\n";
    m_done = true;
}
