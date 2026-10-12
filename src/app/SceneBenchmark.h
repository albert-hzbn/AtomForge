#pragma once
// ATOMFORGE_SCENE_BENCHMARK=REPORT.json: measures drawing the structure
// opened at start (views from its fitted distance down to near its surface,
// and an orbit) with vertical sync off, and writes the frame times, the GPU
// time and what the GPU culling drew (Vulkan) as JSON; then quits.

#include <memory>
#include <string>
#include <vector>

struct Camera;
struct SceneBuffers;

class SceneBenchmark
{
public:
    static std::unique_ptr<SceneBenchmark> fromEnvironment();
    // Before the frame's view is built: the camera of the current view.
    void beforeFrame(Camera& camera, const SceneBuffers& buffers);
    void afterFrame(const SceneBuffers& buffers);
    bool done() const { return m_done; }

private:
    void finish();

    std::string m_report;
    std::vector<std::string> m_entries;
    int m_view = 0;
    int m_frame = 0;       // frames of the current view
    int m_wait = 0;        // frames before the first view (loading, fitting)
    float m_distance = 0;  // fitted camera distance
    float m_yaw = 0, m_pitch = 0;
    double m_start = 0;
    double m_gpuSum = 0;
    bool m_done = false;
};
