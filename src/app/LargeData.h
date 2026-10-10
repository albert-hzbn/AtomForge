#pragma once
// The desktop's large-data path: atom clouds (cloud/AtomCloud.h) shown with
// CloudRenderer instead of the editable structure renderer. Huge LAMMPS dumps
// and XYZ files are converted to a cloud in the background on opening.
#include "cloud/AtomCloud.h"

#include <atomic>
#include <filesystem>
#include <memory>
#include <mutex>
#include <string>
#include <thread>

class CloudRenderer;
struct Camera;
struct FrameView;

// Structures above this many atoms open as atom clouds.
constexpr std::uint64_t kLargeDataAtoms = 5'000'000;

// Whether a file opens through the large-data path: an .afcloud, or a dump or
// XYZ file with more than kLargeDataAtoms atoms (or a format the structure
// loaders do not read: LAMMPS dumps).
bool opensAsCloud(const std::filesystem::path& path);

// Background conversion of a large file to an atom cloud.
class CloudConversion
{
public:
    ~CloudConversion();
    // Starts converting `input` (an existing up-to-date cloud is reused at once).
    void start(const std::filesystem::path& input);
    bool running() const { return m_thread.joinable() && !m_finished; }
    // When the conversion ended: the cloud to open (empty on failure, see error()).
    bool finished(std::filesystem::path& cloud, std::string& error);
    void cancel() { m_cancel = true; }
    double progress() const { return m_progress; }
    std::string message() const;
    const std::filesystem::path& input() const { return m_input; }
    // Where the cloud of `input` is kept: next to it, or in the user's cache folder.
    static std::filesystem::path cachePath(const std::filesystem::path& input);

private:
    std::filesystem::path m_input, m_output;
    std::thread m_thread;
    std::atomic<bool> m_finished{false}, m_cancel{false}, m_reported{false};
    std::atomic<double> m_progress{0};
    mutable std::mutex m_mutex;
    std::string m_message, m_error;
};

// Fits the camera to a cloud (iso view).
void fitCameraToCloud(Camera& camera, const CloudRenderer& cloud, int width, int height);
// Perspective clip planes for a cloud (a sensible near plane inside huge scenes).
void adjustCloudProjection(FrameView& frame, const Camera& camera, const CloudRenderer& cloud, bool orthographic);
// The "Large dataset" panel: statistics and level-of-detail settings.
void drawCloudPanel(CloudRenderer& cloud, float frameMilliseconds);
// Progress window of a running conversion.
void drawConversionProgress(CloudConversion& conversion);

// ATOMFORGE_CLOUD_BENCHMARK=REPORT.json: renders a fixed camera path over the
// open cloud, waits for streaming at each view, times frames and writes the
// report (and screenshots with ATOMFORGE_CLOUD_SHOTS=FOLDER), then quits.
class CloudBenchmark
{
public:
    static std::unique_ptr<CloudBenchmark> fromEnvironment();
    // Called once per frame before drawing (sets the camera), and after drawing.
    void beforeFrame(Camera& camera, const CloudRenderer& cloud, int width, int height);
    void afterFrame(const CloudRenderer& cloud, int width, int height);
    bool done() const { return m_done; }

private:
    std::string m_report, m_shots;
    int m_view = 0, m_frames = 0;
    double m_viewStart = 0, m_measureStart = 0, m_streamSeconds = 0;
    bool m_measuring = false, m_done = false;
    std::string m_results = "[";
};
