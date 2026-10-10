#include "app/LargeData.h"
#include "app/SceneView.h"

#include "Camera.h"
#include "graphics/CloudRenderer.h"
#include "third_party/stb_image_write.h"

#include "imgui.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace cloud = atomforge::cloud;

namespace
{
std::string lower(std::string text)
{
    for (char& c : text) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return text;
}

bool isLammpsDump(const std::filesystem::path& path)
{
    const std::string extension = lower(path.extension().string());
    if (extension == ".dump" || extension == ".lammpstrj") return true;
    std::ifstream in(path, std::ios::binary);
    std::string first;
    std::getline(in, first);
    return first.rfind("ITEM: TIMESTEP", 0) == 0;
}

std::string count(std::uint64_t n)
{
    char text[64];
    if (n >= 1000000000ull) std::snprintf(text, sizeof(text), "%.2f billion", static_cast<double>(n) / 1e9);
    else if (n >= 1000000ull) std::snprintf(text, sizeof(text), "%.1f million", static_cast<double>(n) / 1e6);
    else std::snprintf(text, sizeof(text), "%llu", static_cast<unsigned long long>(n));
    return text;
}

void screenshot(const std::string& file, int width, int height)
{
    std::vector<unsigned char> pixels(static_cast<std::size_t>(width) * height * 3), flipped(pixels.size());
    glReadBuffer(GL_BACK);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
    for (int row = 0; row < height; ++row)
        std::copy_n(pixels.begin() + static_cast<std::ptrdiff_t>(row) * width * 3, width * 3,
                    flipped.begin() + static_cast<std::ptrdiff_t>(height - 1 - row) * width * 3);
    stbi_write_png(file.c_str(), width, height, 3, flipped.data(), width * 3);
}
}

bool opensAsCloud(const std::filesystem::path& path)
{
    std::error_code error;
    if (!std::filesystem::is_regular_file(path, error)) return false;
    if (cloud::isCloudFile(path)) return true;
    if (isLammpsDump(path)) return true;
    const std::string extension = lower(path.extension().string());
    if (extension == ".xyz" || extension == ".extxyz") return cloud::peekAtomCount(path) > kLargeDataAtoms;
    return false;
}

// ---------------------------------------------------------------- conversion
CloudConversion::~CloudConversion()
{
    m_cancel = true;
    if (m_thread.joinable()) m_thread.join();
}

std::filesystem::path CloudConversion::cachePath(const std::filesystem::path& input)
{
    // Next to the input when that folder is writable, else in the user's cache.
    std::filesystem::path beside = input;
    beside += ".afcloud";
    {
        std::ofstream probe(beside.string() + ".probe");
        if (probe) {
            probe.close();
            std::error_code ignored;
            std::filesystem::remove(beside.string() + ".probe", ignored);
            return beside;
        }
    }
    std::filesystem::path root = std::filesystem::temp_directory_path();
    if (const char* local = std::getenv("LOCALAPPDATA")) root = std::filesystem::u8path(local);
    root /= "AtomForge";
    root /= "clouds";
    std::filesystem::create_directories(root);
    return root / (input.filename().u8string() + "-" + std::to_string(std::hash<std::string>{}(input.u8string())) + ".afcloud");
}

void CloudConversion::start(const std::filesystem::path& input)
{
    // A conversion still running is cancelled (its partial file is discarded).
    m_cancel = true;
    if (m_thread.joinable()) m_thread.join();
    m_input = input;
    m_finished = false;
    m_cancel = false;
    m_reported = false;
    m_progress = 0;
    m_error.clear();
    m_message = "Preparing";
    m_output = cachePath(input);
    // A cloud built from this file before (and newer than it) is reused.
    std::error_code ignored;
    if (cloud::isCloudFile(m_output) &&
        std::filesystem::last_write_time(m_output, ignored) >= std::filesystem::last_write_time(input, ignored)) {
        try {
            cloud::CloudFile check(m_output);
            m_finished = true;
            m_thread = std::thread([] {});
            return;
        } catch (const std::exception&) {}
    }
    m_thread = std::thread([this] {
        try {
            const cloud::Progress progress = [this](double fraction, const std::string& message) {
                if (fraction >= 0) m_progress = fraction;
                std::lock_guard<std::mutex> lock(m_mutex);
                m_message = message;
                return !m_cancel.load();
            };
            const auto temporary = std::filesystem::path(m_output).replace_extension(".afcloud-part");
            try {
                if (isLammpsDump(m_input)) cloud::buildFromLammpsDump(m_input, temporary, {}, 0, progress);
                else cloud::buildFromXyz(m_input, temporary, 0, progress);
            } catch (...) {
                std::error_code removeError;
                std::filesystem::remove(temporary, removeError);
                throw;
            }
            std::filesystem::rename(temporary, m_output);
        } catch (const std::exception& error) {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_error = error.what();
        }
        m_finished = true;
    });
}

bool CloudConversion::finished(std::filesystem::path& result, std::string& error)
{
    if (!m_thread.joinable() || !m_finished || m_reported) return false;
    m_thread.join();
    m_reported = true;
    std::lock_guard<std::mutex> lock(m_mutex);
    error = m_error;
    result = m_error.empty() ? m_output : std::filesystem::path();
    return true;
}

std::string CloudConversion::message() const
{
    std::lock_guard<std::mutex> lock(m_mutex);
    return m_message;
}

// ---------------------------------------------------------------- view
void fitCameraToCloud(Camera& camera, const CloudRenderer& cloudRenderer, int width, int height)
{
    camera.yaw = 45.0f;
    camera.pitch = 35.2643897f;
    camera.roll = 0.0f;
    camera.panOffset = glm::vec3(0.0f);
    const float aspect = height > 0 ? static_cast<float>(width) / static_cast<float>(height) : 1.0f;
    const float vertical = glm::radians(45.0f);
    const float horizontal = 2.0f * std::atan(std::tan(vertical * 0.5f) * aspect);
    const float half = std::max(glm::radians(10.0f), 0.5f * std::min(vertical, horizontal));
    camera.distance = std::clamp(cloudRenderer.radius() / std::sin(half) * 1.1f, Camera::kMinDistance, Camera::kMaxDistance);
}

void adjustCloudProjection(FrameView& frame, const Camera& camera, const CloudRenderer& cloudRenderer, bool orthographic)
{
    if (orthographic || frame.framebufferHeight <= 0) return;
    const float aspect = static_cast<float>(frame.framebufferWidth) / static_cast<float>(frame.framebufferHeight);
    const float radius = cloudRenderer.radius();
    const float toCentre = glm::length(frame.cameraPosition - cloudRenderer.center());
    // Near plane: well in front of the nearest atoms, never so close that depth precision is lost.
    const float nearest = toCentre - radius;
    const float nearPlane = nearest > 0 ? std::max(0.5f, nearest * 0.5f) : std::max(0.5f, camera.distance * 0.002f);
    const float farPlane = std::max(nearPlane * 2.0f, toCentre + radius * 1.05f);
    frame.projection = glm::perspective(glm::radians(45.0f), aspect, nearPlane, farPlane);
}

void drawCloudPanel(CloudRenderer& cloudRenderer, float frameMilliseconds)
{
    const auto& stats = cloudRenderer.stats();
    const auto& info = cloudRenderer.info();
    ImGui::SetNextWindowPos(ImVec2(12, ImGui::GetIO().DisplaySize.y - 12), ImGuiCond_FirstUseEver, ImVec2(0, 1));
    ImGui::SetNextWindowSize(ImVec2(400, 0), ImGuiCond_FirstUseEver);
    if (ImGui::Begin("Large dataset", nullptr, ImGuiWindowFlags_NoCollapse)) {
        ImGui::TextUnformatted(cloudRenderer.path().filename().u8string().c_str());
        if (!info.source.empty()) ImGui::TextDisabled("%s", info.source.c_str());
        ImGui::Text("%s atoms in %zu chunks", count(stats.totalAtoms).c_str(), stats.chunks);
        ImGui::Text("Drawn: %s atoms, %zu chunks visible, %zu hidden", count(stats.drawnAtoms).c_str(), stats.visibleChunks, stats.occludedChunks);
        ImGui::Text("GPU: %.1f ms (frame %.1f ms), %.0f MB resident", stats.gpuMilliseconds, frameMilliseconds, stats.residentMB);
        if (!stats.complete) ImGui::TextColored(ImVec4(0.95f, 0.75f, 0.3f, 1), "Streaming (%zu loads pending)", stats.loadingChunks);
        else ImGui::TextColored(ImVec4(0.45f, 0.8f, 0.45f, 1), "Up to date");
        ImGui::Separator();
        auto& s = cloudRenderer.settings;
        float budget = static_cast<float>(s.pointBudget / 1e6);
        if (ImGui::SliderFloat("Atoms per frame (M)", &budget, 1.0f, 500.0f, "%.0f", ImGuiSliderFlags_Logarithmic)) s.pointBudget = budget * 1e6;
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("The most atoms drawn in one frame; lower it for a faster view.");
        float memory = static_cast<float>(s.memoryBudgetMB);
        if (ImGui::SliderFloat("GPU memory (MB)", &memory, 256.0f, 16384.0f, "%.0f", ImGuiSliderFlags_Logarithmic)) s.memoryBudgetMB = memory;
        ImGui::Checkbox("Level of detail", &s.levelOfDetail);
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Draw distant regions with a random subsample of their atoms (as many as the screen resolves).");
        ImGui::BeginDisabled(!s.levelOfDetail);
        ImGui::SliderFloat("Detail", &s.density, 0.5f, 20.0f, "%.1f", ImGuiSliderFlags_Logarithmic);
        ImGui::EndDisabled();
        ImGui::SliderFloat("Atom size", &s.radiusScale, 0.2f, 3.0f, "%.2f");
        const char* colours[] = {"Species", "Height along x", "Height along y", "Height along z"};
        ImGui::Combo("Colour", &cloudRenderer.colourAxis, colours, 4);
        ImGui::Checkbox("Show box", &s.showBox);
        ImGui::SameLine();
        ImGui::Checkbox("Occlusion culling", &s.occlusionCulling);
        ImGui::TextDisabled("Large datasets are view-only; File > Open a smaller file to edit.");
    }
    ImGui::End();
}

void drawConversionProgress(CloudConversion& conversion)
{
    if (!conversion.running()) return;
    ImGui::SetNextWindowPos(ImVec2(ImGui::GetIO().DisplaySize.x * 0.5f, ImGui::GetIO().DisplaySize.y * 0.5f), ImGuiCond_Always, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(480, 0));
    if (ImGui::Begin("Preparing a large dataset", nullptr, ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoSavedSettings)) {
        ImGui::TextWrapped("%s is too large to edit; AtomForge converts it to an atom cloud (once) and shows it with the large-data renderer.",
                           conversion.input().filename().u8string().c_str());
        ImGui::ProgressBar(static_cast<float>(conversion.progress()), ImVec2(-1, 0));
        ImGui::TextDisabled("%s", conversion.message().c_str());
        if (ImGui::Button("Cancel")) conversion.cancel();
    }
    ImGui::End();
}

// ---------------------------------------------------------------- benchmark
namespace
{
struct BenchmarkView { const char* name; float distanceScale; glm::vec3 pan; bool orbit; };
// Distances relative to the fitted view; pans as fractions of the cloud's half extent.
const BenchmarkView kViews[] = {
    {"overview", 1.0f, {0, 0, 0}, false},
    {"half distance", 0.5f, {0, 0, 0}, false},
    {"orbiting overview", 1.0f, {0, 0, 0}, true},
    {"region", 0.12f, {0.6f, 0.6f, 0.6f}, false},
    {"close-up", 0.0f, {0.9f, 0.9f, 0.9f}, false},  // 60 Angstrom from a corner region
};
constexpr int kMeasuredFrames = 120;
constexpr double kStreamTimeout = 60.0;
}

std::unique_ptr<CloudBenchmark> CloudBenchmark::fromEnvironment()
{
    const char* report = std::getenv("ATOMFORGE_CLOUD_BENCHMARK");
    if (!report || !*report) return nullptr;
    auto benchmark = std::make_unique<CloudBenchmark>();
    benchmark->m_report = report;
    if (const char* shots = std::getenv("ATOMFORGE_CLOUD_SHOTS")) benchmark->m_shots = shots;
    return benchmark;
}

void CloudBenchmark::beforeFrame(Camera& camera, const CloudRenderer& cloudRenderer, int width, int height)
{
    if (m_done || !cloudRenderer.active()) return;
    const double now = glfwGetTime();
    if (m_viewStart == 0) { glfwSwapInterval(0); m_viewStart = now; }
    const auto& view = kViews[m_view];
    fitCameraToCloud(camera, cloudRenderer, width, height);
    const auto& info = cloudRenderer.info();
    const glm::vec3 half(0.5 * (info.upper[0] - info.lower[0]), 0.5 * (info.upper[1] - info.lower[1]), 0.5 * (info.upper[2] - info.lower[2]));
    camera.panOffset = view.pan * half;
    camera.distance = view.distanceScale > 0 ? camera.distance * view.distanceScale : 60.0f;
    if (view.orbit) camera.yaw = 45.0f + static_cast<float>(m_frames) * 3.0f;
}

void CloudBenchmark::afterFrame(const CloudRenderer& cloudRenderer, int width, int height)
{
    if (m_done || !cloudRenderer.active()) return;
    const auto& stats = cloudRenderer.stats();
    const double now = glfwGetTime();
    const auto& view = kViews[m_view];
    if (!m_measuring) {
        // Wait for the view's atoms to stream in (an orbit is measured while streaming).
        if (stats.complete || view.orbit || now - m_viewStart > kStreamTimeout) {
            m_streamSeconds = now - m_viewStart;
            m_measuring = true;
            m_measureStart = now;
            m_frames = 0;
        }
        return;
    }
    glFinish();
    static double gpuSum = 0;
    static std::uint64_t drawnSum = 0;
    if (m_frames == 0) { gpuSum = 0; drawnSum = 0; }
    gpuSum += stats.gpuMilliseconds;
    drawnSum += stats.drawnAtoms;
    if (++m_frames < kMeasuredFrames) return;
    const double seconds = glfwGetTime() - m_measureStart;
    char entry[768];
    std::snprintf(entry, sizeof(entry),
                  "%s\n  {\"view\": \"%s\", \"stream_seconds\": %.2f, \"frames\": %d, \"frame_ms\": %.2f, \"fps\": %.1f, \"gpu_ms\": %.2f, "
                  "\"drawn_atoms\": %.0f, \"visible_chunks\": %zu, \"occluded_chunks\": %zu, \"resident_atoms\": %llu, \"resident_mb\": %.0f, \"complete\": %s, \"total_atoms\": %llu}",
                  m_view ? "," : "", view.name, m_streamSeconds, m_frames, 1000.0 * seconds / m_frames, m_frames / seconds, gpuSum / m_frames,
                  static_cast<double>(drawnSum) / m_frames, stats.visibleChunks, stats.occludedChunks, static_cast<unsigned long long>(stats.residentAtoms), stats.residentMB,
                  stats.complete ? "true" : "false", static_cast<unsigned long long>(stats.totalAtoms));
    m_results += entry;
    if (!m_shots.empty()) {
        std::string name = view.name;
        std::replace(name.begin(), name.end(), ' ', '-');
        screenshot(m_shots + "/" + std::to_string(m_view + 1) + "-" + name + ".png", width, height);
    }
    m_measuring = false;
    m_viewStart = now;
    if (++m_view >= static_cast<int>(sizeof(kViews) / sizeof(kViews[0]))) {
        std::ofstream(m_report) << m_results << "\n]\n";
        m_done = true;
    }
}
