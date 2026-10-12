#pragma once
// Renders out-of-core atom clouds (cloud/AtomCloud.h) of up to billions of
// atoms. Atoms are point-sprite sphere impostors (lit, with correct depth),
// 8 bytes each on the GPU. Every frame each chunk is culled against the view
// and given a number of atoms from its size on screen; as chunks are stored
// in random order, drawing a prefix draws a uniform subsample. A total point
// budget keeps the frame rate, a memory budget bounds the GPU memory, and a
// background thread streams the needed prefixes from the file. It draws with
// OpenGL or Vulkan (CloudRendererVulkan.cpp), whichever the desktop uses.
#include "cloud/AtomCloud.h"

#include <GL/glew.h>
#include <glm/glm.hpp>

#include <atomic>
#include <condition_variable>
#include <deque>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

typedef struct VkCommandBuffer_T* VkCommandBuffer;

class CloudRenderer
{
public:
    struct Settings
    {
        double pointBudget = 60e6;            // atoms drawn per frame at most
        double memoryBudgetMB = 2560;         // GPU memory for atom data
        float density = 3.0f;                 // atoms per atom-sized screen area of a chunk (level of detail)
        float radiusScale = 1.0f;
        bool levelOfDetail = true;            // false: every atom of visible chunks (within the budgets)
        bool showBox = true;
        bool occlusionCulling = true;         // skip chunks hidden behind drawn atoms
    };
    struct Stats
    {
        std::uint64_t totalAtoms = 0, drawnAtoms = 0, residentAtoms = 0;
        std::size_t chunks = 0, visibleChunks = 0, occludedChunks = 0, loadingChunks = 0;
        double residentMB = 0, uploadedMBThisFrame = 0;
        bool complete = false;                // every visible chunk has the atoms it wants
        double gpuMilliseconds = 0;           // GPU time of the last measured draw
    };

    CloudRenderer() = default;
    ~CloudRenderer();
    CloudRenderer(const CloudRenderer&) = delete;
    CloudRenderer& operator=(const CloudRenderer&) = delete;

    // Opens a cloud for display (GL context current, or the Vulkan device
    // ready); throws on errors.
    void open(const std::filesystem::path& path);
    void close();
    bool active() const { return static_cast<bool>(m_file); }
    const atomforge::cloud::CloudInfo& info() const;
    const std::filesystem::path& path() const { return m_path; }

    // OpenGL: draws into the current framebuffer (depth test on); `height` in pixels.
    void draw(const glm::mat4& projection, const glm::mat4& view, int width, int height, bool orthographic, bool lightTheme);
    // Vulkan: prepare() records the uploads before the render pass, and
    // drawVulkan() the drawing inside it, of the same view.
    void prepareVulkan(VkCommandBuffer cmd, const glm::mat4& projection, const glm::mat4& view, int width, int height, bool orthographic,
                       bool immediate = false);
    void drawVulkan(VkCommandBuffer cmd, int width, int height, bool lightTheme);
    // Vulkan image export: draws off-screen until the atoms have streamed in
    // (at most 30 s) and returns RGBA pixels, rows top first.
    bool captureVulkan(const glm::mat4& projection, const glm::mat4& view, int width, int height, bool orthographic, bool lightTheme,
                       const glm::vec3& background, std::vector<unsigned char>& rgba, std::string& error);

    Settings settings;
    const Stats& stats() const { return m_stats; }
    glm::vec3 center() const;
    float radius() const;
    // Colours by species, or by height along an axis (0 none, 1 x, 2 y, 3 z).
    int colourAxis = 0;

private:
    struct ChunkState
    {
        GLuint buffer = 0;
        std::uint64_t capacity = 0;   // atoms the buffer holds
        std::uint64_t resident = 0;   // leading atoms uploaded
        std::uint64_t wanted = 0;     // atoms wanted this frame
        std::uint64_t requested = 0;  // prefix requested from the loader
        std::uint64_t lastUsed = 0;   // frame it was last drawn
        float priority = 0;
    };
    struct Load { std::size_t chunk; std::uint64_t first, count; std::vector<atomforge::cloud::PackedAtom> atoms; };
    // A chunk to draw this frame: `cover` atoms fill its screen area; `want`
    // is what it gets (all atoms once they are about a pixel or larger on screen).
    struct Visible { std::size_t chunk; double want; double cover; float distance; float area; };
    struct VulkanState;

    void ensureProgram();
    void loaderLoop();
    void consumeLoads();
    void upload(Load& load);
    void evict(std::uint64_t frame);
    void releaseChunk(std::size_t chunk);
    void request(std::size_t chunk, std::uint64_t upTo);
    // Culling, level of detail, budgets and streaming of a frame (both APIs).
    void plan(const glm::mat4& projection, const glm::mat4& view, int width, int height, bool orthographic);
    void finishStats(std::uint64_t drawn, std::size_t occluded);
    // Vulkan parts (CloudRendererVulkan.cpp).
    void openVulkan(const std::vector<glm::vec3>& outline);
    void closeVulkan();
    bool uploadVulkan(Load& load, std::uint64_t capacity);
    void releaseVulkanChunk(std::size_t chunk);
    void readVulkanQueries();

    std::filesystem::path m_path;
    std::unique_ptr<atomforge::cloud::CloudFile> m_file;
    std::vector<ChunkState> m_chunks;
    GLuint m_program = 0, m_lineProgram = 0, m_vao = 0, m_lineVao = 0, m_lineBuffer = 0, m_speciesTexture = 0;
    GLint m_maxPointSize = 256;
    GLuint m_timers[2] = {0, 0};   // GPU timer queries, alternating frames
    // Occlusion culling: one query per chunk on its box, read a frame later.
    GLuint m_boxVao = 0, m_boxBuffer = 0;
    std::vector<GLuint> m_queries;
    std::vector<char> m_queryPending, m_occluded;
    bool m_timerUsed[2] = {false, false};
    std::uint64_t m_frame = 0;
    std::uint64_t m_residentBytes = 0;
    Stats m_stats;
    std::shared_ptr<VulkanState> m_vk;

    // The plan of the current frame.
    std::vector<Visible> m_visible;
    std::vector<std::size_t> m_tested;  // in the frustum and outside the camera: occlusion-tested
    glm::mat4 m_projection{1.0f}, m_view{1.0f};
    float m_pixelsPerUnit = 1.0f, m_maxRadius = 0.3f;
    bool m_orthographic = false, m_complete = false;

    // Background loading.
    std::thread m_loader;
    std::mutex m_mutex;
    std::condition_variable m_wake;
    std::deque<Load> m_queue, m_done;
    bool m_stop = false;
    std::atomic<std::size_t> m_inFlight{0};
};
