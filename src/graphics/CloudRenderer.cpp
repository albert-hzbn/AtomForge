#include "graphics/CloudRenderer.h"
#include "graphics/Shader.h"

#include <glm/gtc/type_ptr.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

using atomforge::cloud::PackedAtom;

namespace
{
const char* kPointVS = R"(#version 330 core
layout(location = 0) in uvec4 aPacked;
uniform mat4 uView;
uniform mat4 uProjection;
uniform vec3 uLower;
uniform vec3 uScale;
uniform float uPixelsPerUnit;   // projection[1][1] * viewport height / 2
uniform bool uOrthographic;
uniform float uRadiusScale;
uniform float uFill;            // point enlargement for subsampled chunks
uniform float uMaxPointSize;
uniform sampler2D uSpecies;     // rgb colour, a radius
uniform int uColourAxis;
uniform vec2 uColourRange;
out vec3 vColour;
out vec3 vCentre;
out float vRadius;
vec3 ramp(float t)
{
    t = clamp(t, 0.0, 1.0);
    return clamp(vec3(1.5 - abs(4.0 * t - 3.0), 1.5 - abs(4.0 * t - 2.0), 1.5 - abs(4.0 * t - 1.0)), 0.0, 1.0);
}
void main()
{
    vec3 p = uLower + vec3(aPacked.xyz) * uScale;
    int s = int(aPacked.w);
    vec4 species = texelFetch(uSpecies, ivec2(s & 255, s >> 8), 0);
    vColour = species.rgb;
    if (uColourAxis > 0) vColour = ramp((p[uColourAxis - 1] - uColourRange.x) / max(uColourRange.y - uColourRange.x, 1e-6));
    vRadius = species.a * uRadiusScale;
    vec4 viewPosition = uView * vec4(p, 1.0);
    vCentre = viewPosition.xyz;
    gl_Position = uProjection * viewPosition;
    // Rasterise at the sphere's front so that the depth written later is never
    // nearer than this (conservative depth keeps early depth rejection).
    vec4 front = uProjection * vec4(viewPosition.xy, viewPosition.z + vRadius, 1.0);
    gl_Position.z = front.z / front.w * gl_Position.w;
    float pixels = vRadius * uPixelsPerUnit / (uOrthographic ? 1.0 : max(-viewPosition.z, 1e-3));
    gl_PointSize = clamp(2.0 * pixels * uFill, 1.0, uMaxPointSize);
}
)";

const char* kPointFS = R"(#version 330 core
#extension GL_ARB_conservative_depth : enable
#ifdef GL_ARB_conservative_depth
layout(depth_greater) out float gl_FragDepth;
#endif
in vec3 vColour;
in vec3 vCentre;
in float vRadius;
uniform mat4 uProjection;
uniform vec3 uLight;            // view space, normalised
uniform float uAmbient;
out vec4 fragColour;
void main()
{
    vec2 c = gl_PointCoord * 2.0 - 1.0;
    c.y = -c.y;
    float d2 = dot(c, c);
    if (d2 > 1.0) discard;
    vec3 n = vec3(c, sqrt(1.0 - d2));
    // Depth of the sphere surface, so that impostors intersect correctly.
    vec4 clip = uProjection * vec4(vCentre + n * vRadius, 1.0);
    gl_FragDepth = clip.z / clip.w * 0.5 + 0.5;
    float diffuse = max(dot(n, uLight), 0.0);
    float specular = pow(max(dot(reflect(-uLight, n), vec3(0.0, 0.0, 1.0)), 0.0), 40.0) * 0.35;
    vec3 colour = vColour * (uAmbient + (1.0 - uAmbient) * diffuse) + vec3(specular);
    colour *= mix(0.72, 1.0, n.z);  // darken the rim for depth
    fragColour = vec4(colour, 1.0);
}
)";

const char* kLineVS = R"(#version 330 core
layout(location = 0) in vec3 aPosition;
uniform mat4 uMVP;
void main() { gl_Position = uMVP * vec4(aPosition, 1.0); }
)";
const char* kLineFS = R"(#version 330 core
uniform vec3 uColour;
out vec4 fragColour;
void main() { fragColour = vec4(uColour, 1.0); }
)";

constexpr std::uint64_t kBlock = 1u << 16;            // atoms per load granule
constexpr std::uint64_t kMaxLoad = 4u << 20;          // atoms per load request
constexpr std::size_t kMaxInFlight = 12;
constexpr double kUploadBytesPerFrame = 384.0 * 1024 * 1024;
}

CloudRenderer::~CloudRenderer() { close(); }

const atomforge::cloud::CloudInfo& CloudRenderer::info() const
{
    if (!m_file) throw std::logic_error("No atom cloud open");
    return m_file->info();
}

glm::vec3 CloudRenderer::center() const
{
    if (!m_file) return glm::vec3(0);
    const auto& i = m_file->info();
    return glm::vec3(0.5 * (i.lower[0] + i.upper[0]), 0.5 * (i.lower[1] + i.upper[1]), 0.5 * (i.lower[2] + i.upper[2]));
}

float CloudRenderer::radius() const
{
    if (!m_file) return 1.0f;
    const auto& i = m_file->info();
    return std::max(1.0f, 0.5f * glm::length(glm::vec3(i.upper[0] - i.lower[0], i.upper[1] - i.lower[1], i.upper[2] - i.lower[2])));
}

void CloudRenderer::ensureProgram()
{
    if (m_program) return;
    m_program = createProgram(kPointVS, kPointFS);
    m_lineProgram = createProgram(kLineVS, kLineFS);
    if (!m_program || !m_lineProgram) throw std::runtime_error("Cannot compile the atom cloud shaders");
    glGenVertexArrays(1, &m_vao);
    glGenQueries(2, m_timers);
    glGenVertexArrays(1, &m_lineVao);
    glGenBuffers(1, &m_lineBuffer);
    // Unit cube (triangles) for the chunk occlusion queries.
    {
        static const float corners[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
        static const int faces[12][3] = {{0, 2, 1}, {1, 2, 3}, {4, 5, 6}, {5, 7, 6}, {0, 1, 4}, {1, 5, 4},
                                         {2, 6, 3}, {3, 6, 7}, {0, 4, 2}, {2, 4, 6}, {1, 3, 5}, {3, 7, 5}};
        std::vector<float> triangles;
        for (const auto& f : faces)
            for (int v : f) triangles.insert(triangles.end(), corners[v], corners[v] + 3);
        glGenVertexArrays(1, &m_boxVao);
        glGenBuffers(1, &m_boxBuffer);
        glBindVertexArray(m_boxVao);
        glBindBuffer(GL_ARRAY_BUFFER, m_boxBuffer);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(triangles.size() * sizeof(float)), triangles.data(), GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
        glBindVertexArray(0);
    }
    GLfloat range[2] = {1, 256};
    glGetFloatv(GL_POINT_SIZE_RANGE, range);
    m_maxPointSize = static_cast<GLint>(std::max(16.0f, range[1]));
}

void CloudRenderer::open(const std::filesystem::path& path)
{
    close();
    ensureProgram();
    m_file = std::make_unique<atomforge::cloud::CloudFile>(path);
    m_path = path;
    const auto& info = m_file->info();
    m_chunks.assign(info.chunks.size(), ChunkState{});
    m_queries.assign(info.chunks.size(), 0);
    glGenQueries(static_cast<GLsizei>(m_queries.size()), m_queries.data());
    m_queryPending.assign(info.chunks.size(), 0);
    m_occluded.assign(info.chunks.size(), 0);
    m_frame = 0;
    m_residentBytes = 0;
    m_stats = {};
    m_stats.totalAtoms = info.atoms;
    m_stats.chunks = info.chunks.size();

    // Species colours and radii.
    const std::size_t rows = std::max<std::size_t>(1, (info.species.size() + 255) / 256);
    std::vector<float> texels(256 * rows * 4, 0.7f);
    for (std::size_t i = 0; i < info.species.size(); ++i) {
        texels[i * 4 + 0] = info.species[i].color[0];
        texels[i * 4 + 1] = info.species[i].color[1];
        texels[i * 4 + 2] = info.species[i].color[2];
        texels[i * 4 + 3] = info.species[i].radius;
    }
    glGenTextures(1, &m_speciesTexture);
    glBindTexture(GL_TEXTURE_2D, m_speciesTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, 256, static_cast<GLsizei>(rows), 0, GL_RGBA, GL_FLOAT, texels.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);

    // Outline: the cell when known, else the bounding box.
    glm::vec3 o, a, b, c;
    if (info.hasCell) {
        o = glm::vec3(info.cellOrigin[0], info.cellOrigin[1], info.cellOrigin[2]);
        a = glm::vec3(info.cell[0][0], info.cell[0][1], info.cell[0][2]);
        b = glm::vec3(info.cell[1][0], info.cell[1][1], info.cell[1][2]);
        c = glm::vec3(info.cell[2][0], info.cell[2][1], info.cell[2][2]);
    } else {
        o = glm::vec3(info.lower[0], info.lower[1], info.lower[2]);
        a = glm::vec3(info.upper[0] - info.lower[0], 0, 0);
        b = glm::vec3(0, info.upper[1] - info.lower[1], 0);
        c = glm::vec3(0, 0, info.upper[2] - info.lower[2]);
    }
    const glm::vec3 corner[8] = {o, o + a, o + b, o + a + b, o + c, o + a + c, o + b + c, o + a + b + c};
    const int edges[12][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 2}, {1, 3}, {4, 6}, {5, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    std::vector<glm::vec3> lines;
    for (const auto& e : edges) { lines.push_back(corner[e[0]]); lines.push_back(corner[e[1]]); }
    glBindVertexArray(m_lineVao);
    glBindBuffer(GL_ARRAY_BUFFER, m_lineBuffer);
    glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(lines.size() * sizeof(glm::vec3)), lines.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), nullptr);
    glBindVertexArray(0);

    m_stop = false;
    m_loader = std::thread([this] { loaderLoop(); });
}

void CloudRenderer::close()
{
    if (m_loader.joinable()) {
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_stop = true;
        }
        m_wake.notify_all();
        m_loader.join();
    }
    m_queue.clear();
    m_done.clear();
    m_inFlight = 0;
    for (auto& c : m_chunks)
        if (c.buffer) glDeleteBuffers(1, &c.buffer);
    m_chunks.clear();
    if (!m_queries.empty()) glDeleteQueries(static_cast<GLsizei>(m_queries.size()), m_queries.data());
    m_queries.clear();
    m_queryPending.clear();
    m_occluded.clear();
    if (m_speciesTexture) { glDeleteTextures(1, &m_speciesTexture); m_speciesTexture = 0; }
    m_file.reset();
    m_residentBytes = 0;
    m_stats = {};
}

void CloudRenderer::loaderLoop()
{
    std::ifstream file(m_path, std::ios::binary);
    for (;;) {
        Load load;
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_wake.wait(lock, [&] { return m_stop || !m_queue.empty(); });
            if (m_stop) return;
            load = std::move(m_queue.front());
            m_queue.pop_front();
        }
        const auto& chunk = m_file->info().chunks[load.chunk];
        load.atoms.resize(load.count);
        file.clear();
        file.seekg(static_cast<std::streamoff>(chunk.offset + load.first * sizeof(PackedAtom)));
        file.read(reinterpret_cast<char*>(load.atoms.data()), static_cast<std::streamsize>(load.count * sizeof(PackedAtom)));
        if (!file) load.atoms.clear();  // read error: dropped, requested again later
        std::lock_guard<std::mutex> lock(m_mutex);
        m_done.push_back(std::move(load));
    }
}

void CloudRenderer::request(std::size_t chunk, std::uint64_t upTo)
{
    ChunkState& c = m_chunks[chunk];
    const std::uint64_t count = m_file->info().chunks[chunk].count;
    upTo = std::min(count, upTo);
    if (upTo <= c.requested) return;
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        // In pieces, so that one large chunk does not hold up the others.
        for (std::uint64_t first = c.requested; first < upTo; first += kMaxLoad)
            m_queue.push_back({chunk, first, std::min(kMaxLoad, upTo - first), {}});
    }
    m_inFlight += static_cast<std::size_t>((upTo - c.requested + kMaxLoad - 1) / kMaxLoad);
    c.requested = upTo;
    m_wake.notify_one();
}

void CloudRenderer::upload(Load& load)
{
    ChunkState& c = m_chunks[load.chunk];
    if (load.atoms.empty() || load.first != c.resident) {
        // Out of order (after an eviction) or failed: forget the request.
        c.requested = c.resident;
        return;
    }
    const std::uint64_t end = load.first + load.count;
    if (end > c.capacity) {
        const std::uint64_t total = m_file->info().chunks[load.chunk].count;
        const std::uint64_t capacity = std::min(total, std::max(end, c.capacity * 2));
        GLuint buffer = 0;
        glGenBuffers(1, &buffer);
        glBindBuffer(GL_COPY_WRITE_BUFFER, buffer);
        glBufferData(GL_COPY_WRITE_BUFFER, static_cast<GLsizeiptr>(capacity * sizeof(PackedAtom)), nullptr, GL_STATIC_DRAW);
        if (c.buffer && c.resident) {
            glBindBuffer(GL_COPY_READ_BUFFER, c.buffer);
            glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0, static_cast<GLsizeiptr>(c.resident * sizeof(PackedAtom)));
        }
        if (c.buffer) glDeleteBuffers(1, &c.buffer);
        m_residentBytes += (capacity - c.capacity) * sizeof(PackedAtom);
        c.buffer = buffer;
        c.capacity = capacity;
    }
    glBindBuffer(GL_COPY_WRITE_BUFFER, c.buffer);
    glBufferSubData(GL_COPY_WRITE_BUFFER, static_cast<GLintptr>(load.first * sizeof(PackedAtom)),
                    static_cast<GLsizeiptr>(load.count * sizeof(PackedAtom)), load.atoms.data());
    glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
    c.resident = end;
    m_stats.uploadedMBThisFrame += static_cast<double>(load.count * sizeof(PackedAtom)) / (1024.0 * 1024.0);
}

void CloudRenderer::evict(std::uint64_t frame)
{
    const double budget = settings.memoryBudgetMB * 1024.0 * 1024.0;
    if (static_cast<double>(m_residentBytes) <= budget) return;
    // Least recently drawn chunks first; never one with loads in flight.
    std::vector<std::size_t> candidates;
    for (std::size_t i = 0; i < m_chunks.size(); ++i)
        if (m_chunks[i].buffer && m_chunks[i].lastUsed < frame && m_chunks[i].requested == m_chunks[i].resident) candidates.push_back(i);
    std::sort(candidates.begin(), candidates.end(), [&](std::size_t a, std::size_t b) { return m_chunks[a].lastUsed < m_chunks[b].lastUsed; });
    for (std::size_t i : candidates) {
        if (static_cast<double>(m_residentBytes) <= budget * 0.9) break;
        ChunkState& c = m_chunks[i];
        glDeleteBuffers(1, &c.buffer);
        m_residentBytes -= c.capacity * sizeof(PackedAtom);
        c = ChunkState{};
    }
}

void CloudRenderer::draw(const glm::mat4& projection, const glm::mat4& view, int width, int height, bool orthographic, bool lightTheme)
{
    if (!m_file || width <= 0 || height <= 0) return;
    ++m_frame;
    m_stats.uploadedMBThisFrame = 0;
    const auto& info = m_file->info();

    // Finished loads, up to a per-frame upload volume.
    {
        double bytes = 0;
        for (;;) {
            Load load;
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                if (m_done.empty() || bytes > kUploadBytesPerFrame) break;
                load = std::move(m_done.front());
                m_done.pop_front();
            }
            --m_inFlight;
            bytes += static_cast<double>(load.count * sizeof(PackedAtom));
            upload(load);
        }
    }

    // Occlusion results of earlier frames (only those already available: no stall).
    for (std::size_t i = 0; i < m_queries.size(); ++i) {
        if (!m_queryPending[i]) continue;
        GLuint available = 0;
        glGetQueryObjectuiv(m_queries[i], GL_QUERY_RESULT_AVAILABLE, &available);
        if (!available) continue;
        GLuint samples = 0;
        glGetQueryObjectuiv(m_queries[i], GL_QUERY_RESULT, &samples);
        m_occluded[i] = samples < 16;
        m_queryPending[i] = 0;
    }

    // Visibility and level of detail per chunk.
    const glm::mat4 viewProjection = projection * view;
    const float pixelsPerUnit = projection[1][1] * static_cast<float>(height) * 0.5f;
    float maxRadius = 0.3f;
    for (const auto& s : info.species) maxRadius = std::max(maxRadius, s.radius);
    maxRadius *= settings.radiusScale;
    const double screenArea = static_cast<double>(width) * height;
    // Per visible chunk: `cover` atoms fill its screen area; `want` is what it
    // should get (all atoms once they are about a pixel or larger on screen).
    struct Visible { std::size_t chunk; double want; double cover; float distance; float area; };
    std::vector<Visible> visible;
    visible.reserve(m_chunks.size());
    std::vector<std::size_t> tested;  // in the frustum and outside the camera: occlusion-tested
    tested.reserve(m_chunks.size());
    double totalWanted = 0;
    for (std::size_t i = 0; i < info.chunks.size(); ++i) {
        const auto& chunk = info.chunks[i];
        const glm::vec3 lo(chunk.lower[0], chunk.lower[1], chunk.lower[2]), hi(chunk.upper[0], chunk.upper[1], chunk.upper[2]);
        const glm::vec3 pad(maxRadius);
        const glm::vec3 a = lo - pad, b = hi + pad;
        // Frustum test on the box corners (clip space), with screen extent.
        int outside[6] = {0, 0, 0, 0, 0, 0};
        bool behind = false;
        float minX = 1e30f, minY = 1e30f, maxX = -1e30f, maxY = -1e30f;
        for (int k = 0; k < 8; ++k) {
            const glm::vec4 corner(k & 1 ? b.x : a.x, k & 2 ? b.y : a.y, k & 4 ? b.z : a.z, 1.0f);
            const glm::vec4 clip = viewProjection * corner;
            if (clip.x < -clip.w) ++outside[0];
            if (clip.x > clip.w) ++outside[1];
            if (clip.y < -clip.w) ++outside[2];
            if (clip.y > clip.w) ++outside[3];
            if (clip.z < -clip.w) ++outside[4];
            if (clip.z > clip.w) ++outside[5];
            if (clip.w <= 1e-6f) { behind = true; continue; }
            const float x = clip.x / clip.w, y = clip.y / clip.w;
            minX = std::min(minX, x); maxX = std::max(maxX, x);
            minY = std::min(minY, y); maxY = std::max(maxY, y);
        }
        if (std::any_of(std::begin(outside), std::end(outside), [](int n) { return n == 8; })) { m_occluded[i] = 0; continue; }
        if (!behind) tested.push_back(i);
        else m_occluded[i] = 0;
        // Hidden behind atoms drawn last frame: neither drawn nor streamed.
        if (settings.occlusionCulling && m_occluded[i]) continue;
        double area = screenArea;
        if (!behind) {
            const double w = (std::clamp(maxX, -1.0f, 1.0f) - std::clamp(minX, -1.0f, 1.0f)) * 0.5 * width;
            const double h = (std::clamp(maxY, -1.0f, 1.0f) - std::clamp(minY, -1.0f, 1.0f)) * 0.5 * height;
            area = std::max(1.0, w * h);
        }
        const glm::vec3 centre = 0.5f * (lo + hi);
        const float distance = std::max(1e-3f, -(view * glm::vec4(centre, 1.0f)).z);
        const double count = static_cast<double>(chunk.count);
        double want = count, cover = count;
        if (settings.levelOfDetail && !behind) {
            const float halfDiagonal = 0.5f * glm::length(hi - lo);
            const float near = std::max(1e-3f, distance - halfDiagonal), far = distance + halfDiagonal;
            const double nearPixels = maxRadius * pixelsPerUnit / (orthographic ? 1.0f : near);
            const double farPixels = maxRadius * pixelsPerUnit / (orthographic ? 1.0f : far);
            // Enough atoms to cover the chunk's screen area `density` times over...
            const double atomArea = 3.14159265 * std::max(0.5, nearPixels) * std::max(0.5, nearPixels);
            cover = std::min(count, std::max(256.0, settings.density * area / atomArea));
            // ...rising to every atom as atoms grow from 0.75 to 1.5 pixels in radius,
            // where the arrangement of atoms becomes visible.
            const double resolved = std::clamp((farPixels - 0.75) / 0.75, 0.0, 1.0);
            want = cover + (count - cover) * resolved;
        }
        visible.push_back({i, want, cover, distance, static_cast<float>(area)});
        totalWanted += want;
    }
    // Point and memory budgets: every chunk keeps its cover; the remaining budget
    // gives the nearest chunks full detail first.
    const double memoryAtoms = settings.memoryBudgetMB * 1024.0 * 1024.0 / sizeof(PackedAtom) * 0.85;
    const double limit = std::min(settings.pointBudget, memoryAtoms);
    if (totalWanted > limit) {
        double covers = 0;
        for (const auto& v : visible) covers += std::min(v.want, v.cover);
        if (covers > limit) {
            const double scale = limit / covers;
            for (auto& v : visible) v.want = std::max(256.0, std::min(v.want, v.cover) * scale);
        } else {
            std::vector<std::size_t> order(visible.size());
            for (std::size_t k = 0; k < order.size(); ++k) order[k] = k;
            std::sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return visible[x].distance < visible[y].distance; });
            double remaining = limit - covers;
            for (std::size_t k : order) {
                auto& v = visible[k];
                const double base = std::min(v.want, v.cover);
                const double extra = std::min(v.want - base, remaining);
                v.want = base + extra;
                remaining -= extra;
            }
        }
    }

    // Stream what is missing: the largest gaps on screen first.
    std::vector<std::pair<double, std::size_t>> missing;
    bool complete = true;
    for (const auto& v : visible) {
        ChunkState& c = m_chunks[v.chunk];
        c.wanted = static_cast<std::uint64_t>(v.want);
        c.lastUsed = m_frame;
        if (c.resident < c.wanted) complete = false;
        if (c.requested < c.wanted) missing.push_back({v.area * (1.0 - static_cast<double>(c.resident) / std::max(1.0, v.want)), v.chunk});
    }
    std::sort(missing.begin(), missing.end(), [](const auto& x, const auto& y) { return x.first > y.first; });
    for (const auto& [priority, chunk] : missing) {
        if (m_inFlight >= kMaxInFlight) break;
        const std::uint64_t target = (m_chunks[chunk].wanted + kBlock - 1) / kBlock * kBlock;
        request(chunk, target);
    }
    evict(m_frame);

    // GPU time of the draw, read a frame later so that it never stalls.
    const int timer = static_cast<int>(m_frame & 1);
    if (m_timerUsed[timer]) {
        GLuint64 nanoseconds = 0;
        glGetQueryObjectui64v(m_timers[timer], GL_QUERY_RESULT, &nanoseconds);
        m_stats.gpuMilliseconds = static_cast<double>(nanoseconds) * 1e-6;
    }
    glBeginQuery(GL_TIME_ELAPSED, m_timers[timer]);
    m_timerUsed[timer] = true;

    // Draw front to back (cheap depth rejection of hidden atoms).
    std::sort(visible.begin(), visible.end(), [](const Visible& x, const Visible& y) { return x.distance < y.distance; });
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_DEPTH_TEST);
    glUseProgram(m_program);
    glUniformMatrix4fv(glGetUniformLocation(m_program, "uView"), 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(m_program, "uProjection"), 1, GL_FALSE, glm::value_ptr(projection));
    glUniform1f(glGetUniformLocation(m_program, "uPixelsPerUnit"), pixelsPerUnit);
    glUniform1i(glGetUniformLocation(m_program, "uOrthographic"), orthographic ? 1 : 0);
    glUniform1f(glGetUniformLocation(m_program, "uRadiusScale"), settings.radiusScale);
    glUniform1f(glGetUniformLocation(m_program, "uMaxPointSize"), static_cast<float>(m_maxPointSize));
    glUniform1f(glGetUniformLocation(m_program, "uAmbient"), 0.38f);
    const glm::vec3 light = glm::normalize(glm::vec3(0.45f, 0.6f, 0.65f));
    glUniform3fv(glGetUniformLocation(m_program, "uLight"), 1, glm::value_ptr(light));
    glUniform1i(glGetUniformLocation(m_program, "uColourAxis"), colourAxis);
    if (colourAxis > 0) {
        const std::size_t k = static_cast<std::size_t>(colourAxis - 1);
        glUniform2f(glGetUniformLocation(m_program, "uColourRange"), static_cast<float>(info.lower[k]), static_cast<float>(info.upper[k]));
    }
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, m_speciesTexture);
    glUniform1i(glGetUniformLocation(m_program, "uSpecies"), 0);
    const GLint lowerLocation = glGetUniformLocation(m_program, "uLower");
    const GLint scaleLocation = glGetUniformLocation(m_program, "uScale");
    const GLint fillLocation = glGetUniformLocation(m_program, "uFill");
    glBindVertexArray(m_vao);
    glEnableVertexAttribArray(0);
    std::uint64_t drawn = 0, resident = 0;
    for (const auto& v : visible) {
        const ChunkState& c = m_chunks[v.chunk];
        const std::uint64_t n = std::min<std::uint64_t>(c.resident, static_cast<std::uint64_t>(std::ceil(v.want)));
        if (!n) continue;
        const auto& chunk = info.chunks[v.chunk];
        glUniform3f(lowerLocation, static_cast<float>(chunk.lower[0]), static_cast<float>(chunk.lower[1]), static_cast<float>(chunk.lower[2]));
        glUniform3f(scaleLocation, static_cast<float>((chunk.upper[0] - chunk.lower[0]) / 65535.0),
                    static_cast<float>((chunk.upper[1] - chunk.lower[1]) / 65535.0), static_cast<float>((chunk.upper[2] - chunk.lower[2]) / 65535.0));
        // A subsample is drawn with larger points so that the chunk still looks solid.
        const double fraction = static_cast<double>(n) / static_cast<double>(chunk.count);
        glUniform1f(fillLocation, static_cast<float>(std::clamp(std::cbrt(1.0 / fraction), 1.0, 2.5)));
        glBindBuffer(GL_ARRAY_BUFFER, c.buffer);
        glVertexAttribIPointer(0, 4, GL_UNSIGNED_SHORT, sizeof(PackedAtom), nullptr);
        glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(n));
        drawn += n;
    }
    for (const auto& c : m_chunks) resident += c.resident;
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // Occlusion queries: each chunk's box against the atoms drawn now.
    std::size_t occluded = 0;
    if (settings.occlusionCulling) {
        glUseProgram(m_lineProgram);
        glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
        glDepthMask(GL_FALSE);
        glBindVertexArray(m_boxVao);
        const GLint mvpLocation = glGetUniformLocation(m_lineProgram, "uMVP");
        for (std::size_t i : tested) {
            if (m_occluded[i]) ++occluded;
            if (m_queryPending[i]) continue;
            const auto& chunk = info.chunks[i];
            const glm::vec3 a = glm::vec3(chunk.lower[0], chunk.lower[1], chunk.lower[2]) - glm::vec3(maxRadius);
            const glm::vec3 b = glm::vec3(chunk.upper[0], chunk.upper[1], chunk.upper[2]) + glm::vec3(maxRadius);
            glm::mat4 model(1.0f);
            model[0][0] = b.x - a.x; model[1][1] = b.y - a.y; model[2][2] = b.z - a.z;
            model[3] = glm::vec4(a, 1.0f);
            const glm::mat4 mvp = viewProjection * model;
            glUniformMatrix4fv(mvpLocation, 1, GL_FALSE, glm::value_ptr(mvp));
            glBeginQuery(GL_SAMPLES_PASSED, m_queries[i]);
            glDrawArrays(GL_TRIANGLES, 0, 36);
            glEndQuery(GL_SAMPLES_PASSED);
            m_queryPending[i] = 1;
        }
        glBindVertexArray(0);
        glDepthMask(GL_TRUE);
        glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    }

    if (settings.showBox) {
        glUseProgram(m_lineProgram);
        const glm::mat4 mvp = projection * view;
        glUniformMatrix4fv(glGetUniformLocation(m_lineProgram, "uMVP"), 1, GL_FALSE, glm::value_ptr(mvp));
        const glm::vec3 colour = lightTheme ? glm::vec3(0.25f) : glm::vec3(0.85f);
        glUniform3fv(glGetUniformLocation(m_lineProgram, "uColour"), 1, glm::value_ptr(colour));
        glBindVertexArray(m_lineVao);
        glDrawArrays(GL_LINES, 0, 24);
        glBindVertexArray(0);
    }
    glUseProgram(0);
    glEndQuery(GL_TIME_ELAPSED);

    m_stats.drawnAtoms = drawn;
    m_stats.residentAtoms = resident;
    m_stats.visibleChunks = visible.size();
    m_stats.occludedChunks = occluded;
    m_stats.loadingChunks = m_inFlight;
    m_stats.residentMB = static_cast<double>(m_residentBytes) / (1024.0 * 1024.0);
    m_stats.complete = complete && m_inFlight == 0;
}
