#include "ui/SubstitutionalSolidSolutionDialog.h"
#include "algorithms/SubstitutionalSolidSolutionBuilder.h"
#include "camera/Camera.h"
#include "graphics/CylinderMesh.h"
#include "graphics/Renderer.h"
#include "graphics/SphereMesh.h"
#include "graphics/StructureInstanceBuilder.h"
#include "io/StructureLoader.h"
#include "ui/PeriodicTableDialog.h"
#include "util/ElementData.h"
#include "imgui.h"

#include <glm/gtc/matrix_transform.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <numeric>
#include <sstream>

// ===========================================================================
// Construction / destruction
// ===========================================================================

SubstitutionalSolidSolutionDialog::SubstitutionalSolidSolutionDialog()
{
}

SubstitutionalSolidSolutionDialog::~SubstitutionalSolidSolutionDialog()
{
    if (m_previewFBO)      glDeleteFramebuffers(1,  &m_previewFBO);
    if (m_previewColorTex) glDeleteTextures(1,      &m_previewColorTex);
    if (m_previewDepthRbo) glDeleteRenderbuffers(1, &m_previewDepthRbo);

    if (m_previewShadow.depthFBO)     glDeleteFramebuffers(1, &m_previewShadow.depthFBO);
    if (m_previewShadow.depthTexture) glDeleteTextures(1,     &m_previewShadow.depthTexture);

    delete m_previewSphere;
    delete m_previewCylinder;
}

void SubstitutionalSolidSolutionDialog::initRenderResources(Renderer& renderer)
{
    m_renderer        = &renderer;
    m_previewSphere   = new SphereMesh(24, 24);
    m_previewCylinder = new CylinderMesh(16);
    m_previewBuffers.init(m_previewSphere->vao, m_previewCylinder->vao);
    m_previewShadow   = createShadowMap(1, 1);
    m_glReady         = true;
}

// ===========================================================================
// Menu item / drop
// ===========================================================================

void SubstitutionalSolidSolutionDialog::drawMenuItem(bool enabled)
{
    if (ImGui::MenuItem("Substitutional Solid Solution", nullptr, false, enabled))
        m_openRequested = true;
}

void SubstitutionalSolidSolutionDialog::feedDroppedFile(const std::string& path)
{
    m_pendingDropPath = path;
}

// ===========================================================================
// Structure loading helpers
// ===========================================================================

bool SubstitutionalSolidSolutionDialog::loadFromPath(
    const std::string& path,
    const std::vector<float>& radii,
    const std::vector<float>& shininess)
{
    Structure loaded;
    std::string err;
    if (!loadStructureFromFile(path, loaded, err))
    {
        m_statusMsg     = std::string("Load failed: ") + (err.empty() ? path : err);
        m_statusIsError = true;
        return false;
    }

    m_source       = std::move(loaded);
    m_sourceLoaded = true;

    const auto slash = path.find_last_of("/\\");
    m_sourceLabel    = (slash == std::string::npos) ? path : path.substr(slash + 1);

    onStructureLoaded(radii, shininess);
    return true;
}

void SubstitutionalSolidSolutionDialog::loadFromScene(
    const Structure& scene,
    const std::vector<float>& radii,
    const std::vector<float>& shininess)
{
    if (scene.atoms.empty())
    {
        m_statusMsg     = "Current scene has no atoms.";
        m_statusIsError = true;
        return;
    }

    m_source       = scene;
    m_sourceLoaded = true;
    m_sourceLabel  = "scene";

    onStructureLoaded(radii, shininess);
}

void SubstitutionalSolidSolutionDialog::onStructureLoaded(
    const std::vector<float>& radii,
    const std::vector<float>& shininess)
{
    initCompositionFromSource();
    m_previewDirty  = true;
    m_statusIsError = false;

    std::ostringstream msg;
    msg << "Loaded \'" << m_sourceLabel << "\'  ("
        << m_source.atoms.size() << " atoms)";
    m_statusMsg = msg.str();

    if (m_glReady)
    {
        rebuildPreviewBuffers(radii, shininess);
        autoFitPreviewCamera();
    }
}

// ===========================================================================
// Composition helpers
// ===========================================================================

void SubstitutionalSolidSolutionDialog::initCompositionFromSource()
{
    // Count occurrences of each atomic number in the source.
    std::map<int, int> counts;
    for (const AtomSite& a : m_source.atoms)
        counts[a.atomicNumber]++;

    const int total = (int)m_source.atoms.size();
    m_entries.clear();
    for (const auto& kv : counts)
    {
        CompositionEntry e;
        e.atomicNumber = kv.first;
        e.count        = kv.second;
        e.percent      = (total > 0) ? (float)kv.second / (float)total * 100.0f : 0.0f;
        m_entries.push_back(e);
    }
}

void SubstitutionalSolidSolutionDialog::recomputeCountsFromPercents()
{
    const int total = (int)m_source.atoms.size();
    const int n     = (int)m_entries.size();
    if (n == 0 || total == 0) return;

    // Normalise against the sum so any scale of input works.
    float sum = 0.0f;
    for (const auto& e : m_entries) sum += std::max(0.0f, e.percent);

    if (sum <= 0.0f)
    {
        // Distribute evenly when no percents are set yet.
        const int base = total / n;
        const int rem  = total % n;
        for (int i = 0; i < n; ++i)
            m_entries[i].count = base + (i < rem ? 1 : 0);
        return;
    }

    // Real-valued target counts.
    std::vector<float> real(n);
    for (int i = 0; i < n; ++i)
        real[i] = std::max(0.0f, m_entries[i].percent) / sum * (float)total;

    // Largest-remainder method guarantees sum == total.
    int assigned = 0;
    std::vector<int>   floored(n);
    std::vector<float> remainders(n);
    for (int i = 0; i < n; ++i)
    {
        floored[i]    = (int)std::floor(real[i]);
        remainders[i] = real[i] - (float)floored[i];
        assigned += floored[i];
    }
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return remainders[a] > remainders[b];
    });
    const int leftover = total - assigned;
    for (int k = 0; k < leftover && k < n; ++k)
        floored[order[k]]++;
    for (int i = 0; i < n; ++i)
        m_entries[i].count = floored[i];
}

void SubstitutionalSolidSolutionDialog::redistributePercent(int changedIdx, float newPct)
{
    const int n = (int)m_entries.size();
    if (n <= 0) return;

    newPct = std::max(0.0f, std::min(100.0f, newPct));

    if (n == 1)
    {
        m_entries[0].percent = 100.0f;
        recomputeCountsFromPercents();
        return;
    }

    // Sum of all other entries' current percents.
    float otherSum = 0.0f;
    for (int i = 0; i < n; ++i)
        if (i != changedIdx) otherSum += std::max(0.0f, m_entries[i].percent);

    const float remaining = 100.0f - newPct;
    m_entries[changedIdx].percent = newPct;

    if (remaining <= 0.0f)
    {
        for (int i = 0; i < n; ++i)
            if (i != changedIdx) m_entries[i].percent = 0.0f;
    }
    else if (otherSum <= 0.0f)
    {
        // No existing weights — distribute remaining evenly.
        const float each = remaining / (float)(n - 1);
        for (int i = 0; i < n; ++i)
            if (i != changedIdx) m_entries[i].percent = each;
    }
    else
    {
        // Scale others proportionally to keep their mutual ratios fixed.
        const float scale = remaining / otherSum;
        for (int i = 0; i < n; ++i)
            if (i != changedIdx)
                m_entries[i].percent = std::max(0.0f, m_entries[i].percent) * scale;
    }

    recomputeCountsFromPercents();
}

void SubstitutionalSolidSolutionDialog::redistributeCount(int changedIdx,
                                                            int newCount)
{
    const int total = (int)m_source.atoms.size();
    const int n     = (int)m_entries.size();
    if (n <= 0 || total <= 0) return;

    newCount = std::max(0, std::min(total, newCount));
    m_entries[changedIdx].count = newCount;

    // Sum of all other entries.
    int otherSum = 0;
    for (int i = 0; i < n; ++i)
        if (i != changedIdx) otherSum += m_entries[i].count;

    const int remaining = total - newCount;

    if (remaining == 0)
    {
        // Set all others to 0.
        for (int i = 0; i < n; ++i)
            if (i != changedIdx) m_entries[i].count = 0;
        return;
    }

    if (otherSum == 0)
    {
        // Distribute evenly.
        const int others   = n - 1;
        const int base     = remaining / others;
        const int leftover = remaining % others;
        int assigned       = 0;
        for (int i = 0; i < n; ++i)
        {
            if (i == changedIdx) continue;
            m_entries[i].count = base + (assigned < leftover ? 1 : 0);
            ++assigned;
        }
        return;
    }

    // Scale other entries proportionally.
    std::vector<float> scaled(n, 0.0f);
    for (int i = 0; i < n; ++i)
        if (i != changedIdx)
            scaled[i] = (float)m_entries[i].count / (float)otherSum * (float)remaining;

    // Floor + largest-remainder rounding.
    int assignedTotal = 0;
    std::vector<int> floored(n, 0);
    std::vector<float> remainders(n, 0.0f);
    for (int i = 0; i < n; ++i)
    {
        if (i == changedIdx) continue;
        floored[i]    = (int)std::floor(scaled[i]);
        remainders[i] = scaled[i] - (float)floored[i];
        assignedTotal += floored[i];
    }

    int leftover = remaining - assignedTotal;
    if (leftover > 0)
    {
        // Sort indices of other entries by remainder descending.
        std::vector<int> order;
        for (int i = 0; i < n; ++i)
            if (i != changedIdx) order.push_back(i);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return remainders[a] > remainders[b];
        });
        for (int k = 0; k < leftover && k < (int)order.size(); ++k)
            floored[order[k]]++;
    }

    for (int i = 0; i < n; ++i)
        if (i != changedIdx) m_entries[i].count = floored[i];
}

void SubstitutionalSolidSolutionDialog::addEntry()
{
    // Find an unused atomic number.
    int z = 1;
    while (z <= 118)
    {
        bool used = false;
        for (const auto& e : m_entries)
            if (e.atomicNumber == z) { used = true; break; }
        if (!used) break;
        ++z;
    }
    if (z > 118) z = 1;
    CompositionEntry e;
    e.atomicNumber = z;
    e.percent      = 0.0f;
    e.count        = 0;
    m_entries.push_back(e);

    // Give the new entry a 10 at% slice, scaled away from all others
    // proportionally, so the total stays at 100.
    redistributePercent((int)m_entries.size() - 1, 10.0f);
}

void SubstitutionalSolidSolutionDialog::removeEntry(int idx)
{
    if (idx < 0 || idx >= (int)m_entries.size()) return;
    m_entries.erase(m_entries.begin() + idx);
    // Recompute counts from remaining entries' percents.
    recomputeCountsFromPercents();
}

// ===========================================================================
// GL preview helpers
// ===========================================================================

void SubstitutionalSolidSolutionDialog::ensurePreviewFBO(int w, int h)
{
    if (w == m_previewW && h == m_previewH && m_previewFBO != 0)
        return;

    if (m_previewFBO)      { glDeleteFramebuffers(1,  &m_previewFBO);      m_previewFBO = 0; }
    if (m_previewColorTex) { glDeleteTextures(1,      &m_previewColorTex); m_previewColorTex = 0; }
    if (m_previewDepthRbo) { glDeleteRenderbuffers(1, &m_previewDepthRbo); m_previewDepthRbo = 0; }

    glGenTextures(1, &m_previewColorTex);
    glBindTexture(GL_TEXTURE_2D, m_previewColorTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenRenderbuffers(1, &m_previewDepthRbo);
    glBindRenderbuffer(GL_RENDERBUFFER, m_previewDepthRbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    glGenFramebuffers(1, &m_previewFBO);
    glBindFramebuffer(GL_FRAMEBUFFER, m_previewFBO);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                           GL_TEXTURE_2D, m_previewColorTex, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER, m_previewDepthRbo);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    m_previewW = w;
    m_previewH = h;
}

void SubstitutionalSolidSolutionDialog::rebuildPreviewBuffers(
    const std::vector<float>& radii,
    const std::vector<float>& shininess)
{
    if (!m_glReady || m_source.atoms.empty())
        return;

    static const int kIdent[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    const StructureInstanceData data =
        buildStructureInstanceData(m_source, false, kIdent, radii, shininess);

    std::array<bool, 119> noFilter = {};
    m_previewBuffers.upload(data, false, noFilter);
    m_previewDirty = false;
}

void SubstitutionalSolidSolutionDialog::autoFitPreviewCamera()
{
    m_camYaw   = 45.0f;
    m_camPitch = 35.0f;

    if (m_previewBuffers.atomCount == 0)
    {
        m_camDistance = 10.0f;
        return;
    }

    float maxR = 0.0f;
    for (size_t i = 0; i < m_previewBuffers.atomPositions.size(); ++i)
    {
        const float r = (i < m_previewBuffers.atomRadii.size())
                       ? m_previewBuffers.atomRadii[i] : 0.0f;
        const float d = glm::length(m_previewBuffers.atomPositions[i]
                                    - m_previewBuffers.orbitCenter) + r;
        maxR = std::max(maxR, d);
    }
    maxR = std::max(maxR, 1.0f);

    const float halfFov = glm::radians(22.5f);
    float dist = maxR / std::sin(halfFov) * 1.15f;
    dist = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, dist));
    m_camDistance = dist;
}

void SubstitutionalSolidSolutionDialog::renderPreviewToFBO(int w, int h)
{
    if (!m_glReady || !m_renderer || m_previewBuffers.atomCount == 0)
        return;

    ensurePreviewFBO(w, h);

    GLint prevFbo = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
    GLint prevVP[4];
    glGetIntegerv(GL_VIEWPORT, prevVP);

    Camera cam;
    cam.yaw      = m_camYaw;
    cam.pitch    = m_camPitch;
    cam.distance = m_camDistance;

    FrameView frame;
    frame.framebufferWidth  = w;
    frame.framebufferHeight = h;
    buildFrameView(cam, m_previewBuffers, true, frame);

    glBindFramebuffer(GL_FRAMEBUFFER, m_previewFBO);
    glViewport(0, 0, w, h);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.09f, 0.11f, 0.15f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    m_renderer->drawBonds(frame.projection, frame.view,
                          frame.lightPosition, frame.cameraPosition,
                          *m_previewCylinder, m_previewBuffers.bondCount);
    m_renderer->drawAtoms(frame.projection, frame.view,
                          frame.lightMVP, frame.lightPosition, frame.cameraPosition,
                          m_previewShadow, *m_previewSphere,
                          m_previewBuffers.atomCount);
    m_renderer->drawBoxLines(frame.projection, frame.view,
                             m_previewBuffers.lineVAO,
                             m_previewBuffers.boxLines.size());

    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
    glViewport(prevVP[0], prevVP[1], prevVP[2], prevVP[3]);
}

// ===========================================================================
// Main dialog draw
// ===========================================================================

void SubstitutionalSolidSolutionDialog::drawDialog(
    Structure& structure,
    const std::vector<glm::vec3>& /*elementColors*/,
    const std::vector<float>& elementRadii,
    const std::vector<float>& elementShininess,
    const std::function<void(Structure&)>& updateBuffers)
{
    // ---- Consume pending drop ----
    if (!m_pendingDropPath.empty())
    {
        loadFromPath(m_pendingDropPath, elementRadii, elementShininess);
        m_pendingDropPath.clear();
    }

    if (m_openRequested)
    {
        ImGui::OpenPopup("Substitutional Solid Solution");
        m_openRequested = false;
    }

    m_isOpen = ImGui::IsPopupOpen("Substitutional Solid Solution");

    ImGui::SetNextWindowSize(ImVec2(1000.0f, 660.0f), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(720.0f, 480.0f),
                                        ImVec2(3000.0f, 3000.0f));

    bool keepOpen = true;
    if (!ImGui::BeginPopupModal("Substitutional Solid Solution",
                                 &keepOpen,
                                 ImGuiWindowFlags_NoCollapse))
    {
        m_isOpen = false;
        return;
    }
    m_isOpen = true;

    // -----------------------------------------------------------------------
    // Description bar
    // -----------------------------------------------------------------------
    ImGui::TextDisabled(
        "Substitute elements on host lattice sites to match a target composition."
        "  Positions and unit cell are preserved.");
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    const float contentH = ImGui::GetContentRegionAvail().y - 50.0f; // reserve bottom bar
    const float leftW    = ImGui::GetContentRegionAvail().x * 0.60f;

    // =======================================================================
    // LEFT PANEL – host structure loader + 3-D preview
    // =======================================================================
    ImGui::BeginChild("##sss_left", ImVec2(leftW, contentH), false);
    {
        // Compact header row: label + button on the same line.
        ImGui::Text("Host Structure");
        ImGui::SameLine();
        if (ImGui::SmallButton("Use Current Scene"))
            loadFromScene(structure, elementRadii, elementShininess);
        ImGui::Spacing();

        if (m_sourceLoaded)
        {
            ImGui::TextColored(ImVec4(0.40f, 0.90f, 0.40f, 1.0f),
                               "%s", m_sourceLabel.c_str());
            ImGui::SameLine();
            ImGui::TextDisabled("(%d atoms)", (int)m_source.atoms.size());
        }
        else
        {
            ImGui::TextDisabled("No structure loaded yet.");
        }

        ImGui::Spacing();

        // ---- 3-D preview area ----
        const float previewAreaH = ImGui::GetContentRegionAvail().y;
        const float previewAreaW = ImGui::GetContentRegionAvail().x;
        ImVec2 previewMin = ImGui::GetCursorScreenPos();
        ImVec2 previewMax = ImVec2(previewMin.x + previewAreaW,
                                   previewMin.y + previewAreaH);

        // Reserve the area so mouse events register on it.
        ImGui::InvisibleButton("##sss_preview", ImVec2(previewAreaW, previewAreaH));
        const bool previewHovered = ImGui::IsItemHovered();
        const bool previewActive  = ImGui::IsItemActive();

        ImDrawList* dl = ImGui::GetWindowDrawList();

        // Dark background.
        dl->AddRectFilled(previewMin, previewMax, IM_COL32(23, 28, 38, 255));
        dl->AddRect(previewMin, previewMax,
                    IM_COL32(80, 80, 100, 180), 3.0f);

        if (!m_sourceLoaded)
        {
            // Drop-zone hint overlay.
            const char* hint1 = "Drop a structure file here";
            const char* hint2 = "or use  \"Use Current Scene\"";
            const float lh    = ImGui::GetTextLineHeight();
            const ImVec2 mid((previewMin.x + previewMax.x) * 0.5f,
                             (previewMin.y + previewMax.y) * 0.5f);
            const float w1 = ImGui::CalcTextSize(hint1).x;
            const float w2 = ImGui::CalcTextSize(hint2).x;
            const ImU32 hintCol = IM_COL32(130, 140, 160, 200);
            dl->AddText(ImVec2(mid.x - w1 * 0.5f, mid.y - lh * 1.1f), hintCol, hint1);
            dl->AddText(ImVec2(mid.x - w2 * 0.5f, mid.y + lh * 0.1f), hintCol, hint2);
        }
        else
        {
            // Render 3-D preview.
            if (m_glReady)
            {
                if (m_previewDirty)
                    rebuildPreviewBuffers(elementRadii, elementShininess);

                const int pw = std::max(1, (int)previewAreaW);
                const int ph = std::max(1, (int)previewAreaH);
                renderPreviewToFBO(pw, ph);

                dl->AddImage(
                    (ImTextureID)(intptr_t)m_previewColorTex,
                    previewMin, previewMax,
                    ImVec2(0, 1), ImVec2(1, 0));
            }

            // Orbit and zoom via the invisible button.
            if (previewActive && ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.0f))
            {
                const ImVec2 delta = ImGui::GetIO().MouseDelta;
                m_camYaw   -= delta.x * 0.5f;
                m_camPitch += delta.y * 0.5f;
            }
            if (previewHovered)
            {
                const float wheel = ImGui::GetIO().MouseWheel;
                if (wheel != 0.0f)
                {
                    m_camDistance -= wheel * m_camDistance * 0.1f;
                    m_camDistance  = std::max(Camera::kMinDistance,
                                              std::min(Camera::kMaxDistance, m_camDistance));
                }
            }

            // Hint text at bottom of preview.
            const char* orbitHint = "Left-drag = orbit   Scroll = zoom";
            const float hw = ImGui::CalcTextSize(orbitHint).x;
            dl->AddText(ImVec2(previewMin.x + (previewAreaW - hw) * 0.5f,
                               previewMax.y - ImGui::GetTextLineHeight() - 6.0f),
                        IM_COL32(160, 160, 180, 180), orbitHint);
        }
    }
    ImGui::EndChild(); // ##sss_left

    ImGui::SameLine();

    // =======================================================================
    // RIGHT PANEL – composition + options (only shown when source is loaded)
    // =======================================================================
    ImGui::BeginChild("##sss_right", ImVec2(0, contentH), false);
    {
        if (!m_sourceLoaded)
        {
            // Placeholder when no structure is loaded yet.
            ImGui::Spacing();
            ImGui::TextDisabled("Load a host structure on the left to");
            ImGui::TextDisabled("configure the target composition.");
        }
        else
        {
            const int totalAtoms = (int)m_source.atoms.size();

            // ------------------------------------------------------------------
            // Composition section
            // ------------------------------------------------------------------
            ImGui::Text("Target Composition");
            ImGui::SameLine();
            ImGui::TextDisabled("(%d total atoms)", totalAtoms);
            ImGui::Spacing();

            // Column widths.
            constexpr float kSymW  = 56.0f;
            constexpr float kPctW  = 150.0f;
            constexpr float kAtomW = 90.0f;

            // Header row.
            ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetColorU32(ImGuiCol_TextDisabled));
            ImGui::Text("El.");
            ImGui::SameLine(0.0f, kSymW - ImGui::CalcTextSize("El.").x + 8.0f);
            ImGui::Text("at%%");
            ImGui::SameLine(0.0f, kPctW - ImGui::CalcTextSize("at%").x + 8.0f);
            ImGui::Text("atoms");
            ImGui::PopStyleColor();
            ImGui::Separator();
            ImGui::Spacing();

            int  removeIdx  = -1;
            bool openPicker = false;

            for (int i = 0; i < (int)m_entries.size(); ++i)
            {
                CompositionEntry& entry = m_entries[i];
                ImGui::PushID(i);

                // -- Symbol button (opens periodic-table picker) --
                const int z = entry.atomicNumber;
                char symBtn[16];
                std::snprintf(symBtn, sizeof(symBtn), "%s##sym", elementSymbol(z));
                if (ImGui::Button(symBtn, ImVec2(kSymW, 0.0f)))
                {
                    m_pickerTarget = i;
                    openPicker = true;
                }
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("%s", elementName(z));
                ImGui::SameLine(0.0f, 8.0f);

                // -- at% drag-and-edit field --
                // Drag left/right to adjust; Ctrl+click to type exact value.
                // Changing this entry scales all other entries proportionally
                // so the total always stays at 100 at%.
                float pct = entry.percent;
                ImGui::SetNextItemWidth(kPctW);
                if (ImGui::DragFloat("##pct", &pct, 0.15f, 0.0f, 100.0f, "%.2f at%%"))
                    redistributePercent(i, pct);
                if (ImGui::IsItemHovered() && !ImGui::IsItemActive())
                    ImGui::SetTooltip("Drag or Ctrl+click to edit.\nOther elements scale proportionally.");

                ImGui::SameLine(0.0f, 8.0f);

                // -- Derived atom count (read-only) --
                ImGui::Text("%d atoms", entry.count);
                ImGui::SameLine(0.0f, 8.0f);

                // -- Remove button --
                if ((int)m_entries.size() > 1)
                {
                    if (ImGui::SmallButton("X##rm"))
                        removeIdx = i;
                }

                ImGui::PopID();
            }

            if (removeIdx >= 0)
                removeEntry(removeIdx);

            ImGui::Spacing();

            if ((int)m_entries.size() < 118)
            {
                if (ImGui::Button("+ Add Element"))
                    addEntry();
            }

            // Handle periodic-table picker.
            // openPeriodicTable() is intentionally called here (outside PushID loop)
            // so its ID stack depth matches drawPeriodicTable()'s BeginPopupModal call.
            if (openPicker)
                openPeriodicTable();
            {
                std::vector<ElementSelection> sel;
                if (drawPeriodicTable(sel))
                {
                    if (m_pickerTarget >= 0 &&
                        m_pickerTarget < (int)m_entries.size() &&
                        !sel.empty())
                    {
                        m_entries[m_pickerTarget].atomicNumber = sel[0].atomicNumber;
                    }
                    m_pickerTarget = -1;
                }
            }

            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Spacing();

            // ------------------------------------------------------------------
            // Options
            // ------------------------------------------------------------------
            ImGui::Text("Options");
            ImGui::Spacing();
            ImGui::SetNextItemWidth(160.0f);
            ImGui::InputInt("RNG Seed", &m_seed);
            if (m_seed < 0) m_seed = 0;
            ImGui::SameLine();
            ImGui::TextDisabled("(0 = random)");

            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Spacing();

            // ------------------------------------------------------------------
            // Build button
            // ------------------------------------------------------------------
            if (ImGui::Button("Build Solid Solution", ImVec2(-1.0f, 0.0f)))
            {
                SSSParams params;
                params.seed = (unsigned int)m_seed;
                for (const auto& e : m_entries)
                {
                    SSSElementFraction ef;
                    ef.atomicNumber = e.atomicNumber;
                    ef.fraction     = (totalAtoms > 0)
                                     ? (float)e.count / (float)totalAtoms
                                     : 0.0f;
                    params.composition.push_back(ef);
                }

                SSSResult res = buildSubstitutionalSolidSolution(m_source, params);
                if (res.success)
                {
                    structure = std::move(res.output);
                    updateBuffers(structure);
                    m_statusMsg     = res.message;
                    m_statusIsError = false;
                    ImGui::CloseCurrentPopup();
                }
                else
                {
                    m_statusMsg     = res.message;
                    m_statusIsError = true;
                }
            }
        } // if m_sourceLoaded
    }
    ImGui::EndChild(); // ##sss_right

    // =======================================================================
    // Bottom bar – status + close
    // =======================================================================
    ImGui::Separator();
    ImGui::Spacing();

    if (!m_statusMsg.empty())
    {
        if (m_statusIsError)
            ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f),
                               "%s", m_statusMsg.c_str());
        else
            ImGui::TextColored(ImVec4(0.5f, 0.9f, 0.5f, 1.0f),
                               "%s", m_statusMsg.c_str());
        ImGui::SameLine();
    }

    // Push close button to the right.
    const float btnW = 100.0f;
    ImGui::SetCursorPosX(ImGui::GetCursorPosX()
                         + ImGui::GetContentRegionAvail().x - btnW);
    if (ImGui::Button("Close", ImVec2(btnW, 0.0f)))
        ImGui::CloseCurrentPopup();

    ImGui::EndPopup();
}
