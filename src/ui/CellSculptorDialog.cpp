#include "ui/CellSculptorDialog.h"
#include "algorithms/CellSculptorAlgo.h"
#include "app/SceneView.h"
#include "camera/Camera.h"
#include "graphics/CellSculptorShaders.h"
#include "graphics/CylinderMesh.h"
#include "graphics/Renderer.h"
#include "graphics/Shader.h"
#include "graphics/ShadowMap.h"
#include "graphics/SphereMesh.h"
#include "graphics/StructureInstanceBuilder.h"
#include "io/StructureLoader.h"
#include "ui/ThemeUtils.h"
#include "util/ElementData.h"
#include "imgui.h"

#include <GL/glew.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <sstream>

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------
namespace
{
std::string scBaseName(const std::string& path)
{
    const auto sep = path.find_last_of("/\\");
    return (sep == std::string::npos) ? path : path.substr(sep + 1);
}

// Cycling default colours for new slabs.
static const float kSlabColors[5][3] = {
    {0.40f, 0.70f, 1.00f},
    {1.00f, 0.50f, 0.30f},
    {0.40f, 1.00f, 0.50f},
    {1.00f, 0.90f, 0.30f},
    {0.90f, 0.40f, 1.00f},
};
} // namespace

// ===========================================================================
// Lifecycle
// ===========================================================================
CellSculptorDialog::CellSculptorDialog() = default;

CellSculptorDialog::~CellSculptorDialog()
{
    if (m_sourceFBO)      glDeleteFramebuffers(1, &m_sourceFBO);
    if (m_sourceColorTex) glDeleteTextures(1, &m_sourceColorTex);
    if (m_sourceDepthRbo) glDeleteRenderbuffers(1, &m_sourceDepthRbo);

    if (m_resultFBO)      glDeleteFramebuffers(1, &m_resultFBO);
    if (m_resultColorTex) glDeleteTextures(1, &m_resultColorTex);
    if (m_resultDepthRbo) glDeleteRenderbuffers(1, &m_resultDepthRbo);

    if (m_previewShadow.depthFBO)
        glDeleteFramebuffers(1, &m_previewShadow.depthFBO);
    if (m_previewShadow.depthTexture)
        glDeleteTextures(1, &m_previewShadow.depthTexture);

    if (m_planeProgram) glDeleteProgram(m_planeProgram);
    if (m_planeVAO)     glDeleteVertexArrays(1, &m_planeVAO);
    if (m_planeVBO)     glDeleteBuffers(1, &m_planeVBO);
    if (m_edgeVAO)      glDeleteVertexArrays(1, &m_edgeVAO);
    if (m_edgeVBO)      glDeleteBuffers(1, &m_edgeVBO);

    delete m_previewSphere;
    delete m_previewCylinder;
}

void CellSculptorDialog::initRenderResources(Renderer& renderer)
{
    m_renderer        = &renderer;
    m_previewSphere   = new SphereMesh(24, 24);
    m_previewCylinder = new CylinderMesh(16);

    m_sourceBuffers.init(m_previewSphere->vbo, m_previewSphere->ebo,
                         m_previewSphere->indexCount,
                         m_previewCylinder->vbo, m_previewCylinder->vertexCount);
    m_ghostBuffers.init(m_previewSphere->vbo, m_previewSphere->ebo,
                         m_previewSphere->indexCount,
                         m_previewCylinder->vbo, m_previewCylinder->vertexCount);
    m_resultBuffers.init(m_previewSphere->vbo, m_previewSphere->ebo,
                         m_previewSphere->indexCount,
                         m_previewCylinder->vbo, m_previewCylinder->vertexCount);

    m_previewShadow = createShadowMap(1, 1);
    // Clear the shadow map depth to 1.0 (far plane = no shadow).
    // Without this the uninitialized 0.0 depth makes every atom appear
    // fully in shadow → rendered black.
    {
        GLint prevFbo = 0;
        glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
        glBindFramebuffer(GL_FRAMEBUFFER, m_previewShadow.depthFBO);
        glClearDepth(1.0);
        glClear(GL_DEPTH_BUFFER_BIT);
        glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
    }

    // Plane face shader + VAO
    m_planeProgram = createProgram(cellSculptorPlaneVertexShader(),
                                   cellSculptorPlaneFragmentShader());
    glGenVertexArrays(1, &m_planeVAO);
    glGenBuffers(1, &m_planeVBO);
    glBindVertexArray(m_planeVAO);
    glBindBuffer(GL_ARRAY_BUFFER, m_planeVBO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glBindVertexArray(0);

    // Edge line VAO
    glGenVertexArrays(1, &m_edgeVAO);
    glGenBuffers(1, &m_edgeVBO);
    glBindVertexArray(m_edgeVAO);
    glBindBuffer(GL_ARRAY_BUFFER, m_edgeVBO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glBindVertexArray(0);

    m_glReady = true;
}

// ===========================================================================
// Menu / file-drop entry points
// ===========================================================================
void CellSculptorDialog::drawMenuItem(bool /*enabled*/)
{
    if (ImGui::MenuItem("Cell Sculptor", nullptr, false, true))
        m_openRequested = true;
}

void CellSculptorDialog::feedDroppedFile(const std::string& path)
{
    m_pendingDropPaths.push_back(path);
}

// ===========================================================================
// Source loading
// ===========================================================================
bool CellSculptorDialog::loadSourceFromPath(const std::string& path)
{
    Structure loaded;
    std::string error;
    if (!loadStructureFromFile(path, loaded, error))
    {
        m_status = "Failed to load: " + (error.empty() ? path : error);
        return false;
    }
    m_source         = std::move(loaded);
    m_sourceName     = scBaseName(path);
    m_hasSource      = true;
    m_sourceBufDirty = true;
    m_facesDirty     = true;
    m_resultDirty    = true;
    m_srcCameraFit   = true;
    m_resCameraFit   = true;

    std::ostringstream oss;
    oss << "Loaded " << m_sourceName
        << " (" << m_source.atoms.size() << " atoms)";
    m_status = oss.str();
    return true;
}

void CellSculptorDialog::loadSourceFromScene(const Structure& s, const std::string& name)
{
    if (s.atoms.empty()) { m_status = "Scene is empty."; return; }
    m_source         = s;
    m_sourceName     = name;
    m_hasSource      = true;
    m_sourceBufDirty = true;
    m_facesDirty     = true;
    m_resultDirty    = true;
    m_srcCameraFit   = true;
    m_resCameraFit   = true;
    m_status         = "Using scene (" + std::to_string(s.atoms.size()) + " atoms)";
}
// Plane face geometry
// ===========================================================================
void CellSculptorDialog::rebuildFaces()
{
    m_faceBatches.clear();
    m_facesDirty = false;
    if (!m_hasSource || m_slabs.empty()) return;

    const Structure sc_ref = cscBuildSupercell(m_source, m_nx, m_ny, m_nz);
    const glm::vec3 centroid = cscCentroid(sc_ref);

    const size_t N = m_slabs.size();

    // -----------------------------------------------------------------------
    // 3-slab bounded case: build face quads using the EXACT same corner +
    // cell-vector math as cscApplySlabs so that the planes are guaranteed to
    // align with the yellow wireframe and the atom positions.
    // -----------------------------------------------------------------------
    if (N == 3)
    {
        struct SD { glm::vec3 n, ref; float d1, d2; };
        std::array<SD, 3> sd;
        for (int i = 0; i < 3; ++i)
            sd[i] = { cscNormal(m_slabs[i], m_source),
                      cscSlabRef(m_slabs[i], m_source, centroid),
                      cscD1(m_slabs[i], m_source),
                      cscD2(m_slabs[i], m_source) };

        // Solve the 3×3 system  A·p = rhs  where row i is slab i's normal.
        const glm::mat3 A = glm::transpose(glm::mat3(sd[0].n, sd[1].n, sd[2].n));
        if (std::abs(glm::determinant(A)) > 1e-6f)
        {
            const glm::vec3 rhs(sd[0].d1 + glm::dot(sd[0].ref, sd[0].n),
                                sd[1].d1 + glm::dot(sd[1].ref, sd[1].n),
                                sd[2].d1 + glm::dot(sd[2].ref, sd[2].n));

            // p0 = corner where all three d1 planes meet.
            const glm::vec3 p0 = glm::inverse(A) * rhs;
            // Edge vectors along each slab normal (same as cscApplySlabs cellVectors).
            const glm::vec3 a = sd[0].n * (sd[0].d2 - sd[0].d1);
            const glm::vec3 b = sd[1].n * (sd[1].d2 - sd[1].d1);
            const glm::vec3 c = sd[2].n * (sd[2].d2 - sd[2].d1);

            // 8 corners of the parallelepiped.
            const glm::vec3 c000 = p0,         c100 = p0+a,     c010 = p0+b,
                             c001 = p0+c,       c110 = p0+a+b,   c101 = p0+a+c,
                             c011 = p0+b+c,     c111 = p0+a+b+c;

            // Build one parallelogram face: 2 triangles + 4 edges.
            // Winding: CCW when viewed from outside (outward direction).
            auto addFace = [&](glm::vec3 v0, glm::vec3 v1,
                                glm::vec3 v2, glm::vec3 v3,
                                glm::vec3 normal, int slabIdx)
            {
                ScFaceBatch batch;
                batch.normal = normal;
                batch.color  = glm::vec3(m_slabs[slabIdx].color[0],
                                         m_slabs[slabIdx].color[1],
                                         m_slabs[slabIdx].color[2]);
                batch.triangleVerts = { v0,v1,v2, v0,v2,v3 };
                batch.edgeVerts     = { v0,v1, v1,v2, v2,v3, v3,v0 };
                m_faceBatches.push_back(std::move(batch));
            };

            // Slab 0 (a-direction): two faces perpendicular to n0.
            addFace(c000, c010, c011, c001, -sd[0].n, 0); // d1 face, outward = -n0
            addFace(c100, c110, c111, c101,  sd[0].n, 0); // d2 face, outward = +n0
            // Slab 1 (b-direction): two faces perpendicular to n1.
            addFace(c000, c001, c101, c100, -sd[1].n, 1);
            addFace(c010, c110, c111, c011,  sd[1].n, 1);
            // Slab 2 (c-direction): two faces perpendicular to n2.
            addFace(c000, c100, c110, c010, -sd[2].n, 2);
            addFace(c001, c101, c111, c011,  sd[2].n, 2);
            return;
        }
    }

    // -----------------------------------------------------------------------
    // Fallback for < 3 slabs (or degenerate 3-slab case): Sutherland-Hodgman
    // clipping of a large initial quad against all slab halfspaces.
    // -----------------------------------------------------------------------
    float maxExtent = 5.0f;
    for (const auto& a : sc_ref.atoms)
        maxExtent = std::max(maxExtent,
                             glm::length(glm::vec3((float)a.x,(float)a.y,(float)a.z) - centroid));
    maxExtent = maxExtent * 1.4f + 3.0f;

    for (size_t si = 0; si < N; ++si)
    {
        const CellSlabPlane& slab = m_slabs[si];

        const glm::vec3 n_i   = cscNormal(slab, m_source);
        const float     d1_i  = cscD1(slab, m_source);
        const float     d2_i  = cscD2(slab, m_source);
        const glm::vec3 ref_i = cscSlabRef(slab, m_source, centroid);
        const glm::vec3 col(slab.color[0], slab.color[1], slab.color[2]);

        auto buildFace = [&](float faceD, const glm::vec3& faceNormal, bool isLowerFace)
        {
            const glm::vec3 up = (std::abs(glm::dot(n_i, glm::vec3(0,1,0))) < 0.9f)
                                 ? glm::vec3(0,1,0) : glm::vec3(1,0,0);
            const glm::vec3 u = glm::normalize(up - glm::dot(up, n_i) * n_i);
            const glm::vec3 v = glm::cross(n_i, u);
            const glm::vec3 origin = ref_i + n_i * faceD;

            std::vector<glm::vec3> poly = {
                origin + ( u + v) * maxExtent,
                origin + ( u - v) * maxExtent,
                origin + (-u - v) * maxExtent,
                origin + (-u + v) * maxExtent,
            };

            for (size_t sj = 0; sj < N && !poly.empty(); ++sj)
            {
                const glm::vec3 nj   = cscNormal(m_slabs[sj], m_source);
                const float     d1j  = cscD1(m_slabs[sj], m_source);
                const float     d2j  = cscD2(m_slabs[sj], m_source);
                const glm::vec3 refj = cscSlabRef(m_slabs[sj], m_source, centroid);

                if (sj == si)
                {
                    if (isLowerFace)
                        poly = cscClipHalfspace(poly, -n_i, -d2_i - glm::dot(ref_i, n_i));
                    else
                        poly = cscClipHalfspace(poly,  n_i,  d1_i + glm::dot(ref_i, n_i));
                }
                else
                {
                    poly = cscClipHalfspace(poly,  nj,  d1j + glm::dot(refj, nj));
                    if (!poly.empty())
                        poly = cscClipHalfspace(poly, -nj, -d2j - glm::dot(refj, nj));
                }
            }

            if (poly.size() < 3) return;

            ScFaceBatch batch;
            batch.normal = faceNormal;
            batch.color  = col;
            batch.triangleVerts.reserve((poly.size() - 2) * 3);
            for (size_t k = 1; k + 1 < poly.size(); ++k)
            {
                batch.triangleVerts.push_back(poly[0]);
                batch.triangleVerts.push_back(poly[k]);
                batch.triangleVerts.push_back(poly[k + 1]);
            }
            batch.edgeVerts.reserve(poly.size() * 2);
            for (size_t k = 0; k < poly.size(); ++k)
            {
                batch.edgeVerts.push_back(poly[k]);
                batch.edgeVerts.push_back(poly[(k + 1) % poly.size()]);
            }
            m_faceBatches.push_back(std::move(batch));
        };

        buildFace(d1_i, -n_i, true);
        buildFace(d2_i,  n_i, false);
    }
}

// ===========================================================================
// Buffer uploads
// ===========================================================================
void CellSculptorDialog::uploadToPreview(const Structure& s, SceneBuffers& buffers)
{
    if (!m_glReady || s.atoms.empty()) return;
    const auto radii = makeLiteratureCovalentRadii();
    std::vector<float> shininess(119, 40.0f);
    static const int kIdent[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    const StructureInstanceData data =
        buildStructureInstanceData(s, false, kIdent, radii, shininess);
    std::array<bool,119> noFilter = {};
    buffers.upload(data, false, noFilter);
}

void CellSculptorDialog::rebuildSourceBuffers()
{
    if (!m_hasSource || m_source.atoms.empty()) return;
    m_supercell = cscBuildSupercell(m_source, m_nx, m_ny, m_nz);

    // Build instance data once and upload to both buffers.
    // Both share the same geometry; applySlabDimming only patches scale VBOs.
    //
    // PBC boundary images ARE included here (m_supercell passed directly, not a
    // hasUnitCell=false copy) so that the source preview shows a visually complete
    // crystal cell with face/edge/vertex continuity.  applySlabDimming classifies
    // each instance via m_insideMask[atomIndices[i]], and appendPbcBoundaryImages
    // copies the source atom's atomIndices entry to its mirror instances (line 110 of
    // StructureInstanceBuilder.cpp), so every PBC copy is classified identically to
    // its canonical atom — no position-based re-test needed.
    const auto radii = makeLiteratureCovalentRadii();
    std::vector<float> shininess(119, 40.0f);
    static const int kIdent[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    const StructureInstanceData data =
        buildStructureInstanceData(m_supercell, false, kIdent, radii, shininess);
    std::array<bool,119> noFilter = {};
    m_sourceBuffers.upload(data, false, noFilter);
    m_ghostBuffers.upload(data, false, noFilter);

    m_sourceBufDirty = false;
    // Do NOT call applySlabDimming here: the inside mask (m_insideMask) is
    // computed and owned by rebuildResult().  Since m_sourceBufDirty=true also
    // implies m_resultDirty=true, rebuildResult (and therefore applySlabDimming)
    // will always run in the same frame immediately after this function.
    // The GPU scale VBOs start at full-radius (all opaque) from upload() above,
    // so no visual artifact occurs before rebuildResult runs.
}

void CellSculptorDialog::rebuildResult()
{
    if (!m_hasSource) return;
    if (m_supercell.atoms.empty())
        m_supercell = cscBuildSupercell(m_source, m_nx, m_ny, m_nz);

    // Compute the inside/outside mask with the EXACT same loop used inside
    // cscApplySlabs.  This is the single source of truth for dimming:
    // applySlabDimming() simply looks up m_insideMask[i] and never re-runs
    // the slab test, so the dimming is guaranteed to match the cut result.
    {
        const glm::vec3 centroid = cscCentroid(m_supercell);
        struct SD { glm::vec3 n, ref; float d1, d2; bool periodic; };
        std::vector<SD> sds;
        sds.reserve(m_slabs.size());
        for (const auto& slab : m_slabs)
            sds.push_back({ cscNormal(slab, m_source),
                             cscSlabRef(slab, m_source, centroid),
                             cscD1(slab, m_source),
                             cscD2(slab, m_source),
                             slab.usePeriodic && m_source.hasUnitCell });

        // Build the mask per GPU instance (canonical atoms + PBC mirror copies).
        //
        // appendPbcBoundaryImages wraps canonical positions to [0,1) in the
        // supercell AND snaps near-1.0 values to 0.0 (wrapAndSnapFractional).
        // This can move a canonical atom's stored position to the OPPOSITE face,
        // causing a different slab projection than what cscApplySlabs sees.
        //
        // Solution: classify canonical instances (i < N_canonical) with the
        // original unmodified position from m_supercell.atoms[i], which is what
        // cscApplySlabs actually uses.  Mirror instances (i >= N_canonical) were
        // appended by appendPbcBoundaryImages at canonicalPos + a/b/c_supercell;
        // their only position record is atomPositions[i], which IS correct.
        const size_t nCanonical = m_supercell.atoms.size();
        const size_t nInstances = (size_t)m_sourceBuffers.atomCount;
        const bool havePositions = (nInstances > 0 &&
                                    m_sourceBuffers.atomPositions.size() == nInstances);

        m_insideMask.assign(nInstances, true);
        for (size_t i = 0; i < nInstances; ++i)
        {
            glm::vec3 p;
            if (i < nCanonical)
            {
                // Canonical atom: use original position (not wrapped by appendPbc).
                const auto& a = m_supercell.atoms[i];
                p = glm::vec3((float)a.x, (float)a.y, (float)a.z);
            }
            else if (havePositions)
            {
                // PBC mirror atom: use rendered position (canonicalPos + lattice shift).
                p = m_sourceBuffers.atomPositions[i];
            }
            else
            {
                // No position data for mirror; inherit canonical atom classification.
                const size_t ci = (m_sourceBuffers.atomIndices.size() == nInstances)
                                  ? (size_t)m_sourceBuffers.atomIndices[i] : i;
                m_insideMask[i] = (ci < m_insideMask.size()) ? m_insideMask[ci] : true;
                continue;
            }
            for (const auto& sd : sds)
            {
                const float proj = glm::dot(p - sd.ref, sd.n);
                if (proj < sd.d1 - kCscSlabTol)
                    { m_insideMask[i] = false; break; }
                if (sd.periodic ? (proj >= sd.d2 - kCscSlabTol)
                                : (proj >  sd.d2 + kCscSlabTol))
                    { m_insideMask[i] = false; break; }
            }
        }
    }

    m_previewResult = cscApplySlabs(m_supercell, m_slabs, m_source);
    m_resultDirty   = false;
    uploadToPreview(m_previewResult, m_resultBuffers);
    applySlabDimming();
}

void CellSculptorDialog::applySlabDimming()
{
    // Render all atoms uniformly semi-transparent: hide the opaque pass
    // (source scale = 0) and show all atoms through the ghost (transparent)
    // pass.  The ghost alpha slider controls the opacity level.
    if (!m_glReady || m_sourceBuffers.atomCount == 0) return;

    const size_t n          = (size_t)m_sourceBuffers.atomCount;
    const auto&  fullScales = m_sourceBuffers.atomRadii;
    if (fullScales.size() != n) return;

    // No slabs: show all atoms opaque, ghost layer empty.
    if (m_slabs.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, m_sourceBuffers.scaleVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0,
                        (GLsizeiptr)(n * sizeof(float)), fullScales.data());
        if (m_ghostBuffers.atomCount == n)
        {
            const std::vector<float> zeros(n, 0.0f);
            glBindBuffer(GL_ARRAY_BUFFER, m_ghostBuffers.scaleVBO);
            glBufferSubData(GL_ARRAY_BUFFER, 0,
                            (GLsizeiptr)(n * sizeof(float)), zeros.data());
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        return;
    }

    // Slabs active: hide the opaque pass, show all atoms as ghost.
    const std::vector<float> zeros(n, 0.0f);
    glBindBuffer(GL_ARRAY_BUFFER, m_sourceBuffers.scaleVBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, (GLsizeiptr)(n * sizeof(float)), zeros.data());
    if (m_ghostBuffers.atomCount == n)
    {
        glBindBuffer(GL_ARRAY_BUFFER, m_ghostBuffers.scaleVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, (GLsizeiptr)(n * sizeof(float)), fullScales.data());
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// ===========================================================================
// Camera fit
// ===========================================================================
void CellSculptorDialog::autoFitCamera(float& yaw, float& pitch, float& dist,
                                       const SceneBuffers& buffers)
{
    yaw   = 45.0f;
    pitch = 35.0f;
    if (buffers.atomCount == 0) { dist = 10.0f; return; }
    float maxR = 0.0f;
    for (size_t i = 0; i < buffers.atomPositions.size(); ++i)
    {
        const float r = (i < buffers.atomRadii.size()) ? buffers.atomRadii[i] : 0.0f;
        maxR = std::max(maxR,
                        glm::length(buffers.atomPositions[i] - buffers.orbitCenter) + r);
    }
    maxR = std::max(maxR, 1.0f);
    const float halfFov = glm::radians(22.5f);
    dist = std::max(Camera::kMinDistance,
                    std::min(Camera::kMaxDistance, maxR / std::sin(halfFov) * 1.15f));
}

// ===========================================================================
// FBO management
// ===========================================================================
void CellSculptorDialog::ensureFBO(GLuint& fbo, GLuint& colorTex, GLuint& depthRbo,
                                   int& storedW, int& storedH, int w, int h)
{
    if (w == storedW && h == storedH && fbo != 0) return;
    if (fbo)      { glDeleteFramebuffers(1,  &fbo);      fbo      = 0; }
    if (colorTex) { glDeleteTextures(1,       &colorTex); colorTex = 0; }
    if (depthRbo) { glDeleteRenderbuffers(1,  &depthRbo); depthRbo = 0; }

    glGenTextures(1, &colorTex);
    glBindTexture(GL_TEXTURE_2D, colorTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    glGenRenderbuffers(1, &depthRbo);
    glBindRenderbuffer(GL_RENDERBUFFER, depthRbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                           GL_TEXTURE_2D, colorTex, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER, depthRbo);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    storedW = w;
    storedH = h;
}

// ===========================================================================
// Render to FBO
// ===========================================================================
void CellSculptorDialog::renderToFBO(int w, int h, GLuint fbo,
                                     SceneBuffers& buffers,
                                     float yaw, float pitch, float dist,
                                     FrameView& lastFrame, bool drawPlanes)
{
    if (!m_glReady || !m_renderer) return;

    GLint prevFbo = 0; glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
    GLint prevVP[4]; glGetIntegerv(GL_VIEWPORT, prevVP);

    Camera cam;
    cam.yaw = yaw; cam.pitch = pitch; cam.distance = dist;

    FrameView frame;
    frame.framebufferWidth  = w;
    frame.framebufferHeight = h;
    if (buffers.atomCount > 0)
        buildFrameView(cam, buffers, true, frame);
    else
    {
        frame.projection     = glm::perspective(glm::radians(45.0f), (float)w/(float)h, 0.1f, 500.0f);
        frame.view           = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, -dist));
        frame.lightMVP       = glm::mat4(1.0f);
        frame.lightPosition  = glm::vec3(10, 10, 10);
        frame.cameraPosition = glm::vec3(0, 0, dist);
    }
    lastFrame = frame;

    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glViewport(0, 0, w, h);
    glEnable(GL_DEPTH_TEST);
    const ImVec4& bg = ImGui::GetStyle().Colors[ImGuiCol_WindowBg];
    glClearColor(bg.x, bg.y, bg.z, bg.w);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (buffers.atomCount > 0)
    {
        m_renderer->drawBonds(frame.projection, frame.view,
                              frame.lightPosition, frame.cameraPosition,
                              buffers.tabCylinderVAO, buffers.tabCylinderVertexCount,
                              buffers.bondCount);
        m_renderer->drawAtoms(frame.projection, frame.view,
                              frame.lightMVP, frame.lightPosition, frame.cameraPosition,
                              m_previewShadow,
                              buffers.tabSphereVAO, buffers.tabSphereIndexCount,
                              buffers.atomCount);
        m_renderer->drawBoxLines(frame.projection, frame.view,
                                 buffers.lineVAO, buffers.boxLines.size());

        // Ghost pass: draw outside atoms semi-transparently so inside atoms
        // show through them.  Only for the source preview (drawPlanes=true).
        // Uses GL_CONSTANT_ALPHA blend so no shader changes are needed:
        //   result = outside_atom * ghostAlpha + framebuffer * (1 - ghostAlpha)
        // With depth-test ON and depth-write OFF, outside atoms that are in front
        // of inside atoms blend over them; outside atoms behind are clipped.
        if (drawPlanes && m_ghostBuffers.atomCount > 0)
        {
            glEnable(GL_BLEND);
            glBlendColor(0.0f, 0.0f, 0.0f, m_ghostAlpha);
            glBlendFunc(GL_CONSTANT_ALPHA, GL_ONE_MINUS_CONSTANT_ALPHA);
            glDepthMask(GL_FALSE); // don't overwrite depth so inside atoms remain "in front"
            m_renderer->drawAtoms(frame.projection, frame.view,
                                  frame.lightMVP, frame.lightPosition, frame.cameraPosition,
                                  m_previewShadow,
                                  m_ghostBuffers.tabSphereVAO, m_ghostBuffers.tabSphereIndexCount,
                                  m_ghostBuffers.atomCount);
            glDepthMask(GL_TRUE);
            glDisable(GL_BLEND);
        }
    }
    // Draw plane faces (semi-transparent) — only for the source preview.
    // pushes face geometry slightly behind the atom surface to prevent z-fighting
    // at atom/plane intersections.
    if (drawPlanes && m_planeProgram && !m_faceBatches.empty())
    {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
        glDisable(GL_CULL_FACE);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0f, 1.0f);

        glUseProgram(m_planeProgram);
        glUniformMatrix4fv(glGetUniformLocation(m_planeProgram, "projection"),
                           1, GL_FALSE, glm::value_ptr(frame.projection));
        glUniformMatrix4fv(glGetUniformLocation(m_planeProgram, "view"),
                           1, GL_FALSE, glm::value_ptr(frame.view));
        glUniform3f(glGetUniformLocation(m_planeProgram, "lightPos"),
                    frame.lightPosition.x, frame.lightPosition.y, frame.lightPosition.z);
        glUniform3f(glGetUniformLocation(m_planeProgram, "viewPos"),
                    frame.cameraPosition.x, frame.cameraPosition.y, frame.cameraPosition.z);

        glBindVertexArray(m_planeVAO);
        glBindBuffer(GL_ARRAY_BUFFER, m_planeVBO);
        for (const ScFaceBatch& batch : m_faceBatches)
        {
            if (batch.triangleVerts.empty()) continue;
            glBufferData(GL_ARRAY_BUFFER,
                         (GLsizeiptr)(batch.triangleVerts.size() * sizeof(glm::vec3)),
                         batch.triangleVerts.data(), GL_DYNAMIC_DRAW);
            glUniform3f(glGetUniformLocation(m_planeProgram, "faceColor"),
                        batch.color.r, batch.color.g, batch.color.b);
            glUniform3f(glGetUniformLocation(m_planeProgram, "faceNormal"),
                        batch.normal.x, batch.normal.y, batch.normal.z);
            glDrawArrays(GL_TRIANGLES, 0, (GLsizei)batch.triangleVerts.size());
        }
        glBindVertexArray(0);
        glUseProgram(0);

        glDepthMask(GL_TRUE);
        glDisable(GL_BLEND);
        glDisable(GL_POLYGON_OFFSET_FILL);

        // Draw wireframe edges — always on top of atoms so cutting planes are
        // clearly visible regardless of which atoms sit on the plane.
        glDisable(GL_DEPTH_TEST);
        glBindVertexArray(m_edgeVAO);
        glBindBuffer(GL_ARRAY_BUFFER, m_edgeVBO);
        for (const ScFaceBatch& batch : m_faceBatches)
        {
            if (batch.edgeVerts.empty()) continue;
            glBufferData(GL_ARRAY_BUFFER,
                         (GLsizeiptr)(batch.edgeVerts.size() * sizeof(glm::vec3)),
                         batch.edgeVerts.data(), GL_DYNAMIC_DRAW);
            m_renderer->drawBoxLines(frame.projection, frame.view,
                                     m_edgeVAO, batch.edgeVerts.size(),
                                     glm::vec3(batch.color) * 0.9f);
        }
        glBindVertexArray(0);
        glEnable(GL_DEPTH_TEST);
    }

    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
    glViewport(prevVP[0], prevVP[1], prevVP[2], prevVP[3]);
}

void CellSculptorDialog::renderSourceToFBO(int w, int h)
{
    ensureFBO(m_sourceFBO, m_sourceColorTex, m_sourceDepthRbo, m_sourceW, m_sourceH, w, h);
    renderToFBO(w, h, m_sourceFBO, m_sourceBuffers,
                m_srcYaw, m_srcPitch, m_srcDist, m_srcLastFrame, true);
}

void CellSculptorDialog::renderResultToFBO(int w, int h)
{
    ensureFBO(m_resultFBO, m_resultColorTex, m_resultDepthRbo, m_resultW, m_resultH, w, h);
    renderToFBO(w, h, m_resultFBO, m_resultBuffers,
                m_resYaw, m_resPitch, m_resDist, m_resLastFrame, false);
}

// ===========================================================================
// Main dialog
// ===========================================================================
void CellSculptorDialog::drawDialog(
    Structure& structure,
    const std::function<void(Structure&)>& updateBuffers)
{
    // Consume pending file drops
    if (!m_pendingDropPaths.empty())
    {
        loadSourceFromPath(m_pendingDropPaths.front());
        m_pendingDropPaths.clear();
        rebuildSourceBuffers();
        rebuildFaces();
        rebuildResult();
    }

    if (m_openRequested)
    {
        ImGui::OpenPopup("Cell Sculptor");
        m_openRequested = false;
    }

    m_isOpen = ImGui::IsPopupOpen("Cell Sculptor");

    ImGui::SetNextWindowSize(ImVec2(1200.0f, 840.0f), ImGuiCond_FirstUseEver);
    bool dialogOpen = true;
    if (!ImGui::BeginPopupModal("Cell Sculptor", &dialogOpen, 0))
    {
        m_isOpen = false;
        return;
    }
    m_isOpen = true;

    // Lazy rebuild
    if (m_sourceBufDirty) rebuildSourceBuffers();
    if (m_facesDirty)     rebuildFaces();
    if (m_resultDirty)    rebuildResult();

    if (m_srcCameraFit && m_sourceBuffers.atomCount > 0)
    { autoFitCamera(m_srcYaw, m_srcPitch, m_srcDist, m_sourceBuffers); m_srcCameraFit = false; }
    if (m_resCameraFit && m_resultBuffers.atomCount > 0)
    { autoFitCamera(m_resYaw, m_resPitch, m_resDist, m_resultBuffers); m_resCameraFit = false; }

    // ------------------------------------------------------------------
    // Layout: Left = two stacked previews, Right = scrollable controls
    // ------------------------------------------------------------------
    constexpr float kLeftW       = 580.0f;
    constexpr float kColumnH     = 760.0f;
    constexpr float kTopPanelH   = 400.0f;
    const float     kBotPanelH   = kColumnH - kTopPanelH - ImGui::GetStyle().ItemSpacing.y;
    constexpr float kViewPad     = 5.0f;

    // ===================================================================
    // LEFT: stacked previews
    // ===================================================================
    ImGui::BeginGroup();

    // -- TOP: Supercell + planes --
    ImGui::BeginChild("##scSrcPanel", ImVec2(kLeftW, kTopPanelH), true);
    ImGui::TextUnformatted("Supercell + Cutting Planes");
    ImGui::SameLine();
    if (ImGui::Button("Use Scene##scuse"))
    {
        loadSourceFromScene(structure, "scene");
        rebuildSourceBuffers(); rebuildFaces(); rebuildResult();
    }
    if (m_hasSource)
    {
        ImGui::SameLine();
        if (ImGui::Button("Clear##scclear"))
        {
            m_source = {}; m_hasSource = false;
            m_sourceBufDirty = true; m_facesDirty = true;
            m_resultDirty = true; m_previewResult = {}; m_status.clear();
        }
        ImGui::SameLine();
        ImGui::TextDisabled("| %s  %dx%dx%d", m_sourceName.c_str(), m_nx, m_ny, m_nz);
    }
    ImGui::Separator();

    const float srcViewH = kTopPanelH - ImGui::GetCursorPosY()
                           + ImGui::GetStyle().WindowPadding.y - 4.0f;
    ImGui::InvisibleButton("##scSrcView", ImVec2(-1.0f, srcViewH > 20.f ? srcViewH : 20.f));
    ImDrawList* dl = ImGui::GetWindowDrawList();
    ImVec2 vMin = ImGui::GetItemRectMin(), vMax = ImGui::GetItemRectMax();
    const bool srcHov = ImGui::IsItemHovered(), srcAct = ImGui::IsItemActive();

    if (!m_hasSource)
    {
        dl->AddRect(vMin, vMax, ImGui::GetColorU32(ImGuiCol_Border), 4.f, 0, 1.5f);
        const char* h1 = "Drop a structure file here";
        const char* h2 = "or use  'Use Scene'  above";
        const ImVec2 mid((vMin.x+vMax.x)*0.5f, (vMin.y+vMax.y)*0.5f);
        const ImU32 dc = ImGui::GetColorU32(ImGuiCol_TextDisabled);
        const float lh = ImGui::GetTextLineHeight();
        dl->AddText(ImVec2(mid.x - ImGui::CalcTextSize(h1).x*0.5f, mid.y - lh*1.1f), dc, h1);
        dl->AddText(ImVec2(mid.x - ImGui::CalcTextSize(h2).x*0.5f, mid.y + lh*0.1f), dc, h2);
    }
    else if (m_glReady)
    {
        const int pw = std::max(1, (int)(vMax.x - vMin.x - kViewPad*2));
        const int ph = std::max(1, (int)(vMax.y - vMin.y - kViewPad*2));
        renderSourceToFBO(pw, ph);
        dl->AddImage((ImTextureID)(intptr_t)m_sourceColorTex,
                     ImVec2(vMin.x+kViewPad, vMin.y+kViewPad),
                     ImVec2(vMin.x+kViewPad+pw, vMin.y+kViewPad+ph), ImVec2(0,1), ImVec2(1,0));
        if (srcAct && ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.f))
        { const ImVec2 d = ImGui::GetIO().MouseDelta; m_srcYaw -= d.x*0.5f; m_srcPitch += d.y*0.5f; }
        if (srcHov)
        {
            const float wh = ImGui::GetIO().MouseWheel;
            if (wh != 0.f) {
                m_srcDist -= wh * m_srcDist * 0.1f;
                m_srcDist = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, m_srcDist));
            }
        }
        dl->AddText(ImVec2(vMin.x+6, vMax.y - ImGui::GetTextLineHeight()-4),
                    ImGui::GetColorU32(ImGuiCol_TextDisabled), "Drag=orbit  Scroll=zoom");
    }
    ImGui::EndChild();

    // -- BOTTOM: Cut result --
    ImGui::BeginChild("##scResPanel", ImVec2(kLeftW, kBotPanelH), true);
    ImGui::TextUnformatted("Cut Result");
    if (!m_previewResult.atoms.empty())
    {
        ImGui::SameLine();
        ImGui::TextDisabled("| %d atoms", (int)m_previewResult.atoms.size());
    }
    ImGui::Separator();

    const float resViewH = kBotPanelH - ImGui::GetCursorPosY()
                           + ImGui::GetStyle().WindowPadding.y - 4.0f;
    ImGui::InvisibleButton("##scResView", ImVec2(-1.0f, resViewH > 20.f ? resViewH : 20.f));
    ImDrawList* rdl = ImGui::GetWindowDrawList();
    ImVec2 rMin = ImGui::GetItemRectMin(), rMax = ImGui::GetItemRectMax();
    const bool resHov = ImGui::IsItemHovered(), resAct = ImGui::IsItemActive();

    if (m_glReady && !m_previewResult.atoms.empty())
    {
        const int pw = std::max(1, (int)(rMax.x - rMin.x - kViewPad*2));
        const int ph = std::max(1, (int)(rMax.y - rMin.y - kViewPad*2));
        renderResultToFBO(pw, ph);
        rdl->AddImage((ImTextureID)(intptr_t)m_resultColorTex,
                      ImVec2(rMin.x+kViewPad, rMin.y+kViewPad),
                      ImVec2(rMin.x+kViewPad+pw, rMin.y+kViewPad+ph), ImVec2(0,1), ImVec2(1,0));
        if (resAct && ImGui::IsMouseDragging(ImGuiMouseButton_Left, 0.f))
        { const ImVec2 d = ImGui::GetIO().MouseDelta; m_resYaw -= d.x*0.5f; m_resPitch += d.y*0.5f; }
        if (resHov)
        {
            const float wh = ImGui::GetIO().MouseWheel;
            if (wh != 0.f) {
                m_resDist -= wh * m_resDist * 0.1f;
                m_resDist = std::max(Camera::kMinDistance, std::min(Camera::kMaxDistance, m_resDist));
            }
        }
        rdl->AddText(ImVec2(rMin.x+6, rMax.y - ImGui::GetTextLineHeight()-4),
                     ImGui::GetColorU32(ImGuiCol_TextDisabled), "Drag=orbit  Scroll=zoom");
    }
    else
    {
        const char* hint = m_hasSource ? "Add cutting slabs to see the result"
                                       : "Load a structure first";
        const ImVec2 mid((rMin.x+rMax.x)*0.5f, (rMin.y+rMax.y)*0.5f);
        rdl->AddText(ImVec2(mid.x - ImGui::CalcTextSize(hint).x*0.5f,
                            mid.y - ImGui::GetTextLineHeight()*0.5f),
                     ImGui::GetColorU32(ImGuiCol_TextDisabled), hint);
    }
    ImGui::EndChild();

    ImGui::EndGroup();

    // ===================================================================
    // RIGHT: Controls
    // ===================================================================
    ImGui::SameLine();
    ImGui::BeginChild("##scCtrlOuter", ImVec2(0.0f, kColumnH), true);
    ImGui::TextUnformatted("Cell Sculptor — Controls");
    ImGui::Separator();

    // Inner scrollable region, leaving room for Generate/Close buttons
    ImGui::BeginChild("##scCtrlScroll", ImVec2(-1.0f, -46.0f), false);

    if (!m_status.empty()) ImGui::TextWrapped("%s", m_status.c_str());

    // ---- Supercell ----
    ImGui::Spacing();
    // Compute exact height: top/bottom window padding + title row + separator +
    // table header row + table drag-int row + one item-spacing gap.
    const float scPadY  = ImGui::GetStyle().WindowPadding.y;
    const float scFHS   = ImGui::GetFrameHeightWithSpacing();
    const float scFH    = ImGui::GetFrameHeight();
    const float scSepH  = 1.0f + ImGui::GetStyle().ItemSpacing.y;
    const float supercellBoxH = scPadY * 2.0f
                              + ImGui::GetTextLineHeightWithSpacing()  // title
                              + scSepH                                 // separator
                              + scFHS                                  // table header labels
                              + scFH;                                  // table drag-int row
    ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, 4.0f);
    ImGui::BeginChild("##scSupercellPanel", ImVec2(-1.0f, supercellBoxH), true,
                      ImGuiWindowFlags_NoScrollbar);
    ImGui::TextUnformatted("Supercell");
    ImGui::SameLine(); ImGui::TextDisabled("(tiles +/- around origin)");
    ImGui::Separator();
    {
        bool sc = false;
        if (ImGui::BeginTable("##sctbl", 3, ImGuiTableFlags_None))
        {
            ImGui::TableSetupColumn("##cx", ImGuiTableColumnFlags_WidthFixed, 50.0f);
            ImGui::TableSetupColumn("##cy", ImGuiTableColumnFlags_WidthFixed, 50.0f);
            ImGui::TableSetupColumn("##cz", ImGuiTableColumnFlags_WidthFixed, 50.0f);
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::TextDisabled("x");
            ImGui::TableNextColumn(); ImGui::TextDisabled("y");
            ImGui::TableNextColumn(); ImGui::TextDisabled("z");
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); sc |= ImGui::DragInt("##scx", &m_nx, 0.1f, 1, 30);
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); sc |= ImGui::DragInt("##scy", &m_ny, 0.1f, 1, 30);
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); sc |= ImGui::DragInt("##scz", &m_nz, 0.1f, 1, 30);
            ImGui::EndTable();
        }
        m_nx = std::max(1, m_nx); m_ny = std::max(1, m_ny); m_nz = std::max(1, m_nz);
        if (sc) { m_sourceBufDirty = true; m_facesDirty = true; m_resultDirty = true; m_srcCameraFit = true; }
    }
    ImGui::EndChild();
    ImGui::PopStyleVar();

    // ---- Cutting Slabs ----
    ImGui::Spacing();
    ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, 4.0f);
    ImGui::BeginChild("##scSlabsPanel", ImVec2(-1.0f, 0.0f), true,
                      ImGuiWindowFlags_NoScrollbar);
    ImGui::TextUnformatted("Cutting Slabs");
    ImGui::SameLine(); ImGui::TextDisabled("(atoms kept within [d1,d2] along hkl)");
    ImGui::Separator();
    {
        const bool hasPeriodic = m_hasSource && m_source.hasUnitCell;

        for (size_t si = 0; si < m_slabs.size(); ++si)
        {
            CellSlabPlane& slab = m_slabs[si];
            ImGui::PushID((int)si);
            bool changed = false;

            // Header row: color swatch | label | Remove
            float colorArr[3] = {slab.color[0], slab.color[1], slab.color[2]};
            if (ImGui::ColorEdit3("##sc_col", colorArr,
                                  ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel))
            {
                slab.color[0] = colorArr[0]; slab.color[1] = colorArr[1]; slab.color[2] = colorArr[2];
                m_facesDirty = true;
            }
            ImGui::SameLine();
            ImGui::TextUnformatted("Slab"); ImGui::SameLine();
            ImGui::Text("%d", (int)si + 1);
            ImGui::SameLine();
            // Right-align Remove without clipping
            {
                const ImVec2 fp  = ImGui::GetStyle().FramePadding;
                const float  removeW = ImGui::CalcTextSize("Remove").x + fp.x * 2.0f;
                ImGui::SetCursorPosX(ImGui::GetCursorPosX()
                    + ImGui::GetContentRegionAvail().x - removeW);
            }
            if (ImGui::Button("Remove##sr"))
            {
                m_slabs.erase(m_slabs.begin() + (int)si);
                m_facesDirty = true; m_resultDirty = true;
                ImGui::PopID(); break;
            }

            // hkl row with labels above
            if (ImGui::BeginTable("##hkl", 3, ImGuiTableFlags_None))
            {
                ImGui::TableSetupColumn("##ch", ImGuiTableColumnFlags_WidthFixed, 50.0f);
                ImGui::TableSetupColumn("##ck", ImGuiTableColumnFlags_WidthFixed, 50.0f);
                ImGui::TableSetupColumn("##cl", ImGuiTableColumnFlags_WidthFixed, 50.0f);
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::TextDisabled("h");
                ImGui::TableNextColumn(); ImGui::TextDisabled("k");
                ImGui::TableNextColumn(); ImGui::TextDisabled("l");
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); changed |= ImGui::DragInt("##h", &slab.h, 0.15f, -12, 12);
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); changed |= ImGui::DragInt("##k", &slab.k, 0.15f, -12, 12);
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); changed |= ImGui::DragInt("##l", &slab.l, 0.15f, -12, 12);
                ImGui::EndTable();
            }

            // Periodic checkbox — always visible; disabled when source has no unit cell.
            if (!hasPeriodic)
            {
                ImGui::BeginDisabled();
                slab.usePeriodic = false;
            }
            if (ImGui::Checkbox("Periodic##scper", &slab.usePeriodic))
                changed = true;
            if (!hasPeriodic)
            {
                ImGui::EndDisabled();
                if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
                    ImGui::SetTooltip("Load a crystal with a unit cell to enable periodic mode.");
            }

            if (slab.usePeriodic && hasPeriodic)
            {
                const float d_rep  = cscStructuralDhkl(slab, m_source);    // A→A repeat
                const float d_atom = cscAtomLayerSpacing(slab, m_source); // A→B layer
                ImGui::SameLine();
                if (std::abs(d_rep - d_atom) > d_rep * 0.01f)
                    ImGui::TextDisabled("d=%g A (layer %g A)", d_rep, d_atom);
                else
                    ImGui::TextDisabled("d=%g A", d_rep);
                ImGui::SetNextItemWidth(65.0f);
                if (ImGui::DragInt("Start##scsp", &slab.startPlane, 0.1f, -200, 200))
                    changed = true;
                ImGui::SameLine();
                ImGui::SetNextItemWidth(50.0f);
                if (ImGui::DragInt("N##scnp", &slab.nPeriods, 0.1f, 1, 200))
                { slab.nPeriods = std::max(1, slab.nPeriods); changed = true; }
                const float eff1 = cscD1(slab, m_source), eff2 = cscD2(slab, m_source);
                ImGui::TextDisabled("d1=%.3f A   d2=%.3f A   t=%.3f A", eff1, eff2, eff2-eff1);
            }
            else
            {
                ImGui::SetNextItemWidth(70.0f);
                changed |= ImGui::DragFloat("d1##sc", &slab.d1, 0.05f, -500.f, 500.f, "%.2f");
                ImGui::SameLine();
                ImGui::SetNextItemWidth(70.0f);
                changed |= ImGui::DragFloat("d2##sc", &slab.d2, 0.05f, -500.f, 500.f, "%.2f");
            }

            if (changed)
            {
                if (!slab.usePeriodic && slab.d2 < slab.d1 + 0.05f) slab.d2 = slab.d1 + 0.05f;
                m_facesDirty = true; m_resultDirty = true;
            }
            ImGui::Separator();
            ImGui::PopID();
        }

        // Add-slab row
        ImGui::Spacing();
        if (ImGui::BeginTable("##nhkl", 3, ImGuiTableFlags_None))
        {
            ImGui::TableSetupColumn("##nch", ImGuiTableColumnFlags_WidthFixed, 50.0f);
            ImGui::TableSetupColumn("##nck", ImGuiTableColumnFlags_WidthFixed, 50.0f);
            ImGui::TableSetupColumn("##ncl", ImGuiTableColumnFlags_WidthFixed, 50.0f);
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::TextDisabled("h");
            ImGui::TableNextColumn(); ImGui::TextDisabled("k");
            ImGui::TableNextColumn(); ImGui::TextDisabled("l");
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); ImGui::DragInt("##nh", &m_newH, 0.15f, -12, 12);
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); ImGui::DragInt("##nk", &m_newK, 0.15f, -12, 12);
            ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-1); ImGui::DragInt("##nl", &m_newL, 0.15f, -12, 12);
            ImGui::EndTable();
        }
        ImGui::SetNextItemWidth(70.0f);
        ImGui::DragFloat("d1##nnd1", &m_newD1, 0.1f, -500.f, 500.f, "%.1f");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(70.0f);
        ImGui::DragFloat("d2##nnd2", &m_newD2, 0.1f, -500.f, 500.f, "%.1f");
        if (ImGui::Button("+ Add Slab##nnadd", ImVec2(-1.0f, 0.0f)))
        {
            if (m_newD2 < m_newD1 + 0.05f) m_newD2 = m_newD1 + 0.05f;
            CellSlabPlane slab;
            slab.h = m_newH; slab.k = m_newK; slab.l = m_newL;
            slab.d1 = m_newD1; slab.d2 = m_newD2;
            slab.usePeriodic = hasPeriodic; slab.startPlane = 0; slab.nPeriods = 1;
            const size_t ci = m_slabs.size() % 5;
            slab.color[0] = kSlabColors[ci][0]; slab.color[1] = kSlabColors[ci][1]; slab.color[2] = kSlabColors[ci][2];
            m_slabs.push_back(slab);
            m_facesDirty = true; m_resultDirty = true;
        }
    } // Cutting Slabs
    ImGui::EndChild();
    ImGui::PopStyleVar();

    // ---- View options ----
    ImGui::Spacing();
    ImGui::Separator();
    ImGui::TextUnformatted("Outside atom opacity");
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::SliderFloat("##scghostalpha", &m_ghostAlpha, 0.0f, 1.0f, "%.2f");

    // ---- Status ----
    ImGui::Spacing();
    ImGui::Separator();
    const bool bounded = cscIsBounded(m_slabs, m_source);
    if (m_slabs.empty())
        ImGui::TextDisabled("Add >=3 non-parallel slabs to close the region.");
    else if (!bounded)
        ImGui::TextColored(themeStatusWarn(), "Region not closed — add more slabs.");
    else
    {
        ImGui::TextColored(themeStatusGood(), "Closed region.");
        if (!m_previewResult.atoms.empty())
            ImGui::TextDisabled("%d atoms in result", (int)m_previewResult.atoms.size());
    }

    ImGui::EndChild(); // ##scCtrlScroll

    // Fixed bottom buttons
    ImGui::Separator();
    if (!bounded) ImGui::BeginDisabled();
    const bool doGenerate = ImGui::Button("Generate##scgen", ImVec2(100.0f, 0.0f));
    if (!bounded)
    {
        ImGui::EndDisabled();
        if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled))
            ImGui::SetTooltip("Close the region first.");
    }
    ImGui::SameLine();
    if (ImGui::Button("Close##scclose", ImVec2(80.0f, 0.0f)))
        ImGui::CloseCurrentPopup();

    ImGui::EndChild(); // ##scCtrlOuter

    if (doGenerate && bounded && !m_previewResult.atoms.empty())
    {
        structure = m_previewResult;
        updateBuffers(structure);
        ImGui::CloseCurrentPopup();
    }

    ImGui::EndPopup();
}
