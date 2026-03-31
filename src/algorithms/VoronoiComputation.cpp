#include "algorithms/VoronoiComputation.h"
#include "math/StructureMath.h"

#include <algorithm>
#include <cmath>

namespace
{
// A convex polyhedron represented as a set of faces (convex polygons).
struct ConvexPoly
{
    std::vector<VoronoiFace> faces;
};

// Order polygon vertices around their centroid so they form a convex loop.
void sortFaceVertices(VoronoiFace& face)
{
    if (face.vertices.size() < 3)
        return;

    glm::vec3 center(0.0f);
    for (const auto& v : face.vertices)
        center += v;
    center /= (float)face.vertices.size();

    // Compute face normal from first three vertices.
    glm::vec3 normal = glm::normalize(
        glm::cross(face.vertices[1] - face.vertices[0],
                    face.vertices[2] - face.vertices[0]));

    glm::vec3 ref = glm::normalize(face.vertices[0] - center);
    glm::vec3 binormal = glm::cross(normal, ref);

    std::sort(face.vertices.begin(), face.vertices.end(),
        [&](const glm::vec3& a, const glm::vec3& b) {
            glm::vec3 da = a - center;
            glm::vec3 db = b - center;
            float angA = std::atan2(glm::dot(da, binormal), glm::dot(da, ref));
            float angB = std::atan2(glm::dot(db, binormal), glm::dot(db, ref));
            return angA < angB;
        });
}

// Add point if not already close to an existing one.
void addUnique(std::vector<glm::vec3>& pts, const glm::vec3& p, float eps2 = 1e-8f)
{
    for (const auto& q : pts)
    {
        if (glm::dot(p - q, p - q) <= eps2)
            return;
    }
    pts.push_back(p);
}

// Clip a convex polygon (list of vertices in order) against a half-plane:
// keeps the side where dot(v - planePoint, planeNormal) <= 0.
std::vector<glm::vec3> clipPolygonByPlane(const std::vector<glm::vec3>& poly,
                                          const glm::vec3& planePoint,
                                          const glm::vec3& planeNormal)
{
    if (poly.size() < 3)
        return {};

    std::vector<glm::vec3> out;
    out.reserve(poly.size() + 2);

    for (size_t i = 0; i < poly.size(); ++i)
    {
        const glm::vec3& curr = poly[i];
        const glm::vec3& next = poly[(i + 1) % poly.size()];

        float dCurr = glm::dot(curr - planePoint, planeNormal);
        float dNext = glm::dot(next - planePoint, planeNormal);

        bool currInside = dCurr <= 1e-6f;
        bool nextInside = dNext <= 1e-6f;

        if (currInside)
            out.push_back(curr);

        if (currInside != nextInside)
        {
            float denom = dCurr - dNext;
            if (std::abs(denom) > 1e-10f)
            {
                float t = dCurr / denom;
                t = std::max(0.0f, std::min(1.0f, t));
                out.push_back(curr + t * (next - curr));
            }
        }
    }

    return (out.size() >= 3) ? out : std::vector<glm::vec3>();
}

// Clip a convex polyhedron against a half-plane (dot(v - planePoint, planeNormal) <= 0).
// Returns the clipped polyhedron. The cut face (new face from the clipping) is added.
ConvexPoly clipPolyhedronByPlane(const ConvexPoly& poly,
                                 const glm::vec3& planePoint,
                                 const glm::vec3& planeNormal)
{
    ConvexPoly result;
    std::vector<glm::vec3> cutEdgePoints;

    for (const auto& face : poly.faces)
    {
        auto clipped = clipPolygonByPlane(face.vertices, planePoint, planeNormal);
        if (clipped.size() >= 3)
        {
            VoronoiFace f;
            f.vertices = std::move(clipped);
            result.faces.push_back(std::move(f));
        }
    }

    // Collect intersection points on the cutting plane from all clipped faces.
    for (const auto& face : result.faces)
    {
        for (const auto& v : face.vertices)
        {
            float d = glm::dot(v - planePoint, planeNormal);
            if (std::abs(d) < 1e-4f)
                addUnique(cutEdgePoints, v);
        }
    }

    // Create the new cap face from the cutting plane.
    if (cutEdgePoints.size() >= 3)
    {
        VoronoiFace capFace;
        capFace.vertices = cutEdgePoints;
        sortFaceVertices(capFace);
        result.faces.push_back(std::move(capFace));
    }

    return result;
}

// Build a box-shaped convex polyhedron (6 faces, axis-aligned).
ConvexPoly makeBox(const glm::vec3& lo, const glm::vec3& hi)
{
    ConvexPoly box;
    // bottom (y = lo.y)
    box.faces.push_back({{
        {lo.x, lo.y, lo.z}, {hi.x, lo.y, lo.z},
        {hi.x, lo.y, hi.z}, {lo.x, lo.y, hi.z}}});
    // top (y = hi.y)
    box.faces.push_back({{
        {lo.x, hi.y, lo.z}, {lo.x, hi.y, hi.z},
        {hi.x, hi.y, hi.z}, {hi.x, hi.y, lo.z}}});
    // front (z = lo.z)
    box.faces.push_back({{
        {lo.x, lo.y, lo.z}, {lo.x, hi.y, lo.z},
        {hi.x, hi.y, lo.z}, {hi.x, lo.y, lo.z}}});
    // back (z = hi.z)
    box.faces.push_back({{
        {lo.x, lo.y, hi.z}, {hi.x, lo.y, hi.z},
        {hi.x, hi.y, hi.z}, {lo.x, hi.y, hi.z}}});
    // left (x = lo.x)
    box.faces.push_back({{
        {lo.x, lo.y, lo.z}, {lo.x, lo.y, hi.z},
        {lo.x, hi.y, hi.z}, {lo.x, hi.y, lo.z}}});
    // right (x = hi.x)
    box.faces.push_back({{
        {hi.x, lo.y, lo.z}, {hi.x, hi.y, lo.z},
        {hi.x, hi.y, hi.z}, {hi.x, lo.y, hi.z}}});
    return box;
}

// Build a box-shaped convex polyhedron from cell vectors (parallelepiped).
ConvexPoly makeCellBox(const glm::vec3& origin, const glm::mat3& cell)
{
    // 8 corners of the parallelepiped
    glm::vec3 c[8];
    for (int i = 0; i < 8; ++i)
    {
        float fi = (i & 1) ? 1.0f : 0.0f;
        float fj = (i & 2) ? 1.0f : 0.0f;
        float fk = (i & 4) ? 1.0f : 0.0f;
        c[i] = origin + cell * glm::vec3(fi, fj, fk);
    }

    // 6 faces of the parallelepiped (order vertices consistently)
    ConvexPoly box;
    // face at frac_a=0: corners 0,2,6,4
    box.faces.push_back({{c[0], c[2], c[6], c[4]}});
    // face at frac_a=1: corners 1,5,7,3
    box.faces.push_back({{c[1], c[5], c[7], c[3]}});
    // face at frac_b=0: corners 0,4,5,1
    box.faces.push_back({{c[0], c[4], c[5], c[1]}});
    // face at frac_b=1: corners 2,3,7,6
    box.faces.push_back({{c[2], c[3], c[7], c[6]}});
    // face at frac_c=0: corners 0,1,3,2
    box.faces.push_back({{c[0], c[1], c[3], c[2]}});
    // face at frac_c=1: corners 4,6,7,5
    box.faces.push_back({{c[4], c[6], c[7], c[5]}});
    return box;
}
}

VoronoiDiagram computeVoronoi(const Structure& structure)
{
    VoronoiDiagram diagram;
    if (structure.atoms.empty())
        return diagram;

    const bool usePbc = structure.hasUnitCell;

    glm::mat3 cell(1.0f);
    glm::mat3 invCell(1.0f);
    glm::vec3 origin(0.0f);

    if (usePbc)
    {
        if (!tryMakeCellMatrices(structure, cell, invCell))
            return diagram;

        origin = glm::vec3(
            (float)structure.cellOffset[0],
            (float)structure.cellOffset[1],
            (float)structure.cellOffset[2]);
    }

    // Collect atom positions.
    const size_t N = structure.atoms.size();
    std::vector<glm::vec3> positions(N);
    for (size_t i = 0; i < N; ++i)
    {
        positions[i] = glm::vec3(
            (float)structure.atoms[i].x,
            (float)structure.atoms[i].y,
            (float)structure.atoms[i].z);
    }

    // Compute bounding box for non-PBC case.
    glm::vec3 bboxLo(1e30f);
    glm::vec3 bboxHi(-1e30f);
    for (const auto& p : positions)
    {
        bboxLo = glm::min(bboxLo, p);
        bboxHi = glm::max(bboxHi, p);
    }
    // Expand by 50% to give cells room.
    glm::vec3 extent = bboxHi - bboxLo;
    float pad = std::max({extent.x, extent.y, extent.z}) * 0.5f + 2.0f;
    bboxLo -= glm::vec3(pad);
    bboxHi += glm::vec3(pad);

    // Periodic image offsets to consider.
    const int imageRange = 1;

    diagram.cells.resize(N);

    for (size_t i = 0; i < N; ++i)
    {
        const glm::vec3& pi = positions[i];

        // Start with the bounding volume.
        ConvexPoly cell_poly;
        if (usePbc)
            cell_poly = makeCellBox(origin, cell);
        else
            cell_poly = makeBox(bboxLo, bboxHi);

        // Clip against bisector planes with all neighbors.
        for (size_t j = 0; j < N; ++j)
        {
            if (usePbc)
            {
                // Consider all periodic images of atom j.
                for (int da = -imageRange; da <= imageRange; ++da)
                for (int db = -imageRange; db <= imageRange; ++db)
                for (int dc = -imageRange; dc <= imageRange; ++dc)
                {
                    if (j == i && da == 0 && db == 0 && dc == 0)
                        continue;

                    glm::vec3 shift = cell * glm::vec3((float)da, (float)db, (float)dc);
                    glm::vec3 pj_img = positions[j] + shift;

                    glm::vec3 midpoint = 0.5f * (pi + pj_img);
                    glm::vec3 normal = pj_img - pi;
                    float len = glm::length(normal);
                    if (len < 1e-10f)
                        continue;
                    normal /= len;

                    cell_poly = clipPolyhedronByPlane(cell_poly, midpoint, normal);
                    if (cell_poly.faces.empty())
                        break;
                }
            }
            else
            {
                if (j == i)
                    continue;

                glm::vec3 midpoint = 0.5f * (pi + positions[j]);
                glm::vec3 normal = positions[j] - pi;
                float len = glm::length(normal);
                if (len < 1e-10f)
                    continue;
                normal /= len;

                cell_poly = clipPolyhedronByPlane(cell_poly, midpoint, normal);
            }

            if (cell_poly.faces.empty())
                break;
        }

        // Sort each face's vertices for proper rendering.
        for (auto& face : cell_poly.faces)
            sortFaceVertices(face);

        VoronoiCell vc;
        vc.faces = std::move(cell_poly.faces);
        diagram.cells[i] = std::move(vc);
    }

    return diagram;
}
