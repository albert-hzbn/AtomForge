#include "Picking.h"

#include <glm/gtc/matrix_transform.hpp>
#include <cmath>
#include <algorithm>

namespace
{
void unprojectCursorRayEndpoints(double mx, double my, int w, int h,
                                 const glm::mat4& projection,
                                 const glm::mat4& view,
                                 glm::vec3& nearPoint,
                                 glm::vec3& farPoint)
{
    float ndcX =  2.0f * (float)mx / (float)w - 1.0f;
    float ndcY = -2.0f * (float)my / (float)h + 1.0f;

    glm::mat4 invVP = glm::inverse(projection * view);

    glm::vec4 nearP = invVP * glm::vec4(ndcX, ndcY, -1.0f, 1.0f);
    glm::vec4 farP  = invVP * glm::vec4(ndcX, ndcY,  1.0f, 1.0f);

    nearP /= nearP.w;
    farP  /= farP.w;

    nearPoint = glm::vec3(nearP);
    farPoint = glm::vec3(farP);
}
}

glm::vec3 pickRayOrigin(double mx, double my, int w, int h,
                        const glm::mat4& projection,
                        const glm::mat4& view)
{
    glm::vec3 nearPoint(0.0f);
    glm::vec3 farPoint(0.0f);
    unprojectCursorRayEndpoints(mx, my, w, h, projection, view, nearPoint, farPoint);
    return nearPoint;
}

glm::vec3 pickRayDir(double mx, double my, int w, int h,
                     const glm::mat4& projection,
                     const glm::mat4& view)
{
    glm::vec3 nearPoint(0.0f);
    glm::vec3 farPoint(0.0f);
    unprojectCursorRayEndpoints(mx, my, w, h, projection, view, nearPoint, farPoint);
    return glm::normalize(farPoint - nearPoint);
}

int pickAtom(const glm::vec3& origin, const glm::vec3& dir,
             const std::vector<glm::vec3>& positions,
             const std::vector<float>& radii,
             float fallbackRadius)
{
    int   best  = -1;
    float bestT = 1e30f;
    const float dirLength = glm::length(dir);
    if (!std::isfinite(dirLength) || dirLength <= 0.0f)
        return -1;
    const glm::vec3 rayDir = dir / dirLength;

    for (int i = 0; i < (int)positions.size(); ++i)
    {
        glm::vec3 oc = positions[i] - origin;
        float t = glm::dot(oc, rayDir);

        float radius = fallbackRadius;
        if (i < (int)radii.size() && radii[i] > 0.0f)
            radius = radii[i];

        float d2 = glm::dot(oc, oc) - t * t;
        if (radius <= 0.0f || d2 > radius * radius) continue;
        const float halfChord = std::sqrt(std::max(0.0f, radius * radius - d2));
        if (t + halfChord < 0.0f) continue;
        // Compare visible surfaces, since larger atoms may cover nearer centers.
        const float hitT = std::max(0.0f, t - halfChord);
        if (hitT < bestT)
        {
            best  = i;
            bestT = hitT;
        }
    }
    return best;
}
