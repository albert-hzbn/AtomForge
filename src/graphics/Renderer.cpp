#include "Renderer.h"
#include "Shader.h"

#include <array>
#include <cmath>

#include <glm/gtc/type_ptr.hpp>

namespace
{
std::array<glm::vec4, 6> extractFrustumPlanes(const glm::mat4& vp)
{
    const glm::vec4 row0(vp[0][0], vp[1][0], vp[2][0], vp[3][0]);
    const glm::vec4 row1(vp[0][1], vp[1][1], vp[2][1], vp[3][1]);
    const glm::vec4 row2(vp[0][2], vp[1][2], vp[2][2], vp[3][2]);
    const glm::vec4 row3(vp[0][3], vp[1][3], vp[2][3], vp[3][3]);

    std::array<glm::vec4, 6> planes = {
        row3 + row0, // left
        row3 - row0, // right
        row3 + row1, // bottom
        row3 - row1, // top
        row3 + row2, // near
        row3 - row2  // far
    };

    for (glm::vec4& plane : planes)
    {
        const float len = std::sqrt(plane.x * plane.x + plane.y * plane.y + plane.z * plane.z);
        if (len > 1e-6f)
            plane /= len;
    }

    return planes;
}
}

// ---------------------------------------------------------------------------
// Shader sources
// ---------------------------------------------------------------------------

static const char* kAtomVS = R"(
    #version 130

    in vec3 position;
    in vec3 instancePos;
    in vec3 instanceColor;
    in float instanceScale;
    in float instanceShininess;

    uniform mat4 projection;
    uniform mat4 view;
    uniform mat4 lightMVP;
    uniform vec4 frustumPlanes[6];
    uniform int enableFrustumCulling;

    out vec3 fragColor;
    out vec3 fragWorldPos;
    out vec3 fragNormal;
    out float fragShininess;
    out vec4 FragPosLight;

    void main()
    {
        if (enableFrustumCulling != 0)
        {
            for (int i = 0; i < 6; ++i)
            {
                float signedDistance = dot(frustumPlanes[i].xyz, instancePos) + frustumPlanes[i].w;
                if (signedDistance < -instanceScale)
                {
                    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
                    FragPosLight = vec4(0.0);
                    fragColor = vec3(0.0);
                    fragWorldPos = vec3(0.0);
                    fragNormal = vec3(0.0, 0.0, 1.0);
                    fragShininess = instanceShininess;
                    return;
                }
            }
        }

        vec3 worldPos = position * instanceScale + instancePos;
        gl_Position   = projection * view * vec4(worldPos, 1.0);
        FragPosLight  = lightMVP * vec4(worldPos, 1.0);
        fragColor     = instanceColor;
        fragWorldPos  = worldPos;
        fragNormal    = normalize(position);
        fragShininess = instanceShininess;
    }
)";

static const char* kAtomFS = R"(
    #version 130

    in vec3 fragColor;
    in vec3 fragWorldPos;
    in vec3 fragNormal;
    in float fragShininess;
    in vec4 FragPosLight;

    uniform sampler2D shadowMap;
    uniform vec3 lightPos;
    uniform vec3 viewPos;

    out vec4 color;

    float computeShadow(vec4 pos)
    {
        vec3 proj = pos.xyz / pos.w;
        proj = proj * 0.5 + 0.5;

        if (proj.x < 0.0 || proj.x > 1.0 ||
            proj.y < 0.0 || proj.y > 1.0 ||
            proj.z < 0.0 || proj.z > 1.0)
            return 0.0;

        float closest = texture(shadowMap, proj.xy).r;
        float current = proj.z;
        float bias    = 0.003;

        return (current - bias > closest) ? 1.0 : 0.0;
    }

    void main()
    {
        float shadow  = 0.0;
        vec3 norm = normalize(fragNormal);
        vec3 lightPosFillA = vec3(-lightPos.x, lightPos.y, lightPos.z);
        vec3 lightPosFillB = vec3(lightPos.x, -lightPos.y, lightPos.z);
        vec3 lightPosFillC = vec3(-lightPos.x, 0.65 * lightPos.y, -0.60 * lightPos.z);

        vec3 lightDir0 = normalize(lightPos - fragWorldPos);
        vec3 lightDir1 = normalize(lightPosFillA - fragWorldPos);
        vec3 lightDir2 = normalize(lightPosFillB - fragWorldPos);
        vec3 lightDir3 = normalize(lightPosFillC - fragWorldPos);

        vec3 viewDir = normalize(viewPos - fragWorldPos);
        vec3 reflectDir0 = reflect(-lightDir0, norm);
        vec3 reflectDir1 = reflect(-lightDir1, norm);
        vec3 reflectDir2 = reflect(-lightDir2, norm);
        vec3 reflectDir3 = reflect(-lightDir3, norm);

        float diff0 = max(dot(norm, lightDir0), 0.0);
        float diff1 = max(dot(norm, lightDir1), 0.0);
        float diff2 = max(dot(norm, lightDir2), 0.0);
        float diff3 = max(dot(norm, lightDir3), 0.0);
        float diff = 0.45 * diff0 + 0.22 * diff1 + 0.18 * diff2 + 0.15 * diff3;
        diff = max(diff, 0.34);

        float specPower = max(fragShininess * 1.20, 16.0);
        float spec0 = pow(max(dot(viewDir, reflectDir0), 0.0), specPower);
        float spec1 = pow(max(dot(viewDir, reflectDir1), 0.0), specPower);
        float spec2 = pow(max(dot(viewDir, reflectDir2), 0.0), specPower);
        float spec3 = pow(max(dot(viewDir, reflectDir3), 0.0), specPower);
        float spec = 0.45 * spec0 + 0.22 * spec1 + 0.18 * spec2 + 0.15 * spec3;

        // Increase visual separation on light backgrounds.
        float ambient = 0.34;
        vec3 baseColor = fragColor;
        float luma = dot(baseColor, vec3(0.2126, 0.7152, 0.0722));
        baseColor = mix(vec3(luma), baseColor, 1.20);
        baseColor = clamp((baseColor - 0.5) * 1.10 + 0.5, 0.0, 1.0);

        float litFactor = ambient + (1.0 - ambient) * diff * (1.0 - shadow);
        vec3 diffuse = baseColor * litFactor;
        vec3 specular = vec3(0.46 * spec * (1.0 - shadow));

        color = vec4(diffuse + specular, 1.0);
    }
)";

static const char* kShadowVS = R"(
    #version 130

    in vec3 position;
    in vec3 instancePos;
    in float instanceScale;

    uniform mat4 lightMVP;

    void main()
    {
        vec3 worldPos    = position * instanceScale + instancePos;
        gl_Position      = lightMVP * vec4(worldPos, 1.0);
    }
)";

static const char* kShadowFS = R"(
    #version 130
    void main() {}
)";

static const char* kBondVS = R"(
    #version 130

    in vec3 position;
    in vec3 bondStart;
    in vec3 bondEnd;
    in vec3 bondColorA;
    in vec3 bondColorB;
    in float bondRadius;
    in float bondShininessA;
    in float bondShininessB;

    uniform mat4 projection;
    uniform mat4 view;

    out vec3 fragColorA;
    out vec3 fragColorB;
    out float fragAxis;
    out vec3 fragWorldPos;
    out vec3 fragNormal;
    out float fragShininessA;
    out float fragShininessB;

    void main()
    {
        vec3 axis = bondEnd - bondStart;
        float lengthAxis = length(axis);
        vec3 dir = (lengthAxis > 1e-5) ? axis / lengthAxis : vec3(0.0, 0.0, 1.0);
        vec3 helper = (abs(dir.y) < 0.999) ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
        vec3 tangent = normalize(cross(helper, dir));
        vec3 bitangent = cross(dir, tangent);

        vec3 local = vec3(position.xy * bondRadius, position.z * lengthAxis);
        vec3 center = 0.5 * (bondStart + bondEnd);
        vec3 worldPos = center + tangent * local.x + bitangent * local.y + dir * local.z;
        vec3 localNormal = normalize(vec3(position.xy, 0.0));
        vec3 worldNormal = normalize(tangent * localNormal.x + bitangent * localNormal.y + dir * localNormal.z);

        gl_Position = projection * view * vec4(worldPos, 1.0);
        fragColorA = bondColorA;
        fragColorB = bondColorB;
        fragAxis = position.z;
        fragWorldPos = worldPos;
        fragNormal = worldNormal;
        fragShininessA = bondShininessA;
        fragShininessB = bondShininessB;
    }
)";

static const char* kBondFS = R"(
    #version 130

    in vec3 fragColorA;
    in vec3 fragColorB;
    in float fragAxis;
    in vec3 fragWorldPos;
    in vec3 fragNormal;
    in float fragShininessA;
    in float fragShininessB;

    uniform vec3 lightPos;
    uniform vec3 viewPos;

    out vec4 color;

    void main()
    {
        vec3 norm = normalize(fragNormal);
        vec3 lightPosFillA = vec3(-lightPos.x, lightPos.y, lightPos.z);
        vec3 lightPosFillB = vec3(lightPos.x, -lightPos.y, lightPos.z);
        vec3 lightPosFillC = vec3(-lightPos.x, 0.65 * lightPos.y, -0.60 * lightPos.z);

        vec3 lightDir0 = normalize(lightPos - fragWorldPos);
        vec3 lightDir1 = normalize(lightPosFillA - fragWorldPos);
        vec3 lightDir2 = normalize(lightPosFillB - fragWorldPos);
        vec3 lightDir3 = normalize(lightPosFillC - fragWorldPos);

        vec3 viewDir = normalize(viewPos - fragWorldPos);
        vec3 reflectDir0 = reflect(-lightDir0, norm);
        vec3 reflectDir1 = reflect(-lightDir1, norm);
        vec3 reflectDir2 = reflect(-lightDir2, norm);
        vec3 reflectDir3 = reflect(-lightDir3, norm);

        float diff0 = max(dot(norm, lightDir0), 0.0);
        float diff1 = max(dot(norm, lightDir1), 0.0);
        float diff2 = max(dot(norm, lightDir2), 0.0);
        float diff3 = max(dot(norm, lightDir3), 0.0);
        float diff = 0.45 * diff0 + 0.22 * diff1 + 0.18 * diff2 + 0.15 * diff3;
        diff = max(diff, 0.34);

        vec3 baseColor = (fragAxis < 0.0) ? fragColorA : fragColorB;
        float luma = dot(baseColor, vec3(0.2126, 0.7152, 0.0722));
        baseColor = mix(vec3(luma), baseColor, 1.20);
        baseColor = clamp((baseColor - 0.5) * 1.10 + 0.5, 0.0, 1.0);

        float shininess = (fragAxis < 0.0) ? fragShininessA : fragShininessB;
        float specPower = max(shininess * 1.20, 16.0);
        float spec0 = pow(max(dot(viewDir, reflectDir0), 0.0), specPower);
        float spec1 = pow(max(dot(viewDir, reflectDir1), 0.0), specPower);
        float spec2 = pow(max(dot(viewDir, reflectDir2), 0.0), specPower);
        float spec3 = pow(max(dot(viewDir, reflectDir3), 0.0), specPower);
        float spec = 0.45 * spec0 + 0.22 * spec1 + 0.18 * spec2 + 0.15 * spec3;

        float ambient = 0.34;
        float litFactor = ambient + (1.0 - ambient) * diff;
        vec3 diffuse = baseColor * litFactor;
        vec3 specular = vec3(0.46 * spec);

        color = vec4(diffuse + specular, 1.0);
    }
)";

static const char* kBondShadowVS = R"(
    #version 130

    in vec3 position;
    in vec3 bondStart;
    in vec3 bondEnd;
    in float bondRadius;

    uniform mat4 lightMVP;

    void main()
    {
        vec3 axis = bondEnd - bondStart;
        float lengthAxis = length(axis);
        vec3 dir = (lengthAxis > 1e-5) ? axis / lengthAxis : vec3(0.0, 0.0, 1.0);
        vec3 helper = (abs(dir.y) < 0.999) ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
        vec3 tangent = normalize(cross(helper, dir));
        vec3 bitangent = cross(dir, tangent);

        vec3 local = vec3(position.xy * bondRadius, position.z * lengthAxis);
        vec3 center = 0.5 * (bondStart + bondEnd);
        vec3 worldPos = center + tangent * local.x + bitangent * local.y + dir * local.z;

        gl_Position = lightMVP * vec4(worldPos, 1.0);
    }
)";

// Low-poly instancing shaders (12 facet icosahedron)
static const char* kAtomLowPolyVS = R"(
    #version 130

    in vec3 position;
    in vec3 instancePos;
    in vec3 instanceColor;
    in float instanceScale;
    in float instanceShininess;

    uniform mat4 projection;
    uniform mat4 view;
    uniform mat4 lightMVP;
    uniform vec4 frustumPlanes[6];
    uniform int enableFrustumCulling;

    out vec3 fragColor;
    out vec3 fragWorldPos;
    out vec3 fragNormal;
    out float fragShininess;
    out vec4 FragPosLight;

    void main()
    {
        if (enableFrustumCulling != 0)
        {
            for (int i = 0; i < 6; ++i)
            {
                float signedDistance = dot(frustumPlanes[i].xyz, instancePos) + frustumPlanes[i].w;
                if (signedDistance < -instanceScale)
                {
                    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
                    FragPosLight = vec4(0.0);
                    fragColor = vec3(0.0);
                    fragWorldPos = vec3(0.0);
                    fragNormal = vec3(0.0, 0.0, 1.0);
                    fragShininess = instanceShininess;
                    return;
                }
            }
        }

        vec3 worldPos = position * instanceScale + instancePos;
        gl_Position   = projection * view * vec4(worldPos, 1.0);
        FragPosLight  = lightMVP * vec4(worldPos, 1.0);
        fragColor     = instanceColor;
        fragWorldPos  = worldPos;
        fragNormal    = normalize(position);
        fragShininess = instanceShininess;
    }
)";

static const char* kAtomLowPolyFS = kAtomFS;  // Reuse the same fragment shader

// Billboard imposters shader (quad rendering for very high atom counts)
static const char* kAtomBillboardVS = R"(
    #version 130

    in vec3 position;
    in vec3 instancePos;
    in vec3 instanceColor;
    in float instanceScale;
    in float instanceShininess;

    uniform mat4 projection;
    uniform mat4 view;
    uniform mat4 lightMVP;
    uniform vec4 frustumPlanes[6];
    uniform int enableFrustumCulling;
    uniform vec3 viewPos;

    out vec3 fragColor;
    out vec3 fragWorldCenter;
    out float fragRadius;
    out float fragShininess;
    out vec2 quad_uv;
    out vec4 FragPosLight;
    out mat4 fragView;
    out mat4 fragProj;

    void main()
    {
        if (enableFrustumCulling != 0)
        {
            for (int i = 0; i < 6; ++i)
            {
                float signedDistance = dot(frustumPlanes[i].xyz, instancePos) + frustumPlanes[i].w;
                if (signedDistance < -instanceScale)
                {
                    gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
                    fragColor = vec3(0.0);
                    fragWorldCenter = vec3(0.0);
                    fragRadius = 0.0;
                    fragShininess = instanceShininess;
                    quad_uv = vec2(0.0);
                    FragPosLight = vec4(0.0);
                    fragView = view;
                    fragProj = projection;
                    return;
                }
            }
        }

        // Expand quad towards camera
        vec3 viewDir = normalize(instancePos - viewPos);
        vec3 right = normalize(cross(vec3(0.0, 1.0, 0.0), viewDir));
        vec3 up = normalize(cross(viewDir, right));

        vec3 quadOffset = right * position.x + up * position.y;
        vec3 worldPos = instancePos + quadOffset * instanceScale;

        gl_Position = projection * view * vec4(worldPos, 1.0);
        FragPosLight = lightMVP * vec4(instancePos, 1.0);
        fragColor = instanceColor;
        fragWorldCenter = instancePos;
        fragRadius = instanceScale;
        fragShininess = instanceShininess;
        quad_uv = position.xy;
        fragView = view;
        fragProj = projection;
    }
)";

static const char* kAtomBillboardFS = R"(
    #version 130

    in vec3 fragColor;
    in vec3 fragWorldCenter;
    in float fragRadius;
    in float fragShininess;
    in vec2 quad_uv;
    in vec4 FragPosLight;
    in mat4 fragView;
    in mat4 fragProj;

    uniform sampler2D shadowMap;
    uniform vec3 lightPos;
    uniform vec3 viewPos;

    out vec4 color;

    float computeShadow(vec4 pos)
    {
        vec3 proj = pos.xyz / pos.w;
        proj = proj * 0.5 + 0.5;

        if (proj.x < 0.0 || proj.x > 1.0 ||
            proj.y < 0.0 || proj.y > 1.0 ||
            proj.z < 0.0 || proj.z > 1.0)
            return 0.0;

        float closest = texture(shadowMap, proj.xy).r;
        float current = proj.z;
        float bias    = 0.003;

        return (current - bias > closest) ? 1.0 : 0.0;
    }

    void main()
    {
        // Discard pixels outside the sphere silhouette.
        float r2 = dot(quad_uv, quad_uv);
        if (r2 > 1.0)
            discard;

        // Reconstruct spherical normal in view-aligned space and map to world.
        float nz = sqrt(1.0 - r2);
        vec3 viewDir = normalize(viewPos - fragWorldCenter);
        vec3 right = normalize(cross(vec3(0.0, 1.0, 0.0), viewDir));
        vec3 up = normalize(cross(viewDir, right));
        vec3 norm = normalize(right * quad_uv.x + up * quad_uv.y + viewDir * nz);

        // Approximate world position on the sphere surface.
        vec3 fragWorldPos = fragWorldCenter + norm * fragRadius;

        // Write correct depth so spheres intersect properly.
        vec4 clipPos = fragProj * fragView * vec4(fragWorldPos, 1.0);
        float ndcDepth = clipPos.z / clipPos.w;
        gl_FragDepth = ndcDepth * 0.5 + 0.5;

        float shadow = computeShadow(FragPosLight);

        // Multi-light Phong matching the standard atom shader.
        vec3 lightPosFillA = vec3(-lightPos.x, lightPos.y, lightPos.z);
        vec3 lightPosFillB = vec3(lightPos.x, -lightPos.y, lightPos.z);
        vec3 lightPosFillC = vec3(-lightPos.x, 0.65 * lightPos.y, -0.60 * lightPos.z);

        vec3 lightDir0 = normalize(lightPos - fragWorldPos);
        vec3 lightDir1 = normalize(lightPosFillA - fragWorldPos);
        vec3 lightDir2 = normalize(lightPosFillB - fragWorldPos);
        vec3 lightDir3 = normalize(lightPosFillC - fragWorldPos);

        vec3 vDir = normalize(viewPos - fragWorldPos);
        vec3 reflectDir0 = reflect(-lightDir0, norm);
        vec3 reflectDir1 = reflect(-lightDir1, norm);
        vec3 reflectDir2 = reflect(-lightDir2, norm);
        vec3 reflectDir3 = reflect(-lightDir3, norm);

        float diff0 = max(dot(norm, lightDir0), 0.0);
        float diff1 = max(dot(norm, lightDir1), 0.0);
        float diff2 = max(dot(norm, lightDir2), 0.0);
        float diff3 = max(dot(norm, lightDir3), 0.0);
        float diff = 0.45 * diff0 + 0.22 * diff1 + 0.18 * diff2 + 0.15 * diff3;
        diff = max(diff, 0.34);

        float specPower = max(fragShininess * 1.20, 16.0);
        float spec0 = pow(max(dot(vDir, reflectDir0), 0.0), specPower);
        float spec1 = pow(max(dot(vDir, reflectDir1), 0.0), specPower);
        float spec2 = pow(max(dot(vDir, reflectDir2), 0.0), specPower);
        float spec3 = pow(max(dot(vDir, reflectDir3), 0.0), specPower);
        float spec = 0.45 * spec0 + 0.22 * spec1 + 0.18 * spec2 + 0.15 * spec3;

        float ambient = 0.34;
        vec3 baseColor = fragColor;
        float luma = dot(baseColor, vec3(0.2126, 0.7152, 0.0722));
        baseColor = mix(vec3(luma), baseColor, 1.20);
        baseColor = clamp((baseColor - 0.5) * 1.10 + 0.5, 0.0, 1.0);

        float litFactor = ambient + (1.0 - ambient) * diff * (1.0 - shadow);
        vec3 diffuse = baseColor * litFactor;
        vec3 specular = vec3(0.46 * spec * (1.0 - shadow));

        color = vec4(diffuse + specular, 1.0);
    }
)";

// Low-poly shadow shader
static const char* kShadowLowPolyVS = R"(
    #version 130

    in vec3 position;
    in vec3 instancePos;
    in float instanceScale;

    uniform mat4 lightMVP;

    void main()
    {
        vec3 worldPos    = position * instanceScale + instancePos;
        gl_Position      = lightMVP * vec4(worldPos, 1.0);
    }
)";

// Billboard shadow shader
static const char* kShadowBillboardVS = R"(
    #version 130

    in vec3 position;
    in vec3 instancePos;
    in float instanceScale;

    uniform mat4 projection;
    uniform mat4 view;
    uniform mat4 lightMVP;
    uniform vec3 viewPos;

    out vec2 quad_uv;

    void main()
    {
        // Project quad towards light source
        vec3 viewDir = normalize(instancePos - viewPos);
        vec3 right = normalize(cross(vec3(0.0, 1.0, 0.0), viewDir));
        vec3 up = normalize(cross(viewDir, right));

        vec3 quadOffset = right * position.x + up * position.y;
        vec3 worldPos = instancePos + quadOffset * instanceScale;

        gl_Position = lightMVP * vec4(worldPos, 1.0);
        quad_uv = position.xy;
    }
)";

static const char* kShadowBillboardFS = R"(
    #version 130

    in vec2 quad_uv;

    void main()
    {
        if (dot(quad_uv, quad_uv) > 1.0)
            discard;
    }
)";

static const char* kLineVS = R"(
    #version 130

    in vec3 position;

    uniform mat4 projection;
    uniform mat4 view;

    void main()
    {
        gl_Position = projection * view * vec4(position, 1.0);
    }
)";

static const char* kLineFS = R"(
    #version 130

    uniform vec3 uColor;
    out vec4 color;

    void main()
    {
        color = vec4(uColor, 1.0);
    }
)";

// ---------------------------------------------------------------------------
// Renderer implementation
// ---------------------------------------------------------------------------

void Renderer::init()
{
    atomProgram   = createProgram(kAtomVS,   kAtomFS);
    atomLowPolyProgram = createProgram(kAtomLowPolyVS, kAtomLowPolyFS);
    atomBillboardProgram = createProgram(kAtomBillboardVS, kAtomBillboardFS);
    bondProgram   = createProgram(kBondVS,   kBondFS);
    shadowProgram = createProgram(kShadowVS, kShadowFS);
    shadowLowPolyProgram = createProgram(kShadowLowPolyVS, kShadowFS);
    shadowBillboardProgram = createProgram(kShadowBillboardVS, kShadowBillboardFS);
    bondShadowProgram = createProgram(kBondShadowVS, kShadowFS);
    lineProgram   = createProgram(kLineVS,   kLineFS);
}

RenderingMode Renderer::selectRenderingMode(size_t atomCount) const
{
    // Selection thresholds:
    // N > 10^7        -> Billboard imposters (quad rendering)
    // 10^5 < N < 10^7 -> Low-poly instancing (12 facet)
    // N <= 10^5       -> Standard instancing (full sphere)
    
    if (atomCount >= 10000000)  // 10 million
        return RenderingMode::BillboardImposters;
    else if (atomCount > 100000)  // 100k
        return RenderingMode::LowPolyInstancing;
    else
        return RenderingMode::StandardInstancing;
}

void Renderer::drawShadowPass(const ShadowMap& shadow,
                               const SphereMesh& sphere,
                               const glm::mat4& lightMVP,
                               size_t atomCount)
{
    beginShadowPass(shadow);

    glUseProgram(shadowProgram);
    glUniformMatrix4fv(glGetUniformLocation(shadowProgram, "lightMVP"),
                       1, GL_FALSE, glm::value_ptr(lightMVP));

    glBindVertexArray(sphere.vao);
    glDrawElementsInstanced(GL_TRIANGLES, sphere.indexCount, GL_UNSIGNED_INT, 0, (GLsizei)atomCount);

    endShadowPass();
}

void Renderer::drawBondShadowPass(const ShadowMap& shadow,
                                  const CylinderMesh& cylinder,
                                  const glm::mat4& lightMVP,
                                  size_t bondCount)
{
    if (bondCount == 0)
        return;

    beginShadowPass(shadow);

    glUseProgram(bondShadowProgram);
    glUniformMatrix4fv(glGetUniformLocation(bondShadowProgram, "lightMVP"),
                       1, GL_FALSE, glm::value_ptr(lightMVP));

    glBindVertexArray(cylinder.vao);
    glDrawArraysInstanced(GL_TRIANGLES, 0, cylinder.vertexCount, (GLsizei)bondCount);

    endShadowPass();
}

void Renderer::drawAtoms(const glm::mat4& projection,
                          const glm::mat4& view,
                          const glm::mat4& lightMVP,
                          const glm::vec3& lightPos,
                          const glm::vec3& viewPos,
                          const ShadowMap& shadow,
                          const SphereMesh& sphere,
                          size_t atomCount)
{
    glUseProgram(atomProgram);

    glUniformMatrix4fv(glGetUniformLocation(atomProgram, "projection"),
                       1, GL_FALSE, glm::value_ptr(projection));
    glUniformMatrix4fv(glGetUniformLocation(atomProgram, "view"),
                       1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(atomProgram, "lightMVP"),
                       1, GL_FALSE, glm::value_ptr(lightMVP));

    const std::array<glm::vec4, 6> frustumPlanes = extractFrustumPlanes(projection * view);
    glUniform4fv(glGetUniformLocation(atomProgram, "frustumPlanes"),
                 6,
                 glm::value_ptr(frustumPlanes[0]));
    glUniform1i(glGetUniformLocation(atomProgram, "enableFrustumCulling"), 1);

    glUniform3fv(glGetUniformLocation(atomProgram, "lightPos"),
                 1, glm::value_ptr(lightPos));
    glUniform3fv(glGetUniformLocation(atomProgram, "viewPos"),
                 1, glm::value_ptr(viewPos));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, shadow.depthTexture);
    glUniform1i(glGetUniformLocation(atomProgram, "shadowMap"), 0);

    glBindVertexArray(sphere.vao);
    glDrawElementsInstanced(GL_TRIANGLES, sphere.indexCount, GL_UNSIGNED_INT, 0, (GLsizei)atomCount);
}

void Renderer::drawShadowPassLowPoly(const ShadowMap& shadow,
                                     const LowPolyMesh& mesh,
                                     const glm::mat4& lightMVP,
                                     size_t atomCount)
{
    beginShadowPass(shadow);

    glUseProgram(shadowLowPolyProgram);
    glUniformMatrix4fv(glGetUniformLocation(shadowLowPolyProgram, "lightMVP"),
                       1, GL_FALSE, glm::value_ptr(lightMVP));

    glBindVertexArray(mesh.vao);
    glDrawElementsInstanced(GL_TRIANGLES, mesh.indexCount, GL_UNSIGNED_INT, 0, (GLsizei)atomCount);

    endShadowPass();
}

void Renderer::drawShadowPassBillboard(const ShadowMap& shadow,
                                       const BillboardMesh& mesh,
                                       const glm::mat4& lightMVP,
                                       const glm::mat4& view,
                                       size_t atomCount)
{
    beginShadowPass(shadow);

    glUseProgram(shadowBillboardProgram);
    glUniformMatrix4fv(glGetUniformLocation(shadowBillboardProgram, "lightMVP"),
                       1, GL_FALSE, glm::value_ptr(lightMVP));
    glUniformMatrix4fv(glGetUniformLocation(shadowBillboardProgram, "projection"),
                       1, GL_FALSE, glm::value_ptr(glm::ortho(-1.0f, 1.0f, -1.0f, 1.0f, 0.1f, 100.0f)));
    glUniformMatrix4fv(glGetUniformLocation(shadowBillboardProgram, "view"),
                       1, GL_FALSE, glm::value_ptr(view));

    glBindVertexArray(mesh.vao);
    glDrawElementsInstanced(GL_TRIANGLES, mesh.indexCount, GL_UNSIGNED_INT, 0, (GLsizei)atomCount);

    endShadowPass();
}

void Renderer::drawAtomsLowPoly(const glm::mat4& projection,
                                const glm::mat4& view,
                                const glm::mat4& lightMVP,
                                const glm::vec3& lightPos,
                                const glm::vec3& viewPos,
                                const ShadowMap& shadow,
                                const LowPolyMesh& mesh,
                                size_t atomCount)
{
    glUseProgram(atomLowPolyProgram);

    glUniformMatrix4fv(glGetUniformLocation(atomLowPolyProgram, "projection"),
                       1, GL_FALSE, glm::value_ptr(projection));
    glUniformMatrix4fv(glGetUniformLocation(atomLowPolyProgram, "view"),
                       1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(atomLowPolyProgram, "lightMVP"),
                       1, GL_FALSE, glm::value_ptr(lightMVP));

    const std::array<glm::vec4, 6> frustumPlanes = extractFrustumPlanes(projection * view);
    glUniform4fv(glGetUniformLocation(atomLowPolyProgram, "frustumPlanes"),
                 6,
                 glm::value_ptr(frustumPlanes[0]));
    glUniform1i(glGetUniformLocation(atomLowPolyProgram, "enableFrustumCulling"), 1);

    glUniform3fv(glGetUniformLocation(atomLowPolyProgram, "lightPos"),
                 1, glm::value_ptr(lightPos));
    glUniform3fv(glGetUniformLocation(atomLowPolyProgram, "viewPos"),
                 1, glm::value_ptr(viewPos));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, shadow.depthTexture);
    glUniform1i(glGetUniformLocation(atomLowPolyProgram, "shadowMap"), 0);

    glBindVertexArray(mesh.vao);
    glDrawElementsInstanced(GL_TRIANGLES, mesh.indexCount, GL_UNSIGNED_INT, 0, (GLsizei)atomCount);
}

void Renderer::drawAtomsBillboard(const glm::mat4& projection,
                                  const glm::mat4& view,
                                  const glm::mat4& lightMVP,
                                  const glm::vec3& lightPos,
                                  const glm::vec3& viewPos,
                                  const ShadowMap& shadow,
                                  const BillboardMesh& mesh,
                                  size_t atomCount)
{
    glUseProgram(atomBillboardProgram);

    glUniformMatrix4fv(glGetUniformLocation(atomBillboardProgram, "projection"),
                       1, GL_FALSE, glm::value_ptr(projection));
    glUniformMatrix4fv(glGetUniformLocation(atomBillboardProgram, "view"),
                       1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(atomBillboardProgram, "lightMVP"),
                       1, GL_FALSE, glm::value_ptr(lightMVP));

    const std::array<glm::vec4, 6> frustumPlanes = extractFrustumPlanes(projection * view);
    glUniform4fv(glGetUniformLocation(atomBillboardProgram, "frustumPlanes"),
                 6,
                 glm::value_ptr(frustumPlanes[0]));
    glUniform1i(glGetUniformLocation(atomBillboardProgram, "enableFrustumCulling"), 1);

    glUniform3fv(glGetUniformLocation(atomBillboardProgram, "lightPos"),
                 1, glm::value_ptr(lightPos));
    glUniform3fv(glGetUniformLocation(atomBillboardProgram, "viewPos"),
                 1, glm::value_ptr(viewPos));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, shadow.depthTexture);
    glUniform1i(glGetUniformLocation(atomBillboardProgram, "shadowMap"), 0);

    glBindVertexArray(mesh.vao);
    glDrawElementsInstanced(GL_TRIANGLES, mesh.indexCount, GL_UNSIGNED_INT, 0, (GLsizei)atomCount);
}

void Renderer::drawBonds(const glm::mat4& projection,
                         const glm::mat4& view,
                         const glm::vec3& lightPos,
                         const glm::vec3& viewPos,
                         const CylinderMesh& cylinder,
                         size_t bondCount)
{
    if (bondCount == 0)
        return;

    glUseProgram(bondProgram);

    glUniformMatrix4fv(glGetUniformLocation(bondProgram, "projection"),
                       1, GL_FALSE, glm::value_ptr(projection));
    glUniformMatrix4fv(glGetUniformLocation(bondProgram, "view"),
                       1, GL_FALSE, glm::value_ptr(view));
    glUniform3fv(glGetUniformLocation(bondProgram, "lightPos"),
                 1, glm::value_ptr(lightPos));
    glUniform3fv(glGetUniformLocation(bondProgram, "viewPos"),
                 1, glm::value_ptr(viewPos));

    glBindVertexArray(cylinder.vao);
    glDrawArraysInstanced(GL_TRIANGLES, 0, cylinder.vertexCount, (GLsizei)bondCount);
}

void Renderer::drawBoxLines(const glm::mat4& projection,
                             const glm::mat4& view,
                             GLuint lineVAO,
                             size_t lineVertexCount)
{
    if (lineVertexCount == 0)
        return;

    glUseProgram(lineProgram);

    glUniformMatrix4fv(glGetUniformLocation(lineProgram, "projection"),
                       1, GL_FALSE, glm::value_ptr(projection));
    glUniformMatrix4fv(glGetUniformLocation(lineProgram, "view"),
                       1, GL_FALSE, glm::value_ptr(view));
    glUniform3f(glGetUniformLocation(lineProgram, "uColor"), 0.85f, 0.85f, 0.85f);

    glLineWidth(2.0f);
    glBindVertexArray(lineVAO);
    glDrawArrays(GL_LINES, 0, (GLsizei)lineVertexCount);
}
