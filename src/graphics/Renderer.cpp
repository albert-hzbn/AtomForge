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

        vec3 lightDir0 = normalize(lightPos - fragWorldPos);
        vec3 lightDir1 = normalize(lightPosFillA - fragWorldPos);
        vec3 lightDir2 = normalize(lightPosFillB - fragWorldPos);

        vec3 viewDir = normalize(viewPos - fragWorldPos);
        vec3 reflectDir0 = reflect(-lightDir0, norm);
        vec3 reflectDir1 = reflect(-lightDir1, norm);

        float diff0 = max(dot(norm, lightDir0), 0.0);
        float diff1 = max(dot(norm, lightDir1), 0.0);
        float diff2 = max(dot(norm, lightDir2), 0.0);
        float diff = 0.55 * diff0 + 0.30 * diff1 + 0.15 * diff2;
        diff = max(diff, 0.30);

        float spec0 = pow(max(dot(viewDir, reflectDir0), 0.0), max(fragShininess, 1.0));
        float spec1 = pow(max(dot(viewDir, reflectDir1), 0.0), max(fragShininess, 1.0));
        float spec = 0.75 * spec0 + 0.25 * spec1;

        float ambient = 0.45;
        float litFactor = ambient + (1.0 - ambient) * diff * (1.0 - shadow);
        vec3 diffuse = fragColor * litFactor;
        vec3 specular = vec3(0.40 * spec * (1.0 - shadow));

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

    uniform mat4 projection;
    uniform mat4 view;

    out vec3 fragColorA;
    out vec3 fragColorB;
    out float fragAxis;

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

        gl_Position = projection * view * vec4(worldPos, 1.0);
        fragColorA = bondColorA;
        fragColorB = bondColorB;
        fragAxis = position.z;
    }
)";

static const char* kBondFS = R"(
    #version 130

    in vec3 fragColorA;
    in vec3 fragColorB;
    in float fragAxis;

    out vec4 color;

    void main()
    {
        vec3 baseColor = (fragAxis < 0.0) ? fragColorA : fragColorB;
        color = vec4(baseColor, 1.0);
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
    bondProgram   = createProgram(kBondVS,   kBondFS);
    shadowProgram = createProgram(kShadowVS, kShadowFS);
    bondShadowProgram = createProgram(kBondShadowVS, kShadowFS);
    lineProgram   = createProgram(kLineVS,   kLineFS);
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
    glDrawArraysInstanced(GL_TRIANGLES, 0, sphere.vertexCount, (GLsizei)atomCount);

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
    glDrawArraysInstanced(GL_TRIANGLES, 0, sphere.vertexCount, (GLsizei)atomCount);
}

void Renderer::drawBonds(const glm::mat4& projection,
                         const glm::mat4& view,
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
