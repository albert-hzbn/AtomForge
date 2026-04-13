#include "graphics/CustomStructureShaders.h"

const char* customStructureMeshVertexShader()
{
    return R"(
    #version 130
    in vec3 position;
    in vec3 normal;
    uniform mat4 projection;
    uniform mat4 view;
    out vec3 vNormal;
    out vec3 vFragPos;
    void main()
    {
        vFragPos = position;
        vNormal  = normal;
        gl_Position = projection * view * vec4(position, 1.0);
    }
)";
}

const char* customStructureMeshFragmentShader()
{
    return R"(
    #version 130
    in vec3 vNormal;
    in vec3 vFragPos;
    uniform vec3 uLightDir;
    uniform vec3 uViewPos;
    uniform vec3 uColor;
    out vec4 fragColor;
    void main()
    {
        vec3 n = normalize(vNormal);
        // Two-sided lighting
        vec3 lightDir = normalize(uLightDir);
        float diff = max(dot(n, lightDir), 0.0);
        float diffBack = max(dot(-n, lightDir), 0.0);
        diff = max(diff, diffBack * 0.7);

        vec3 viewDir = normalize(uViewPos - vFragPos);
        vec3 halfDir = normalize(lightDir + viewDir);
        float spec = pow(max(dot(n, halfDir), 0.0), 32.0) * 0.3;
        float specBack = pow(max(dot(-n, halfDir), 0.0), 32.0) * 0.15;
        spec = max(spec, specBack);

        float ambient = 0.15;
        vec3 color = uColor * (ambient + diff * 0.75) + vec3(spec);
        fragColor = vec4(color, 1.0);
    }
)";
}
