#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "Camera.h"
#include "SphereMesh.h"
#include "Shader.h"
#include "ShadowMap.h"
#include "StructureLoader.h"
#include "ui/FileBrowser.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <vector>
#include <cstring>
#include <iostream>

int main()
{
    if(!glfwInit())
        return -1;

    GLFWwindow* window =
        glfwCreateWindow(1280,800,"Atoms Editor",nullptr,nullptr);

    if(!window)
        return -1;

    glfwMakeContextCurrent(window);

    glewExperimental = GL_TRUE;

    if(glewInit()!=GLEW_OK)
        return -1;

    glEnable(GL_DEPTH_TEST);

    // ------------------------------------------------
    // Camera
    // ------------------------------------------------

    Camera camera;
    Camera::instance=&camera;

    glfwSetMouseButtonCallback(window,Camera::mouseButton);
    glfwSetCursorPosCallback(window,Camera::cursor);
    glfwSetScrollCallback(window,Camera::scroll);

    // ------------------------------------------------
    // ImGui
    // ------------------------------------------------

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();

    // Improve default appearance and readability
    ImGui::StyleColorsDark();
    ImGuiIO& io = ImGui::GetIO();

    // Use ImGui's default font but scale it up for better readability
    io.FontGlobalScale = 1.0f;

    // Better UX when dragging windows (only from title bar)
    io.ConfigWindowsMoveFromTitleBarOnly = true;

    // Slightly round corners and increase padding for a cleaner look
    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 6.0f;
    style.FrameRounding = 4.0f;
    style.WindowPadding = ImVec2(10, 10);
    style.FramePadding = ImVec2(6, 4);
    style.ItemSpacing = ImVec2(8, 6);
    style.ScrollbarSize = 18.0f;
    style.GrabMinSize = 10.0f;

    // Subtle color tweaks
    style.Colors[ImGuiCol_Header] = ImVec4(0.12f, 0.55f, 0.76f, 0.52f);
    style.Colors[ImGuiCol_HeaderHovered] = ImVec4(0.26f, 0.59f, 0.98f, 0.8f);
    style.Colors[ImGuiCol_HeaderActive] = ImVec4(0.26f, 0.59f, 0.98f, 1.0f);
    style.Colors[ImGuiCol_Button] = ImVec4(0.20f, 0.50f, 0.80f, 0.40f);
    style.Colors[ImGuiCol_ButtonHovered] = ImVec4(0.24f, 0.60f, 0.90f, 0.60f);
    style.Colors[ImGuiCol_ButtonActive] = ImVec4(0.06f, 0.53f, 0.83f, 0.80f);

    ImGui_ImplGlfw_InitForOpenGL(window,true);
    ImGui_ImplOpenGL3_Init("#version 130");

    // ------------------------------------------------
    // Geometry
    // ------------------------------------------------

    SphereMesh sphere(40,40);

    // ------------------------------------------------
    // Load structure
    // ------------------------------------------------

    std::string filename="../data/Cu6Sn5.cif";

    Structure structure = loadStructure(filename);

    FileBrowser fileBrowser;
    fileBrowser.initFromPath(filename);

    // ------------------------------------------------
    // Instance buffers
    // ------------------------------------------------

    std::vector<glm::vec3> positions;
    std::vector<glm::vec3> colors;
    std::vector<glm::vec3> boxLines;
    glm::vec3 orbitCenter(0.0f, 0.0f, 0.0f);
    size_t atomCount = 0;

    GLuint instanceVBO, colorVBO;
    GLuint lineVAO, lineVBO;

    glBindVertexArray(sphere.vao);

    glGenBuffers(1, &instanceVBO);
    glGenBuffers(1, &colorVBO);
    glGenVertexArrays(1, &lineVAO);
    glGenBuffers(1, &lineVBO);

    // Setup a simple VAO for line rendering (unit cell bounding box)
    glBindVertexArray(lineVAO);
    glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
    glBindVertexArray(0);


    auto updateBuffers = [&](const Structure& s) {
        positions.clear();
        colors.clear();
        boxLines.clear();

        std::vector<glm::vec3> basePositions;
        std::vector<glm::vec3> baseColors;
        basePositions.reserve(s.atoms.size());
        baseColors.reserve(s.atoms.size());

        for (const auto& atom : s.atoms)
        {
            basePositions.emplace_back((float)atom.x, (float)atom.y, (float)atom.z);
            baseColors.emplace_back(atom.r, atom.g, atom.b);
        }

        bool supercell = fileBrowser.isSupercellEnabled() && s.hasUnitCell;
        bool transform = fileBrowser.isTransformMatrixEnabled() && s.hasUnitCell;

        if (transform)
        {
            constexpr float eps = 1e-5f;

            glm::vec3 origin((float)s.cellOffset[0], (float)s.cellOffset[1], (float)s.cellOffset[2]);
            glm::vec3 a((float)s.cellVectors[0][0], (float)s.cellVectors[0][1], (float)s.cellVectors[0][2]);
            glm::vec3 b((float)s.cellVectors[1][0], (float)s.cellVectors[1][1], (float)s.cellVectors[1][2]);
            glm::vec3 c((float)s.cellVectors[2][0], (float)s.cellVectors[2][1], (float)s.cellVectors[2][2]);

            const int (&matrix)[3][3] = fileBrowser.getTransformMatrix();

            glm::mat3 cellMat(a, b, c);
            glm::mat3 invCellMat = glm::inverse(cellMat);

            // GLM matrices are column-major; fill M so rows match matrix[i][j].
            glm::mat3 M(1.0f);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    M[j][i] = (float)matrix[i][j];

            float detM = glm::determinant(M);
            if (std::abs(detM) <= eps)
            {
                // Degenerate matrix: keep original atoms instead of producing invalid output.
                positions = basePositions;
                colors = baseColors;
            }
            else
            {
                glm::mat3 invM = glm::inverse(M);

                int expectedCopies = std::max(1, (int)std::round(std::abs(detM)));
                positions.reserve(basePositions.size() * (size_t)expectedCopies);
                colors.reserve(baseColors.size() * (size_t)expectedCopies);

                int nMin[3] = {0, 0, 0};
                int nMax[3] = {0, 0, 0};
                for (int i = 0; i < 3; ++i)
                {
                    int rowMin = 0;
                    int rowMax = 0;
                    for (int j = 0; j < 3; ++j)
                    {
                        if (matrix[i][j] < 0) rowMin += matrix[i][j];
                        if (matrix[i][j] > 0) rowMax += matrix[i][j];
                    }
                    // Pad by one to be robust around boundaries and negative entries.
                    nMin[i] = rowMin - 1;
                    nMax[i] = rowMax + 1;
                }

                for (size_t ai = 0; ai < basePositions.size(); ++ai)
                {
                    glm::vec3 cart = basePositions[ai];
                    glm::vec3 frac = invCellMat * (cart - origin);

                    // Fold into [0,1) so boundary filtering yields one unique image.
                    frac.x -= std::floor(frac.x);
                    frac.y -= std::floor(frac.y);
                    frac.z -= std::floor(frac.z);

                    for (int ix = nMin[0]; ix <= nMax[0]; ++ix)
                    {
                        for (int iy = nMin[1]; iy <= nMax[1]; ++iy)
                        {
                            for (int iz = nMin[2]; iz <= nMax[2]; ++iz)
                            {
                                glm::vec3 shiftedFrac = frac + glm::vec3((float)ix, (float)iy, (float)iz);
                                glm::vec3 newFrac = invM * shiftedFrac;

                                if (newFrac.x < -eps || newFrac.y < -eps || newFrac.z < -eps)
                                    continue;
                                if (newFrac.x >= 1.0f - eps || newFrac.y >= 1.0f - eps || newFrac.z >= 1.0f - eps)
                                    continue;

                                glm::vec3 newCart = origin + shiftedFrac.x * a + shiftedFrac.y * b + shiftedFrac.z * c;
                                positions.push_back(newCart);
                                colors.push_back(baseColors[ai]);
                            }
                        }
                    }
                }
            }
        }
        else if (supercell)
        {
            glm::vec3 a((float)s.cellVectors[0][0], (float)s.cellVectors[0][1], (float)s.cellVectors[0][2]);
            glm::vec3 b((float)s.cellVectors[1][0], (float)s.cellVectors[1][1], (float)s.cellVectors[1][2]);
            glm::vec3 c((float)s.cellVectors[2][0], (float)s.cellVectors[2][1], (float)s.cellVectors[2][2]);

            for (int ix = -1; ix <= 1; ++ix)
            {
                for (int iy = -1; iy <= 1; ++iy)
                {
                    for (int iz = -1; iz <= 1; ++iz)
                    {
                        glm::vec3 shift = (float)ix * a + (float)iy * b + (float)iz * c;
                        for (size_t ai = 0; ai < basePositions.size(); ++ai)
                        {
                            positions.push_back(basePositions[ai] + shift);
                            colors.push_back(baseColors[ai]);
                        }
                    }
                }
            }
        }
        else
        {
            positions = std::move(basePositions);
            colors = std::move(baseColors);
        }

        // Compute bounds from the final position set
        glm::vec3 minPos(1e9f), maxPos(-1e9f);
        for (const auto& p : positions)
        {
            minPos = glm::min(minPos, p);
            maxPos = glm::max(maxPos, p);
        }

        atomCount = positions.size();

        // Build bounding box edges from min/max or unit cell lattice vectors
        if (atomCount > 0)
        {
            glm::vec3 corners[8];

            if (s.hasUnitCell)
            {
                glm::vec3 origin(
                    (float)s.cellOffset[0],
                    (float)s.cellOffset[1],
                    (float)s.cellOffset[2]);

                glm::vec3 a(
                    (float)s.cellVectors[0][0],
                    (float)s.cellVectors[0][1],
                    (float)s.cellVectors[0][2]);

                glm::vec3 b(
                    (float)s.cellVectors[1][0],
                    (float)s.cellVectors[1][1],
                    (float)s.cellVectors[1][2]);

                glm::vec3 c(
                    (float)s.cellVectors[2][0],
                    (float)s.cellVectors[2][1],
                    (float)s.cellVectors[2][2]);

                if (fileBrowser.isTransformMatrixEnabled())
                {
                    // Compute transformed lattice vectors
                    const int (&matrix)[3][3] = fileBrowser.getTransformMatrix();
                    glm::vec3 aT = (float)matrix[0][0]*a + (float)matrix[0][1]*b + (float)matrix[0][2]*c;
                    glm::vec3 bT = (float)matrix[1][0]*a + (float)matrix[1][1]*b + (float)matrix[1][2]*c;
                    glm::vec3 cT = (float)matrix[2][0]*a + (float)matrix[2][1]*b + (float)matrix[2][2]*c;
                    // Build box corners from transformed basis
                    corners[0] = origin;
                    corners[1] = origin + aT;
                    corners[2] = origin + aT + bT;
                    corners[3] = origin + bT;
                    corners[4] = origin + cT;
                    corners[5] = origin + aT + cT;
                    corners[6] = origin + aT + bT + cT;
                    corners[7] = origin + bT + cT;
                }
                else if (fileBrowser.isSupercellEnabled())
                {
                    glm::vec3 origin3 = origin - a - b - c;
                    glm::vec3 a3 = a * 3.0f;
                    glm::vec3 b3 = b * 3.0f;
                    glm::vec3 c3 = c * 3.0f;

                    corners[0] = origin3;
                    corners[1] = origin3 + a3;
                    corners[2] = origin3 + a3 + b3;
                    corners[3] = origin3 + b3;
                    corners[4] = origin3 + c3;
                    corners[5] = origin3 + a3 + c3;
                    corners[6] = origin3 + a3 + b3 + c3;
                    corners[7] = origin3 + b3 + c3;
                }
                else
                {
                    corners[0] = origin;
                    corners[1] = origin + a;
                    corners[2] = origin + a + b;
                    corners[3] = origin + b;
                    corners[4] = origin + c;
                    corners[5] = origin + a + c;
                    corners[6] = origin + a + b + c;
                    corners[7] = origin + b + c;
                }
            }
            else
            {
                corners[0] = {minPos.x, minPos.y, minPos.z};
                corners[1] = {maxPos.x, minPos.y, minPos.z};
                corners[2] = {maxPos.x, maxPos.y, minPos.z};
                corners[3] = {minPos.x, maxPos.y, minPos.z};
                corners[4] = {minPos.x, minPos.y, maxPos.z};
                corners[5] = {maxPos.x, minPos.y, maxPos.z};
                corners[6] = {maxPos.x, maxPos.y, maxPos.z};
                corners[7] = {minPos.x, maxPos.y, maxPos.z};
            }

            // Orbit around the geometric center of the displayed cell/box.
            orbitCenter = 0.5f * (corners[0] + corners[6]);

            const int edges[12][2] = {
                {0,1},{1,2},{2,3},{3,0},
                {4,5},{5,6},{6,7},{7,4},
                {0,4},{1,5},{2,6},{3,7}
            };

            for (auto& e : edges)
            {
                boxLines.push_back(corners[e[0]]);
                boxLines.push_back(corners[e[1]]);
            }

            glBindBuffer(GL_ARRAY_BUFFER, lineVBO);
            glBufferData(GL_ARRAY_BUFFER,
                         boxLines.size() * sizeof(glm::vec3),
                         boxLines.data(),
                         GL_STATIC_DRAW);
        }
        else
        {
            orbitCenter = glm::vec3(0.0f, 0.0f, 0.0f);
        }

        glBindBuffer(GL_ARRAY_BUFFER, instanceVBO);
        glBufferData(GL_ARRAY_BUFFER,
                     positions.size() * sizeof(glm::vec3),
                     positions.data(),
                     GL_STATIC_DRAW);

        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
        glVertexAttribDivisor(1, 1);

        glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
        glBufferData(GL_ARRAY_BUFFER,
                     colors.size() * sizeof(glm::vec3),
                     colors.data(),
                     GL_STATIC_DRAW);

        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
        glVertexAttribDivisor(2, 1);
    };

    updateBuffers(structure);

    // ------------------------------------------------
    // Shader
    // ------------------------------------------------

    const char* vs = R"(

    #version 130

    in vec3 position;
    in vec3 instancePos;
    in vec3 instanceColor;

    uniform mat4 projection;
    uniform mat4 view;
    uniform mat4 lightMVP;

    out vec3 fragColor;
    out vec4 FragPosLight;

    void main()
    {
        vec3 worldPos = position + instancePos;

        gl_Position = projection * view * vec4(worldPos,1.0);

        FragPosLight = lightMVP * vec4(worldPos,1.0);

        fragColor = instanceColor;
    }

    )";

    const char* fs = R"(

    #version 130

    in vec3 fragColor;
    in vec4 FragPosLight;

    uniform sampler2D shadowMap;

    out vec4 color;

    float computeShadow(vec4 pos)
    {
        vec3 proj = pos.xyz / pos.w;
        proj = proj * 0.5 + 0.5;

        // Outside the shadow map frustum should not be treated as shadowed.
        if (proj.x < 0.0 || proj.x > 1.0 ||
            proj.y < 0.0 || proj.y > 1.0 ||
            proj.z < 0.0 || proj.z > 1.0)
            return 0.0;

        float closest = texture(shadowMap,proj.xy).r;
        float current = proj.z;

        float bias = 0.003;

        if(current - bias > closest)
            return 1.0;

        return 0.0;
    }

    void main()
    {
        float shadow = computeShadow(FragPosLight);

        float ambient = 0.25;
        vec3 lighting = (ambient + (1.0 - ambient) * (1.0 - shadow)) * fragColor;

        color = vec4(lighting,1.0);
    }

    )";

    GLuint program = createProgram(vs,fs);

    // ------------------------------------------------
    // Shadow shader
    // ------------------------------------------------

    const char* shadowVS = R"(

    #version 130
    in vec3 position;
    in vec3 instancePos;

    uniform mat4 lightMVP;

    void main()
    {
        vec3 worldPos = position + instancePos;
        gl_Position = lightMVP * vec4(worldPos,1.0);
    }

    )";

    const char* shadowFS = R"(

    #version 130
    void main(){}

    )";

    GLuint shadowProgram = createProgram(shadowVS,shadowFS);

    // ------------------------------------------------
    // Line shader (for bounding box / lattice vectors)
    // ------------------------------------------------
    const char* lineVS = R"(

    #version 130

    in vec3 position;

    uniform mat4 projection;
    uniform mat4 view;

    void main()
    {
        gl_Position = projection * view * vec4(position, 1.0);
    }

    )";

    const char* lineFS = R"(

    #version 130

    uniform vec3 uColor;
    out vec4 color;

    void main()
    {
        color = vec4(uColor, 1.0);
    }

    )";

    GLuint lineProgram = createProgram(lineVS, lineFS);

    // ------------------------------------------------
    // Shadow map
    // ------------------------------------------------

    ShadowMap shadow = createShadowMap(1024,1024);

    // ------------------------------------------------
    // Render loop
    // ------------------------------------------------

    while(!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        fileBrowser.draw(structure, updateBuffers);

        ImGui::Render();

        int w,h;
        glfwGetFramebufferSize(window,&w,&h);

        // --------------------------------------------
        // Light matrices
        // --------------------------------------------

        glm::mat4 lightProj =
            glm::perspective(glm::radians(45.0f),1.0f,0.1f,1000.0f);

        glm::vec3 lightPos = orbitCenter + glm::vec3(40.0f, 40.0f, 40.0f);
        glm::mat4 lightView =
            glm::lookAt(lightPos,
                        orbitCenter,
                        glm::vec3(0,1,0));

        glm::mat4 lightMVP = lightProj * lightView;

        // --------------------------------------------
        // Shadow pass
        // --------------------------------------------

        beginShadowPass(shadow);

        glUseProgram(shadowProgram);

        GLuint lightLoc =
            glGetUniformLocation(shadowProgram,"lightMVP");

        glUniformMatrix4fv(lightLoc,1,GL_FALSE,
                           glm::value_ptr(lightMVP));

        glBindVertexArray(sphere.vao);

        glDrawArraysInstanced(
            GL_TRIANGLES,
            0,
            sphere.vertexCount,
            atomCount
        );

        endShadowPass();

        // --------------------------------------------
        // Normal pass
        // --------------------------------------------

        glViewport(0,0,w,h);

        glClearColor(0.2f,0.2f,0.2f,1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glm::mat4 projection =
            glm::perspective(glm::radians(45.0f),
                             (float)w/h,
                             0.1f,
                             1000.0f);

        float yaw   = glm::radians(camera.yaw);
        float pitch = glm::radians(camera.pitch);

        glm::vec3 camOffset(
            camera.distance*cos(pitch)*sin(yaw),
            camera.distance*sin(pitch),
            camera.distance*cos(pitch)*cos(yaw)
        );

        glm::vec3 camPos = orbitCenter + camOffset;

        glm::mat4 view =
            glm::lookAt(camPos,
                        orbitCenter,
                        glm::vec3(0,1,0));

        glUseProgram(program);

        GLuint projLoc = glGetUniformLocation(program,"projection");
        GLuint viewLoc = glGetUniformLocation(program,"view");
        GLuint lightLoc2 = glGetUniformLocation(program,"lightMVP");

        glUniformMatrix4fv(projLoc,1,GL_FALSE,
                           glm::value_ptr(projection));

        glUniformMatrix4fv(viewLoc,1,GL_FALSE,
                           glm::value_ptr(view));

        glUniformMatrix4fv(lightLoc2,1,GL_FALSE,
                           glm::value_ptr(lightMVP));

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D,shadow.depthTexture);

        GLuint shadowLoc =
            glGetUniformLocation(program,"shadowMap");

        glUniform1i(shadowLoc,0);

        glBindVertexArray(sphere.vao);

        glDrawArraysInstanced(
            GL_TRIANGLES,
            0,
            sphere.vertexCount,
            atomCount
        );

        // Render bounding box / lattice lines
        if (!boxLines.empty()) {
            glUseProgram(lineProgram);

            GLuint projLoc2 = glGetUniformLocation(lineProgram, "projection");
            GLuint viewLoc2 = glGetUniformLocation(lineProgram, "view");
            GLuint colorLoc = glGetUniformLocation(lineProgram, "uColor");

            glUniformMatrix4fv(projLoc2, 1, GL_FALSE, glm::value_ptr(projection));
            glUniformMatrix4fv(viewLoc2, 1, GL_FALSE, glm::value_ptr(view));
            glUniform3f(colorLoc, 0.85f, 0.85f, 0.85f);

            glLineWidth(2.0f);
            glBindVertexArray(lineVAO);
            glDrawArrays(GL_LINES, 0, (GLsizei)boxLines.size());
        }

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    glfwTerminate();
}