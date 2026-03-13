#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "Camera.h"
#include "SphereMesh.h"
#include "Shader.h"
#include "Controls.h"
#include "ShadowMap.h"
#include "StructureLoader.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <vector>
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

    // ------------------------------------------------
    // Instance buffers
    // ------------------------------------------------

    std::vector<glm::vec3> positions;
    std::vector<glm::vec3> colors;

    for(const auto& atom : structure.atoms)
    {
        positions.emplace_back(atom.x,atom.y,atom.z);
        colors.emplace_back(atom.r,atom.g,atom.b);
    }

    GLuint instanceVBO,colorVBO;

    glBindVertexArray(sphere.vao);

    glGenBuffers(1,&instanceVBO);
    glBindBuffer(GL_ARRAY_BUFFER,instanceVBO);

    glBufferData(GL_ARRAY_BUFFER,
                 positions.size()*sizeof(glm::vec3),
                 positions.data(),
                 GL_STATIC_DRAW);

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(glm::vec3),(void*)0);
    glVertexAttribDivisor(1,1);

    glGenBuffers(1,&colorVBO);
    glBindBuffer(GL_ARRAY_BUFFER,colorVBO);

    glBufferData(GL_ARRAY_BUFFER,
                 colors.size()*sizeof(glm::vec3),
                 colors.data(),
                 GL_STATIC_DRAW);

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2,3,GL_FLOAT,GL_FALSE,sizeof(glm::vec3),(void*)0);
    glVertexAttribDivisor(2,1);

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

        vec3 lighting = (1.0-shadow) * fragColor;

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

    uniform mat4 lightMVP;

    void main()
    {
        gl_Position = lightMVP * vec4(position,1.0);
    }

    )";

    const char* shadowFS = R"(

    #version 130
    void main(){}

    )";

    GLuint shadowProgram = createProgram(shadowVS,shadowFS);

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

        drawControls(camera);

        ImGui::Render();

        int w,h;
        glfwGetFramebufferSize(window,&w,&h);

        // --------------------------------------------
        // Light matrices
        // --------------------------------------------

        glm::mat4 lightProj =
            glm::perspective(glm::radians(45.0f),1.0f,0.1f,100.0f);

        glm::mat4 lightView =
            glm::lookAt(glm::vec3(10,10,10),
                        glm::vec3(0,0,0),
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
            structure.atoms.size()
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

        glm::vec3 camPos(
            camera.distance*cos(pitch)*sin(yaw),
            camera.distance*sin(pitch),
            camera.distance*cos(pitch)*cos(yaw)
        );

        glm::mat4 view =
            glm::lookAt(camPos,
                        glm::vec3(0,0,0),
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
            structure.atoms.size()
        );

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    glfwTerminate();
}