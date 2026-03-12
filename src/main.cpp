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

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <iostream>

int main()
{
    // ------------------------------------------------------------
    // GLFW INITIALIZATION
    // ------------------------------------------------------------

    if(!glfwInit())
    {
        std::cerr<<"GLFW init failed\n";
        return -1;
    }

    GLFWwindow* window =
        glfwCreateWindow(1280,800,"Shadow Mapping Demo",nullptr,nullptr);

    if(!window)
    {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);

    // ------------------------------------------------------------
    // GLEW INITIALIZATION
    // ------------------------------------------------------------

    glewExperimental = GL_TRUE;

    if(glewInit()!=GLEW_OK)
    {
        std::cerr<<"GLEW init failed\n";
        return -1;
    }

    std::cout<<"OpenGL "<<glGetString(GL_VERSION)<<std::endl;

    glEnable(GL_DEPTH_TEST);

    // ------------------------------------------------------------
    // CAMERA
    // ------------------------------------------------------------

    Camera camera;
    Camera::instance=&camera;

    glfwSetMouseButtonCallback(window,Camera::mouseButton);
    glfwSetCursorPosCallback(window,Camera::cursor);
    glfwSetScrollCallback(window,Camera::scroll);

    // ------------------------------------------------------------
    // IMGUI
    // ------------------------------------------------------------

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();

    ImGui_ImplGlfw_InitForOpenGL(window,true);
    ImGui_ImplOpenGL3_Init("#version 130");

    // ------------------------------------------------------------
    // GEOMETRY
    // ------------------------------------------------------------

    SphereMesh sphere(50,50); // sphere mesh

    // simple ground plane
    float planeVertices[] =
    {
        -10,-1,-10,
         10,-1,-10,
         10,-1, 10,
        -10,-1,-10,
         10,-1, 10,
        -10,-1, 10
    };

    GLuint planeVAO,planeVBO;

    glGenVertexArrays(1,&planeVAO);
    glGenBuffers(1,&planeVBO);

    glBindVertexArray(planeVAO);

    glBindBuffer(GL_ARRAY_BUFFER,planeVBO);
    glBufferData(GL_ARRAY_BUFFER,sizeof(planeVertices),
                 planeVertices,GL_STATIC_DRAW);

    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,(void*)0);
    glEnableVertexAttribArray(0);

    // ------------------------------------------------------------
    // MAIN SHADER (with shadow support)
    // ------------------------------------------------------------

    const char* vs = R"(

    #version 130

    in vec3 position;

    uniform mat4 MVP;
    uniform mat4 lightMVP;

    out vec4 FragPosLight;

    void main()
    {
        gl_Position = MVP * vec4(position,1.0);
        FragPosLight = lightMVP * vec4(position,1.0);
    }

    )";

    const char* fs = R"(

    #version 130

    in vec4 FragPosLight;

    uniform sampler2D shadowMap;

    out vec4 color;

    float computeShadow(vec4 pos)
    {
        vec3 proj = pos.xyz / pos.w;
        proj = proj * 0.5 + 0.5;

        float closest = texture(shadowMap,proj.xy).r;
        float current = proj.z;

        float bias = 0.005;

        if(current - bias > closest)
            return 1.0;

        return 0.0;
    }

    void main()
    {
        float shadow = computeShadow(FragPosLight);

        vec3 base = vec3(0.8,0.3,0.3);

        vec3 lighting = (1.0-shadow)*base;

        color = vec4(lighting,1.0);
    }

    )";

    GLuint program = createProgram(vs,fs);

    // ------------------------------------------------------------
    // SHADOW SHADER (depth only)
    // ------------------------------------------------------------

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

    // ------------------------------------------------------------
    // SHADOW MAP
    // ------------------------------------------------------------

    ShadowMap shadow = createShadowMap(1024,1024);

    // ------------------------------------------------------------
    // MAIN LOOP
    // ------------------------------------------------------------

    while(!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        // ----------------------------------------------------
        // IMGUI
        // ----------------------------------------------------

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        drawControls(camera);

        ImGui::Render();

        int w,h;
        glfwGetFramebufferSize(window,&w,&h);

        // ----------------------------------------------------
        // LIGHT MATRICES
        // ----------------------------------------------------

        glm::mat4 lightProj =
            glm::perspective(glm::radians(45.0f),1.0f,0.1f,100.0f);

        glm::mat4 lightView =
            glm::lookAt(glm::vec3(10,10,10),
                        glm::vec3(0,0,0),
                        glm::vec3(0,1,0));

        glm::mat4 lightMVP = lightProj * lightView;

        // ----------------------------------------------------
        // SHADOW PASS
        // ----------------------------------------------------

        beginShadowPass(shadow);

        glUseProgram(shadowProgram);

        GLuint loc = glGetUniformLocation(shadowProgram,"lightMVP");

        glUniformMatrix4fv(loc,1,GL_FALSE,glm::value_ptr(lightMVP));

        glBindVertexArray(sphere.vao);
        sphere.draw();

        glBindVertexArray(planeVAO);
        glDrawArrays(GL_TRIANGLES,0,6);

        endShadowPass();

        // ----------------------------------------------------
        // NORMAL PASS
        // ----------------------------------------------------

        glViewport(0,0,w,h);

        glClearColor(0.2f,0.2f,0.2f,1.0f);
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

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

        glm::mat4 model = glm::mat4(1.0f);

        glm::mat4 MVP = projection * view * model;

        glUseProgram(program);

        GLuint mvpLoc = glGetUniformLocation(program,"MVP");
        GLuint lightLoc = glGetUniformLocation(program,"lightMVP");

        glUniformMatrix4fv(mvpLoc,1,GL_FALSE,glm::value_ptr(MVP));
        glUniformMatrix4fv(lightLoc,1,GL_FALSE,glm::value_ptr(lightMVP));

        // bind shadow map texture
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D,shadow.depthTexture);

        GLuint shadowLoc = glGetUniformLocation(program,"shadowMap");
        glUniform1i(shadowLoc,0);

        // draw sphere
        glBindVertexArray(sphere.vao);
        sphere.draw();

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    // ------------------------------------------------------------
    // CLEANUP
    // ------------------------------------------------------------

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwTerminate();

    return 0;
}