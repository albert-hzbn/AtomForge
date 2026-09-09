#include "graphics/ElectronicViewport.h"
#include <GLFW/glfw3.h>
#include <array>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{
std::vector<unsigned char> pixels(GLuint texture)
{
    std::vector<unsigned char> result(128*128*4);
    glBindTexture(GL_TEXTURE_2D,texture);
    glGetTexImage(GL_TEXTURE_2D,0,GL_RGBA,GL_UNSIGNED_BYTE,result.data());
    return result;
}
}

int main()
{
    if(!glfwInit()) return 77;
    glfwWindowHint(GLFW_VISIBLE,GLFW_FALSE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR,3);
    GLFWwindow* window=glfwCreateWindow(128,128,"Electronic viewport test",nullptr,nullptr);
    if(!window) { glfwTerminate(); return 77; }
    glfwMakeContextCurrent(window);
    glewExperimental=GL_TRUE;
    if(glewInit()!=GLEW_OK) { glfwDestroyWindow(window); glfwTerminate(); return 77; }
    while(glGetError()!=GL_NO_ERROR) {}
    int result=0;
    try
    {
        ElectronicViewport viewport;
        atomforge::electronic::Mesh mesh;
        const glm::dvec3 p[]={{1,0,0},{-1,0,0},{0,.7,0},{0,-.7,0},{0,0,.45},{0,0,-.45}};
        const int faces[][3]={{0,2,4},{2,1,4},{1,3,4},{3,0,4},{2,0,5},{1,2,5},{3,1,5},{0,3,5}};
        for(const auto& face:faces) for(int i:face) { mesh.vertices.push_back(p[i]); mesh.colors.push_back(.25); }
        viewport.setMesh(mesh);
        glViewport(3,4,50,60); glEnable(GL_SCISSOR_TEST); glDisable(GL_DEPTH_TEST);
        glClearColor(.1f,.2f,.3f,.4f);
        const auto texture=viewport.render(128,128,0,0,1,{0,0},1,0,1,0,.3f,40);
        GLint rectangle[4]; glGetIntegerv(GL_VIEWPORT,rectangle);
        if(rectangle[0]!=3 || rectangle[1]!=4 || rectangle[2]!=50 || rectangle[3]!=60 || !glIsEnabled(GL_SCISSOR_TEST) || glIsEnabled(GL_DEPTH_TEST))
            throw std::runtime_error("Offscreen render leaked editor GL state");
        const auto a=pixels(texture);
        if(a[0]<240 || a[1]<240 || a[2]<240) throw std::runtime_error("Viewport background is not light");
        const auto center=(64*128+64)*4;
        if(a[center+2]<a[center] || a[center]>200) throw std::runtime_error("Density color mapping was not rendered");
        const auto b=pixels(viewport.render(128,128,.8f,.3f,1,{0,0},1,0,1,0,.3f,40));
        if(a==b) throw std::runtime_error("Orbit did not change rendered geometry");
        const auto c=pixels(viewport.render(128,128,0,0,1,{0,0},1,0,.3f,0,0,40));
        if(c[center]<=a[center]) throw std::runtime_error("Changing the density range did not change the surface color");
        const auto d=pixels(viewport.render(128,128,0,0,1,{0,0},1,0,1,0,0,40));
        if(a==d) throw std::runtime_error("Specular lighting did not affect the surface");
        atomforge::electronic::Grid volume;
        volume.shape={16,16,16}; volume.values.resize(16*16*16);
        for(int z=0;z<16;++z) for(int y=0;y<16;++y) for(int x=0;x<16;++x)
        {
            const auto p=glm::dvec3(x,y,z)/15.0-glm::dvec3(.5);
            volume.values[volume.index(x,y,z)]=std::exp(-24*glm::dot(p,p));
        }
        viewport.setVolume(volume);
        const auto clear=pixels(viewport.renderVolume(128,128,0,0,1,{0,0},0,.05f,0,1,0));
        const auto filled=pixels(viewport.renderVolume(128,128,0,0,1,{0,0},.5f,.05f,0,1,0));
        if(clear[center]<240 || clear==filled) throw std::runtime_error("Volume transparency did not reveal interior density");
        const auto empty=pixels(viewport.renderVolume(128,128,0,0,1,{0,0},.5f,2,0,1,0));
        if(clear!=empty) throw std::runtime_error("Volume threshold did not exclude out-of-range data");
        volume.values.assign(16*16*16,.4);
        viewport.setVolume(volume);
        const auto constantVolume=pixels(viewport.renderVolume(128,128,0,0,1,{0,0},.5f,.4f,.4f,.4f,0));
        if(clear==constantVolume) throw std::runtime_error("Constant field volume was not rendered");
        if(glGetError()!=GL_NO_ERROR) throw std::runtime_error("OpenGL error during viewport rendering");
        std::cout << "Viewport pixels, orbit, lighting, palette and GL state passed\n";
    }
    catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; result=1; }
    glfwDestroyWindow(window); glfwTerminate();
    return result;
}
