#include "graphics/ElectronicViewport.h"
#include "graphics/Shader.h"
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <stdexcept>

namespace
{
const char* vertexShader = R"(#version 330 core
layout(location=0) in vec3 position;
layout(location=1) in vec3 normal;
layout(location=2) in float scalar;
uniform mat4 view;
uniform mat4 projection;
out vec3 eyePosition;
out vec3 eyeNormal;
out float value;
void main() {
    vec4 p = view * vec4(position,1.0);
    eyePosition = p.xyz;
    eyeNormal = mat3(view) * normal;
    value = scalar;
    gl_Position = projection * p;
})";
const char* fragmentShader = R"(#version 330 core
in vec3 eyePosition;
in vec3 eyeNormal;
in float value;
uniform float low;
uniform float high;
uniform float opacity;
uniform float specular;
uniform float shininess;
uniform int palette;
out vec4 outputColor;
vec3 ramp(float t) {
    if (palette == 1) return t < 0.5 ? mix(vec3(0.10,0.28,0.80),vec3(0.96),2.0*t)
                                                  : mix(vec3(0.96),vec3(0.80,0.10,0.12),2.0*t-1.0);
    if (palette == 2) return mix(vec3(0.91,0.95,1.0),vec3(0.08,0.24,0.70),t);
    if (t < 0.25) return mix(vec3(0.10,0.20,0.85),vec3(0.0,0.75,0.95),4.0*t);
    if (t < 0.5) return mix(vec3(0.0,0.75,0.95),vec3(0.12,0.8,0.3),4.0*t-1.0);
    if (t < 0.75) return mix(vec3(0.12,0.8,0.3),vec3(1.0,0.85,0.08),4.0*t-2.0);
    return mix(vec3(1.0,0.85,0.08),vec3(0.85,0.10,0.08),4.0*t-3.0);
}
void main() {
    float t = high > low ? clamp((value-low)/(high-low),0.0,1.0) : 0.5;
    vec3 n = normalize(eyeNormal);
    vec3 v = normalize(-eyePosition);
    if (dot(n,v) < 0.0) n = -n;
    vec3 key = normalize(vec3(-0.4,0.6,1.0));
    vec3 fill = normalize(vec3(0.7,0.1,0.6));
    float diffuse = 0.27 + 0.57*max(dot(n,key),0.0) + 0.16*max(dot(n,fill),0.0);
    float highlight = specular * pow(max(dot(n,normalize(key+v)),0.0),shininess);
    outputColor = vec4(clamp(ramp(t)*diffuse + vec3(highlight),0.0,1.0),opacity);
})";

// Offscreen drawing must not change the editor or ImGui renderer's GL state.
struct State
{
    GLint drawFbo, readFbo, viewport[4], program, vao, buffer, texture, texture3D, renderbuffer;
    GLint srcRgb, dstRgb, srcAlpha, dstAlpha, equationRgb, equationAlpha, depthFunc;
    GLboolean depth, blend, cull, scissor, depthMask, colorMask[4];
    GLfloat clearColor[4];
    GLdouble clearDepth;
    State()
    {
        glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING,&drawFbo); glGetIntegerv(GL_READ_FRAMEBUFFER_BINDING,&readFbo);
        glGetIntegerv(GL_VIEWPORT,viewport); glGetIntegerv(GL_CURRENT_PROGRAM,&program);
        glGetIntegerv(GL_VERTEX_ARRAY_BINDING,&vao); glGetIntegerv(GL_ARRAY_BUFFER_BINDING,&buffer);
        glGetIntegerv(GL_TEXTURE_BINDING_2D,&texture); glGetIntegerv(GL_RENDERBUFFER_BINDING,&renderbuffer);
        glGetIntegerv(GL_TEXTURE_BINDING_3D,&texture3D);
        glGetIntegerv(GL_BLEND_SRC_RGB,&srcRgb); glGetIntegerv(GL_BLEND_DST_RGB,&dstRgb);
        glGetIntegerv(GL_BLEND_SRC_ALPHA,&srcAlpha); glGetIntegerv(GL_BLEND_DST_ALPHA,&dstAlpha);
        glGetIntegerv(GL_BLEND_EQUATION_RGB,&equationRgb); glGetIntegerv(GL_BLEND_EQUATION_ALPHA,&equationAlpha);
        glGetIntegerv(GL_DEPTH_FUNC,&depthFunc); glGetBooleanv(GL_DEPTH_WRITEMASK,&depthMask);
        glGetBooleanv(GL_COLOR_WRITEMASK,colorMask); glGetFloatv(GL_COLOR_CLEAR_VALUE,clearColor);
        glGetDoublev(GL_DEPTH_CLEAR_VALUE,&clearDepth);
        depth=glIsEnabled(GL_DEPTH_TEST); blend=glIsEnabled(GL_BLEND); cull=glIsEnabled(GL_CULL_FACE); scissor=glIsEnabled(GL_SCISSOR_TEST);
    }
    ~State()
    {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER,drawFbo); glBindFramebuffer(GL_READ_FRAMEBUFFER,readFbo);
        glViewport(viewport[0],viewport[1],viewport[2],viewport[3]); glUseProgram(program);
        glBindVertexArray(vao); glBindBuffer(GL_ARRAY_BUFFER,buffer);
        glBindTexture(GL_TEXTURE_2D,texture); glBindRenderbuffer(GL_RENDERBUFFER,renderbuffer);
        glBindTexture(GL_TEXTURE_3D,texture3D);
        glBlendFuncSeparate(srcRgb,dstRgb,srcAlpha,dstAlpha); glBlendEquationSeparate(equationRgb,equationAlpha);
        glDepthFunc(depthFunc); glDepthMask(depthMask); glColorMask(colorMask[0],colorMask[1],colorMask[2],colorMask[3]);
        glClearColor(clearColor[0],clearColor[1],clearColor[2],clearColor[3]); glClearDepth(clearDepth);
        auto restore=[](GLenum flag,GLboolean enabled) { if(enabled) glEnable(flag); else glDisable(flag); };
        restore(GL_DEPTH_TEST,depth); restore(GL_BLEND,blend); restore(GL_CULL_FACE,cull); restore(GL_SCISSOR_TEST,scissor);
    }
};
}

glm::vec3 ElectronicViewport::color(float t, int palette)
{
    t=std::clamp(t,0.0f,1.0f);
    if(palette==1) return t<0.5f ? glm::mix(glm::vec3(.10f,.28f,.80f),glm::vec3(.96f),2*t)
                                      : glm::mix(glm::vec3(.96f),glm::vec3(.80f,.10f,.12f),2*t-1);
    if(palette==2) return glm::mix(glm::vec3(.91f,.95f,1),glm::vec3(.08f,.24f,.70f),t);
    const glm::vec3 stops[]={{.10f,.20f,.85f},{0,.75f,.95f},{.12f,.8f,.3f},{1,.85f,.08f},{.85f,.10f,.08f}};
    const int i=std::min(3,static_cast<int>(t*4));
    return glm::mix(stops[i],stops[i+1],t*4-i);
}

ElectronicViewport::~ElectronicViewport() { if(glfwGetCurrentContext()) release(); }
void ElectronicViewport::release()
{
    glDeleteFramebuffers(1,&m_fbo); glDeleteTextures(1,&m_texture); glDeleteRenderbuffers(1,&m_depth);
    glDeleteBuffers(1,&m_vbo); glDeleteBuffers(1,&m_ebo); glDeleteVertexArrays(1,&m_vao);
    if(m_program) glDeleteProgram(m_program);
    glDeleteTextures(1,&m_volumeTexture);
    if(m_volumeProgram) glDeleteProgram(m_volumeProgram);
}

void ElectronicViewport::setMesh(const atomforge::electronic::Mesh& mesh)
{
    m_vertices.clear(); m_indices.clear(); m_dirty=true; m_sortedYaw=100;
    if(mesh.vertices.empty()) return;
    auto low=mesh.vertices[0], high=low;
    for(auto p:mesh.vertices) { low=glm::min(low,p); high=glm::max(high,p); }
    const auto center=(low+high)*0.5;
    const double radius=std::max(1e-12,glm::length(high-low)*0.5);
    using Key=std::array<long long,3>;
    auto key=[](glm::vec3 p) { return Key{std::llround(p.x*1e6),std::llround(p.y*1e6),std::llround(p.z*1e6)}; };
    std::map<Key,glm::vec3> normals;
    m_vertices.reserve(mesh.vertices.size());
    for(std::size_t i=0;i<mesh.vertices.size();++i)
        m_vertices.push_back({glm::vec3((mesh.vertices[i]-center)/radius),glm::vec3(0),static_cast<float>(mesh.colors.at(i))});
    for(std::size_t i=0;i+2<m_vertices.size();i+=3)
    {
        auto n=glm::cross(m_vertices[i+1].position-m_vertices[i].position,m_vertices[i+2].position-m_vertices[i].position);
        for(int j=0;j<3;++j)
        {
            auto k=key(m_vertices[i+j].position);
            auto found=normals.emplace(k,glm::vec3(0)).first;
            found->second+=n;
        }
    }
    for(auto& v:m_vertices)
    {
        const auto n=normals.at(key(v.position));
        v.normal=glm::length(n)>1e-12f ? glm::normalize(n) : glm::vec3(0,0,1);
    }
    m_indices.resize(m_vertices.size());
    std::iota(m_indices.begin(),m_indices.end(),0u);
}

GLuint ElectronicViewport::render(int width,int height,float yaw,float pitch,float zoom,glm::vec2 pan,
    float opacity,float low,float high,int palette,float specular,float shininess)
{
    State state;
    width=std::clamp(width,1,4096); height=std::clamp(height,1,4096);
    if(!m_program)
    {
        m_program=createProgram(vertexShader,fragmentShader);
        if(!m_program) throw std::runtime_error("Electronic viewport shader failed to compile");
        glGenFramebuffers(1,&m_fbo); glGenTextures(1,&m_texture); glGenRenderbuffers(1,&m_depth);
        glGenVertexArrays(1,&m_vao); glGenBuffers(1,&m_vbo); glGenBuffers(1,&m_ebo);
    }
    glBindFramebuffer(GL_FRAMEBUFFER,m_fbo);
    if(width!=m_width || height!=m_height)
    {
        glBindTexture(GL_TEXTURE_2D,m_texture);
        glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA8,width,height,0,GL_RGBA,GL_UNSIGNED_BYTE,nullptr);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR); glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
        glBindRenderbuffer(GL_RENDERBUFFER,m_depth); glRenderbufferStorage(GL_RENDERBUFFER,GL_DEPTH_COMPONENT24,width,height);
        glFramebufferTexture2D(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0,GL_TEXTURE_2D,m_texture,0);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER,GL_DEPTH_ATTACHMENT,GL_RENDERBUFFER,m_depth);
        if(glCheckFramebufferStatus(GL_FRAMEBUFFER)!=GL_FRAMEBUFFER_COMPLETE) throw std::runtime_error("Electronic viewport framebuffer is incomplete");
        m_width=width; m_height=height;
    }
    glViewport(0,0,width,height); glDisable(GL_SCISSOR_TEST); glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST); glDepthFunc(GL_LESS); glDepthMask(GL_TRUE); glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
    glClearColor(.96f,.97f,.985f,1); glClearDepth(1); glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    glBindVertexArray(m_vao); glBindBuffer(GL_ARRAY_BUFFER,m_vbo); glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,m_ebo);
    const auto rotation=glm::rotate(glm::rotate(glm::mat4(1),pitch,glm::vec3(1,0,0)),yaw,glm::vec3(0,1,0));
    if(m_dirty)
    {
        glBufferData(GL_ARRAY_BUFFER,m_vertices.size()*sizeof(Vertex),m_vertices.data(),GL_STATIC_DRAW);
        for(int a=0;a<3;++a) glEnableVertexAttribArray(a);
        glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(Vertex),reinterpret_cast<void*>(offsetof(Vertex,position)));
        glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,sizeof(Vertex),reinterpret_cast<void*>(offsetof(Vertex,normal)));
        glVertexAttribPointer(2,1,GL_FLOAT,GL_FALSE,sizeof(Vertex),reinterpret_cast<void*>(offsetof(Vertex,scalar)));
    }
    if(m_dirty || yaw!=m_sortedYaw || pitch!=m_sortedPitch)
    {
        std::vector<unsigned int> triangles(m_vertices.size()/3);
        std::iota(triangles.begin(),triangles.end(),0u);
        std::vector<float> depths(triangles.size());
        for(auto i:triangles) depths[i]=(rotation*glm::vec4((m_vertices[i*3].position+m_vertices[i*3+1].position+m_vertices[i*3+2].position)/3.0f,1)).z;
        std::sort(triangles.begin(),triangles.end(),[&](unsigned int a,unsigned int b) { return depths[a]<depths[b]; });
        for(std::size_t i=0;i<triangles.size();++i) for(unsigned int j=0;j<3;++j) m_indices[i*3+j]=triangles[i]*3+j;
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,m_indices.size()*sizeof(unsigned int),m_indices.data(),GL_DYNAMIC_DRAW);
        m_sortedYaw=yaw; m_sortedPitch=pitch; m_dirty=false;
    }
    const float aspect=static_cast<float>(width)/height, extent=1.25f/std::max(.1f,zoom);
    const auto projection=glm::ortho(-extent*aspect,extent*aspect,-extent,extent,.1f,10.0f);
    const auto view=glm::translate(glm::mat4(1),glm::vec3(pan.x,pan.y,-3))*rotation;
    glUseProgram(m_program);
    auto uniform=[&](const char* name,float value) { glUniform1f(glGetUniformLocation(m_program,name),value); };
    glUniformMatrix4fv(glGetUniformLocation(m_program,"view"),1,GL_FALSE,glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(m_program,"projection"),1,GL_FALSE,glm::value_ptr(projection));
    uniform("low",low); uniform("high",high); uniform("opacity",opacity); uniform("specular",specular); uniform("shininess",shininess);
    glUniform1i(glGetUniformLocation(m_program,"palette"),palette);
    glEnable(GL_BLEND); glBlendEquation(GL_FUNC_ADD); glBlendFuncSeparate(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA,GL_ONE,GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(opacity>=.999f ? GL_TRUE : GL_FALSE);
    glDrawElements(GL_TRIANGLES,static_cast<GLsizei>(m_indices.size()),GL_UNSIGNED_INT,nullptr);
    return m_texture;
}

void ElectronicViewport::setVolume(const atomforge::electronic::Grid& grid)
{
    grid.validate();
    State state;
    // Bound GPU memory and upload cost. Numerical grids and exports remain full resolution.
    for(int i=0;i<3;++i) m_volumeShape[i]=std::clamp(grid.shape[i]+(grid.periodic ? 1 : 0),2,96);
    std::vector<float> values(static_cast<std::size_t>(m_volumeShape.x)*m_volumeShape.y*m_volumeShape.z);
    for(int z=0;z<m_volumeShape.z;++z) for(int y=0;y<m_volumeShape.y;++y) for(int x=0;x<m_volumeShape.x;++x)
    {
        const glm::dvec3 f(static_cast<double>(x)/(m_volumeShape.x-1),static_cast<double>(y)/(m_volumeShape.y-1),static_cast<double>(z)/(m_volumeShape.z-1));
        values[(static_cast<std::size_t>(z)*m_volumeShape.y+y)*m_volumeShape.x+x]=static_cast<float>(grid.sample(grid.origin+grid.cell*f));
    }
    double radius=0;
    for(int z : {-1,1}) for(int y : {-1,1}) for(int x : {-1,1})
        radius=std::max(radius,glm::length(grid.cell*glm::dvec3(x,y,z)*.5));
    m_inverseCell=glm::mat3(glm::inverse(grid.cell/radius));
    if(!m_volumeTexture) glGenTextures(1,&m_volumeTexture);
    glBindTexture(GL_TEXTURE_3D,m_volumeTexture);
    glTexImage3D(GL_TEXTURE_3D,0,GL_R32F,m_volumeShape.x,m_volumeShape.y,m_volumeShape.z,0,GL_RED,GL_FLOAT,values.data());
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    for(auto axis : {GL_TEXTURE_WRAP_S,GL_TEXTURE_WRAP_T,GL_TEXTURE_WRAP_R}) glTexParameteri(GL_TEXTURE_3D,axis,GL_CLAMP_TO_EDGE);
}

GLuint ElectronicViewport::renderVolume(int width,int height,float yaw,float pitch,float zoom,glm::vec2 pan,
    float opacity,float threshold,float low,float high,int palette)
{
    // Reuse the surface framebuffer allocation and state restoration. A fully
    // transparent surface contributes nothing to the subsequent volume image.
    render(width,height,yaw,pitch,zoom,pan,0,low,high,palette,0,4);
    State state;
    if(!m_volumeProgram)
    {
        const char* vertex=R"(#version 330 core
out vec2 screen;
void main() {
    vec2 p=vec2((gl_VertexID<<1)&2,gl_VertexID&2);
    screen=p*2.0-1.0;
    gl_Position=vec4(screen,0,1);
})";
        // Share the exact transfer palette with the isosurface renderer.
        const std::string surface(fragmentShader);
        const auto begin=surface.find("vec3 ramp"), end=surface.find("void main",begin);
        const std::string fragment=std::string(R"(#version 330 core
in vec2 screen;
out vec4 outputColor;
uniform sampler3D density;
uniform mat3 inverseCell;
uniform mat3 inverseRotation;
uniform vec3 dimensions;
uniform vec2 pan;
uniform float aspect, extent, opacity, threshold, low, high;
uniform int palette;
)")+surface.substr(begin,end-begin)+R"(
void main() {
    vec3 background=vec3(.96,.97,.985);
    vec3 origin=inverseCell*(inverseRotation*vec3(screen.x*extent*aspect-pan.x,screen.y*extent-pan.y,2.0))+vec3(.5);
    vec3 direction=inverseCell*(inverseRotation*vec3(0,0,-1));
    // Slab intersection handles rays parallel to a cell face without NaNs.
    float entry=0.0, leave=1e20;
    for(int axis=0;axis<3;++axis) {
        if(abs(direction[axis])<1e-8) {
            if(origin[axis]<0.0 || origin[axis]>1.0) { outputColor=vec4(background,1); return; }
        } else {
            float a=-origin[axis]/direction[axis], b=(1.0-origin[axis])/direction[axis];
            entry=max(entry,min(a,b)); leave=min(leave,max(a,b));
        }
    }
    if(leave<=entry) { outputColor=vec4(background,1); return; }
    float stepLength=(leave-entry)/256.0;
    vec4 accumulated=vec4(0);
    for(int i=0;i<256;++i) {
        vec3 p=origin+direction*(entry+(float(i)+.5)*stepLength);
        float value=texture(density,(p*(dimensions-1.0)+.5)/dimensions).r;
        if(value<threshold) continue;
        float t=high>low ? clamp((value-low)/(high-low),0.0,1.0) : .5;
        float alpha=1.0-exp(-opacity*12.0*max(.05,t)*stepLength);
        accumulated.rgb+=(1.0-accumulated.a)*alpha*ramp(t);
        accumulated.a+=(1.0-accumulated.a)*alpha;
        if(accumulated.a>.995) break;
    }
    outputColor=vec4(accumulated.rgb+(1.0-accumulated.a)*background,1);
})";
        m_volumeProgram=createProgram(vertex,fragment.c_str());
        if(!m_volumeProgram) throw std::runtime_error("Volume shader failed to compile");
    }
    glBindFramebuffer(GL_FRAMEBUFFER,m_fbo);
    glViewport(0,0,m_width,m_height);
    glDisable(GL_DEPTH_TEST); glDisable(GL_BLEND); glDisable(GL_CULL_FACE); glDisable(GL_SCISSOR_TEST);
    glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
    glBindVertexArray(m_vao); glUseProgram(m_volumeProgram);
    GLint active; glGetIntegerv(GL_ACTIVE_TEXTURE,&active);
    glBindTexture(GL_TEXTURE_3D,m_volumeTexture);
    glUniform1i(glGetUniformLocation(m_volumeProgram,"density"),active-GL_TEXTURE0);
    const auto inverseRotation=glm::transpose(glm::mat3(glm::rotate(glm::rotate(glm::mat4(1),pitch,glm::vec3(1,0,0)),yaw,glm::vec3(0,1,0))));
    glUniformMatrix3fv(glGetUniformLocation(m_volumeProgram,"inverseCell"),1,GL_FALSE,glm::value_ptr(m_inverseCell));
    glUniformMatrix3fv(glGetUniformLocation(m_volumeProgram,"inverseRotation"),1,GL_FALSE,glm::value_ptr(inverseRotation));
    const glm::vec3 dimensions(m_volumeShape);
    glUniform3fv(glGetUniformLocation(m_volumeProgram,"dimensions"),1,glm::value_ptr(dimensions));
    glUniform2fv(glGetUniformLocation(m_volumeProgram,"pan"),1,glm::value_ptr(pan));
    auto uniform=[&](const char* name,float value) { glUniform1f(glGetUniformLocation(m_volumeProgram,name),value); };
    uniform("aspect",static_cast<float>(m_width)/m_height); uniform("extent",1.25f/std::max(.1f,zoom));
    uniform("opacity",opacity); uniform("threshold",threshold); uniform("low",low); uniform("high",high);
    glUniform1i(glGetUniformLocation(m_volumeProgram,"palette"),palette);
    glDrawArrays(GL_TRIANGLES,0,3);
    return m_texture;
}
