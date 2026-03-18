#include "Shader.h"

#include <algorithm>
#include <iostream>
#include <vector>

static GLuint compile(GLenum type,const char* src)
{
    GLuint s=glCreateShader(type);
    glShaderSource(s,1,&src,nullptr);
    glCompileShader(s);

    GLint ok = GL_FALSE;
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (ok != GL_TRUE)
    {
        GLint logLength = 0;
        glGetShaderiv(s, GL_INFO_LOG_LENGTH, &logLength);
        std::vector<char> log((size_t)std::max(logLength, 1), '\0');
        glGetShaderInfoLog(s, logLength, nullptr, log.data());
        const char* shaderName = (type == GL_VERTEX_SHADER) ? "vertex" : "fragment";
        std::cerr << "Failed to compile " << shaderName << " shader:\n" << log.data() << "\n";
        glDeleteShader(s);
        return 0;
    }

    return s;
}

GLuint createProgram(const char* vs,const char* fs)
{
    GLuint v=compile(GL_VERTEX_SHADER,vs);
    GLuint f=compile(GL_FRAGMENT_SHADER,fs);

    if (v == 0 || f == 0)
    {
        if (v != 0) glDeleteShader(v);
        if (f != 0) glDeleteShader(f);
        return 0;
    }

    GLuint p=glCreateProgram();

    glAttachShader(p,v);
    glAttachShader(p,f);

    glBindAttribLocation(p,0,"position");

    glLinkProgram(p);

    GLint ok = GL_FALSE;
    glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (ok != GL_TRUE)
    {
        GLint logLength = 0;
        glGetProgramiv(p, GL_INFO_LOG_LENGTH, &logLength);
        std::vector<char> log((size_t)std::max(logLength, 1), '\0');
        glGetProgramInfoLog(p, logLength, nullptr, log.data());
        std::cerr << "Failed to link shader program:\n" << log.data() << "\n";

        glDeleteProgram(p);
        glDeleteShader(v);
        glDeleteShader(f);
        return 0;
    }

    glDeleteShader(v);
    glDeleteShader(f);

    return p;
}