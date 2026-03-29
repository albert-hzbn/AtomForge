#pragma once
#include <GL/glew.h>

GLuint createProgram(const char* vs,const char* fs);
GLuint createComputeProgram(const char* cs);