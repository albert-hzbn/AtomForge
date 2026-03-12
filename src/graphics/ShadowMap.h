#pragma once

#include <GL/glew.h>

struct ShadowMap
{
    GLuint depthFBO;
    GLuint depthTexture;

    int width;
    int height;
};

ShadowMap createShadowMap(int width, int height);

void beginShadowPass(const ShadowMap& shadow);

void endShadowPass();

void bindShadowTexture(const ShadowMap& shadow, int textureUnit);