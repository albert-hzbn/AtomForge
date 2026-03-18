#include "ShadowMap.h"
#include <iostream>

ShadowMap createShadowMap(int width, int height)
{
    ShadowMap shadow;
    shadow.width = width;
    shadow.height = height;

    // --------------------------------------
    // Create framebuffer
    // --------------------------------------

    glGenFramebuffers(1, &shadow.depthFBO);

    // --------------------------------------
    // Create depth texture
    // --------------------------------------

    glGenTextures(1, &shadow.depthTexture);
    glBindTexture(GL_TEXTURE_2D, shadow.depthTexture);

    glTexImage2D(GL_TEXTURE_2D,
                 0,
                 GL_DEPTH_COMPONENT24,
                 width,
                 height,
                 0,
                 GL_DEPTH_COMPONENT,
                 GL_FLOAT,
                 nullptr);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

    float borderColor[] = {1.0f,1.0f,1.0f,1.0f};
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

    // --------------------------------------
    // Attach depth texture to framebuffer
    // --------------------------------------

    glBindFramebuffer(GL_FRAMEBUFFER, shadow.depthFBO);

    glFramebufferTexture2D(GL_FRAMEBUFFER,
                           GL_DEPTH_ATTACHMENT,
                           GL_TEXTURE_2D,
                           shadow.depthTexture,
                           0);

    // Disable color output
    glDrawBuffer(GL_NONE);
    glReadBuffer(GL_NONE);

    // --------------------------------------
    // Check framebuffer completeness
    // --------------------------------------

    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        std::cerr<<"Shadow framebuffer not complete!\n";
    }

    glBindFramebuffer(GL_FRAMEBUFFER,0);

    return shadow;
}

void beginShadowPass(const ShadowMap& shadow)
{
    glViewport(0,0,shadow.width,shadow.height);

    glBindFramebuffer(GL_FRAMEBUFFER,shadow.depthFBO);

    glEnable(GL_DEPTH_TEST);

    glClear(GL_DEPTH_BUFFER_BIT);

    // disable color writes (important for shadow pass)
    glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);
}

void endShadowPass()
{
    glBindFramebuffer(GL_FRAMEBUFFER,0);

    // re-enable color writes
    glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
}