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

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_COMPARE_REF_TO_TEXTURE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);

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

    // Depth clamping must be OFF during the shadow pass: clamped vertices would
    // write incorrect near/far-boundary depths into the depth texture, producing
    // phantom shadow samples in the colour pass.
    glDisable(GL_DEPTH_CLAMP);

    // Only write back faces into the shadow map.  Front faces can then never
    // appear to be shadowed by their own surface, eliminating self-shadow acne
    // even where spheres overlap each other.
    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);

    // Push shadow depths slightly away from the light so sphere front-faces
    // never self-intersect their own shadow map entry (eliminates acne/moire).
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(3.0f, 6.0f);

    glClear(GL_DEPTH_BUFFER_BIT);

    // disable color writes (important for shadow pass)
    glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);
}

void endShadowPass()
{
    glDisable(GL_POLYGON_OFFSET_FILL);
    glDisable(GL_CULL_FACE);

    // Restore depth clamping for the colour pass (prevents near-plane clipping
    // of sphere geometry when the camera zooms in close).
    glEnable(GL_DEPTH_CLAMP);

    glBindFramebuffer(GL_FRAMEBUFFER,0);

    // re-enable color writes
    glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
}