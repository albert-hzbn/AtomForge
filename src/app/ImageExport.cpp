#include "app/ImageExport.h"

#include "graphics/CylinderMesh.h"
#include "graphics/Renderer.h"
#include "graphics/SceneBuffers.h"
#include "graphics/ShadowMap.h"
#include "graphics/SphereMesh.h"
#include "graphics/LowPolyMesh.h"
#include "graphics/BillboardMesh.h"

#include <GL/glew.h>
#include <glm/glm.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#if defined(__clang__) || defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif
#include "third_party/stb_image_write.h"
#if defined(__clang__) || defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace
{
int clampToByte(float value)
{
    const float clamped = std::max(0.0f, std::min(1.0f, value));
    return (int)std::lround(clamped * 255.0f);
}

std::string svgFloat(float value)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.3f", value);
    return std::string(buf);
}

std::string svgRgb(const glm::vec3& color)
{
    char buf[64];
    std::snprintf(buf,
                  sizeof(buf),
                  "rgb(%d,%d,%d)",
                  clampToByte(color.r),
                  clampToByte(color.g),
                  clampToByte(color.b));
    return std::string(buf);
}

bool isOrthographicProjection(const glm::mat4& projection)
{
    return std::abs(projection[3][3] - 1.0f) <= 1e-4f;
}

float pixelsPerWorldUnit(const ImageExportView& view, float viewDepth)
{
    const float scale = std::abs(view.projection[1][1]);
    if (isOrthographicProjection(view.projection))
        return 0.5f * (float)view.height * scale;

    const float safeDepth = std::max(viewDepth, 1e-4f);
    return 0.5f * (float)view.height * scale / safeDepth;
}

bool projectPoint(const glm::vec3& world,
                  const ImageExportView& view,
                  float& screenX,
                  float& screenY,
                  float& ndcDepth,
                  float& viewDepth)
{
    const glm::vec4 clip = view.projection * view.view * glm::vec4(world, 1.0f);
    if (clip.w <= 1e-6f)
        return false;

    const glm::vec3 ndc = glm::vec3(clip) / clip.w;
    if (ndc.z < -1.0f || ndc.z > 1.0f)
        return false;

    screenX = (ndc.x * 0.5f + 0.5f) * (float)view.width;
    screenY = (1.0f - (ndc.y * 0.5f + 0.5f)) * (float)view.height;
    ndcDepth = ndc.z;

    const glm::vec4 viewPos = view.view * glm::vec4(world, 1.0f);
    viewDepth = -viewPos.z;
    return true;
}

void drawSceneToCurrentFramebuffer(const ImageExportView& view,
                                   const glm::vec4& clearColor,
                                   bool showBonds,
                                   bool showAtoms,
                                   bool showBoundingBox,
                                   const SceneBuffers& sceneBuffers,
                                   Renderer& renderer,
                                   const ShadowMap& shadow)
{
    // Shadow passes first (render into shadow FBO)
    // Save the currently-bound framebuffer so we can restore it after the shadow
    // passes (endShadowPass() binds FBO 0, which is wrong when we're rendering
    // into an off-screen export FBO).
    GLint exportFBO = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &exportFBO);

    if (showAtoms)
    {
        switch (sceneBuffers.renderMode)
        {
            case RenderingMode::StandardInstancing:
                renderer.drawShadowPass(shadow, sceneBuffers.tabSphereVAO, sceneBuffers.tabSphereIndexCount, view.lightMVP, sceneBuffers.atomCount);
                break;
            case RenderingMode::LowPolyInstancing:
                renderer.drawShadowPassLowPoly(shadow, sceneBuffers.tabLowPolyVAO, sceneBuffers.tabLowPolyIndexCount, view.lightMVP, sceneBuffers.atomCount);
                break;
            case RenderingMode::BillboardImposters:
                renderer.drawShadowPassBillboard(shadow, sceneBuffers.tabBillboardVAO, sceneBuffers.tabBillboardIndexCount, view.lightMVP, view.view, sceneBuffers.atomCount);
                break;
        }
    }

    if (showBonds && showAtoms)
        renderer.drawBondShadowPass(shadow, sceneBuffers.tabCylinderVAO, sceneBuffers.tabCylinderVertexCount, view.lightMVP, sceneBuffers.bondCount);

    // Restore the export FBO and set viewport to export dimensions
    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)exportFBO);
    glViewport(0, 0, view.width, view.height);
    glEnable(GL_DEPTH_TEST);
    glClearColor(clearColor.r, clearColor.g, clearColor.b, clearColor.a);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (showBonds && showAtoms)
    {
        renderer.drawBonds(
            view.projection,
            view.view,
            view.lightPosition,
            view.cameraPosition,
            sceneBuffers.tabCylinderVAO,
            sceneBuffers.tabCylinderVertexCount,
            sceneBuffers.bondCount);
    }

    // Draw atoms based on rendering mode
    if (showAtoms)
    {
        switch (sceneBuffers.renderMode)
        {
            case RenderingMode::StandardInstancing:
                renderer.drawAtoms(
                    view.projection,
                    view.view,
                    view.lightMVP,
                    view.lightPosition,
                    view.cameraPosition,
                    shadow,
                    sceneBuffers.tabSphereVAO,
                    sceneBuffers.tabSphereIndexCount,
                    sceneBuffers.atomCount);
                break;
            case RenderingMode::LowPolyInstancing:
                renderer.drawAtomsLowPoly(
                    view.projection,
                    view.view,
                    view.lightMVP,
                    view.lightPosition,
                    view.cameraPosition,
                    shadow,
                    sceneBuffers.tabLowPolyVAO,
                    sceneBuffers.tabLowPolyIndexCount,
                    sceneBuffers.atomCount);
                break;
            case RenderingMode::BillboardImposters:
                renderer.drawAtomsBillboard(
                    view.projection,
                    view.view,
                    view.lightMVP,
                    view.lightPosition,
                    view.cameraPosition,
                    shadow,
                    sceneBuffers.tabBillboardVAO,
                    sceneBuffers.tabBillboardIndexCount,
                    sceneBuffers.atomCount);
                break;
        }
    }

    if (showBoundingBox)
    {
        renderer.drawBoxLines(
            view.projection,
            view.view,
            sceneBuffers.lineVAO,
            sceneBuffers.boxLines.size(),
            clearColor.r > 0.5f ? glm::vec3(0.25f) : glm::vec3(0.85f));
    }
}

void flipImageRows(std::vector<unsigned char>& pixels, int width, int height, int channels)
{
    const size_t rowSize = (size_t)width * (size_t)channels;
    std::vector<unsigned char> tempRow(rowSize);

    for (int y = 0; y < height / 2; ++y)
    {
        unsigned char* top = pixels.data() + (size_t)y * rowSize;
        unsigned char* bottom = pixels.data() + (size_t)(height - 1 - y) * rowSize;

        std::copy(top, top + rowSize, tempRow.data());
        std::copy(bottom, bottom + rowSize, top);
        std::copy(tempRow.data(), tempRow.data() + rowSize, bottom);
    }
}

bool captureSceneToRgba(const ImageExportView& view,
                        const glm::vec4& clearColor,
                        bool showBonds,
                        bool showAtoms,
                        bool showBoundingBox,
                        const SceneBuffers& sceneBuffers,
                        Renderer& renderer,
                        const ShadowMap& shadow,
                        std::vector<unsigned char>& outPixels,
                        std::string& errorMessage)
{
    if (view.width <= 0 || view.height <= 0)
    {
        errorMessage = "Invalid viewport size for image export.";
        return false;
    }

    GLint previousFramebuffer = 0;
    GLint previousViewport[4] = {0, 0, 0, 0};
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &previousFramebuffer);
    glGetIntegerv(GL_VIEWPORT, previousViewport);

    GLuint fbo = 0;
    GLuint colorTexture = 0;
    GLuint depthBuffer = 0;

    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    glGenTextures(1, &colorTexture);
    glBindTexture(GL_TEXTURE_2D, colorTexture);
    glTexImage2D(GL_TEXTURE_2D,
                 0,
                 GL_RGBA8,
                 view.width,
                 view.height,
                 0,
                 GL_RGBA,
                 GL_UNSIGNED_BYTE,
                 nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture2D(GL_FRAMEBUFFER,
                           GL_COLOR_ATTACHMENT0,
                           GL_TEXTURE_2D,
                           colorTexture,
                           0);

    glGenRenderbuffers(1, &depthBuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, depthBuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, view.width, view.height);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER,
                              GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER,
                              depthBuffer);

    const GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    if (status != GL_FRAMEBUFFER_COMPLETE)
    {
        errorMessage = "Failed to initialize framebuffer for image export.";

        glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)previousFramebuffer);
        glViewport(previousViewport[0], previousViewport[1], previousViewport[2], previousViewport[3]);

        if (depthBuffer != 0)
            glDeleteRenderbuffers(1, &depthBuffer);
        if (colorTexture != 0)
            glDeleteTextures(1, &colorTexture);
        if (fbo != 0)
            glDeleteFramebuffers(1, &fbo);
        return false;
    }

    drawSceneToCurrentFramebuffer(view,
                                  clearColor,
                                  showBonds,
                                  showAtoms,
                                  showBoundingBox,
                                  sceneBuffers,
                                  renderer,
                                  shadow);

    outPixels.resize((size_t)view.width * (size_t)view.height * 4);
    glReadPixels(0,
                 0,
                 view.width,
                 view.height,
                 GL_RGBA,
                 GL_UNSIGNED_BYTE,
                 outPixels.data());
    flipImageRows(outPixels, view.width, view.height, 4);

    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)previousFramebuffer);
    glViewport(previousViewport[0], previousViewport[1], previousViewport[2], previousViewport[3]);

    glDeleteRenderbuffers(1, &depthBuffer);
    glDeleteTextures(1, &colorTexture);
    glDeleteFramebuffers(1, &fbo);

    return true;
}

// ---------------------------------------------------------------------------
// Software gizmo rasterizer
// ---------------------------------------------------------------------------

void blendPixelRgba(std::vector<unsigned char>& pixels, int width, int height,
                    int px, int py,
                    unsigned char r, unsigned char g, unsigned char b, unsigned char alpha)
{
    if (px < 0 || py < 0 || px >= width || py >= height)
        return;
    const size_t idx = ((size_t)py * (size_t)width + (size_t)px) * 4;
    const float fa = alpha / 255.0f;
    const float fb = 1.0f - fa;
    pixels[idx + 0] = (unsigned char)(pixels[idx + 0] * fb + r * fa);
    pixels[idx + 1] = (unsigned char)(pixels[idx + 1] * fb + g * fa);
    pixels[idx + 2] = (unsigned char)(pixels[idx + 2] * fb + b * fa);
    pixels[idx + 3] = 255;
}

void drawSoftCircleFilled(std::vector<unsigned char>& pixels, int width, int height,
                          float cx, float cy, float radius,
                          unsigned char r, unsigned char g, unsigned char b, unsigned char a)
{
    const int x0 = std::max(0, (int)(cx - radius - 1.5f));
    const int x1 = std::min(width - 1, (int)(cx + radius + 1.5f));
    const int y0 = std::max(0, (int)(cy - radius - 1.5f));
    const int y1 = std::min(height - 1, (int)(cy + radius + 1.5f));
    for (int y = y0; y <= y1; ++y)
    {
        for (int x = x0; x <= x1; ++x)
        {
            const float dx = (float)x - cx;
            const float dy = (float)y - cy;
            const float dist = std::sqrt(dx * dx + dy * dy);
            if (dist <= radius + 0.5f)
            {
                const float blend = std::min(1.0f, radius + 0.5f - dist);
                blendPixelRgba(pixels, width, height, x, y, r, g, b,
                               (unsigned char)((float)a * blend));
            }
        }
    }
}

void drawSoftThickLine(std::vector<unsigned char>& pixels, int width, int height,
                       float x0, float y0, float x1, float y1, float thickness,
                       unsigned char r, unsigned char g, unsigned char b, unsigned char a)
{
    const float dx = x1 - x0;
    const float dy = y1 - y0;
    const float len = std::sqrt(dx * dx + dy * dy);
    if (len < 1e-3f)
        return;
    const float nx = -dy / len;
    const float ny =  dx / len;
    const float half = thickness * 0.5f;

    const float xMin = std::min(x0, x1);
    const float xMax = std::max(x0, x1);
    const float yMin = std::min(y0, y1);
    const float yMax = std::max(y0, y1);

    const int bx0 = std::max(0, (int)(xMin - half - 1.5f));
    const int bx1 = std::min(width  - 1, (int)(xMax + half + 1.5f));
    const int by0 = std::max(0, (int)(yMin - half - 1.5f));
    const int by1 = std::min(height - 1, (int)(yMax + half + 1.5f));

    const float invLen2 = 1.0f / (len * len);
    for (int y = by0; y <= by1; ++y)
    {
        for (int x = bx0; x <= bx1; ++x)
        {
            const float px = (float)x - x0;
            const float py = (float)y - y0;
            const float t = (px * dx + py * dy) * invLen2;
            const float perp = std::abs(px * nx + py * ny);
            if (t >= -0.01f && t <= 1.01f && perp <= half + 0.5f)
            {
                const float blend = std::min(1.0f, half + 0.5f - perp);
                blendPixelRgba(pixels, width, height, x, y, r, g, b,
                               (unsigned char)((float)a * blend));
            }
        }
    }
}

void drawSoftRing(std::vector<unsigned char>& pixels, int width, int height,
                  float cx, float cy, float radius, float ringThickness,
                  unsigned char r, unsigned char g, unsigned char b, unsigned char a)
{
    const int x0 = std::max(0, (int)(cx - radius - ringThickness - 1.5f));
    const int x1 = std::min(width  - 1, (int)(cx + radius + ringThickness + 1.5f));
    const int y0 = std::max(0, (int)(cy - radius - ringThickness - 1.5f));
    const int y1 = std::min(height - 1, (int)(cy + radius + ringThickness + 1.5f));
    for (int y = y0; y <= y1; ++y)
    {
        for (int x = x0; x <= x1; ++x)
        {
            const float dx = (float)x - cx;
            const float dy = (float)y - cy;
            const float dist = std::sqrt(dx * dx + dy * dy);
            const float distFromRing = std::abs(dist - radius);
            if (distFromRing <= ringThickness * 0.5f + 0.5f)
            {
                const float blend = std::min(1.0f, ringThickness * 0.5f + 0.5f - distFromRing);
                blendPixelRgba(pixels, width, height, x, y, r, g, b, (unsigned char)((float)a * blend));
            }
        }
    }
}

void drawGizmoOnPixels(std::vector<unsigned char>& pixels,
                       int width, int height,
                       const glm::mat4& view,
                       bool lightTheme,
                       int scale)
{
    // 5x7 bitmap glyphs for X, Y, Z (row-major, MSB = left pixel)
    static const unsigned char kGlyphX[7] = { 0b10001, 0b10001, 0b01010, 0b00100, 0b01010, 0b10001, 0b10001 };
    static const unsigned char kGlyphY[7] = { 0b10001, 0b10001, 0b01010, 0b00100, 0b00100, 0b00100, 0b00100 };
    static const unsigned char kGlyphZ[7] = { 0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b10000, 0b11111 };
    static const unsigned char* kGlyphs[3] = { kGlyphX, kGlyphY, kGlyphZ };

    // Draw a 5x7 glyph at pixel position (ox, oy), each pixel rendered as a filled cellSize square
    auto drawGlyph = [&](float ox, float oy, const unsigned char* glyph,
                         unsigned char r, unsigned char g, unsigned char b, float ps)
    {
        const float cellSize = std::max(1.0f, ps * 1.5f);
        for (int row = 0; row < 7; ++row)
        {
            for (int col = 0; col < 5; ++col)
            {
                if (glyph[row] & (1 << (4 - col)))
                {
                    const float px = ox + (float)col * cellSize;
                    const float py = oy + (float)row * cellSize;
                    const int ix0 = (int)px;
                    const int iy0 = (int)py;
                    const int ix1 = std::min((int)(px + cellSize), ix0 + (int)cellSize + 1);
                    const int iy1 = std::min((int)(py + cellSize), iy0 + (int)cellSize + 1);
                    for (int y = iy0; y <= iy1; ++y)
                        for (int x = ix0; x <= ix1; ++x)
                            blendPixelRgba(pixels, width, height, x, y, r, g, b, 255);
                }
            }
        }
    };

    const float s         = (float)std::max(1, scale);
    const float bgRadius  = 53.0f * s;
    const float cx        = 73.0f * s;
    const float cy        = (float)height - 73.0f * s;
    const float projScale = 53.0f * s;
    const float axisLength = 0.77f;

    // --- Backdrop: filled disk + 2 ring outlines (matching drawGizmoBackdrop) ---
    const unsigned char bgR = lightTheme ? 235 : 18;
    const unsigned char bgG = lightTheme ? 238 : 22;
    const unsigned char bgB = lightTheme ? 242 : 28;
    const unsigned char bgA = lightTheme ? 200 : 190;
    drawSoftCircleFilled(pixels, width, height, cx, cy, bgRadius, bgR, bgG, bgB, bgA);

    // Outer ring
    const unsigned char ringA1 = lightTheme ? 55 : 55;
    const unsigned char ringR1 = lightTheme ? 0 : 255;
    drawSoftRing(pixels, width, height, cx, cy, bgRadius, s, ringR1, ringR1, ringR1, ringA1);

    // Inner ring (radius - 7)
    const unsigned char ringA2 = lightTheme ? 20 : 20;
    const unsigned char ringR2 = lightTheme ? 0 : 255;
    drawSoftRing(pixels, width, height, cx, cy, bgRadius - 7.0f * s, s, ringR2, ringR2, ringR2, ringA2);

    // --- Axis colors ---
    struct AxisColor { unsigned char r, g, b; };
    const AxisColor axisColors[3] = {
        lightTheme ? AxisColor{200,  50,  50} : AxisColor{235,  92,  92},  // X
        lightTheme ? AxisColor{ 40, 160,  50} : AxisColor{110, 220, 120},  // Y
        lightTheme ? AxisColor{ 40, 110, 220} : AxisColor{110, 175, 255}   // Z
    };
    // Ring outline color around endpoint dot (matching IM_COL32(30,30,40,200) / IM_COL32(245,245,250,200))
    const unsigned char outlineR = lightTheme ?  30 : 245;
    const unsigned char outlineG = lightTheme ?  30 : 245;
    const unsigned char outlineB = lightTheme ?  40 : 250;

    const glm::vec3 worldAxes[3] = {
        {1.0f, 0.0f, 0.0f},
        {0.0f, 1.0f, 0.0f},
        {0.0f, 0.0f, 1.0f}
    };

    struct Axis2D { float sx, sy, depth; int ci; };
    Axis2D axes2D[3];
    const glm::mat3 viewRot(view);
    for (int i = 0; i < 3; ++i)
    {
        const glm::vec3 dir = glm::normalize(viewRot * worldAxes[i]);
        axes2D[i].sx    = cx + dir.x * axisLength * projScale;
        axes2D[i].sy    = cy - dir.y * axisLength * projScale;
        axes2D[i].depth = -dir.z;   // toGizmoDepth: negative Z = toward viewer
        axes2D[i].ci    = i;
    }

    // Sort: most-positive depth (deepest) drawn first — painter's algorithm
    int order[3] = {0, 1, 2};
    std::sort(order, order + 3, [&](int a, int b) {
        return axes2D[a].depth < axes2D[b].depth;
    });

    const float glyphScale = std::max(1.0f, s * 1.6f);
    const float glyphW = 5.0f * glyphScale * 1.5f;
    const float glyphH = 7.0f * glyphScale * 1.5f;

    for (int o = 0; o < 3; ++o)
    {
        const int i = order[o];
        const AxisColor& c = axisColors[axes2D[i].ci];
        const float depth01    = glm::clamp((axes2D[i].depth + 1.0f) * 0.5f, 0.0f, 1.0f);
        const float thickness  = (1.8f + 1.2f * depth01) * s;
        const float dotRadius  = (2.5f + 2.0f * depth01) * s;

        // Line from origin to endpoint
        drawSoftThickLine(pixels, width, height, cx, cy, axes2D[i].sx, axes2D[i].sy,
                          thickness, c.r, c.g, c.b, 255);

        // Filled endpoint dot
        drawSoftCircleFilled(pixels, width, height, axes2D[i].sx, axes2D[i].sy, dotRadius,
                             c.r, c.g, c.b, 255);

        // Ring outline around endpoint dot (matching draw3DAxis AddCircle at dotRadius+1)
        drawSoftRing(pixels, width, height, axes2D[i].sx, axes2D[i].sy,
                     dotRadius + 1.0f * s, s, outlineR, outlineG, outlineB, 200);

        // Label at endpoint + (5.5, 4.0) screen pixels (matching draw3DAxis AddText offset)
        const float lx = axes2D[i].sx + 5.5f * s - glyphW * 0.5f;
        const float ly = axes2D[i].sy + 4.0f * s - glyphH * 0.5f;
        drawGlyph(lx, ly, kGlyphs[axes2D[i].ci], c.r, c.g, c.b, glyphScale);
    }

    // Center dot (matching drawGizmoCenter)
    const unsigned char cdR = lightTheme ?  20 : 240;
    const unsigned char cdG = lightTheme ?  20 : 240;
    const unsigned char cdB = lightTheme ?  25 : 245;
    drawSoftCircleFilled(pixels, width, height, cx, cy, 3.2f * s, cdR, cdG, cdB, 240);
}

bool writeRasterImage(const ImageExportRequest& request,
                      const ImageExportView& view,
                      const glm::vec4& backgroundColor,
                      bool showBonds,
                      bool showAtoms,
                      bool showBoundingBox,
                      const SceneBuffers& sceneBuffers,
                      Renderer& renderer,
                      const ShadowMap& shadow,
                      std::string& errorMessage)
{
    const glm::vec4 clearColor = request.includeBackground
        ? backgroundColor
        : glm::vec4(0.0f, 0.0f, 0.0f, 0.0f);

    std::vector<unsigned char> rgbaPixels;
    if (!captureSceneToRgba(view,
                            clearColor,
                            showBonds,
                            showAtoms,
                            showBoundingBox,
                            sceneBuffers,
                            renderer,
                            shadow,
                            rgbaPixels,
                            errorMessage))
    {
        return false;
    }

    if (request.includeGizmo)
    {
        const bool lightTheme = backgroundColor.r > 0.5f;
        drawGizmoOnPixels(rgbaPixels, view.width, view.height,
                          view.view, lightTheme, request.resolutionScale);
    }

    if (request.format == ImageExportFormat::Png)
    {
        const int ok = stbi_write_png(request.outputPath.c_str(),
                                      view.width,
                                      view.height,
                                      4,
                                      rgbaPixels.data(),
                                      view.width * 4);
        if (ok == 0)
        {
            errorMessage = "Failed to write PNG file.";
            return false;
        }
        return true;
    }

    std::vector<unsigned char> rgbPixels((size_t)view.width * (size_t)view.height * 3);
    for (size_t i = 0, j = 0; i < rgbaPixels.size(); i += 4, j += 3)
    {
        const float a = (float)rgbaPixels[i + 3] / 255.0f;
        if (request.includeBackground)
        {
            rgbPixels[j + 0] = rgbaPixels[i + 0];
            rgbPixels[j + 1] = rgbaPixels[i + 1];
            rgbPixels[j + 2] = rgbaPixels[i + 2];
        }
        else
        {
            rgbPixels[j + 0] = (unsigned char)std::lround((float)rgbaPixels[i + 0] * a + 255.0f * (1.0f - a));
            rgbPixels[j + 1] = (unsigned char)std::lround((float)rgbaPixels[i + 1] * a + 255.0f * (1.0f - a));
            rgbPixels[j + 2] = (unsigned char)std::lround((float)rgbaPixels[i + 2] * a + 255.0f * (1.0f - a));
        }
    }

    const int ok = stbi_write_jpg(request.outputPath.c_str(),
                                  view.width,
                                  view.height,
                                  3,
                                  rgbPixels.data(),
                                  95);
    if (ok == 0)
    {
        errorMessage = "Failed to write JPEG file.";
        return false;
    }

    return true;
}

struct SvgBond
{
    float depth;
    float x1;
    float y1;
    float x2;
    float y2;
    float strokeWidth;
    glm::vec3 color;
};

struct SvgAtom
{
    float depth;
    float x;
    float y;
    float radius;
    glm::vec3 color;
};

bool writeSvgImage(const ImageExportRequest& request,
                   const ImageExportView& view,
                   const glm::vec4& backgroundColor,
                   bool showBonds,
                   bool showAtoms,
                   const SceneBuffers& sceneBuffers,
                   std::string& errorMessage)
{
    std::ofstream output(request.outputPath.c_str(), std::ios::out | std::ios::trunc);
    if (!output)
    {
        errorMessage = "Failed to open SVG output file.";
        return false;
    }

    std::vector<SvgBond> bonds;
    std::vector<SvgAtom> atoms;

    if (showBonds && showAtoms)
    {
        const size_t bondCount = std::min(
            std::min(sceneBuffers.bondStarts.size(), sceneBuffers.bondEnds.size()),
            std::min(sceneBuffers.bondColorsA.size(), sceneBuffers.bondColorsB.size()));

        for (size_t i = 0; i < bondCount; ++i)
        {
            float x1 = 0.0f;
            float y1 = 0.0f;
            float z1 = 0.0f;
            float depth1 = 0.0f;
            float x2 = 0.0f;
            float y2 = 0.0f;
            float z2 = 0.0f;
            float depth2 = 0.0f;
            if (!projectPoint(sceneBuffers.bondStarts[i], view, x1, y1, z1, depth1))
                continue;
            if (!projectPoint(sceneBuffers.bondEnds[i], view, x2, y2, z2, depth2))
                continue;

            const glm::vec3 mid = 0.5f * (sceneBuffers.bondStarts[i] + sceneBuffers.bondEnds[i]);
            float xMid = 0.0f;
            float yMid = 0.0f;
            float zMid = 0.0f;
            float depthMid = 0.0f;
            if (!projectPoint(mid, view, xMid, yMid, zMid, depthMid))
                continue;

            const float bondRadius = (i < sceneBuffers.bondRadiiCpu.size()) ? sceneBuffers.bondRadiiCpu[i] : 0.08f;
            const float strokeWidth = std::max(1.0f, bondRadius * pixelsPerWorldUnit(view, depthMid) * 2.0f);

            SvgBond bond;
            bond.depth = zMid;
            bond.x1 = x1;
            bond.y1 = y1;
            bond.x2 = x2;
            bond.y2 = y2;
            bond.strokeWidth = strokeWidth;
            bond.color = 0.5f * (sceneBuffers.bondColorsA[i] + sceneBuffers.bondColorsB[i]);
            bonds.push_back(bond);
        }

        std::sort(bonds.begin(), bonds.end(),
                  [](const SvgBond& a, const SvgBond& b) { return a.depth > b.depth; });
    }

    const size_t atomCount = std::min(
        std::min(sceneBuffers.atomPositions.size(), sceneBuffers.atomRadii.size()),
        sceneBuffers.atomColors.size());

    if (showAtoms)
    {
    for (size_t i = 0; i < atomCount; ++i)
    {
        float x = 0.0f;
        float y = 0.0f;
        float z = 0.0f;
        float depth = 0.0f;
        if (!projectPoint(sceneBuffers.atomPositions[i], view, x, y, z, depth))
            continue;

        const float radiusPx = std::max(1.0f, sceneBuffers.atomRadii[i] * pixelsPerWorldUnit(view, depth));

        SvgAtom atom;
        atom.depth = z;
        atom.x = x;
        atom.y = y;
        atom.radius = radiusPx;
        atom.color = sceneBuffers.atomColors[i];
        atoms.push_back(atom);
    }

    std::sort(atoms.begin(), atoms.end(),
              [](const SvgAtom& a, const SvgAtom& b) { return a.depth > b.depth; });
    } // if (showAtoms)

    output << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    output << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\""
           << view.width << "\" height=\"" << view.height << "\" viewBox=\"0 0 "
           << view.width << " " << view.height << "\">\n";

    if (request.includeBackground)
    {
        const glm::vec3 bg(backgroundColor.r, backgroundColor.g, backgroundColor.b);
        output << "  <rect x=\"0\" y=\"0\" width=\"" << view.width
               << "\" height=\"" << view.height << "\" fill=\""
               << svgRgb(bg) << "\" />\n";
    }

    for (size_t i = 0; i < bonds.size(); ++i)
    {
        output << "  <line x1=\"" << svgFloat(bonds[i].x1)
               << "\" y1=\"" << svgFloat(bonds[i].y1)
               << "\" x2=\"" << svgFloat(bonds[i].x2)
               << "\" y2=\"" << svgFloat(bonds[i].y2)
               << "\" stroke=\"" << svgRgb(bonds[i].color)
               << "\" stroke-width=\"" << svgFloat(bonds[i].strokeWidth)
               << "\" stroke-linecap=\"round\" />\n";
    }

    for (size_t i = 0; i < atoms.size(); ++i)
    {
        const glm::vec3 stroke = glm::clamp(atoms[i].color * 0.65f, glm::vec3(0.0f), glm::vec3(1.0f));
        output << "  <circle cx=\"" << svgFloat(atoms[i].x)
               << "\" cy=\"" << svgFloat(atoms[i].y)
               << "\" r=\"" << svgFloat(atoms[i].radius)
               << "\" fill=\"" << svgRgb(atoms[i].color)
               << "\" stroke=\"" << svgRgb(stroke)
               << "\" stroke-width=\"1\" />\n";
    }

    output << "</svg>\n";

    if (!output.good())
    {
        errorMessage = "Failed while writing SVG file.";
        return false;
    }

    return true;
}

} // namespace

bool exportStructureImage(const ImageExportRequest& request,
                          const ImageExportView& view,
                          const glm::vec4& backgroundColor,
                          bool showBonds,
                          bool showAtoms,
                          bool showBoundingBox,
                          const SceneBuffers& sceneBuffers,
                          Renderer& renderer,
                          const ShadowMap& shadow,
                          std::string& errorMessage)
{
    if (request.format == ImageExportFormat::Svg)
    {
        return writeSvgImage(request,
                             view,
                             backgroundColor,
                             showBonds,
                             showAtoms,
                             sceneBuffers,
                             errorMessage);
    }

    return writeRasterImage(request,
                            view,
                            backgroundColor,
                            showBonds,
                            showAtoms,
                            showBoundingBox,
                            sceneBuffers,
                            renderer,
                            shadow,
                            errorMessage);
}
