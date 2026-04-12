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
#include "third_party/stb_image_write.h"

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
                                   const SceneBuffers& sceneBuffers,
                                   Renderer& renderer,
                                   const ShadowMap& shadow,
                                   const SphereMesh& sphere,
                                   const LowPolyMesh& lowPolyMesh,
                                   const BillboardMesh& billboardMesh,
                                   const CylinderMesh& cylinder)
{
    // Shadow passes first (render into shadow FBO)
    // Save the currently-bound framebuffer so we can restore it after the shadow
    // passes (endShadowPass() binds FBO 0, which is wrong when we're rendering
    // into an off-screen export FBO).
    GLint exportFBO = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &exportFBO);

    switch (sceneBuffers.renderMode)
    {
        case RenderingMode::StandardInstancing:
            renderer.drawShadowPass(shadow, sphere, view.lightMVP, sceneBuffers.atomCount);
            break;
        case RenderingMode::LowPolyInstancing:
            renderer.drawShadowPassLowPoly(shadow, lowPolyMesh, view.lightMVP, sceneBuffers.atomCount);
            break;
        case RenderingMode::BillboardImposters:
            renderer.drawShadowPassBillboard(shadow, billboardMesh, view.lightMVP, view.view, sceneBuffers.atomCount);
            break;
    }

    if (showBonds)
        renderer.drawBondShadowPass(shadow, cylinder, view.lightMVP, sceneBuffers.bondCount);

    // Restore the export FBO and set viewport to export dimensions
    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)exportFBO);
    glViewport(0, 0, view.width, view.height);
    glEnable(GL_DEPTH_TEST);
    glClearColor(clearColor.r, clearColor.g, clearColor.b, clearColor.a);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (showBonds)
    {
        renderer.drawBonds(
            view.projection,
            view.view,
            view.lightPosition,
            view.cameraPosition,
            cylinder,
            sceneBuffers.bondCount);
    }

    // Draw atoms based on rendering mode
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
                sphere,
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
                lowPolyMesh,
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
                billboardMesh,
                sceneBuffers.atomCount);
            break;
    }

    renderer.drawBoxLines(
        view.projection,
        view.view,
        sceneBuffers.lineVAO,
        sceneBuffers.boxLines.size(),
        clearColor.r > 0.5f ? glm::vec3(0.25f) : glm::vec3(0.85f));
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
                        const SceneBuffers& sceneBuffers,
                        Renderer& renderer,
                        const ShadowMap& shadow,
                        const SphereMesh& sphere,
                        const LowPolyMesh& lowPolyMesh,
                        const BillboardMesh& billboardMesh,
                        const CylinderMesh& cylinder,
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
                                  sceneBuffers,
                                  renderer,
                                  shadow,
                                  sphere,
                                  lowPolyMesh,
                                  billboardMesh,
                                  cylinder);

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

bool writeRasterImage(const ImageExportRequest& request,
                      const ImageExportView& view,
                      const glm::vec4& backgroundColor,
                      bool showBonds,
                      const SceneBuffers& sceneBuffers,
                      Renderer& renderer,
                      const ShadowMap& shadow,
                      const SphereMesh& sphere,
                      const LowPolyMesh& lowPolyMesh,
                      const BillboardMesh& billboardMesh,
                      const CylinderMesh& cylinder,
                      std::string& errorMessage)
{
    const glm::vec4 clearColor = request.includeBackground
        ? backgroundColor
        : glm::vec4(0.0f, 0.0f, 0.0f, 0.0f);

    std::vector<unsigned char> rgbaPixels;
    if (!captureSceneToRgba(view,
                            clearColor,
                            showBonds,
                            sceneBuffers,
                            renderer,
                            shadow,
                            sphere,
                            lowPolyMesh,
                            billboardMesh,
                            cylinder,
                            rgbaPixels,
                            errorMessage))
    {
        return false;
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

    if (showBonds)
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
                          const SceneBuffers& sceneBuffers,
                          Renderer& renderer,
                          const ShadowMap& shadow,
                          const SphereMesh& sphere,
                          const LowPolyMesh& lowPolyMesh,
                          const BillboardMesh& billboardMesh,
                          const CylinderMesh& cylinder,
                          std::string& errorMessage)
{
    if (request.format == ImageExportFormat::Svg)
    {
        return writeSvgImage(request,
                             view,
                             backgroundColor,
                             showBonds,
                             sceneBuffers,
                             errorMessage);
    }

    return writeRasterImage(request,
                            view,
                            backgroundColor,
                            showBonds,
                            sceneBuffers,
                            renderer,
                            shadow,
                            sphere,
                            lowPolyMesh,
                            billboardMesh,
                            cylinder,
                            errorMessage);
}
