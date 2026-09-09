#pragma once

#include "electronic/Volume.h"
#include "util/BackgroundTask.h"
#include "ui/PathPicker.h"
#include "ui/ElectronicSliceViewport.h"
#include "graphics/ElectronicViewport.h"
#include <string>
#include <deque>
#include <utility>

struct ElectronicPostProcessingDialog
{
    void drawMenuItem();
    void drawDialog();
    bool isOpen() const { return m_open; }
    void feedDroppedFile(const std::string& path);

private:
    struct Output
    {
        atomforge::electronic::Volume volume;
        std::vector<double> table;
        int columns = 0;
        std::string heading;
        atomforge::electronic::Mesh mesh;
        bool loaded = false;
        bool referenceLoaded = false;
        bool appendMesh = false;
        bool surfaceReady = false;
        double colorLow = 0, colorHigh = 1;
        std::string colorUnit;
        std::string sourcePath;
    };
    void drawPreview();
    void draw3DPreview();
    void load(const std::string& path, bool reference);
    void save(const std::string& path);
    void resetCamera();
    bool m_open = false;
    bool m_dropReference = false;
    std::deque<std::pair<std::string,bool>> m_pendingDrops;
    char m_output[2048] = "electronic-result.xsf";
    char m_reflections[8192] = "0 0 0\n1 0 0\n-1 0 0";
    char m_charges[8192]{};
    int m_quantity = 0;
    int m_cubeUnits = 0;
    int m_selected = 0;
    int m_reference = 0;
    int m_operation = 22;
    int m_axis = 2;
    int m_count = 100;
    int m_window = 3;
    int m_radius = 2;
    int m_exportFormat = 0;
    float m_scalar = 0.1f;
    float m_surfaceLevel = .1f;
    float m_suggestedLevel = .1f;
    int m_renderMode = 0;
    bool m_volumeDirty = true;
    bool m_generateSurface = false;
    int m_viewLayout = 0;
    ElectronicSliceViewport m_sliceViewport;
    bool m_showSlicePlane = true;
    float m_sigma = 0.5f;
    float m_start[3]{};
    float m_end[3] = {1, 0, 0};
    float m_u[3] = {1, 0, 0};
    float m_v[3] = {0, 1, 0};
    float m_alpha = 1;
    float m_realCutoff = 8;
    float m_reciprocalCutoff = 12;
    float m_yaw = 0.6f;
    float m_pitch = 0.4f;
    float m_opacity = .35f;
    float m_zoom = 1.0f;
    glm::vec2 m_pan{0.0f};
    float m_specular = 0.28f;
    float m_shininess = 40.0f;
    float m_colorLow = 0, m_colorHigh = 1;
    float m_autoLow = 0, m_autoHigh = 1;
    int m_palette = 0;
    bool m_autoRange = true;
    std::string m_colorUnit = "raw";
    std::string m_loadedPath;
    std::string m_referencePath;
    int m_sliceField = -1;
    float m_sliceLow = 0, m_sliceHigh = 1;
    int m_pickerAction = 0;
    PathPicker m_picker;
    ElectronicViewport m_viewport;
    atomforge::electronic::Mesh m_surface;
    bool m_colorSurface = false;
    bool m_appendSurface = false;
    std::string m_error;
    atomforge::electronic::Volume m_volume;
    atomforge::electronic::Volume m_referenceVolume;
    Output m_result;
    atomforge::BackgroundTask<Output> m_task;
};
