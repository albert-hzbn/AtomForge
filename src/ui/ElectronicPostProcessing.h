#pragma once

#include "electronic/Volume.h"
#include "util/BackgroundTask.h"
#include <string>

struct ElectronicPostProcessingDialog
{
    void drawMenuItem();
    void drawDialog();

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
    };
    void drawPreview();
    bool m_open = false;
    char m_input[2048]{};
    char m_output[2048] = "electronic-result.xsf";
    char m_reflections[8192] = "0 0 0\n1 0 0\n-1 0 0";
    char m_charges[8192]{};
    int m_quantity = 0;
    int m_cubeUnits = 0;
    int m_selected = 0;
    int m_reference = 0;
    int m_operation = 0;
    int m_axis = 2;
    int m_count = 100;
    int m_window = 3;
    int m_radius = 2;
    int m_exportFormat = 0;
    float m_scalar = 0.1f;
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
    float m_opacity = 0.75f;
    bool m_colorSurface = false;
    bool m_appendSurface = false;
    std::string m_error;
    atomforge::electronic::Volume m_volume;
    atomforge::electronic::Volume m_referenceVolume;
    Output m_result;
    atomforge::BackgroundTask<Output> m_task;
};
