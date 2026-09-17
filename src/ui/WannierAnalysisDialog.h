#pragma once

#include "electronic/Wannier.h"
#include "ui/PathPicker.h"
#include "util/BackgroundTask.h"

#include <string>
#include <vector>

// Bands, Berry curvature and Chern numbers from a Wannier90 seedname_hr.dat
// file. Unlike ElectronicPostProcessingDialog this has no Grid/3D-viewport
// concept, so it is a lightweight modal analysis dialog (RDF/ADF pattern):
// compute and plot, no export or workspace persistence.
struct WannierAnalysisDialog
{
    void drawMenuItem();
    void drawDialog();

private:
    struct BandsResult { std::vector<double> energies; int numWann = 0; int numK = 0; };
    void drawBandsPlot(const BandsResult& result);
    void load(const std::string& path);

    bool m_openRequested = false;
    bool m_open = false;
    PathPicker m_picker;
    std::string m_loadedPath;
    std::string m_loadError;
    atomforge::electronic::WannierHamiltonian m_model;  // numWann==0: nothing loaded

    char m_kpath[4096] = "0 0 0\n0.5 0 0\n0.5 0.5 0\n0 0 0";
    int m_samplesPerSegment = 40;
    atomforge::BackgroundTask<BandsResult> m_bandsTask;
    BandsResult m_bands;

    float m_curvatureK[3]{};
    int m_curvaturePlane = 0;  // 0: kx-ky, 1: kx-kz, 2: ky-kz
    atomforge::BackgroundTask<std::vector<double>> m_curvatureTask;
    std::vector<double> m_curvature;

    int m_chernBand = 0;
    int m_chernPlane = 0;
    int m_chernGrid = 30;
    float m_chernFixed = 0;
    atomforge::BackgroundTask<double> m_chernTask;
    bool m_chernComputed = false;
    double m_chernNumber = 0;
};
