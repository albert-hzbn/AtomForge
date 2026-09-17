#pragma once

#include "electronic/Lobster.h"
#include "ui/PathPicker.h"

#include <string>

// COHP/COOP/COBI bonding curves and their per-bond Fermi-level ICOHP/ICOOP/
// ICOBI summary, from LOBSTER output files. Like WannierAnalysisDialog this
// has no Grid/3D-viewport concept, so it is a lightweight modal analysis
// dialog (RDF/ADF pattern): load, browse, plot -- no export or workspace
// persistence. File parsing is fast enough to run synchronously (no
// BackgroundTask), matching WannierAnalysisDialog::load().
struct LobsterAnalysisDialog
{
    void drawMenuItem();
    void drawDialog();

private:
    void loadCohpcar(const std::string& path);
    void loadIcohplist(const std::string& path);

    bool m_openRequested = false;
    bool m_open = false;

    PathPicker m_cohpcarPicker;
    std::string m_cohpcarPath;
    std::string m_cohpcarError;
    atomforge::electronic::CohpData m_cohpcar;  // energies.empty(): nothing loaded
    int m_selectedBond = -1;  // -1: the "average" entry
    int m_spinChannel = 0;

    PathPicker m_icohplistPicker;
    std::string m_icohplistPath;
    std::string m_icohplistError;
    atomforge::electronic::IcohpList m_icohplist;
};
