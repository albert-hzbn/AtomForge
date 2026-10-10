#pragma once

#include <functional>
#include <string>
#include <vector>

#include <glm/glm.hpp>

#include "algorithms/BulkCrystalBuilder.h"
#include "model/Structure.h"
#include "ui/StepEditing.h"

struct BulkCrystalBuilderDialog
{
    BulkCrystalBuilderDialog();

    void drawMenuItem(bool enabled);
    void drawDialog(Structure& structure,
                    const std::vector<glm::vec3>& elementColors,
                    const std::function<void(Structure&)>& updateBuffers);

    // Opens the dialog on a pipeline step: filled from the step's settings, with
    // "Update step" writing them back instead of building.
    bool editStep(StepEdit edit);
    // The dialog's current settings as the step's options (command-line flags of
    // AtomForge --build bulk).
    std::string stepOptions() const;

private:
    bool m_openRequested = false;
    StepEdit m_step;

    // Builder settings.
    int m_crystalSystemIndex = (int)CrystalSystem::Cubic;
    int m_selectedSpaceGroup = 225;
    int m_lastCrystalSystemIndex = (int)CrystalSystem::Cubic;
    LatticeParameters m_latticeParams;
    std::vector<AtomSite> m_asymmetricAtoms;
    // Atoms read from a step get their colours once the element colours are known.
    bool m_recolorAtoms = false;
};
