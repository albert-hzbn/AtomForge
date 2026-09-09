#pragma once

#include <functional>
#include <vector>

#include <glm/glm.hpp>

#include "model/Structure.h"

struct BulkCrystalBuilderDialog
{
    BulkCrystalBuilderDialog();

    void drawMenuItem(bool enabled);
    void drawDialog(Structure& structure,
                    const std::vector<glm::vec3>& elementColors,
                    const std::function<void(Structure&)>& updateBuffers);

private:
    bool m_openRequested = false;
};