#include "algorithms/NanoCrystalBuilder.h"
#include <iostream>

int main()
{
    Structure reference;
    reference.hasUnitCell = true;
    reference.cellVectors = {{{3.61,0,0},{0,3.61,0},{0,0,3.61}}};
    for (const auto& p : {glm::dvec3(0,0,0), glm::dvec3(0,1.805,1.805),
                          glm::dvec3(1.805,0,1.805), glm::dvec3(1.805,1.805,0)})
        reference.atoms.push_back({"Cu",29,p.x,p.y,p.z});
    for (const auto mode : {NanoGenerationMode::Shape, NanoGenerationMode::WulffConstruction})
    {
#ifndef ATOMS_ENABLE_SPGLIB
        if (mode == NanoGenerationMode::WulffConstruction) continue;
#endif
        Structure output = reference;
        output.dislocationLoopPoints = {{1,2,3},{3,2,1},{1,1,1}};
        output.dislocationDetectionDone = true;
        NanoParams params;
        params.generationMode = mode;
        params.sphereRadius = 5;
        params.wulffMaxRadius = 5;
        params.wulffPlanes = {{1,0,0,1}, {1,1,1,1}};
        const auto result = buildNanocrystal(output, reference, params, {}, {}, {});
        if (!result.success || !output.dislocationLoopPoints.empty())
        {
            std::cerr << "Nanocrystal retained a stale dislocation overlay: " << result.message << '\n';
            return 1;
        }
    }
}
