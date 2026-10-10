#include "ui/MenuParity.h"

const std::vector<MenuLocation>& menuLocations()
{
    static const std::vector<MenuLocation> table = {
        // Build menu: builders with their own dialogs...
        {"build-bulk", "Build", "", "Bulk Crystal", true},
        {"build-sss", "Build", "", "Substitutional Solid Solution", true},
        {"build-stacking-fault", "Build", "", "Stacking Fault", true},
        {"build-gb", "Build", "", "CSL Grain Boundary", true},
        {"build-interface", "Build", "", "Interface Builder", true},
        {"build-nano", "Build", "", "Nanocrystal", true},
        {"build-custom", "Build", "", "Custom Structure", true},
        {"build-poly", "Build", "", "Polycrystal", true},
        {"build-amorphous", "Build", "", "Amorphous Structure", true},
        // ...and through the operation dialog.
        {"build-vacancy", "Build", "", "Vacancies", false},
        {"build-strain", "Build", "", "Strained Cell", false},
        {"build-primitive", "Build", "", "Primitive Cell", false},
        {"build-surface", "Build", "", "Surface Slab", false},
        {"build-sqs", "Build", "", "Special Quasirandom Structure (SQS)", false},
        {"build-nanowire", "Build", "", "Nanowire", false},
        {"build-core-shell", "Build", "", "Core-Shell Particle", false},
        // Edit menu: dialogs that do the same as these steps.
        {"add-atom", "Edit", "", "Edit Structure (add atom)", true},
        {"set-cell", "Edit", "", "Edit Structure (lattice vectors)", true},
        {"build-dislocation", "Edit", "", "Insert Dislocation", true},
        {"insert-interstitials", "Edit", "", "Add Interstitial Atoms", true},
        {"supercell", "Edit", "", "Transform Structure", true},
        {"merge", "Edit", "", "Merge Structures", true},
        {"build-sculpt", "Edit", "", "Cell Sculptor", true},
        // Edit > Structure Operations (operation dialog).
        {"delete-selected", "Edit", "Structure Operations", "Delete Atoms", false},
        {"assign-element", "Edit", "Structure Operations", "Assign Element", false},
        {"displace", "Edit", "Structure Operations", "Displace Atoms", false},
        {"random-displacement", "Edit", "Structure Operations", "Random Displacement", false},
        {"slice", "Edit", "Structure Operations", "Slice", false},
        {"replicate", "Edit", "Structure Operations", "Replicate", false},
        {"transform", "Edit", "Structure Operations", "Affine Transformation", false},
        {"strain", "Edit", "Structure Operations", "Strain", false},
        {"wrap", "Edit", "Structure Operations", "Wrap into Cell", false},
        {"center", "Edit", "Structure Operations", "Center in Cell", false},
        {"add-vacuum", "Edit", "Structure Operations", "Add Vacuum", false},
        {"compute-property", "Edit", "Structure Operations", "Compute Per-Atom Property", false},
        // Edit > Select Atoms (sets the selection in the view).
        {"select-element", "Edit", "Select Atoms", "By Element", false},
        {"select-expression", "Edit", "Select Atoms", "By Expression", false},
        {"select-slab", "Edit", "Select Atoms", "In a Slab", false},
        {"select-sphere", "Edit", "Select Atoms", "In a Sphere", false},
        {"select-property", "Edit", "Select Atoms", "By Property Range", false},
        {"select-random", "Edit", "Select Atoms", "Random Fraction", false},
        {"expand-selection", "Edit", "Select Atoms", "Expand to Neighbours", false},
        {"invert-selection", "Edit", "Select Atoms", "Invert Selection", false},
        {"clear-selection", "Edit", "Select Atoms", "Clear Selection", false},
        // Analysis menu: scientific tools whose result is a structure (the tool's window).
        {"relax", "Analysis", "Simulation", "Structure relaxation", true},
        {"nvt-dynamics", "Analysis", "Simulation", "NVT dynamics", true},
        {"npt-dynamics", "Analysis", "Simulation", "NPT dynamics", true},
        {"standardize-cell", "Analysis", "Reciprocal space", "Symmetry and Wyckoff positions", true},
    };
    return table;
}

std::string menuPathOf(const std::string& step)
{
    for (const auto& location : menuLocations())
        if (step == location.step)
            return std::string(location.menu) + (*location.submenu ? std::string(" > ") + location.submenu : "") + " > " + location.label;
    return "";
}
