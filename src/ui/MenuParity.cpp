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
        // Build > Low-Dimensional and Build > Clusters and Molecules (operation dialog).
        {"build-nanotube", "Build", "Low-Dimensional", "Nanotube", false},
        {"build-monolayer", "Build", "Low-Dimensional", "2D Monolayer (Graphene, h-BN, MX2)", false},
        {"build-nanoribbon", "Build", "Low-Dimensional", "Nanoribbon", false},
        {"build-twisted-bilayer", "Build", "Low-Dimensional", "Twisted Bilayer (Moire)", false},
        {"build-cluster", "Build", "Clusters and Molecules", "Magic-Number Cluster", false},
        {"build-molecule", "Build", "Clusters and Molecules", "Molecule", false},
        {"build-packing", "Build", "Clusters and Molecules", "Molecular Packing (Liquid, Gas)", false},
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
        {"assign-color", "Edit", "Structure Operations", "Assign Colour", false},
        {"ambient-occlusion", "Edit", "Structure Operations", "Ambient Occlusion", false},
        {"clear-property", "Edit", "Structure Operations", "Remove Per-Atom Property", false},
        {"freeze-property", "Edit", "Structure Operations", "Freeze Per-Atom Property", false},
        {"python-script", "Edit", "Structure Operations", "Python Script", false},
        {"displace", "Edit", "Structure Operations", "Displace Atoms", false},
        {"fix-atoms", "Edit", "Structure Operations", "Fix Atoms", false},
        {"free-atoms", "Edit", "Structure Operations", "Free Atoms", false},
        {"random-displacement", "Edit", "Structure Operations", "Random Displacement", false},
        {"slice", "Edit", "Structure Operations", "Slice", false},
        {"replicate", "Edit", "Structure Operations", "Replicate", false},
        {"transform", "Edit", "Structure Operations", "Affine Transformation", false},
        {"strain", "Edit", "Structure Operations", "Strain", false},
        {"wrap", "Edit", "Structure Operations", "Wrap into Cell", false},
        {"set-pbc", "Edit", "Structure Operations", "Periodic Boundaries", false},
        {"center", "Edit", "Structure Operations", "Center in Cell", false},
        {"add-vacuum", "Edit", "Structure Operations", "Add Vacuum", false},
        {"compute-property", "Edit", "Structure Operations", "Compute Per-Atom Property", false},
        {"compute-expression", "Edit", "Structure Operations", "Compute Property from Formula", false},
        {"segment-grains", "Edit", "Structure Operations", "Segment Grains", false},
        {"elastic-strain", "Edit", "Structure Operations", "Elastic Strain", false},
        {"displacements", "Edit", "Structure Operations", "Displacement Vectors", false},
        {"bin-property", "Edit", "Structure Operations", "Spatial Binning", false},
        {"reduce-property", "Edit", "Structure Operations", "Property Statistics", false},
        {"shrink-wrap", "Edit", "Structure Operations", "Shrink-Wrap Cell", false},
        // Edit > Geometry, Edit > Cell, Edit > Chemistry and Defects (operation dialog).
        {"mirror", "Edit", "Geometry", "Mirror", false},
        {"invert", "Edit", "Geometry", "Invert through a Point", false},
        {"rotate", "Edit", "Geometry", "Rotate", false},
        {"set-geometry", "Edit", "Geometry", "Set Distance, Angle or Dihedral", false},
        {"measure-geometry", "Edit", "Geometry", "Measure Distance, Angle or Dihedral", false},
        {"duplicate-selected", "Edit", "Geometry", "Duplicate Selection", false},
        {"interpolate", "Edit", "Geometry", "Interpolate to a Structure", false},
        {"unwrap-molecules", "Edit", "Geometry", "Make Molecules Whole", false},
        {"align-direction", "Edit", "Cell", "Align a Direction with an Axis", false},
        {"near-cubic-supercell", "Edit", "Cell", "Near-Cubic Supercell", false},
        {"niggli-reduce", "Edit", "Cell", "Niggli-Reduced Cell", false},
        {"scale-cell", "Edit", "Cell", "Scale Cell (Volume, Density)", false},
        {"sort-atoms", "Edit", "Cell", "Sort Atoms", false},
        {"add-adsorbate", "Edit", "Chemistry and Defects", "Add Adsorbate", false},
        {"saturate-hydrogen", "Edit", "Chemistry and Defects", "Saturate with Hydrogen", false},
        {"swap-elements", "Edit", "Chemistry and Defects", "Swap Elements", false},
        {"create-antisite", "Edit", "Chemistry and Defects", "Antisite Defect", false},
        {"create-frenkel-pair", "Edit", "Chemistry and Defects", "Frenkel Pair", false},
        // Edit > Select Atoms (sets the selection in the view).
        {"select-element", "Edit", "Select Atoms", "By Element", false},
        {"select-expression", "Edit", "Select Atoms", "By Expression", false},
        {"select-slab", "Edit", "Select Atoms", "In a Slab", false},
        {"select-sphere", "Edit", "Select Atoms", "In a Sphere", false},
        {"select-property", "Edit", "Select Atoms", "By Property Range", false},
        {"select-random", "Edit", "Select Atoms", "Random Fraction", false},
        {"select-overlapping", "Edit", "Select Atoms", "Overlapping Atoms", false},
        {"select-surface", "Edit", "Select Atoms", "Surface Atoms", false},
        {"select-region", "Edit", "Select Atoms", "In a LAMMPS Region", false},
        {"expand-selection", "Edit", "Select Atoms", "Expand to Neighbours", false},
        {"invert-selection", "Edit", "Select Atoms", "Invert Selection", false},
        {"clear-selection", "Edit", "Select Atoms", "Clear Selection", false},
        {"select-fragment", "Edit", "Select Atoms", "Bonded Fragments", false},
        // Analysis menu: scientific tools whose result is a structure (the tool's window).
        {"relax", "Analysis", "Simulation", "Structure Relaxation", true},
        {"nvt-dynamics", "Analysis", "Simulation", "NVT Dynamics", true},
        {"npt-dynamics", "Analysis", "Simulation", "NPT Dynamics", true},
        {"standardize-cell", "Analysis", "Reciprocal Space", "Symmetry and Wyckoff Positions", true},
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
