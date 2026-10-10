#pragma once

// Registers the Build and Edit operations of AtomForge (every `--build` mode
// and the Cell Sculptor) as structure-pipeline modifiers, e.g.
//     build-bulk | build-vacancy --count 4 | replicate 2 2 2
// Each step runs the same code as the menus and the command line: modes that
// take a structure receive the pipeline's current structure as their input.
// Safe to call more than once.
void registerBuildModifiers();
