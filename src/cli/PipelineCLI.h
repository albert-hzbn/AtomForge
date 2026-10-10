#pragma once

// AtomForge --pipe "replicate 2 2 1 | select-element O | delete-selected" --input IN --output OUT
// AtomForge --pipeline STEPS.json|STEPS.txt --input IN --output OUT
// "-" reads or writes extended XYZ on stdin/stdout, so runs chain with shell pipes.
int runPipelineCLI(int argc, char* argv[]);
