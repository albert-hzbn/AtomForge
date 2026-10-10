#pragma once

// AtomForge --cloud build --input DUMP|XYZ|FILE --output BIG.afcloud
// AtomForge --cloud generate --lattice fcc --a 3.615 --cells 630 630 630 --element Cu --output BIG.afcloud
// AtomForge --cloud info BIG.afcloud
// Out-of-core atom clouds for the desktop's large-data viewer (billions of atoms).
int runCloudCLI(int argc, char* argv[]);
