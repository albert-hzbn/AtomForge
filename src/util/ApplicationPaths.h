#pragma once

#include <filesystem>

// Resolve assets relative to the executable, not the launch working directory.
std::filesystem::path applicationDirectory();
std::filesystem::path findManualPdf(const std::filesystem::path& executableDirectory);
// The example structures (installed next to the manual; the source tree's
// docs/manual/examples in development builds). Empty when there are none.
std::filesystem::path findExamplesDirectory(const std::filesystem::path& executableDirectory);
