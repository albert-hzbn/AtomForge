#pragma once

#include <filesystem>

// Resolve assets relative to the executable, not the launch working directory.
std::filesystem::path applicationDirectory();
std::filesystem::path findManualPdf(const std::filesystem::path& executableDirectory);
