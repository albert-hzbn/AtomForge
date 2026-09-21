#pragma once

#include <filesystem>
#include <string>
#include <vector>

// Execute an argument vector without a shell. Output is captured to a file;
// TaskControl cancellation terminates this process and its child processes.
int runScientificProcess(const std::vector<std::string>& arguments, const std::filesystem::path& logPath);
