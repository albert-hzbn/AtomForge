#include "app/EditorApplication.h"
#include "cli/BuildModifiers.h"
#include "cli/CLIMode.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <glob.h>
#endif

// ---------------------------------------------------------------------------
// Expand a single argument that may contain wildcard characters (* or ?) into
// a list of matching file paths.  Arguments without wildcards are returned
// unchanged (even if they don't exist yet) so that the normal missing-file
// error path still fires.
// ---------------------------------------------------------------------------
static void expandWildcard(const std::string& pattern, std::vector<std::string>& out)
{
    // If no wildcard characters are present just pass through unchanged.
    if (pattern.find('*') == std::string::npos &&
        pattern.find('?') == std::string::npos)
    {
        out.push_back(pattern);
        return;
    }

#ifdef _WIN32
    // Split the pattern into a directory prefix and the filename pattern.
    std::string dir;
    const std::size_t sep = pattern.find_last_of("/\\");
    if (sep != std::string::npos)
        dir = pattern.substr(0, sep + 1); // includes trailing separator

    WIN32_FIND_DATAA fd{};
    HANDLE h = FindFirstFileA(pattern.c_str(), &fd);
    if (h == INVALID_HANDLE_VALUE)
    {
        // No match – pass through so the caller can emit an appropriate error.
        out.push_back(pattern);
        return;
    }
    std::vector<std::string> matches;
    do {
        // Skip virtual directory entries.
        if (std::strcmp(fd.cFileName, ".") == 0 || std::strcmp(fd.cFileName, "..") == 0)
            continue;
        matches.push_back(dir + fd.cFileName);
    } while (FindNextFileA(h, &fd));
    FindClose(h);

    std::sort(matches.begin(), matches.end());
    for (auto& m : matches)
        out.push_back(std::move(m));
#else
    glob_t g{};
    const int rc = glob(pattern.c_str(), GLOB_TILDE, nullptr, &g);
    if (rc == 0)
    {
        for (std::size_t i = 0; i < g.gl_pathc; ++i)
            out.emplace_back(g.gl_pathv[i]);
    }
    else
    {
        // No match – pass through.
        out.push_back(pattern);
    }
    globfree(&g);
#endif
}

int main(int argc, char* argv[])
{
    // Build and Edit operations are also structure-pipeline steps (desktop and CLI).
    registerBuildModifiers();
    if (isCLIMode(argc, argv))
        return runCLI(argc, argv);

    std::vector<std::string> startupPaths;
    for (int i = 1; i < argc; ++i)
        expandWildcard(argv[i], startupPaths);
    return runAtomsEditor(startupPaths);
}

