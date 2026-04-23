#include "app/EditorApplication.h"

#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[])
{
    std::vector<std::string> startupPaths;
    for (int i = 1; i < argc; ++i)
        startupPaths.emplace_back(argv[i]);
    return runAtomsEditor(startupPaths);
}

