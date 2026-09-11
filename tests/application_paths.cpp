#include "util/ApplicationPaths.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

int main()
{
    namespace fs = std::filesystem;
    const auto root = fs::temp_directory_path() / ("atomforge manual paths " +
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(root);
    try
    {
        for (const auto& layout : {"portable", "installed", "development"})
        {
            const auto base = root / layout;
            const auto executable = base / "bin";
            fs::create_directories(executable);
            if (!findManualPdf(executable).empty()) throw std::runtime_error("unexpected manual");
            const auto manual = std::string(layout) == "portable" ? executable / "AtomForge-manual.pdf"
                : std::string(layout) == "installed" ? base / "share/doc/atomforge/AtomForge-manual.pdf"
                : base / "docs/manual/AtomForge-manual.pdf";
            fs::create_directories(manual.parent_path());
            std::ofstream(manual) << "%PDF-1.4";
            if (findManualPdf(executable) != manual) throw std::runtime_error("manual layout not resolved");
        }
        if (applicationDirectory().empty()) throw std::runtime_error("missing executable directory");
        fs::remove_all(root);
    }
    catch (const std::exception& error)
    {
        fs::remove_all(root);
        std::cerr << error.what() << '\n';
        return 1;
    }
}
