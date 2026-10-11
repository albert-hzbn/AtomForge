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
        for (const auto& layout : {"portable", "portable-bin", "installed", "development"})
        {
            const auto base = root / layout;
            const auto executable = base / "bin";
            fs::create_directories(executable);
            if (!findManualPdf(executable).empty()) throw std::runtime_error("unexpected manual");
            const auto manual = std::string(layout) == "portable" ? executable / "AtomForge-manual.pdf"
                : std::string(layout) == "portable-bin" ? base / "AtomForge-manual.pdf"
                : std::string(layout) == "installed" ? base / "share/doc/atomforge/AtomForge-manual.pdf"
                : base / "docs/manual/AtomForge-manual.pdf";
            fs::create_directories(manual.parent_path());
            std::ofstream(manual) << "%PDF-1.4";
            if (findManualPdf(executable) != manual) throw std::runtime_error("manual layout not resolved");
            if (std::string(layout) == "development")
            {
                // A stale copy beside the executable must not shadow the source manual.
                std::ofstream(executable / "AtomForge-manual.pdf") << "old manual";
                if (findManualPdf(executable) != manual)
                    throw std::runtime_error("development manual shadowed by stale copy");
                const auto configuration = executable / "Release";
                fs::create_directories(configuration);
                if (findManualPdf(configuration) != manual)
                    throw std::runtime_error("multi-configuration manual not resolved");
            }
        }
        // Example structures in each package layout.
        for (const auto& layout : {"portable", "portable-bin", "installed"})
        {
            const auto base = root / (std::string("examples-") + layout);
            const auto executable = base / "bin";
            fs::create_directories(executable);
            const auto examples = std::string(layout) == "portable" ? executable / "examples"
                : std::string(layout) == "portable-bin" ? base / "examples"
                : base / "share/doc/atomforge/examples";
            fs::create_directories(examples);
            if (findExamplesDirectory(executable) != examples) throw std::runtime_error("examples layout not resolved");
        }
        if (findExamplesDirectory(applicationDirectory()) != fs::path(ATOMFORGE_SOURCE_MANUAL).parent_path() / "examples")
            throw std::runtime_error("built application cannot find the source examples");
        if (applicationDirectory().empty()) throw std::runtime_error("missing executable directory");
        if (findManualPdf(applicationDirectory()) != fs::path(ATOMFORGE_SOURCE_MANUAL))
            throw std::runtime_error("built application cannot find its published manual");
        fs::remove_all(root);
    }
    catch (const std::exception& error)
    {
        fs::remove_all(root);
        std::cerr << error.what() << '\n';
        return 1;
    }
}
