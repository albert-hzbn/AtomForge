#include "util/ScientificProcess.h"
#include "util/TaskControl.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <thread>

int main(int argc, char** argv)
{
    if (argc != 2) return 2;
    const auto folder = std::filesystem::temp_directory_path() /
        ("AtomForge process test " + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        std::filesystem::create_directories(folder);
        const auto script = folder / "argument test.py", log = folder / "result.log";
        { std::ofstream output(script); output << "import sys\nfor arg in sys.argv[1:]: print(arg)\n"; }
        const std::vector<std::string> data = {"two words", "quote\"inside", "trailing\\", "$HOME;echo nope", ""};
        std::vector<std::string> arguments{argv[1], script.u8string()};
        arguments.insert(arguments.end(), data.begin(), data.end());
        if (runScientificProcess(arguments, log) != 0) throw std::runtime_error("Argument process failed");
        std::ifstream input(log);
        for (const auto& expected : data) {
            std::string actual;
            if (!std::getline(input, actual)) throw std::runtime_error("Missing argument");
            if (!actual.empty() && actual.back() == '\r') actual.pop_back();
            if (actual != expected) throw std::runtime_error("Arguments were not preserved");
        }
        input.close();
        if (runScientificProcess({argv[1], "-c", "raise SystemExit(17)"}, log) != 17)
            throw std::runtime_error("Exit status was lost");
        atomforge::TaskControl control;
        auto task = std::async(std::launch::async, [&] {
            atomforge::TaskControlScope scope(&control);
            return runScientificProcess({argv[1], "-c", "import time; time.sleep(30)"}, log);
        });
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        const auto start = std::chrono::steady_clock::now();
        control.cancelled = true;
        bool cancelled = false;
        try { task.get(); } catch (const std::runtime_error&) { cancelled = true; }
        if (!cancelled || std::chrono::steady_clock::now() - start > std::chrono::seconds(5))
            throw std::runtime_error("Process cancellation failed");
        std::filesystem::remove_all(folder);
        std::cout << "Scientific process argument, exit-status and cancellation tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
