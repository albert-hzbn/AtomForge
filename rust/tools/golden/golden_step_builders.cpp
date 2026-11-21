// The C++ Nanocrystal, Polycrystal and Amorphous builder dialogs' editStep()
// followed by stepOptions() for a list of step option texts, as the golden
// rust/tests/golden/step_builders.tsv (step, options, stepOptions per line).
// Linked against the desktop app's object files (see build.sh).
#include "ui/AmorphousBuilderDialog.h"
#include "ui/NanoCrystalBuilderDialog.h"
#include "ui/PolyCrystalBuilderDialog.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace
{
StepEdit edit(const std::string& step, const std::string& options)
{
    StepEdit e;
    e.step = step;
    e.parameters["options"] = options;
    e.commit = [](const atomforge::pipeline::Json&) {};
    return e;
}
}  // namespace

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "usage: golden_step_builders CASES.tsv OUT.tsv\n";
        return 2;
    }
    std::ifstream in(argv[1]);
    std::ofstream out(argv[2], std::ios::binary);
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        const auto tab = line.find('\t');
        const std::string step = line.substr(0, tab), options = tab == std::string::npos ? "" : line.substr(tab + 1);
        std::string result;
        if (step == "build-nano") {
            NanoCrystalBuilderDialog d;
            d.editStep(edit(step, options));
            result = d.stepOptions();
        } else if (step == "build-poly") {
            PolyCrystalBuilderDialog d;
            d.editStep(edit(step, options));
            result = d.stepOptions();
        } else if (step == "build-amorphous") {
            AmorphousBuilderDialog d;
            d.editStep(edit(step, options));
            result = d.stepOptions();
        } else {
            std::cerr << "unknown step " << step << "\n";
            return 1;
        }
        out << step << '\t' << options << '\t' << result << '\n';
    }
    return 0;
}
