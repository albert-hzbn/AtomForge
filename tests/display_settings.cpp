// Display settings: light directions, presets, and the saved text form
// (round trip, clamping, unknown keys and colours).
#include "model/DisplaySettings.h"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
using namespace atomforge::display;

void check(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

bool close(const glm::vec3& a, const glm::vec3& b, float tolerance = 1e-5f)
{
    return glm::length(a - b) <= tolerance;
}
}

int main()
{
    try {
        // The default light is along (1, 1, 1), as before the setting existed.
        const DisplaySettings defaults;
        check(close(lightDirection(defaults), glm::normalize(glm::vec3(1.0f))), "default light along (1, 1, 1)");
        DisplaySettings s;
        s.lightAzimuth = 0;
        s.lightElevation = 0;
        check(close(lightDirection(s), {1, 0, 0}), "azimuth 0 lights along +x");
        s.lightAzimuth = 90;
        check(close(lightDirection(s), {0, 0, 1}), "azimuth 90 lights along +z");
        s.lightElevation = 89;  // clamped to 85: the shadow map keeps a horizontal component
        check(std::abs(lightDirection(s).y - std::sin(85.0f * 3.14159265f / 180.0f)) < 1e-5f, "elevation clamped to 85");
        s.lightFollowsCamera = true;
        s.lightAzimuth = 0;
        s.lightElevation = 0;
        check(close(lightDirection(s), {0, 0, 1}), "a headlight at azimuth 0 comes from the viewer");
        s.lightAzimuth = 90;
        check(close(lightDirection(s), {1, 0, 0}), "a headlight at azimuth 90 comes from the right");

        // Presets: the defaults are "Default"; each preset is recognised after applying it.
        check(matchingPreset(defaults) == 0, "the defaults match the Default preset");
        for (std::size_t i = 0; i < presets().size(); ++i) {
            DisplaySettings look;
            look.atomScale = 1.4f;  // sizes are not part of a look
            applyPreset(look, presets()[i]);
            check(matchingPreset(look) == static_cast<int>(i), std::string("preset recognised: ") + presets()[i].name);
            check(look.atomScale == 1.4f, "a preset keeps the atom size");
        }
        DisplaySettings publication;
        applyPreset(publication, presets().back());
        check(publication.customBackground && publication.background == glm::vec3(1.0f) &&
              publication.customCellColour && publication.cellColour == glm::vec3(0.0f), "Publication: white background, black cell");
        applyPreset(publication, presets()[1]);
        check(!publication.customBackground && !publication.customCellColour, "other presets return to the theme's colours");
        DisplaySettings customised;
        customised.ambient = 0.41f;
        check(matchingPreset(customised) == -1, "a changed slider is a custom look");
        customised = DisplaySettings();
        customised.customBackground = true;
        check(matchingPreset(customised) == -1, "a custom background is not the Default look");

        // Text: every field round-trips exactly.
        DisplaySettings all;
        all.ambient = 0.123456f;
        all.shadowStrength = 0.9f;
        all.contrast = 2.1f;
        all.saturation = 0.3f;
        all.lightAzimuth = 271.5f;
        all.lightElevation = -20.25f;
        all.lightFollowsCamera = true;
        all.specularIntensity = 1.7f;
        all.shininessScale = 0.2f;
        all.shininessFloor = 200.0f;
        all.atomScale = 0.35f;
        all.bondScale = 2.5f;
        all.uniformBondColour = true;
        all.bondColour = {0.1f, 0.2f, 0.3f};
        all.customBackground = true;
        all.background = {0.0f, 0.0f, 0.0f};
        all.customCellColour = true;
        all.cellColour = {0.9f, 0.5f, 0.25f};
        check(fromText(toText(all)) == all, "every setting round-trips through text");
        check(fromText(toText(defaults)) == defaults, "the defaults round-trip");
        check(toText(defaults).find("background=theme") != std::string::npos, "the theme background is written as 'theme'");

        // Reading: clamping, wrapping, bad values, unknown keys, partial files.
        const DisplaySettings read = fromText(
            "# comment\n"
            "ambient=7\n"
            "atom_scale=-1\n"
            "light_azimuth=-30\n"
            "light_elevation=120\n"
            "contrast=nonsense\n"
            "unknown_key=3\n"
            "bond_colour=1 0 2\n"
            "background = 0.5 0.5\n"
            " cell_colour = theme \n");
        check(read.ambient == kAmbient.high, "ambient clamped");
        check(read.atomScale == kAtomScale.low, "atom size clamped");
        check(std::abs(read.lightAzimuth - 330.0f) < 1e-4f, "azimuth wrapped into 0-360");
        check(read.lightElevation == kElevation.high, "elevation clamped");
        check(read.contrast == defaults.contrast, "an unreadable number keeps the default");
        check(read.uniformBondColour && read.bondColour == glm::vec3(1, 0, 1), "colour channels clamped to 0-1");
        check(!read.customBackground, "an incomplete colour is ignored");
        check(!read.customCellColour, "'theme' keeps the theme's cell colour");
        DisplaySettings base;
        base.bondScale = 2.0f;
        check(fromText("ambient=0.5\n", base).bondScale == 2.0f, "missing keys keep the base value");
        check(fromText("ambient=0.5\n", base).ambient == 0.5f, "present keys are read");

        // Saving and loading a file; a missing file gives the defaults.
        const auto path = std::filesystem::temp_directory_path() / "atomforge_display_settings_test.txt";
        check(saveDefaults(all, path), "settings saved");
        check(loadDefaults(path) == all, "saved settings load back");
        std::filesystem::remove(path);
        check(loadDefaults(path) == DisplaySettings(), "a missing file gives the defaults");
    } catch (const std::exception& e) {
        std::cerr << "display_settings failed: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "display_settings: all checks passed" << std::endl;
    return 0;
}
