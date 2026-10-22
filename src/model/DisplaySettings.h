#pragma once

// How the 3D view draws a structure: lighting, surface finish, atom and bond
// sizes, bond colouring, background and cell colours. One per window tab; the
// Display Settings dialog edits it and the renderers read it each frame.

#include <glm/glm.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace atomforge::display
{

struct DisplaySettings
{
    // Light: brightness of the shadow side, darkness of shadows, colour.
    float ambient = 0.40f;
    float shadowStrength = 0.20f;
    float contrast = 1.25f;
    float saturation = 1.55f;

    // Where the light comes from: azimuth around the vertical axis and
    // elevation above the horizontal, in degrees. In the scene's frame, or
    // relative to the viewer when the light follows the camera.
    float lightAzimuth = 45.0f;
    float lightElevation = 35.26439f;
    bool lightFollowsCamera = false;

    // Surface finish: highlight brightness and sharpness.
    float specularIntensity = 0.65f;
    float shininessScale = 1.5f;
    float shininessFloor = 32.0f;

    // Sizes, multiplying the element radii and the automatic bond radius.
    float atomScale = 1.0f;
    float bondScale = 1.0f;

    // Bonds take half the colour of each atom unless one colour is chosen.
    bool uniformBondColour = false;
    glm::vec3 bondColour = glm::vec3(0.62f);

    // Background and cell outline: the theme's colours unless chosen.
    bool customBackground = false;
    glm::vec3 background = glm::vec3(1.0f);
    bool customCellColour = false;
    glm::vec3 cellColour = glm::vec3(0.15f);

};

bool operator==(const DisplaySettings& a, const DisplaySettings& b);
inline bool operator!=(const DisplaySettings& a, const DisplaySettings& b) { return !(a == b); }

// The allowed range of every number (values read from files are clamped to it).
struct Range { float low, high; };
inline constexpr Range kAmbient{0.0f, 1.0f}, kShadowStrength{0.0f, 1.0f}, kContrast{0.5f, 2.5f},
    kSaturation{0.0f, 3.0f}, kAzimuth{0.0f, 360.0f}, kElevation{-85.0f, 85.0f},
    kSpecularIntensity{0.0f, 2.0f}, kShininessScale{0.1f, 5.0f}, kShininessFloor{1.0f, 256.0f},
    kAtomScale{0.2f, 3.0f}, kBondScale{0.2f, 4.0f};

// Unit vector towards the light. In the scene's frame (y up) for a fixed
// light; for one following the camera, in the camera's frame: x right, y up
// and z towards the viewer (azimuth 0, elevation 0 lights from the viewer).
glm::vec3 lightDirection(const DisplaySettings& settings);

// A named look: lighting, surface finish, background and cell colours.
// Sizes, bond colours and the light direction are left as they are.
struct Preset
{
    const char* name;
    const char* description;
    float ambient, shadowStrength, contrast, saturation;
    float specularIntensity, shininessScale, shininessFloor;
    bool whiteBackground;  // a white background with a black cell (else the theme's)
};
const std::vector<Preset>& presets();
void applyPreset(DisplaySettings& settings, const Preset& preset);
// The preset the settings currently match, or -1 for a customised look.
int matchingPreset(const DisplaySettings& settings);

// Text form, one "key=value" per line; colours are "r g b" (0-1) or "theme"
// (background and cell) / "atoms" (bonds). Reading keeps `base` for missing
// or unreadable keys and clamps numbers to their ranges.
std::string toText(const DisplaySettings& settings);
DisplaySettings fromText(const std::string& text, const DisplaySettings& base = DisplaySettings());

// The user's saved default (~/.atomforge_display.txt), used for new windows.
std::filesystem::path savedDefaultsPath();
bool saveDefaults(const DisplaySettings& settings, const std::filesystem::path& path = savedDefaultsPath());
DisplaySettings loadDefaults(const std::filesystem::path& path = savedDefaultsPath());

}  // namespace atomforge::display
