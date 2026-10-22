#include "model/DisplaySettings.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace atomforge::display
{
namespace
{
constexpr float kDegrees = 3.14159265358979f / 180.0f;

float clampTo(float value, Range range) { return std::clamp(value, range.low, range.high); }

bool near(float a, float b) { return std::abs(a - b) <= 1e-4f * std::max(1.0f, std::abs(b)); }

std::string colourText(const glm::vec3& c)
{
    std::ostringstream out;
    out.precision(9);
    out << c.r << ' ' << c.g << ' ' << c.b;
    return out.str();
}

// A colour "r g b" in 0-1; false when unreadable.
bool readColour(const std::string& text, glm::vec3& colour)
{
    std::istringstream in(text);
    glm::vec3 c;
    if (!(in >> c.r >> c.g >> c.b)) return false;
    for (int k = 0; k < 3; ++k)
        if (!std::isfinite(c[k])) return false;
    colour = glm::clamp(c, glm::vec3(0.0f), glm::vec3(1.0f));
    return true;
}

bool readNumber(const std::string& text, float& value, Range range)
{
    char* end = nullptr;
    const float number = std::strtof(text.c_str(), &end);
    if (end == text.c_str() || !std::isfinite(number)) return false;
    value = clampTo(number, range);
    return true;
}

bool readSwitch(const std::string& text, bool& value)
{
    if (text == "true" || text == "1" || text == "yes") { value = true; return true; }
    if (text == "false" || text == "0" || text == "no") { value = false; return true; }
    return false;
}

std::string trimmed(const std::string& text)
{
    const auto first = text.find_first_not_of(" \t\r");
    if (first == std::string::npos) return {};
    return text.substr(first, text.find_last_not_of(" \t\r") - first + 1);
}
}  // namespace

bool operator==(const DisplaySettings& a, const DisplaySettings& b)
{
    return a.ambient == b.ambient && a.shadowStrength == b.shadowStrength && a.contrast == b.contrast &&
           a.saturation == b.saturation && a.lightAzimuth == b.lightAzimuth && a.lightElevation == b.lightElevation &&
           a.lightFollowsCamera == b.lightFollowsCamera && a.specularIntensity == b.specularIntensity &&
           a.shininessScale == b.shininessScale && a.shininessFloor == b.shininessFloor && a.atomScale == b.atomScale &&
           a.bondScale == b.bondScale && a.uniformBondColour == b.uniformBondColour && a.bondColour == b.bondColour &&
           a.customBackground == b.customBackground && a.background == b.background &&
           a.customCellColour == b.customCellColour && a.cellColour == b.cellColour;
}

glm::vec3 lightDirection(const DisplaySettings& s)
{
    const float azimuth = s.lightAzimuth * kDegrees;
    const float elevation = clampTo(s.lightElevation, kElevation) * kDegrees;
    const float horizontal = std::cos(elevation);
    if (s.lightFollowsCamera)
        return {horizontal * std::sin(azimuth), std::sin(elevation), horizontal * std::cos(azimuth)};
    return {horizontal * std::cos(azimuth), std::sin(elevation), horizontal * std::sin(azimuth)};
}

const std::vector<Preset>& presets()
{
    static const std::vector<Preset> all = {
        {"Default", "Balanced light with soft shadows and moderate highlights.",
         0.40f, 0.20f, 1.25f, 1.55f, 0.65f, 1.5f, 32.0f, false},
        {"Soft", "Even light, faint shadows and gentle highlights.",
         0.55f, 0.10f, 1.05f, 1.25f, 0.35f, 1.0f, 16.0f, false},
        {"Glossy", "Polished spheres with small, bright highlights.",
         0.35f, 0.30f, 1.30f, 1.60f, 1.25f, 3.0f, 96.0f, false},
        {"Matte", "Chalky surfaces without highlights.",
         0.45f, 0.20f, 1.15f, 1.30f, 0.05f, 0.6f, 8.0f, false},
        {"Dramatic", "Low ambient light and deep shadows for strong depth.",
         0.15f, 0.75f, 1.60f, 1.70f, 0.90f, 2.0f, 48.0f, false},
        {"Flat", "No shadows or highlights: plain colours, as in a schematic.",
         0.75f, 0.0f, 1.0f, 1.30f, 0.0f, 1.0f, 32.0f, false},
        {"Publication", "White background, black cell and soft light, for figures.",
         0.50f, 0.15f, 1.10f, 1.35f, 0.50f, 1.5f, 32.0f, true},
    };
    return all;
}

void applyPreset(DisplaySettings& s, const Preset& p)
{
    s.ambient = p.ambient;
    s.shadowStrength = p.shadowStrength;
    s.contrast = p.contrast;
    s.saturation = p.saturation;
    s.specularIntensity = p.specularIntensity;
    s.shininessScale = p.shininessScale;
    s.shininessFloor = p.shininessFloor;
    s.customBackground = p.whiteBackground;
    s.background = glm::vec3(1.0f);
    s.customCellColour = p.whiteBackground;
    s.cellColour = glm::vec3(0.0f);
}

int matchingPreset(const DisplaySettings& s)
{
    const auto& all = presets();
    for (std::size_t i = 0; i < all.size(); ++i) {
        const Preset& p = all[i];
        const bool look = near(s.ambient, p.ambient) && near(s.shadowStrength, p.shadowStrength) &&
                          near(s.contrast, p.contrast) && near(s.saturation, p.saturation) &&
                          near(s.specularIntensity, p.specularIntensity) && near(s.shininessScale, p.shininessScale) &&
                          near(s.shininessFloor, p.shininessFloor);
        const bool colours = p.whiteBackground
            ? s.customBackground && s.background == glm::vec3(1.0f) && s.customCellColour && s.cellColour == glm::vec3(0.0f)
            : !s.customBackground && !s.customCellColour;
        if (look && colours) return static_cast<int>(i);
    }
    return -1;
}

std::string toText(const DisplaySettings& s)
{
    std::ostringstream out;
    out.precision(9);
    out << "ambient=" << s.ambient << '\n'
        << "shadow_strength=" << s.shadowStrength << '\n'
        << "contrast=" << s.contrast << '\n'
        << "saturation=" << s.saturation << '\n'
        << "light_azimuth=" << s.lightAzimuth << '\n'
        << "light_elevation=" << s.lightElevation << '\n'
        << "light_follows_camera=" << (s.lightFollowsCamera ? "true" : "false") << '\n'
        << "specular_intensity=" << s.specularIntensity << '\n'
        << "shininess_scale=" << s.shininessScale << '\n'
        << "shininess_floor=" << s.shininessFloor << '\n'
        << "atom_scale=" << s.atomScale << '\n'
        << "bond_scale=" << s.bondScale << '\n'
        << "bond_colour=" << (s.uniformBondColour ? colourText(s.bondColour) : "atoms") << '\n'
        << "background=" << (s.customBackground ? colourText(s.background) : "theme") << '\n'
        << "cell_colour=" << (s.customCellColour ? colourText(s.cellColour) : "theme") << '\n';
    return out.str();
}

DisplaySettings fromText(const std::string& text, const DisplaySettings& base)
{
    DisplaySettings s = base;
    std::istringstream in(text);
    for (std::string line; std::getline(in, line);) {
        const auto equals = line.find('=');
        if (equals == std::string::npos) continue;
        const std::string key = trimmed(line.substr(0, equals)), value = trimmed(line.substr(equals + 1));
        if (key == "ambient") readNumber(value, s.ambient, kAmbient);
        else if (key == "shadow_strength") readNumber(value, s.shadowStrength, kShadowStrength);
        else if (key == "contrast") readNumber(value, s.contrast, kContrast);
        else if (key == "saturation") readNumber(value, s.saturation, kSaturation);
        else if (key == "light_azimuth") {
            float azimuth = 0;
            if (readNumber(value, azimuth, {-1e6f, 1e6f})) s.lightAzimuth = std::fmod(std::fmod(azimuth, 360.0f) + 360.0f, 360.0f);
        }
        else if (key == "light_elevation") readNumber(value, s.lightElevation, kElevation);
        else if (key == "light_follows_camera") readSwitch(value, s.lightFollowsCamera);
        else if (key == "specular_intensity") readNumber(value, s.specularIntensity, kSpecularIntensity);
        else if (key == "shininess_scale") readNumber(value, s.shininessScale, kShininessScale);
        else if (key == "shininess_floor") readNumber(value, s.shininessFloor, kShininessFloor);
        else if (key == "atom_scale") readNumber(value, s.atomScale, kAtomScale);
        else if (key == "bond_scale") readNumber(value, s.bondScale, kBondScale);
        else if (key == "bond_colour") {
            if (value == "atoms") s.uniformBondColour = false;
            else if (readColour(value, s.bondColour)) s.uniformBondColour = true;
        }
        else if (key == "background") {
            if (value == "theme") s.customBackground = false;
            else if (readColour(value, s.background)) s.customBackground = true;
        }
        else if (key == "cell_colour") {
            if (value == "theme") s.customCellColour = false;
            else if (readColour(value, s.cellColour)) s.customCellColour = true;
        }
    }
    return s;
}

std::filesystem::path savedDefaultsPath()
{
#ifdef _WIN32
    const char* home = std::getenv("USERPROFILE");
#else
    const char* home = std::getenv("HOME");
#endif
    return std::filesystem::u8path(home && *home ? home : ".") / ".atomforge_display.txt";
}

bool saveDefaults(const DisplaySettings& settings, const std::filesystem::path& path)
{
    std::ofstream out(path, std::ios::trunc);
    if (!out) return false;
    out << "# AtomForge display settings (View > Display Settings > Save as My Default)\n" << toText(settings);
    return static_cast<bool>(out);
}

DisplaySettings loadDefaults(const std::filesystem::path& path)
{
    std::ifstream in(path);
    if (!in) return DisplaySettings();
    std::ostringstream text;
    text << in.rdbuf();
    return fromText(text.str());
}

}  // namespace atomforge::display
