#include "ui/ModifierEditor.h"
#include "ui/ResponsiveLayout.h"
#include "imgui.h"

#include <algorithm>
#include <cstdio>
#include <sstream>
#include <vector>

using atomforge::pipeline::Json;

void drawModifierParameters(const atomforge::pipeline::ModifierType& type, const Json& parameters,
                            const std::function<void(const std::string&, const Json&)>& changed)
{
    const float width = std::min(ImGui::GetContentRegionAvail().x * 0.55f, responsive::dp(320));
    for (const auto& p : type.parameters) {
        ImGui::PushID(p.name);
        const Json value = atomforge::pipeline::parameter(parameters, type, p.name);
        const std::string kind = p.kind;
        ImGui::SetNextItemWidth(width);
        if (kind == "bool") {
            bool v = value.boolean();
            if (ImGui::Checkbox(p.label, &v)) changed(p.name, Json(v));
        } else if (kind == "int") {
            int v = static_cast<int>(value.number());
            if (ImGui::InputInt(p.label, &v)) changed(p.name, Json(v));
        } else if (kind == "float") {
            double v = value.number();
            if (ImGui::InputDouble(p.label, &v, 0, 0, "%.6g")) changed(p.name, Json(v));
        } else if (kind == "vector" || kind == "ivec3") {
            double v[3] = {value.items()[0].number(), value.items()[1].number(), value.items()[2].number()};
            bool edited = false;
            if (kind == "ivec3") {
                int n[3] = {static_cast<int>(v[0]), static_cast<int>(v[1]), static_cast<int>(v[2])};
                edited = ImGui::InputInt3(p.label, n);
                for (int k = 0; k < 3; ++k) v[k] = n[k];
            } else edited = ImGui::InputScalarN(p.label, ImGuiDataType_Double, v, 3, nullptr, nullptr, "%.6g");
            if (edited) changed(p.name, Json::array({v[0], v[1], v[2]}));
        } else if (kind == "choice") {
            std::vector<std::string> options;
            std::stringstream stream(p.options);
            for (std::string option; std::getline(stream, option, '|');) options.push_back(option);
            const std::string current = value.string();
            if (ImGui::BeginCombo(p.label, current.c_str())) {
                for (const auto& option : options)
                    if (ImGui::Selectable(option.c_str(), option == current)) changed(p.name, Json(option));
                ImGui::EndCombo();
            }
        } else {
            // string, element(s), expression or command-line options
            char buffer[2048];
            std::snprintf(buffer, sizeof(buffer), "%s", value.string().c_str());
            // Expressions and option text get the full width, with their label above.
            const bool wide = kind == "expression" || kind == "options";
            if (wide) { ImGui::TextUnformatted(p.label); ImGui::SetNextItemWidth(-1); }
            if (ImGui::InputText(wide ? "##text" : p.label, buffer, sizeof(buffer))) changed(p.name, Json(std::string(buffer)));
        }
        ImGui::PopID();
    }
}
