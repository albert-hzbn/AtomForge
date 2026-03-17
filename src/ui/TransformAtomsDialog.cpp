#include "TransformAtomsDialog.h"

#include "imgui.h"

void TransformAtomsDialog::clearTransform()
{
    useTransformMatrix = false;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            int value = (i == j) ? 1 : 0;
            transformMatrix[i][j] = value;
            pendingMatrix[i][j] = value;
        }
    }
}

void TransformAtomsDialog::drawMenuItem(bool hasUnitCell)
{
    if (ImGui::MenuItem("Transform Structure...", NULL, false, hasUnitCell))
        showDialog = true;
}

void TransformAtomsDialog::drawDialog(const std::function<void()>& onApply)
{
    if (showDialog)
    {
        ImGui::OpenPopup("Transform Structure");
        showDialog = false;
    }

    bool transformAtomsOpen = true;
    if (ImGui::BeginPopupModal("Transform Structure", &transformAtomsOpen, ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::Text("Enter 3x3 integer transformation matrix:");

        for (int i = 0; i < 3; ++i)
        {
            ImGui::PushID(i);
            ImGui::InputInt3("", pendingMatrix[i]);
            ImGui::PopID();
        }

        if (ImGui::Button("Apply"))
        {
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    transformMatrix[i][j] = pendingMatrix[i][j];

            useTransformMatrix = true;
            ImGui::CloseCurrentPopup();
            onApply();
        }
        ImGui::EndPopup();
    }
    if (!transformAtomsOpen)
        ImGui::CloseCurrentPopup();
}
