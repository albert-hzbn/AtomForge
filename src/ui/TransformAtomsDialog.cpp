#include "TransformAtomsDialog.h"

#include "imgui.h"
#include "ui/DialogLayout.h"

namespace
{
constexpr int kMatrixSize = 3;
const char* kTransformPopupTitle = "Transform Structure";

void setIdentityMatrix(int (&matrix)[kMatrixSize][kMatrixSize])
{
    for (int row = 0; row < kMatrixSize; ++row)
    {
        for (int col = 0; col < kMatrixSize; ++col)
            matrix[row][col] = (row == col) ? 1 : 0;
    }
}

void copyMatrix(const int (&source)[kMatrixSize][kMatrixSize],
                int (&target)[kMatrixSize][kMatrixSize])
{
    for (int row = 0; row < kMatrixSize; ++row)
    {
        for (int col = 0; col < kMatrixSize; ++col)
            target[row][col] = source[row][col];
    }
}
} // namespace

void TransformAtomsDialog::clearTransform()
{
    useTransformMatrix = false;
    setIdentityMatrix(transformMatrix);
    setIdentityMatrix(pendingMatrix);
}

void TransformAtomsDialog::drawMenuItem(bool hasUnitCell)
{
    if (ImGui::MenuItem("Transform Structure", NULL, false, hasUnitCell))
        showDialog = true;
}

void TransformAtomsDialog::drawDialog(const std::function<void()>& onApply)
{
    if (showDialog)
    {
        ImGui::OpenPopup(kTransformPopupTitle);
        showDialog = false;
    }

    bool transformAtomsOpen = true;
    if (ImGui::BeginPopupModal(kTransformPopupTitle, &transformAtomsOpen, ImGuiWindowFlags_AlwaysAutoResize))
    {
        dialogLayout::section("Transformation matrix");
        ImGui::TextUnformatted("Enter three rows of integer coefficients.");

        for (int row = 0; row < kMatrixSize; ++row)
        {
            ImGui::PushID(row);
            ImGui::InputInt3("", pendingMatrix[row]);
            ImGui::PopID();
        }

        ImGui::Spacing();
        ImGui::Separator();
        if (dialogLayout::primaryButton("Apply",dialogLayout::actionSize()))
        {
            copyMatrix(pendingMatrix, transformMatrix);

            useTransformMatrix = true;
            ImGui::CloseCurrentPopup();
            onApply();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel",dialogLayout::actionSize())) ImGui::CloseCurrentPopup();
        ImGui::EndPopup();
    }
    if (!transformAtomsOpen)
        ImGui::CloseCurrentPopup();
}
