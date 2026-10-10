#include "ui/ResponsiveLayout.h"
#include "TransformAtomsDialog.h"

#include "imgui.h"
#include "ui/DialogLayout.h"

#include <sstream>
#include <string>

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
    {
        m_step.finish();
        showDialog = true;
    }
}

bool TransformAtomsDialog::editStep(StepEdit edit)
{
    m_step = std::move(edit);
    setIdentityMatrix(pendingMatrix);
    // "matrix" holds nine integers row by row; a missing or short matrix keeps the identity.
    const atomforge::pipeline::Json* matrix = m_step.parameters.find("matrix");
    if (matrix && matrix->isString())
    {
        std::istringstream in(matrix->string());
        int values[kMatrixSize * kMatrixSize];
        int count = 0;
        while (count < kMatrixSize * kMatrixSize && in >> values[count]) ++count;
        if (count == kMatrixSize * kMatrixSize)
            for (int k = 0; k < count; ++k) pendingMatrix[k / kMatrixSize][k % kMatrixSize] = values[k];
    }
    showDialog = true;
    return true;
}

atomforge::pipeline::Json TransformAtomsDialog::stepParameters() const
{
    atomforge::pipeline::Json parameters = m_step.parameters;
    std::string text;
    for (int row = 0; row < kMatrixSize; ++row)
        for (int col = 0; col < kMatrixSize; ++col)
            text += (text.empty() ? "" : " ") + std::to_string(pendingMatrix[row][col]);
    parameters["matrix"] = text;
    return parameters;
}

void TransformAtomsDialog::drawDialog(const std::function<void()>& onApply)
{
    if (showDialog)
    {
        ImGui::OpenPopup(kTransformPopupTitle);
        showDialog = false;
    }

    bool transformAtomsOpen = true;
    if (responsive::beginModal(kTransformPopupTitle, &transformAtomsOpen, ImGuiWindowFlags_AlwaysAutoResize))
    {
        if (m_step.active())
            ImGui::TextColored(ImVec4(0.45f, 0.75f, 1.0f, 1.0f), "Editing pipeline step %s: Update step writes these settings into the pipeline.", m_step.step.c_str());
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
        if (m_step.active())
        {
            if (dialogLayout::primaryButton("Update step", dialogLayout::actionSize()))
            {
                m_step.commitParameters(stepParameters());
                ImGui::CloseCurrentPopup();
                m_step.finish();
            }
        }
        else if (dialogLayout::primaryButton("Apply",dialogLayout::actionSize()))
        {
            copyMatrix(pendingMatrix, transformMatrix);

            useTransformMatrix = true;
            ImGui::CloseCurrentPopup();
            onApply();
        }
        ImGui::SameLine();
        if (responsive::button("Cancel",dialogLayout::actionSize()))
        {
            ImGui::CloseCurrentPopup();
            m_step.finish();
        }
        ImGui::EndPopup();
    }
    if (!transformAtomsOpen)
    {
        ImGui::CloseCurrentPopup();
        m_step.finish();
    }
}
