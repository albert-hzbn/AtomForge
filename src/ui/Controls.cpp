#include "Controls.h"
#include "imgui.h"

void drawControls(Camera& camera, std::function<void()> onDeleteSelected)
{
    ImGui::Begin("Controls");

    ImGui::SliderFloat("Yaw",   &camera.yaw,   -180.0f, 180.0f);
    ImGui::SliderFloat("Pitch", &camera.pitch, -89.0f,  89.0f);
    ImGui::SliderFloat("Distance", &camera.distance, 1.0f, 50.0f);

    ImGui::Separator();

    if(ImGui::Button("Reset Camera"))
    {
        camera.yaw = 0.0f;
        camera.pitch = 20.0f;
        camera.distance = 10.0f;
    }

    ImGui::Separator();

    // Delete button for selected atoms
    if (onDeleteSelected)
    {
        if (ImGui::Button("Delete Selected", ImVec2(-1, 0)))
        {
            onDeleteSelected();
        }
    }

    ImGui::Separator();

    ImGui::Text("Mouse:");
    ImGui::Text("Left Drag : Rotate");
    ImGui::Text("Scroll    : Zoom");
    ImGui::Text("Ctrl+Click: Multi-select");

    ImGui::End();
}