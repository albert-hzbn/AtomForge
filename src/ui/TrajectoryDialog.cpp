#include "ui/TrajectoryDialog.h"
#include "ui/ResponsiveLayout.h"
#include "util/ElementData.h"
#include "imgui.h"
#include <algorithm>

void TrajectoryDialog::drawMenuItem()
{
    if (ImGui::MenuItem("Trajectory playback...")) open=true;
}

void TrajectoryDialog::draw(Structure& structure,const std::function<void(Structure&)>& update)
{
    bool changed=false;
    if (task.poll()) {
        error=task.error();
        if (task.result()) { frames=std::move(*task.result()); frame=0; changed=true; }
        task.clearResult();
    }
    if (!open) {
        const bool wasPlaying=playing; playing=false;
        if (wasPlaying && !frames.empty()) { structure=frames[frame]; update(structure); }
        return;
    }
    responsive::windowSize(ImVec2(560,390),ImGuiCond_FirstUseEver);
    if (responsive::begin("Trajectory playback",&open,ImGuiWindowFlags_NoCollapse)) {
        ImGui::BeginDisabled(task.running());
        if (responsive::button("Open XYZ trajectory...")) picker.open("Open XYZ trajectory",false,"trajectory.xyz");
        ImGui::EndDisabled();
        if (task.running()) { ImGui::TextUnformatted("Loading frames..."); if (responsive::button("Cancel")) task.cancel(); }
        if (!frames.empty()) {
            if (responsive::button(playing ? "Pause" : "Play")) { playing=!playing; previous=ImGui::GetTime(); changed=!playing; }
            ImGui::SameLine();
            if (responsive::button("First frame")) { frame=0; changed=true; playing=false; }
            ImGui::TextUnformatted("Frame (zero-based)"); ImGui::SetNextItemWidth(-FLT_MIN);
            changed |= ImGui::SliderInt("##frame",&frame,0,static_cast<int>(frames.size())-1);
            ImGui::TextUnformatted("Frames per second"); ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##fps",&fps,1,60,"%.0f");
            if (playing && ImGui::GetTime()-previous>=1.0/fps) {
                frame=(frame+1)%static_cast<int>(frames.size()); previous=ImGui::GetTime(); changed=true;
            }
            ImGui::Text("%zu frames, %zu atoms in current frame",frames.size(),frames[frame].atoms.size());
        }
        ImGui::TextWrapped("Playback replaces the active structure. Camera controls remain available in the main view. GIF export and additional trajectory formats are available through the Python science API.");
        if (!error.empty()) ImGui::TextWrapped("%s",error.c_str());
        if (auto path=picker.draw()) { playing=false; task.start([path=*path] { return loadXyzTrajectory(path); }); }
        if (changed && !frames.empty()) { structure=frames[frame]; update(structure); }
    }
    ImGui::End();
}
