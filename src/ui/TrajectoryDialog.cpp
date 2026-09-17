#include "ui/TrajectoryDialog.h"
#include "ui/ResponsiveLayout.h"
#include "util/ElementData.h"
#include "imgui.h"
#include <algorithm>

void TrajectoryDialog::drawMenuItem()
{
    if (ImGui::MenuItem("Trajectory playback")) open=true;
}

void TrajectoryDialog::draw(Structure& structure,const std::function<void(Structure&)>& update)
{
    bool changed=false;
    if (task.poll()) {
        error=task.error();
        if (task.result()) { frames=std::move(*task.result()); playback.frame=0; playback.playing=false; changed=true; }
        task.clearResult();
    }
    if (!open) {
        const bool wasPlaying=playback.playing; playback.playing=false;
        if (wasPlaying && !frames.empty()) { structure=frames[static_cast<std::size_t>(playback.frame)]; update(structure); }
        return;
    }
    responsive::windowSize(ImVec2(560,390),ImGuiCond_FirstUseEver);
    if (responsive::begin("Trajectory playback",&open,ImGuiWindowFlags_NoCollapse)) {
        ImGui::BeginDisabled(task.running());
        if (responsive::button("Open trajectory...")) picker.open("Open trajectory (XYZ/extXYZ or LAMMPS dump)",false,"trajectory.xyz");
        ImGui::EndDisabled();
        if (task.running()) { ImGui::TextUnformatted("Loading frames..."); if (responsive::button("Cancel")) task.cancel(); }
        if (!frames.empty()) {
            const bool togglePlay = responsive::button(playback.playing ? "Pause" : "Play");
            ImGui::SameLine();
            const bool jumpToFirst = responsive::button("First frame");
            ImGui::TextUnformatted("Frame (zero-based)"); ImGui::SetNextItemWidth(-FLT_MIN);
            const bool sliderChanged = ImGui::SliderInt("##frame",&playback.frame,0,static_cast<int>(frames.size())-1);
            ImGui::TextUnformatted("Frames per second"); ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##fps",&fps,1,60,"%.0f");
            changed |= playback.tick(togglePlay,jumpToFirst,sliderChanged,static_cast<int>(frames.size()),fps,ImGui::GetTime());
            ImGui::Text("%zu frames, %zu atoms in current frame",frames.size(),frames[static_cast<std::size_t>(playback.frame)].atoms.size());
        }
        ImGui::TextWrapped("Playback replaces the active structure. Camera controls remain available in the main view. GIF export and additional trajectory formats are available through the Python science API.");
        if (!error.empty()) ImGui::TextWrapped("%s",error.c_str());
        if (auto path=picker.draw()) { playback.playing=false; task.start([path=*path] { return loadTrajectory(path); }); }
        if (changed && !frames.empty()) { structure=frames[static_cast<std::size_t>(playback.frame)]; update(structure); }
    }
    ImGui::End();
}
