#include "ui/TrajectoryDialog.h"
#include "science/ScienceData.h"
#include "ui/ResponsiveLayout.h"
#include "util/ElementData.h"
#include "imgui.h"
#include <algorithm>
#include <cmath>
#include <filesystem>

void TrajectoryDialog::drawMenuItem()
{
    if (ImGui::MenuItem("Trajectory playback")) open=true;
}

void TrajectoryDialog::addGifFrame(const std::vector<unsigned char>& rgba,int width,int height)
{
    if (gifFrame<0) return;
    if (width<=0 || height<=0 || rgba.size()<static_cast<std::size_t>(width)*static_cast<std::size_t>(height)*4) { finishGif("GIF export failed: the view could not be captured"); return; }
    try {
        // Box-downsample large views so the GIF stays a sensible size.
        const int factor=std::max(1,(std::max(width,height)+959)/960);
        const int w=std::max(1,width/factor), h=std::max(1,height/factor);
        std::vector<unsigned char> small(static_cast<std::size_t>(w)*static_cast<std::size_t>(h)*3);
        for (int y=0;y<h;++y)
            for (int x=0;x<w;++x)
                for (int c=0;c<3;++c) {
                    int sum=0;
                    for (int dy=0;dy<factor;++dy)
                        for (int dx=0;dx<factor;++dx)
                            sum+=rgba[(static_cast<std::size_t>(y*factor+dy)*static_cast<std::size_t>(width)+static_cast<std::size_t>(x*factor+dx))*4+static_cast<std::size_t>(c)];
                    small[(static_cast<std::size_t>(y)*static_cast<std::size_t>(w)+static_cast<std::size_t>(x))*3+static_cast<std::size_t>(c)]=static_cast<unsigned char>(sum/(factor*factor));
                }
        if (!gif) gif=std::make_unique<atomforge::science::GifWriter>(std::filesystem::u8path(gifPath),w,h,std::max(2,static_cast<int>(std::lround(100.0/fps))));
        gif->addFrame(small,3);
        gifFrame+=gifStride;
        gifShown=false;
        if (gifFrame>=static_cast<int>(frames.size())) finishGif("Saved "+std::to_string(gif->frames())+" frames to "+gifPath);
    } catch (const std::exception& e) {
        finishGif(std::string("GIF export failed: ")+e.what());
    }
}

void TrajectoryDialog::finishGif(const std::string& message)
{
    try { if (gif) gif->close(); } catch (const std::exception&) {}
    gif.reset();
    gifFrame=-1;
    gifShown=false;
    gifMessage=message;
    gifRestorePending=true;
}

void TrajectoryDialog::draw(Structure& structure,const std::function<void(Structure&)>& update)
{
    bool changed=false;
    if (task.poll()) {
        error=task.error();
        if (task.result()) { frames=std::move(*task.result()); playback.frame=0; playback.playing=false; changed=true; loaded=pending; }
        task.clearResult();
    }
    if (gifFrame>=0 && gifFrame>=static_cast<int>(frames.size())) finishGif("GIF export stopped: the trajectory changed");
    if (gifFrame>=0 && !gifShown) {
        structure=frames[static_cast<std::size_t>(gifFrame)]; update(structure);
        playback.frame=gifFrame;
        gifShown=true;
    }
    if (gifRestorePending) {
        gifRestorePending=false;
        if (!frames.empty()) {
            playback.frame=std::min(gifRestoreFrame,static_cast<int>(frames.size())-1);
            structure=frames[static_cast<std::size_t>(playback.frame)]; update(structure);
        }
    }
    if (!open) {
        const bool wasPlaying=playback.playing; playback.playing=false;
        if (wasPlaying && !frames.empty()) { structure=frames[static_cast<std::size_t>(playback.frame)]; update(structure); }
        return;
    }
    responsive::windowSize(ImVec2(560,390),ImGuiCond_FirstUseEver);
    if (responsive::begin("Trajectory playback",&open,ImGuiWindowFlags_NoCollapse)) {
        const bool exporting=gifFrame>=0;
        ImGui::BeginDisabled(task.running() || exporting);
        if (responsive::button("Open trajectory...")) picker.open("Open trajectory (XYZ/extXYZ, XDATCAR, POSCAR or LAMMPS dump)",false,"trajectory.xyz");
        ImGui::EndDisabled();
        if (task.running()) { ImGui::TextUnformatted("Loading frames..."); if (responsive::button("Cancel")) task.cancel(); }
        if (!frames.empty()) {
            ImGui::BeginDisabled(exporting);
            const bool togglePlay = responsive::button(playback.playing ? "Pause" : "Play");
            ImGui::SameLine();
            const bool jumpToFirst = responsive::button("First frame");
            ImGui::TextUnformatted("Frame (zero-based)"); ImGui::SetNextItemWidth(-FLT_MIN);
            const bool sliderChanged = ImGui::SliderInt("##frame",&playback.frame,0,static_cast<int>(frames.size())-1);
            ImGui::TextUnformatted("Frames per second"); ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##fps",&fps,1,60,"%.0f");
            if (!exporting) changed |= playback.tick(togglePlay,jumpToFirst,sliderChanged,static_cast<int>(frames.size()),fps,ImGui::GetTime());
            ImGui::Text("%zu frames, %zu atoms in current frame",frames.size(),frames[static_cast<std::size_t>(playback.frame)].atoms.size());
            ImGui::TextUnformatted("GIF: use every Nth frame"); ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderInt("##gifstride",&gifStride,1,std::max(1,static_cast<int>(frames.size())/2));
            gifStride=std::max(1,gifStride);
            if (responsive::button("Export GIF...")) gifPicker.open("Export trajectory animation as GIF",true,"trajectory.gif");
            ImGui::EndDisabled();
            if (exporting) {
                ImGui::Text("Writing GIF: frame %d of %zu",gifFrame+1,frames.size());
                if (responsive::button("Stop export")) finishGif("GIF export stopped; "+std::to_string(gif ? gif->frames() : 0)+" frames written to "+gifPath);
            }
            if (!gifMessage.empty()) ImGui::TextWrapped("%s",gifMessage.c_str());
        }
        ImGui::TextWrapped("Playback replaces the active structure. Camera controls remain available in the main view. Loaded trajectories can be analysed directly with the trajectory and transport tools. GIF export renders every chosen frame with the current camera and frame rate, at most 960 pixels wide.");
        if (!error.empty()) ImGui::TextWrapped("%s",error.c_str());
        if (auto path=picker.draw()) {
            playback.playing=false;
            pending=*path;
            task.start([path=*path] {
                std::vector<Structure> result;
                for (auto& frame : atomforge::science::readFrames(std::filesystem::u8path(path))) result.push_back(std::move(frame.structure));
                return result;
            });
        }
        if (auto path=gifPicker.draw()) {
            gifPath=*path;
            if (std::filesystem::u8path(gifPath).extension().empty()) gifPath+=".gif";
            playback.playing=false;
            gifRestoreFrame=playback.frame;
            gifFrame=0;
            gifShown=false;
            gifMessage.clear();
        }
        if (changed && !frames.empty()) { structure=frames[static_cast<std::size_t>(playback.frame)]; update(structure); }
    }
    ImGui::End();
}
