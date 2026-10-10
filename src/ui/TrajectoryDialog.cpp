#include "ui/TrajectoryDialog.h"
#include "science/ScienceTools.h"
#include "ui/ResponsiveLayout.h"
#include "util/ElementData.h"
#include "imgui.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>


void TrajectoryDialog::drawMenuItem()
{
    if (ImGui::MenuItem("Trajectory playback")) open=true;
}

const std::vector<std::pair<std::string, std::string>>& TrajectoryDialog::frameAnalyses()
{
    static const std::vector<std::pair<std::string, std::string>> tools = {
        {"structure-type", "Structure type (Ackland-Jones)"},
        {"centrosymmetry", "Centrosymmetry parameter"},
        {"bond-order", "Steinhardt q6 bond order"},
    };
    return tools;
}

atomforge::science::AtomProperty TrajectoryDialog::analyseFrame(const Structure& frame, const std::string& tool, double cutoff)
{
    return atomforge::science::analysePerAtom(frame,tool,cutoff);
}

void TrajectoryDialog::show(int index, Structure& structure, const std::function<void(Structure&)>& update, const std::function<void()>& showAtomProperty)
{
    if (!stream || index<0 || static_cast<std::size_t>(index)>=stream->size()) return;
    try {
        structure=stream->frame(static_cast<std::size_t>(index)).structure;
        shownAtoms=structure.atoms.size();
        if (analyse) {
            const auto start=std::chrono::steady_clock::now();
            const auto& tool=frameAnalyses()[static_cast<std::size_t>(std::clamp(analysis,0,static_cast<int>(frameAnalyses().size())-1))];
            const auto property=analyseFrame(structure,tool.first,analysisCutoff);
            structure.atomProperty=property.values;
            structure.atomPropertyName=property.name;
            if (showAtomProperty) showAtomProperty();
            const double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
            char text[96];
            std::snprintf(text,sizeof(text),"%s: %.0f ms per frame",tool.second.c_str(),ms);
            analysisMessage=text;
        }
        update(structure);
    } catch (const std::exception& e) {
        error=e.what();
        if (analyse) { analyse=false; analysisMessage="Per-frame analysis stopped: "+error; }
    }
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
        if (gifFrame>=static_cast<int>(frameCount())) finishGif("Saved "+std::to_string(gif->frames())+" frames to "+gifPath);
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

void TrajectoryDialog::draw(Structure& structure,const std::function<void(Structure&)>& update,const std::function<void()>& showAtomProperty)
{
    bool changed=false;
    if (task.poll()) {
        error=task.error();
        if (task.result()) { stream=*task.result(); playback.frame=0; playback.playing=false; changed=true; loaded=pending; }
        task.clearResult();
    }
    const int count=static_cast<int>(frameCount());
    if (gifFrame>=0 && gifFrame>=count) finishGif("GIF export stopped: the trajectory changed");
    if (gifFrame>=0 && !gifShown) {
        show(gifFrame,structure,update,showAtomProperty);
        playback.frame=gifFrame;
        gifShown=true;
    }
    if (gifRestorePending) {
        gifRestorePending=false;
        if (count>0) {
            playback.frame=std::min(gifRestoreFrame,count-1);
            show(playback.frame,structure,update,showAtomProperty);
        }
    }
    if (!open) {
        const bool wasPlaying=playback.playing; playback.playing=false;
        if (wasPlaying && count>0) show(playback.frame,structure,update,showAtomProperty);
        return;
    }
    responsive::windowSize(ImVec2(560,470),ImGuiCond_FirstUseEver);
    if (responsive::begin("Trajectory playback",&open,ImGuiWindowFlags_NoCollapse)) {
        const bool exporting=gifFrame>=0;
        ImGui::BeginDisabled(task.running() || exporting);
        if (responsive::button("Open trajectory")) picker.open("Open trajectory (XYZ/extXYZ, XDATCAR, POSCAR or LAMMPS dump)",false,"trajectory.xyz");
        ImGui::EndDisabled();
        if (task.running()) { ImGui::TextUnformatted("Indexing frames..."); if (responsive::button("Cancel")) task.cancel(); }
        if (count>0) {
            ImGui::BeginDisabled(exporting);
            const bool togglePlay = responsive::button(playback.playing ? "Pause" : "Play");
            ImGui::SameLine();
            const bool jumpToFirst = responsive::button("First frame");
            ImGui::TextUnformatted("Frame (zero-based)"); ImGui::SetNextItemWidth(-FLT_MIN);
            const bool sliderChanged = ImGui::SliderInt("##frame",&playback.frame,0,count-1);
            ImGui::TextUnformatted("Frames per second"); ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##fps",&fps,1,60,"%.0f");
            if (!exporting) changed |= playback.tick(togglePlay,jumpToFirst,sliderChanged,count,fps,ImGui::GetTime());
            ImGui::Text("%d frames%s, %zu atoms in current frame",count,stream->streamed() ? " (read on demand)" : "",shownAtoms);
            // Per-frame analysis: recomputed and shown as atom colours for every frame.
            ImGui::SeparatorText("Colour each frame");
            if (ImGui::Checkbox("Analyse every displayed frame",&analyse)) { changed=true; analysisMessage.clear(); }
            ImGui::SetNextItemWidth(std::min(ImGui::GetContentRegionAvail().x,responsive::dp(320)));
            if (ImGui::BeginCombo("##frameanalysis",frameAnalyses()[static_cast<std::size_t>(analysis)].second.c_str())) {
                for (std::size_t t=0;t<frameAnalyses().size();++t)
                    if (ImGui::Selectable(frameAnalyses()[t].second.c_str(),static_cast<int>(t)==analysis)) { analysis=static_cast<int>(t); changed=analyse; }
                ImGui::EndCombo();
            }
            if (frameAnalyses()[static_cast<std::size_t>(analysis)].first!="structure-type" || !structure.hasUnitCell) {
                ImGui::SetNextItemWidth(std::min(ImGui::GetContentRegionAvail().x,responsive::dp(200)));
                if (ImGui::InputFloat("Neighbour cutoff (Angstrom)",&analysisCutoff,0.1f,0.5f,"%.2f")) { analysisCutoff=std::max(0.5f,analysisCutoff); changed=analyse; }
            }
            if (!analysisMessage.empty()) ImGui::TextDisabled("%s",analysisMessage.c_str());
            ImGui::SeparatorText("Animation");
            ImGui::TextUnformatted("GIF: use every Nth frame"); ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderInt("##gifstride",&gifStride,1,std::max(1,count/2));
            gifStride=std::max(1,gifStride);
            if (responsive::button("Export GIF")) gifPicker.open("Export trajectory animation as GIF",true,"trajectory.gif");
            ImGui::EndDisabled();
            if (exporting) {
                ImGui::Text("Writing GIF: frame %d of %d",gifFrame+1,count);
                if (responsive::button("Stop export")) finishGif("GIF export stopped; "+std::to_string(gif ? gif->frames() : 0)+" frames written to "+gifPath);
            }
            if (!gifMessage.empty()) ImGui::TextWrapped("%s",gifMessage.c_str());
        }
        ImGui::TextWrapped("Playback replaces the active structure. Frames are read from the file as they are shown, so long trajectories need little memory. Loaded trajectories can be analysed directly with the trajectory and transport tools. GIF export renders every chosen frame with the current camera and frame rate, at most 960 pixels wide.");
        if (!error.empty()) ImGui::TextWrapped("%s",error.c_str());
        if (auto path=picker.draw()) {
            playback.playing=false;
            pending=*path;
            error.clear();
            task.start([path=*path] {
                return std::shared_ptr<const atomforge::science::TrajectoryStream>(
                    std::make_shared<atomforge::science::TrajectoryStream>(std::filesystem::u8path(path)));
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
        if (changed && count>0) show(playback.frame,structure,update,showAtomProperty);
    }
    ImGui::End();
}
