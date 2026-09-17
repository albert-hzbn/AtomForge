#include "ui/ResponsiveLayout.h"
#include "ui/ElectronicPostProcessing.h"
#include "electronic/DisplayRange.h"
#include "electronic/ChargeAnalysis.h"
#include "electronic/Topology.h"
#include "electronic/BaderPartition.h"
#include "ui/DialogLayout.h"
#include "imgui.h"
#include "third_party/stb_image_write.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

using namespace atomforge::electronic;

namespace
{
glm::dvec3 vector(const float* p) { return {p[0],p[1],p[2]}; }

// Match the padded label/value rows used by the other analysis dialogs.
template<class DrawControl>
bool field(const char* label, DrawControl drawControl)
{
    bool changed = false;
    ImGui::PushID(label);
    if (ImGui::BeginTable("Field",2,ImGuiTableFlags_SizingStretchProp))
    {
        ImGui::TableSetupColumn("Label",ImGuiTableColumnFlags_WidthStretch,.40f);
        ImGui::TableSetupColumn("Value",ImGuiTableColumnFlags_WidthStretch,.60f);
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::AlignTextToFramePadding();
        ImGui::TextWrapped("%s",label);
        ImGui::TableNextColumn();
        ImGui::SetNextItemWidth(-FLT_MIN);
        changed = drawControl();
        ImGui::EndTable();
    }
    ImGui::PopID();
    return changed;
}

bool combo(const char* label, int* selected, const char* items)
{
    return field(label,[&] { return ImGui::Combo("##value",selected,items); });
}

void inputFloat(const char* label, float* value, float step=0, float fastStep=0, const char* format="%.3f")
{
    field(label,[&] { return ImGui::InputFloat("##value",value,step,fastStep,format); });
}

void inputInt(const char* label, int* value)
{
    field(label,[&] { return ImGui::InputInt("##value",value); });
}

void inputVector(const char* label, float* value)
{
    field(label,[&] { return ImGui::InputFloat3("##value",value); });
}

void slider(const char* label, float* value, float low, float high, const char* format)
{
    field(label,[&] { return ImGui::SliderFloat("##value",value,low,high,format); });
}

void filename(const char* label, const std::string& path)
{
    if (path.empty()) return;
    ImGui::TextWrapped("%s: %s",label,std::filesystem::u8path(path).filename().u8string().c_str());
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s",path.c_str());
}
}

atomforge::Workspace ElectronicPostProcessingDialog::snapshot(bool includeData) const
{
    atomforge::Workspace out;
    if (includeData) {
        out.volume=m_volume; out.reference=m_referenceVolume; out.surface=m_surface;
        out.table=m_result.table; out.columns=m_result.columns; out.heading=m_result.heading;
    }
    out.sourcePath=m_loadedPath; out.referencePath=m_referencePath; out.history=m_history;
    out.sliceSettings=m_sliceViewport.settings();
    out.settings["m_quantity"]=m_quantity;
    out.settings["m_cubeUnits"]=m_cubeUnits;
    out.settings["m_selected"]=m_selected;
    out.settings["m_reference"]=m_reference;
    out.settings["m_localReference"]=m_localReference;
    out.settings["m_operation"]=m_operation;
    out.settings["m_toolGroup"]=m_toolGroup;
    out.settings["m_axis"]=m_axis;
    out.settings["m_count"]=m_count;
    out.settings["m_window"]=m_window;
    out.settings["m_radius"]=m_radius;
    out.settings["m_exportFormat"]=m_exportFormat;
    out.settings["m_scalar"]=m_scalar;
    out.settings["m_maskLow"]=m_maskLow;
    out.settings["m_maskHigh"]=m_maskHigh;
    out.settings["m_referenceWeight"]=m_referenceWeight;
    out.settings["m_booleanOperation"]=m_booleanOperation;
    out.settings["m_surfaceLevel"]=m_surfaceLevel;
    out.settings["m_suggestedLevel"]=m_suggestedLevel;
    out.settings["m_renderMode"]=m_renderMode;
    out.settings["m_viewLayout"]=m_viewLayout;
    out.settings["m_showSlicePlane"]=m_showSlicePlane;
    out.settings["m_sigma"]=m_sigma;
    out.settings["m_bondingFloor"]=m_bondingFloor;
    out.settings["m_alpha"]=m_alpha;
    out.settings["m_realCutoff"]=m_realCutoff;
    out.settings["m_reciprocalCutoff"]=m_reciprocalCutoff;
    out.settings["m_yaw"]=m_yaw;
    out.settings["m_pitch"]=m_pitch;
    out.settings["m_opacity"]=m_opacity;
    out.settings["m_zoom"]=m_zoom;
    out.settings["m_specular"]=m_specular;
    out.settings["m_shininess"]=m_shininess;
    out.settings["m_colorLow"]=m_colorLow;
    out.settings["m_colorHigh"]=m_colorHigh;
    out.settings["m_autoLow"]=m_autoLow;
    out.settings["m_autoHigh"]=m_autoHigh;
    out.settings["m_palette"]=m_palette;
    out.settings["m_quality"]=m_quality;
    out.settings["m_autoRange"]=m_autoRange;
    out.settings["m_sliceLow"]=m_sliceLow;
    out.settings["m_sliceHigh"]=m_sliceHigh;
    out.settings["m_colorSurface"]=m_colorSurface;
    out.settings["m_appendSurface"]=m_appendSurface;
    out.settings["m_pan.x"]=m_pan.x;
    out.settings["m_pan.y"]=m_pan.y;
    out.settings["m_start[0]"]=m_start[0];
    out.settings["m_start[1]"]=m_start[1];
    out.settings["m_start[2]"]=m_start[2];
    out.settings["m_end[0]"]=m_end[0];
    out.settings["m_end[1]"]=m_end[1];
    out.settings["m_end[2]"]=m_end[2];
    out.settings["m_u[0]"]=m_u[0];
    out.settings["m_u[1]"]=m_u[1];
    out.settings["m_u[2]"]=m_u[2];
    out.settings["m_v[0]"]=m_v[0];
    out.settings["m_v[1]"]=m_v[1];
    out.settings["m_v[2]"]=m_v[2];
    return out;
}

void ElectronicPostProcessingDialog::restore(const atomforge::Workspace& saved)
{
    if (m_task.running()) throw std::runtime_error("Wait for or cancel the calculation before restoring a workspace");
    // Validate before changing any live state or narrowing a stored double.
    for (const auto* settings:{&saved.settings,&saved.sliceSettings})
        for (const auto& entry:*settings)
            if (!std::isfinite(entry.second) || std::abs(entry.second)>std::numeric_limits<int>::max())
                throw std::runtime_error("Workspace setting outside the supported range: "+entry.first);
    m_volume=saved.volume; m_referenceVolume=saved.reference; m_surface=saved.surface;
    m_result={}; m_result.table=saved.table; m_result.columns=saved.columns; m_result.heading=saved.heading;
    m_loadedPath=saved.sourcePath; m_referencePath=saved.referencePath; m_history=saved.history;
    const auto setting=[&](const char* name,auto& target) {
        const auto found=saved.settings.find(name);
        if (found!=saved.settings.end()) target=static_cast<std::decay_t<decltype(target)>>(found->second);
    };
    setting("m_quantity",m_quantity);
    setting("m_cubeUnits",m_cubeUnits);
    setting("m_selected",m_selected);
    setting("m_reference",m_reference);
    setting("m_localReference",m_localReference);
    setting("m_operation",m_operation);
    setting("m_toolGroup",m_toolGroup);
    setting("m_axis",m_axis);
    setting("m_count",m_count);
    setting("m_window",m_window);
    setting("m_radius",m_radius);
    setting("m_exportFormat",m_exportFormat);
    setting("m_scalar",m_scalar);
    setting("m_maskLow",m_maskLow);
    setting("m_maskHigh",m_maskHigh);
    setting("m_referenceWeight",m_referenceWeight);
    setting("m_booleanOperation",m_booleanOperation);
    setting("m_surfaceLevel",m_surfaceLevel);
    setting("m_suggestedLevel",m_suggestedLevel);
    setting("m_renderMode",m_renderMode);
    setting("m_viewLayout",m_viewLayout);
    setting("m_showSlicePlane",m_showSlicePlane);
    setting("m_sigma",m_sigma);
    setting("m_bondingFloor",m_bondingFloor);
    setting("m_alpha",m_alpha);
    setting("m_realCutoff",m_realCutoff);
    setting("m_reciprocalCutoff",m_reciprocalCutoff);
    setting("m_yaw",m_yaw);
    setting("m_pitch",m_pitch);
    setting("m_opacity",m_opacity);
    setting("m_zoom",m_zoom);
    setting("m_specular",m_specular);
    setting("m_shininess",m_shininess);
    setting("m_colorLow",m_colorLow);
    setting("m_colorHigh",m_colorHigh);
    setting("m_autoLow",m_autoLow);
    setting("m_autoHigh",m_autoHigh);
    setting("m_palette",m_palette);
    setting("m_quality",m_quality);
    setting("m_autoRange",m_autoRange);
    setting("m_sliceLow",m_sliceLow);
    setting("m_sliceHigh",m_sliceHigh);
    setting("m_colorSurface",m_colorSurface);
    setting("m_appendSurface",m_appendSurface);
    setting("m_pan.x",m_pan.x);
    setting("m_pan.y",m_pan.y);
    setting("m_start[0]",m_start[0]);
    setting("m_start[1]",m_start[1]);
    setting("m_start[2]",m_start[2]);
    setting("m_end[0]",m_end[0]);
    setting("m_end[1]",m_end[1]);
    setting("m_end[2]",m_end[2]);
    setting("m_u[0]",m_u[0]);
    setting("m_u[1]",m_u[1]);
    setting("m_u[2]",m_u[2]);
    setting("m_v[0]",m_v[0]);
    setting("m_v[1]",m_v[1]);
    setting("m_v[2]",m_v[2]);
    m_selected=std::clamp(m_selected,0,std::max(0,static_cast<int>(m_volume.fields.size())-1));
    const auto& references=m_localReference || m_referenceVolume.fields.empty() ? m_volume : m_referenceVolume;
    m_reference=std::clamp(m_reference,0,std::max(0,static_cast<int>(references.fields.size())-1));
    m_operation=std::clamp(m_operation,0,37); m_axis=std::clamp(m_axis,0,2);
    m_quantity=std::clamp(m_quantity,0,4); m_cubeUnits=std::clamp(m_cubeUnits,0,1);
    m_exportFormat=std::clamp(m_exportFormat,0,7); m_booleanOperation=std::clamp(m_booleanOperation,0,3);
    m_quality=std::clamp(m_quality,0,3); m_palette=std::clamp(m_palette,0,2);
    m_renderMode=std::clamp(m_renderMode,0,1); m_viewLayout=std::clamp(m_viewLayout,0,2);
    m_opacity=std::clamp(m_opacity,0.0f,1.0f); m_zoom=std::clamp(m_zoom,.1f,100.0f);
    m_viewport.setQuality(m_quality); m_viewport.setMesh(m_surface); m_volumeDirty=true;
    m_sliceViewport.setResolution(std::array<int,4>{65,129,257,513}[m_quality]);
    m_sliceViewport.restoreSettings(saved.sliceSettings); m_sliceField=m_selected;
    m_open=!m_volume.fields.empty() || !m_referenceVolume.fields.empty();
}

void ElectronicPostProcessingDialog::remember()
{
    std::size_t bytes=m_surface.vertices.size()*sizeof(glm::dvec3)+m_surface.colors.size()*sizeof(double);
    for (const auto& grid:m_volume.fields) bytes+=grid.values.size()*sizeof(double);
    for (const auto& grid:m_referenceVolume.fields) bytes+=grid.values.size()*sizeof(double);
    m_redo.clear();
    // Keep interactive undo bounded; disk projects retain full-resolution data.
    if (bytes>32*1024*1024) { m_undo.clear(); return; }
    m_undo.push_back(snapshot());
    while (m_undo.size()>4) m_undo.pop_front();
}

void ElectronicPostProcessingDialog::drawMenuItem()
{
    if (ImGui::MenuItem("Electronic Post-processing")) m_open = true;
}

void ElectronicPostProcessingDialog::feedDroppedFile(const std::string& path)
{
    if (m_open && !path.empty()) m_pendingDrops.emplace_back(path,m_dropReference);
}

void ElectronicPostProcessingDialog::drawDialog()
{
    if (m_task.poll())
    {
        m_error = m_task.error();
        if (m_task.result())
        {
            remember();
            if (!m_pendingOperation.empty()) m_history.push_back(m_pendingOperation);
            m_result = *m_task.result();
            if (m_result.surfaceReady)
            {
                if (!m_result.appendMesh) m_surface = {};
                if (m_surface.vertices.size() + m_result.mesh.vertices.size() <= 12000000)
                {
                    m_surface.vertices.insert(m_surface.vertices.end(),m_result.mesh.vertices.begin(),m_result.mesh.vertices.end());
                    m_surface.colors.insert(m_surface.colors.end(),m_result.mesh.colors.begin(),m_result.mesh.colors.end());
                    m_viewport.setMesh(m_surface);
                    m_autoLow = m_result.appendMesh ? std::min(m_autoLow,static_cast<float>(m_result.colorLow)) : static_cast<float>(m_result.colorLow);
                    m_autoHigh = m_result.appendMesh ? std::max(m_autoHigh,static_cast<float>(m_result.colorHigh)) : static_cast<float>(m_result.colorHigh);
                    m_colorUnit = m_result.colorUnit;
                    if (m_autoRange) { m_colorLow=m_autoLow; m_colorHigh=m_autoHigh; }
                }
                else m_error = "Combined surfaces exceed the mesh limit; clear or replace the surface.";
                m_result.mesh = {};
            }
            if (m_result.referenceLoaded) { m_referenceVolume = m_result.volume; m_reference = 0; m_referencePath=m_result.sourcePath; }
            else if (m_result.loaded)
            {
                m_volume = m_result.volume; m_selected = m_reference = 0;
                m_loadedPath=m_result.sourcePath; m_sliceField=-1; m_sliceViewport.invalidate(); m_sliceViewport.reset();
                m_surface = {}; m_viewport.setMesh(m_surface); resetCamera();
            }
            else if (!m_result.volume.fields.empty())
            {
                m_selected = static_cast<int>(m_volume.fields.size());
                m_volume.fields.insert(m_volume.fields.end(),m_result.volume.fields.begin(),m_result.volume.fields.end());
            }
        }
        m_task.clearResult();
    }
    if (!m_open) { m_pendingDrops.clear(); return; }
    if (!m_task.running() && !m_pendingDrops.empty())
    {
        const auto dropped = std::move(m_pendingDrops.front());
        m_pendingDrops.pop_front();
        load(dropped.first,dropped.second);
    }
    responsive::windowSize(ImVec2(1380,860),ImGuiCond_FirstUseEver);
    responsive::windowConstraints(ImVec2(980,680),ImVec2(FLT_MAX,FLT_MAX));
    if (!responsive::begin("Electronic Post-processing",&m_open,ImGuiWindowFlags_NoScrollbar|ImGuiWindowFlags_NoCollapse)) { ImGui::End(); return; }
    const bool stackPanels = responsive::stacked(1000);
    const float sidebar = stackPanels ? ImGui::GetContentRegionAvail().x : std::clamp(ImGui::GetContentRegionAvail().x * .38f,responsive::dp(340),responsive::dp(440));
    responsive::beginChild("Electronic controls",ImVec2(sidebar,stackPanels ? responsive::dp(420) : 0),true);
    ImGui::BeginDisabled(m_task.running());
    if (responsive::button("Open workspace...")) { m_pickerAction=3; m_picker.open("Open electronic workspace",false,"workspace.afproject"); }
    ImGui::SameLine();
    if (responsive::button("Save workspace...")) { m_pickerAction=4; m_picker.open("Save electronic workspace",true,"workspace.afproject"); }
    ImGui::BeginDisabled(m_undo.empty());
    if (responsive::button("Undo calculation")) {
        m_redo.push_back(snapshot()); auto previous=std::move(m_undo.back()); m_undo.pop_back(); restore(previous);
    }
    ImGui::EndDisabled(); ImGui::SameLine(); ImGui::BeginDisabled(m_redo.empty());
    if (responsive::button("Redo")) {
        m_undo.push_back(snapshot()); auto next=std::move(m_redo.back()); m_redo.pop_back(); restore(next);
    }
    ImGui::EndDisabled();
    ImGui::TextWrapped("Undo retains 4 states up to 32 MiB each. Save larger workspaces to disk.");
    if (ImGui::CollapsingHeader("Operation history"))
        for (const auto& entry:m_history) ImGui::TextWrapped("%s",entry.c_str());
    dialogLayout::section("Tools");
    if (ImGui::RadioButton("Charge transfer",m_toolGroup==1)) { m_toolGroup=1; m_operation=25; }
    ImGui::SameLine();
    if (ImGui::RadioButton("General analysis",m_toolGroup==0)) { m_toolGroup=0; m_operation=0; }
    ImGui::SameLine();
    if (ImGui::RadioButton("Bonding & topology",m_toolGroup==2)) { m_toolGroup=2; m_operation=33; }
    if (m_toolGroup==1)
    {
        const char* labels[]={"Density difference","Threshold mask","Boolean masks","Apply mask",
                              "Accumulation / depletion","Charge summary","Cumulative profile","Invert mask"};
        if (ImGui::BeginTable("Charge tools",2,ImGuiTableFlags_SizingStretchSame))
        {
            for (int i=0;i<8;++i)
            {
                ImGui::TableNextColumn();
                ImGui::PushID(i);
                if (ImGui::Selectable(labels[i],m_operation==25+i)) m_operation=25+i;
                ImGui::PopID();
            }
            ImGui::EndTable();
        }
        if (m_volume.fields.empty()) ImGui::TextWrapped("Load a density file below to use the selected tool.");
    }
    else if (m_toolGroup==2)
    {
        const char* labels[]={"Reduced density gradient","Signed density (NCI)","DORI","Betti curve","Bader partition"};
        if (ImGui::BeginTable("Bonding tools",2,ImGuiTableFlags_SizingStretchSame))
        {
            for (int i=0;i<5;++i)
            {
                ImGui::TableNextColumn();
                ImGui::PushID(i);
                if (ImGui::Selectable(labels[i],m_operation==33+i)) m_operation=33+i;
                ImGui::PopID();
            }
            ImGui::EndTable();
        }
        if (m_volume.fields.empty()) ImGui::TextWrapped("Load a density file below to use the selected tool.");
    }
    else combo("Tool",&m_operation,
            "Total integral\0Add reference\0Subtract reference (difference density)\0Multiply reference\0Divide by reference\0Scale\0Gaussian smoothing\0Cartesian gradient\0Laplacian\0Energy-density conversion\0Line profile\0Planar average\0Macroscopic average\0Plane section\0Contour segments (z=0)\0Peak search\0Voronoi site integration\0Sphere integration\0Structure factors\0Fourier synthesis\0Patterson density\0Ewald site potentials\0Isosurface\0Resample onto reference\0Verify periodic endpoint planes\0");
    dialogLayout::section("Input data");
    const float buttonWidth = (ImGui::GetContentRegionAvail().x-ImGui::GetStyle().ItemSpacing.x)*.5f;
    if (dialogLayout::primaryButton("Open volume...",ImVec2(buttonWidth,0))) { m_pickerAction=0; m_picker.open("Open electronic volume",false,m_loadedPath); }
    ImGui::SameLine();
    if (responsive::button("Reference...",ImVec2(buttonWidth,0))) { m_pickerAction=1; m_picker.open("Open reference volume",false,m_referencePath.empty() ? m_loadedPath : m_referencePath); }
    ImGui::BeginDisabled(m_loadedPath.empty());
    if (responsive::button("Reload volume",ImVec2(buttonWidth,0))) load(m_loadedPath,false);
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (responsive::button("Clear",ImVec2(buttonWidth,0)))
    {
        m_volume={}; m_referenceVolume={}; m_result={}; m_surface={}; m_pendingDrops.clear();
        m_viewport.setMesh(m_surface); m_selected=m_reference=0; m_sliceField=-1; m_loadedPath.clear(); m_referencePath.clear(); resetCamera();
    }
    filename("Volume",m_loadedPath);
    filename("Reference",m_referencePath);
    ImGui::Checkbox("Drop files as reference",&m_dropReference);
    ImGui::TextWrapped("Drop a VASP, Cube or XSF file onto the app to load it here.");
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Import options"))
    {
        combo("Quantity",&m_quantity,"Auto (Cube/XSF: raw)\0Density\0Potential\0ELF\0Raw\0");
        combo("Cube units",&m_cubeUnits,"Bohr\0Angstrom\0");
        ImGui::TextWrapped("Choose import settings before browsing, or reload the selected volume to apply changes.");
    }
    if (!m_volume.fields.empty())
    {
        auto selectField = [&](const char* label, int& selected, const Volume& volume)
        {
            selected = std::clamp(selected,0,static_cast<int>(volume.fields.size())-1);
            field(label,[&]
            {
                if (ImGui::BeginCombo("##value",volume.fields[selected].name.c_str()))
                {
                    for (int i = 0; i < static_cast<int>(volume.fields.size()); ++i)
                    {
                        ImGui::PushID(i);
                        if (ImGui::Selectable(volume.fields[i].name.c_str(),i == selected)) selected = i;
                        ImGui::PopID();
                    }
                    ImGui::EndCombo();
                }
                return false;
            });
        };
        dialogLayout::section("Analysis");
        selectField("Field",m_selected,m_volume);
        if (!m_referenceVolume.fields.empty()) ImGui::Checkbox("Use loaded/result field as reference",&m_localReference);
        const auto& referenceVolume = m_referenceVolume.fields.empty() || m_localReference ? m_volume : m_referenceVolume;
        m_reference = std::clamp(m_reference,0,static_cast<int>(referenceVolume.fields.size())-1);
        const Grid& g = m_volume.fields[m_selected];
        ImGui::TextDisabled("%d x %d x %d | %s",g.shape[0],g.shape[1],g.shape[2],g.unit.c_str());
        if ((m_operation >= 1 && m_operation <= 4) || m_operation == 23 || m_operation == 22) selectField("Reference",m_reference,referenceVolume);
        if (m_operation == 5) inputFloat("Scale factor",&m_scalar);
        if (m_operation == 6) { inputFloat("Sigma (A)",&m_sigma); inputInt("Radius (steps)",&m_radius); }
        if (m_operation == 9) ImGui::TextWrapped("Gradient-expansion kinetic and local-virial energy densities (eV/A^3). Requires nonnegative e/A^3. Samples <= 1e-12 e/A^3 are masked to zero.");
        if (m_operation == 10 || m_operation == 13 || m_operation == 17) inputVector("Origin (A)",m_start);
        if (m_operation == 10) inputVector("End (A)",m_end);
        if (m_operation == 10 || m_operation == 13) inputInt("Sample count",&m_count);
        if (m_operation == 11 || m_operation == 12) combo("Plane normal",&m_axis,"a*\0b*\0c*\0");
        if (m_operation == 12) inputInt("Window (odd)",&m_window);
        if (m_operation == 13) { inputVector("Plane span u (A)",m_u); inputVector("Plane span v (A)",m_v); }
        if (m_operation == 14) inputFloat("Isovalue",&m_scalar);
        if (m_operation == 17) inputFloat("Sphere radius (A)",&m_scalar);
        if (m_operation == 15) ImGui::TextWrapped("26-neighbor grid maxima; adjacent equal values use a deterministic tie break. Refine the grid to check peak positions.");
        if (m_operation == 16) ImGui::TextWrapped("Integrates nearest-site Voronoi cells around imported atoms. Boundary samples share weight. These are geometric partitions, not Bader basins.");
        if (m_operation == 18 || m_operation == 19)
        {
            ImGui::TextWrapped("%s",m_operation == 18 ? "One h k l per line" : "One h k l real imag per line; include conjugate pairs");
            ImGui::InputTextMultiline("##reflections",m_reflections,sizeof(m_reflections),ImVec2(-1,100));
        }
        if (m_operation == 21)
        {
            ImGui::TextWrapped("Explicit charges in electron units, one per imported atom in file order. Requires a neutral cell; conducting boundary. Increase both cutoffs to check convergence.");
            ImGui::TextUnformatted("Charges");
            ImGui::InputTextMultiline("##charges",m_charges,sizeof(m_charges),ImVec2(-1,90));
            inputFloat("Ewald alpha (1/A)",&m_alpha);
            inputFloat("Real cutoff (A)",&m_realCutoff);
            inputFloat("Reciprocal cutoff (1/A, includes 2*pi)",&m_reciprocalCutoff);
        }
        if (m_operation == 22)
        {
            ImGui::Checkbox("Color from reference",&m_colorSurface);
            ImGui::Checkbox("Append surface level",&m_appendSurface);
        }
        if (m_operation == 24) ImGui::TextWrapped("For Cube/XSF covering a complete periodic cell: checks all opposite endpoint values, then removes duplicate planes. Matching values alone do not establish physical periodicity.");
        if (m_operation==25 || m_operation==27 || m_operation==28) selectField("Operand",m_reference,referenceVolume);
        if (m_operation==25)
        {
            inputFloat("Reference weight",&m_referenceWeight);
            ImGui::TextWrapped("Result = field - weight * reference. Use aligned fragment densities from the same cell and geometry. Select the result and repeat for additional fragments.");
        }
        if (m_operation==26)
        {
            inputFloat("Lower bound",&m_maskLow,0,0,"%.6g"); inputFloat("Upper bound",&m_maskHigh,0,0,"%.6g");
            ImGui::TextWrapped("Creates a binary mask for the inclusive value range.");
        }
        if (m_operation==27) combo("Boolean operation",&m_booleanOperation,"Union\0Intersection\0Difference (field minus operand)\0XOR\0");
        if (m_operation==28) ImGui::TextWrapped("The operand must be a binary mask. Keeps field values inside the mask; integrate the result to measure regional charge.");
        if (m_operation==29 || m_operation==30) ImGui::TextWrapped("Apply to a difference density in e/A^3. Accumulation and depletion are positive electron counts, not atom-assigned transfer charges.");
        if (m_operation==31)
        {
            combo("Plane normal",&m_axis,"a*\0b*\0c*\0");
            ImGui::TextWrapped("Cumulative electrons from the lower cell face. The final row includes the whole cell; periodic profiles depend on the cell origin.");
        }
        if (m_operation>=33 && m_operation<=35)
        {
            inputFloat("Density floor (e/A^3)",&m_bondingFloor,0,0,"%.3g");
            ImGui::TextWrapped("Requires nonnegative electron density in e/A^3. Samples at or below the floor are masked to zero, avoiding the divergent gradient/density ratio near vacuum.");
            if (m_operation==34) ImGui::TextWrapped("Sign(lambda_2)*rho: negative values mark bonding/hydrogen-bonding accumulation, positive values mark steric repulsion.");
            if (m_operation==35) ImGui::TextWrapped("Density Overlap Regions Indicator, in [0,1); values approaching 1 mark boundaries between density basins.");
        }
        if (m_operation==36)
        {
            ImGui::TextWrapped("One density threshold per line (e/A^3). Finite (non-periodic) grids only.");
            ImGui::InputTextMultiline("##thresholds",m_thresholds,sizeof(m_thresholds),ImVec2(-1,80));
        }
        if (m_operation==37) ImGui::TextWrapped("Native QTAIM/Bader on-grid steepest-ascent partitioning (requires nonnegative e/A^3). Adds a basin-index field and a per-basin charge/volume/maximum table.");
        ImGui::Spacing();
        if (dialogLayout::primaryButton("Calculate",ImVec2(-FLT_MIN,0)) || (m_generateSurface && !m_task.running()))
        {
            const Grid source = g, reference = referenceVolume.fields[m_reference];
            const auto sites = m_volume.sites;
            const int op = m_generateSurface ? 22 : m_operation;
            m_generateSurface = false;
            const int count = m_count, axis = m_axis, window = m_window, radius = m_radius;
            const double maskLow=m_maskLow, maskHigh=m_maskHigh, referenceWeight=m_referenceWeight;
            const int booleanOperation=m_booleanOperation;
            const double scalar = op == 22 ? m_surfaceLevel : m_scalar, sigma = m_sigma, alpha = m_alpha, realCutoff = m_realCutoff, reciprocalCutoff = m_reciprocalCutoff;
            const auto start = vector(m_start), end = vector(m_end), u = vector(m_u), v = vector(m_v);
            const std::string reflectionText(m_reflections), chargeText(m_charges), thresholdText(m_thresholds);
            const double bondingFloor = m_bondingFloor;
            const bool colored = m_colorSurface;
            const bool appendSurface = m_appendSurface;
            const bool existingSurface = !m_surface.vertices.empty();
            const std::string existingColorUnit = m_colorUnit;
            m_error.clear();
            std::ostringstream provenance;
            provenance << "Operation " << op << ", input " << source.name << ", reference " << reference.name;
            for (const auto& parameter:snapshot(false).settings) provenance << "; " << parameter.first << '=' << parameter.second;
            m_pendingOperation=provenance.str();
            m_task.start([=]
            {
                if (op==22 && !displayRange(source.values).contains(scalar))
                    throw std::invalid_argument("3D isovalue is outside the selected field's range.");
                Output out;
                out.volume.sites = sites;
                auto add = [&](Grid f) { out.volume.fields.push_back(std::move(f)); };
                auto profile = [&](const std::vector<glm::dvec2>& rows)
                { out.columns = 2; out.heading = "distance_A,value"; for (auto row : rows) out.table.insert(out.table.end(),{row.x,row.y}); };
                if (op == 0) { out.columns = 1; out.heading = "integral (field unit * A^3)"; out.table = {integrate(source)}; }
                else if (op >= 1 && op <= 4) add(arithmetic(source,reference,std::array<const char*,4>{"add","subtract","multiply","divide"}[op-1]));
                else if (op == 5) add(scale(source,scalar));
                else if (op == 25) add(densityDifference(source,reference,referenceWeight));
                else if (op == 26) add(thresholdMask(source,maskLow,maskHigh));
                else if (op == 27) add(booleanMask(source,reference,std::array<const char*,4>{"union","intersection","difference","xor"}[booleanOperation]));
                else if (op == 28) add(applyMask(source,reference));
                else if (op == 29) for (auto& f : splitDensity(source)) add(std::move(f));
                else if (op == 30)
                {
                    out.columns=3; out.heading="accumulation_e,depletion_e,net_e";
                    const auto q=chargeSummary(source); out.table.assign(q.begin(),q.end());
                }
                else if (op == 31) { profile(cumulativeCharge(source,axis)); out.heading="distance_A,cumulative_e"; }
                else if (op == 32) add(booleanMask(thresholdMask(source,0,1),source,"difference"));
                else if (op == 6) add(smooth(source,sigma,radius));
                else if (op == 7) for (auto& f : gradient(source)) add(std::move(f));
                else if (op == 8) add(laplacian(source));
                else if (op == 9) for (auto& f : energyDensity(source)) add(std::move(f));
                else if (op == 10) profile(lineProfile(source,start,end,count));
                else if (op == 11) profile(planarAverage(source,axis));
                else if (op == 12) profile(macroscopicAverage(source,axis,window));
                else if (op == 13) add(section(source,start,u,v,count,count));
                else if (op == 14)
                {
                    out.columns = 3; out.heading = "x_A,y_A,z_A (successive pairs are segment endpoints)";
                    for (auto p : contours(source,scalar)) out.table.insert(out.table.end(),{p.x,p.y,p.z});
                }
                else if (op == 15)
                {
                    out.columns = 4; out.heading = "x_A,y_A,z_A,value";
                    for (auto p : peaks(source)) out.table.insert(out.table.end(),{p.position.x,p.position.y,p.position.z,p.value});
                }
                else if (op == 16)
                {
                    std::vector<glm::dvec3> positions;
                    for (auto s : sites) positions.push_back(s.position);
                    out.columns = 2; out.heading = "integral,volume_A^3 (atom order)";
                    for (auto b : voronoiIntegrate(source,positions)) out.table.insert(out.table.end(),{b.integral,b.volume});
                }
                else if (op == 17) { out.columns = 1; out.heading = "sphere_integral"; out.table = {integrateSphere(source,start,scalar)}; }
                else if (op == 18 || op == 19)
                {
                    std::istringstream input(reflectionText);
                    std::vector<glm::ivec3> indices;
                    std::vector<Reflection> reflections;
                    for (std::string row; std::getline(input,row);)
                    {
                        if (row.find_first_not_of(" \t\r") == std::string::npos) continue;
                        std::istringstream fields(row);
                        glm::ivec3 h;
                        double real = 0, imaginary = 0;
                        if (!(fields >> h.x >> h.y >> h.z) || (op == 19 && !(fields >> real >> imaginary))) throw std::invalid_argument("Malformed reflection row");
                        std::string extra;
                        if (fields >> extra) throw std::invalid_argument("Extra reflection columns");
                        indices.push_back(h); reflections.push_back({h,{real,imaginary}});
                    }
                    if (op == 19) add(fourierSynthesis(source,reflections));
                    else
                    {
                        out.columns = 5; out.heading = "h,k,l,real,imag";
                        for (const auto& r : structureFactors(source,indices)) out.table.insert(out.table.end(),{static_cast<double>(r.hkl.x),static_cast<double>(r.hkl.y),static_cast<double>(r.hkl.z),r.value.real(),r.value.imag()});
                    }
                }
                else if (op == 20) add(patterson(source));
                else if (op == 21)
                {
                    if (!source.periodic) throw std::invalid_argument("Ewald requires periodic geometry");
                    std::istringstream input(chargeText);
                    std::vector<double> charges;
                    for (double q; input >> q;) charges.push_back(q);
                    if (!input.eof()) throw std::invalid_argument("Invalid charge list");
                    std::vector<glm::dvec3> positions;
                    for (auto s : sites) positions.push_back(s.position);
                    const auto e = ewald(source.cell,positions,charges,alpha,realCutoff,reciprocalCutoff);
                    out.columns = 1; out.heading = "cell_energy_eV_then_site_potentials_V"; out.table = {e.energy};
                    out.table.insert(out.table.end(),e.potentials.begin(),e.potentials.end());
                }
                else if (op == 22)
                {
                    const auto& colorGrid = colored ? reference : source;
                    if (appendSurface && existingSurface && existingColorUnit!=colorGrid.unit)
                        throw std::invalid_argument("Appended surfaces must use the same coloring units.");
                    out.mesh = isosurface(source,scalar,colored ? &reference : nullptr);
                    out.appendMesh = appendSurface; out.surfaceReady = true;
                    const auto range = std::minmax_element(colorGrid.values.begin(),colorGrid.values.end());
                    out.colorLow=*range.first; out.colorHigh=*range.second; out.colorUnit=colorGrid.unit;
                }
                else if (op == 23) add(resample(source,reference));
                else if (op == 24)
                {
                    add(periodicEndpoints(source));
                }
                else if (op == 33) add(reducedDensityGradient(source,bondingFloor));
                else if (op == 34) add(signedDensity(source,bondingFloor));
                else if (op == 35) add(dori(source,bondingFloor));
                else if (op == 36)
                {
                    std::istringstream input(thresholdText);
                    std::vector<double> thresholds;
                    for (double t; input >> t;) thresholds.push_back(t);
                    if (thresholds.empty()) throw std::invalid_argument("Supply at least one threshold");
                    out.columns = 4; out.heading = "threshold,betti0,betti1,betti2";
                    for (auto row : bettiCurve(source,thresholds)) out.table.insert(out.table.end(),{row.x,row.y,row.z,row.w});
                }
                else if (op == 37)
                {
                    const auto partition = baderOnGrid(source);
                    Grid basins = source;
                    basins.name = "Bader basin index"; basins.unit = "basin index";
                    for (std::size_t i=0;i<basins.values.size();++i) basins.values[i] = partition.basin[i];
                    basins.validate();
                    add(std::move(basins));
                    out.columns = 6; out.heading = "basin,charge_e,volume_A3,max_x_A,max_y_A,max_z_A";
                    for (std::size_t i=0;i<partition.basins.size();++i)
                    {
                        const auto& basin = partition.basins[i];
                        out.table.insert(out.table.end(),{static_cast<double>(i),basin.charge,basin.volume,basin.maximum.x,basin.maximum.y,basin.maximum.z});
                    }
                }
                return out;
            });
        }
        ImGui::Spacing();
        if (ImGui::CollapsingHeader("Export",ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (combo("Format",&m_exportFormat,"Grid: XSF\0Grid: Cube\0Grid: VASP\0Table: CSV\0Surface: OBJ\0Surface: PLY\0View: 3D PNG\0View: 2D SVG\0"))
            {
                auto output = std::filesystem::u8path(m_output);
                output.replace_extension(std::array<const char*,8>{".xsf",".cube",".vasp",".csv",".obj",".ply",".png",".svg"}[m_exportFormat]);
                std::snprintf(m_output,sizeof(m_output),"%s",output.u8string().c_str());
            }
            if (responsive::button("Save as...",ImVec2(-FLT_MIN,0)))
            {
                auto output = std::filesystem::u8path(m_output);
                if (!output.has_parent_path() && !m_loadedPath.empty())
                    output = std::filesystem::u8path(m_loadedPath).parent_path()/output;
                m_pickerAction=2;
                m_picker.open("Save electronic result",true,output.u8string());
            }
        }
    }
    ImGui::EndDisabled();
    if (m_task.running())
    {
        ImGui::TextUnformatted("Calculating...");
        if (m_task.progress()>=0) ImGui::ProgressBar(m_task.progress(),ImVec2(-1,0),"Current calculation stage");
        if (responsive::button("Cancel calculation")) m_task.cancel();
    }
    if (!m_pendingDrops.empty()) ImGui::Text("Queued files: %zu",m_pendingDrops.size());
    if (!m_error.empty()) ImGui::TextWrapped("%s",m_error.c_str());
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Appearance"))
    {
        if (combo("Preview quality",&m_quality,"Fast\0Balanced\0High\0Very high\0"))
        {
            m_viewport.setQuality(m_quality);
            m_sliceViewport.setResolution(std::array<int,4>{65,129,257,513}[m_quality]);
            m_volumeDirty=true;
        }
        combo("Colors",&m_palette,"Spectrum\0Blue-white-red\0Sequential blue\0");
        if (ImGui::Checkbox("Automatic range",&m_autoRange) && m_autoRange) { m_colorLow=m_autoLow; m_colorHigh=m_autoHigh; }
        if (!m_autoRange)
        {
            inputFloat("Minimum",&m_colorLow,0,0,"%.5g");
            inputFloat("Maximum",&m_colorHigh,0,0,"%.5g");
        }
        if (ImGui::TreeNode("Lighting"))
        {
            slider("Specular",&m_specular,0,1,"%.2f");
            slider("Shininess",&m_shininess,4,128,"%.0f");
            ImGui::TreePop();
        }
    }
    ImGui::EndChild();
    responsive::nextPanel(stackPanels);
    responsive::beginChild("Electronic results",responsive::size(0,stackPanels ? 620 : 0),false,ImGuiWindowFlags_NoScrollWithMouse|ImGuiWindowFlags_NoScrollbar);
    if (!m_result.table.empty() && ImGui::CollapsingHeader("Results table",ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::TextWrapped("%s",m_result.heading.c_str());
        if (m_result.heading == "distance_A,value")
        {
            std::vector<float> values;
            for (std::size_t i = 1; i < m_result.table.size(); i += 2) values.push_back(static_cast<float>(m_result.table[i]));
            ImGui::PlotLines("Profile",values.data(),static_cast<int>(values.size()),0,nullptr,FLT_MAX,FLT_MAX,ImVec2(-1,140));
        }
        else if (m_result.heading == "threshold,betti0,betti1,betti2")
        {
            std::vector<float> betti0,betti1,betti2;
            for (std::size_t i = 0; i+3 < m_result.table.size(); i += 4)
            {
                betti0.push_back(static_cast<float>(m_result.table[i+1]));
                betti1.push_back(static_cast<float>(m_result.table[i+2]));
                betti2.push_back(static_cast<float>(m_result.table[i+3]));
            }
            ImGui::PlotLines("Connected components (betti0)",betti0.data(),static_cast<int>(betti0.size()),0,nullptr,0,FLT_MAX,ImVec2(-1,80));
            ImGui::PlotLines("Loops (betti1)",betti1.data(),static_cast<int>(betti1.size()),0,nullptr,0,FLT_MAX,ImVec2(-1,80));
            ImGui::PlotLines("Enclosed cavities (betti2)",betti2.data(),static_cast<int>(betti2.size()),0,nullptr,0,FLT_MAX,ImVec2(-1,80));
        }
        responsive::beginChild("Electronic table",responsive::size(0,100),true);
        for (std::size_t i = 0; i < std::min<std::size_t>(m_result.table.size(),m_result.columns * 200); i += m_result.columns)
        {
            std::ostringstream row;
            row << std::setprecision(9);
            for (int j = 0; j < m_result.columns; ++j) row << m_result.table[i+j] << "   ";
            ImGui::TextUnformatted(row.str().c_str());
        }
        ImGui::EndChild();
    }
    drawPreview();
    ImGui::EndChild();
    if (auto path=m_picker.draw())
    {
        try {
            if (m_pickerAction==4) atomforge::saveWorkspace({snapshot()},*path);
            else if (m_pickerAction==3) { auto saved=atomforge::loadWorkspace(*path); remember(); restore(saved.front()); }
            else if (m_pickerAction==2) save(*path);
            else load(*path,m_pickerAction==1);
        } catch (const std::exception& error) { m_error=error.what(); }
    }
    ImGui::End();
}

void ElectronicPostProcessingDialog::load(const std::string& path, bool reference)
{
    const std::string q=std::array<const char*,5>{"auto","density","potential","elf","raw"}[m_quantity];
    const std::string units=m_cubeUnits ? "angstrom" : "bohr";
    m_error.clear();
    m_pendingOperation=std::string(reference ? "Load reference: " : "Load: ")+path;
    if (!m_task.start([path,q,units,reference]
    {
        Output out; out.volume=loadVolume(path,q,units);
        out.loaded=!reference; out.referenceLoaded=reference; out.sourcePath=path; return out;
    })) m_error=m_task.error();
}

void ElectronicPostProcessingDialog::save(const std::string& path)
{
    try
    {
        if (m_volume.fields.empty()) throw std::invalid_argument("Load a volume first.");
        if (m_exportFormat<3) saveVolume({m_volume.sites,{m_volume.fields.at(m_selected)}},path,std::array<const char*,3>{"xsf","cube","vasp"}[m_exportFormat]);
        else if (m_exportFormat==6)
        {
            constexpr int width=1600,height=1200;
            if (m_renderMode==0) {
                m_viewport.setVolume(m_volume.fields.at(m_selected));
                m_viewport.renderVolume(width,height,m_yaw,m_pitch,m_zoom,m_pan,m_opacity,m_surfaceLevel,m_sliceLow,m_sliceHigh,m_palette);
            } else {
                if (m_surface.vertices.empty()) throw std::invalid_argument("Calculate an isosurface first");
                m_viewport.render(width,height,m_yaw,m_pitch,m_zoom,m_pan,m_opacity,m_colorLow,m_colorHigh,m_palette,m_specular,m_shininess);
            }
            const auto pixels=m_viewport.pixels();
            std::ofstream stream(std::filesystem::u8path(path),std::ios::binary);
            const auto write=[](void* context,void* data,int size) {
                static_cast<std::ostream*>(context)->write(static_cast<const char*>(data),size);
            };
            const bool written=stbi_write_png_to_func(write,&stream,width,height,4,pixels.data(),width*4)!=0;
            stream.close();
            if (!written || !stream) throw std::runtime_error("Failed to write viewport PNG");
        }
        else if (m_exportFormat==7) m_sliceViewport.saveSvg(m_volume.fields.at(m_selected),path,m_sliceLow,m_sliceHigh,m_palette);
        else if (m_exportFormat>=4)
        {
            if (m_surface.vertices.empty()) throw std::invalid_argument("Calculate an isosurface first.");
            saveMesh(m_surface,path);
        }
        else
        {
            if (!m_result.columns) throw std::invalid_argument("No table to export");
            std::ofstream stream{std::filesystem::u8path(path)};
            if (!stream) throw std::runtime_error("Cannot open CSV output");
            stream << "# " << m_result.heading << '\n' << std::setprecision(17);
            for (std::size_t i=0;i<m_result.table.size();++i) stream << m_result.table[i] << ((i+1)%m_result.columns ? ',' : '\n');
            if (!stream) throw std::runtime_error("Failed writing CSV");
        }
        std::snprintf(m_output,sizeof(m_output),"%s",path.c_str());
        m_error="Saved "+std::filesystem::u8path(path).filename().u8string();
    }
    catch (const std::exception& e) { m_error=e.what(); }
}

void ElectronicPostProcessingDialog::resetCamera()
{
    m_yaw=.6f; m_pitch=.4f; m_zoom=1; m_pan=glm::vec2(0);
}

void ElectronicPostProcessingDialog::drawPreview()
{
    ImGui::SetNextItemWidth(responsive::dp(-1));
    ImGui::Combo("##view layout",&m_viewLayout,"3D and 2D\0Only 3D\0Only 2D\0");
    if (m_volume.fields.empty())
    {
        ImGui::TextWrapped("Open or drop a VASP, Cube or XSF file to display both views.");
        return;
    }
    const auto& grid=m_volume.fields[m_selected];
    if (m_sliceField!=m_selected)
    {
        const auto range=displayRange(grid.values);
        m_sliceLow=static_cast<float>(range.low); m_sliceHigh=static_cast<float>(range.high);
        m_surfaceLevel=m_suggestedLevel=static_cast<float>(range.suggested);
        m_volumeDirty=true;
        m_sliceField=m_selected;
        m_sliceViewport.invalidate();
        if (m_surface.vertices.empty())
        {
            m_autoLow=m_sliceLow; m_autoHigh=m_sliceHigh; m_colorUnit=grid.unit;
            if (m_autoRange) { m_colorLow=m_autoLow; m_colorHigh=m_autoHigh; }
        }
    }
    const float width=ImGui::GetContentRegionAvail().x;
    const float viewWidth=m_viewLayout==0 ? (width-ImGui::GetStyle().ItemSpacing.x)*.5f : width;
    const auto flags=ImGuiWindowFlags_NoScrollbar|ImGuiWindowFlags_NoScrollWithMouse;
    if (m_viewLayout!=2)
    {
        responsive::beginChild("3D view",ImVec2(viewWidth,0),true,flags);
        ImGui::TextUnformatted("3D view");
        ImGui::SetNextItemWidth(responsive::dp(-1));
        ImGui::Combo("##render mode",&m_renderMode,"Volume (interior density)\0Isosurface (boundary)\0");
        ImGui::BeginDisabled(m_task.running());
        ImGui::TextUnformatted(m_renderMode==0 ? "Density threshold" : "3D isovalue");
        ImGui::SetNextItemWidth(responsive::dp(-1));
        ImGui::InputFloat("##surface level",&m_surfaceLevel,0,0,"%.5g");
        if (responsive::button("Estimate level")) m_surfaceLevel=m_suggestedLevel;
        ImGui::TextWrapped("Field range: %.5g to %.5g %s",m_sliceLow,m_sliceHigh,grid.unit.c_str());
        const bool invalid=!std::isfinite(m_surfaceLevel) || m_surfaceLevel<m_sliceLow || m_surfaceLevel>m_sliceHigh;
        if (invalid)
        {
            ImGui::PushStyleColor(ImGuiCol_Text,ImVec4(.85f,.32f,.08f,1));
            ImGui::TextWrapped("Warning: value is outside the field range. Choose a value within the range above.");
            ImGui::PopStyleColor();
        }
        else if (m_sliceLow==m_sliceHigh) ImGui::TextWrapped("Constant field: no distinct isosurface exists.");
        if (m_renderMode==1)
        {
            ImGui::BeginDisabled(invalid || m_sliceLow==m_sliceHigh);
            if (responsive::button("Update surface",ImVec2(-FLT_MIN,0))) m_generateSurface=true;
            ImGui::EndDisabled();
        }
        ImGui::EndDisabled();
        ImGui::TextUnformatted("Transparency");
        ImGui::SetNextItemWidth(responsive::dp(-1));
        float transparency=1-m_opacity;
        if (ImGui::SliderFloat("##transparency",&transparency,0,1,"%.2f")) m_opacity=1-transparency;
        ImGui::Checkbox("Show slice plane",&m_showSlicePlane);
        draw3DPreview();
        ImGui::EndChild();
    }
    if (m_viewLayout==0) ImGui::SameLine();
    if (m_viewLayout!=1)
    {
        responsive::beginChild("2D view",ImVec2(viewWidth,0),true,ImGuiWindowFlags_NoScrollWithMouse);
        try
        {
            // The section's density scale remains independent of a potential
            // mapped onto the 3D surface.
            m_sliceViewport.draw(grid,m_sliceLow,m_sliceHigh,m_palette);
        }
        catch (const std::exception& error) { ImGui::TextWrapped("%s",error.what()); }
        ImGui::EndChild();
    }
}

void ElectronicPostProcessingDialog::draw3DPreview()
{
    if (responsive::button("Fit view")) resetCamera();
    ImGui::TextWrapped("Drag: orbit | Right drag: pan | Wheel: zoom");
    const auto pos=ImGui::GetCursorScreenPos();
    const ImVec2 size(std::max(100.0f,ImGui::GetContentRegionAvail().x),std::max(150.0f,ImGui::GetContentRegionAvail().y-52));
    ImGui::InvisibleButton("Electronic viewport",size,ImGuiButtonFlags_MouseButtonLeft|ImGuiButtonFlags_MouseButtonRight);
    if (ImGui::IsItemHovered())
    {
        m_zoom=std::clamp(m_zoom*std::exp(ImGui::GetIO().MouseWheel*.12f),.15f,15.0f);
        if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) resetCamera();
    }
    if (ImGui::IsItemActive())
    {
        const auto delta=ImGui::GetIO().MouseDelta;
        if (ImGui::IsMouseDown(ImGuiMouseButton_Left))
        {
            m_yaw=std::remainder(m_yaw+delta.x*.008f,6.2831853f);
            m_pitch=std::clamp(m_pitch+delta.y*.008f,-1.55f,1.55f);
        }
        if (ImGui::IsMouseDown(ImGuiMouseButton_Right)) m_pan+=glm::vec2(delta.x,-delta.y)*(2.5f/(size.y*m_zoom));
    }
    auto* draw=ImGui::GetWindowDrawList();
    draw->PushClipRect(pos,ImVec2(pos.x+size.x,pos.y+size.y),true);
    draw->AddRectFilled(pos,ImVec2(pos.x+size.x,pos.y+size.y),IM_COL32(245,247,251,255));
    if (m_renderMode==0 || !m_surface.vertices.empty())
    {
        try
        {
            if (!std::isfinite(m_colorLow) || !std::isfinite(m_colorHigh) || m_colorHigh<m_colorLow)
                throw std::invalid_argument("Color range must be finite with maximum >= minimum.");
            const auto scale=ImGui::GetIO().DisplayFramebufferScale;
            GLuint texture;
            if (m_renderMode==0)
            {
                if (m_volumeDirty) { m_viewport.setVolume(m_volume.fields[m_selected]); m_volumeDirty=false; }
                const float threshold=std::isfinite(m_surfaceLevel) ? m_surfaceLevel : m_sliceHigh;
                texture=m_viewport.renderVolume(static_cast<int>(size.x*scale.x),static_cast<int>(size.y*scale.y),m_yaw,m_pitch,m_zoom,m_pan,m_opacity,threshold,m_sliceLow,m_sliceHigh,m_palette);
            }
            else texture=m_viewport.render(static_cast<int>(size.x*scale.x),static_cast<int>(size.y*scale.y),m_yaw,m_pitch,m_zoom,m_pan,m_opacity,m_colorLow,m_colorHigh,m_palette,m_specular,m_shininess);
            draw->AddImage((ImTextureID)(intptr_t)texture,pos,ImVec2(pos.x+size.x,pos.y+size.y),ImVec2(0,1),ImVec2(1,0));
            if (m_showSlicePlane)
            {
                // A translucent overlay keeps the entire guide readable through
                // dense volumes and surfaces. Its geometry follows the 3D camera.
                try
                {
                    const auto plane=m_sliceViewport.plane(m_volume.fields[m_selected]);
                    std::vector<ImVec2> points;
                    for (auto p : plane.boundary)
                    {
                        const auto q=m_viewport.project(p,m_renderMode==0,size.x/size.y,m_yaw,m_pitch,m_zoom,m_pan);
                        points.emplace_back(pos.x+(q.x+1)*size.x*.5f,pos.y+(1-q.y)*size.y*.5f);
                    }
                    if (points.size()>=3)
                    {
                        // ImGui convex fills require clockwise screen winding.
                        float area=0;
                        for (std::size_t i=0;i<points.size();++i)
                        {
                            const auto a=points[i], b=points[(i+1)%points.size()];
                            area+=a.x*b.y-b.x*a.y;
                        }
                        if (area<0) std::reverse(points.begin(),points.end());
                        draw->AddConvexPolyFilled(points.data(),static_cast<int>(points.size()),IM_COL32(65,135,200,48));
                        draw->AddPolyline(points.data(),static_cast<int>(points.size()),IM_COL32(35,95,155,220),ImDrawFlags_Closed,1.5f);
                    }
                }
                catch (const std::invalid_argument&) { /* Invalid plane inputs are reported in the 2D controls. */ }
            }
            const std::string caption=m_renderMode==0 ? "Volume density" : std::to_string(m_surface.vertices.size()/3)+" triangles";
            draw->AddText(ImVec2(pos.x+12,pos.y+12),IM_COL32(70,80,95,255),caption.c_str());
        }
        catch (const std::exception& e) { m_error=e.what(); }
    }
    else draw->AddText(ImVec2(pos.x+20,pos.y+25),IM_COL32(90,100,115,255),"Set the 3D isovalue and click Update surface.");
    draw->PopClipRect();
    const auto bar=ImGui::GetCursorScreenPos();
    const float width=std::max(100.0f,size.x-4);
    for(int i=0;i<128;++i)
    {
        const auto c=ElectronicViewport::color(static_cast<float>(i)/127,m_palette);
        draw->AddRectFilled(ImVec2(bar.x+width*i/128,bar.y),ImVec2(bar.x+width*(i+1)/128,bar.y+10),ImGui::ColorConvertFloat4ToU32(ImVec4(c.x,c.y,c.z,1)));
    }
    ImGui::Dummy(ImVec2(width,12));
    ImGui::Text("%.5g",m_renderMode==0 ? m_sliceLow : m_colorLow); ImGui::SameLine(width*.4f);
    ImGui::TextUnformatted(m_renderMode==0 ? m_volume.fields[m_selected].unit.c_str() : m_colorUnit.c_str());
    ImGui::SameLine(std::max(width-85.0f,150.0f)); ImGui::Text("%.5g",m_renderMode==0 ? m_sliceHigh : m_colorHigh);
}
