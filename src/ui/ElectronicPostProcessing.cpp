#include "ui/ElectronicPostProcessing.h"
#include "imgui.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
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

void ElectronicPostProcessingDialog::drawMenuItem()
{
    if (ImGui::MenuItem("Electronic Post-processing...")) m_open = true;
}

void ElectronicPostProcessingDialog::drawDialog()
{
    if (m_task.poll())
    {
        m_error = m_task.error();
        if (m_task.result())
        {
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
                m_loadedPath=m_result.sourcePath; m_sliceField=-1;
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
    if (!m_open) return;
    ImGui::SetNextWindowSize(ImVec2(1080,700),ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(780,480),ImVec2(FLT_MAX,FLT_MAX));
    if (!ImGui::Begin("Electronic Post-processing",&m_open,ImGuiWindowFlags_NoScrollbar|ImGuiWindowFlags_NoCollapse)) { ImGui::End(); return; }
    const float sidebar = std::clamp(ImGui::GetContentRegionAvail().x * .38f,340.0f,440.0f);
    ImGui::BeginChild("Electronic controls",ImVec2(sidebar,0),true);
    ImGui::BeginDisabled(m_task.running());
    const float buttonWidth = (ImGui::GetContentRegionAvail().x-ImGui::GetStyle().ItemSpacing.x)*.5f;
    if (ImGui::Button("Open volume...",ImVec2(buttonWidth,0))) { m_pickerAction=0; m_picker.open("Open electronic volume",false,m_loadedPath); }
    ImGui::SameLine();
    if (ImGui::Button("Reference...",ImVec2(buttonWidth,0))) { m_pickerAction=1; m_picker.open("Open reference volume",false,m_referencePath.empty() ? m_loadedPath : m_referencePath); }
    ImGui::BeginDisabled(m_loadedPath.empty());
    if (ImGui::Button("Reload volume",ImVec2(buttonWidth,0))) load(m_loadedPath,false);
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (ImGui::Button("Clear",ImVec2(buttonWidth,0)))
    {
        m_volume={}; m_referenceVolume={}; m_result={}; m_surface={};
        m_viewport.setMesh(m_surface); m_selected=m_reference=0; m_sliceField=-1; m_loadedPath.clear(); m_referencePath.clear(); resetCamera();
    }
    filename("Volume",m_loadedPath);
    filename("Reference",m_referencePath);
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
        selectField("Field",m_selected,m_volume);
        const auto& referenceVolume = m_referenceVolume.fields.empty() ? m_volume : m_referenceVolume;
        m_reference = std::clamp(m_reference,0,static_cast<int>(referenceVolume.fields.size())-1);
        const Grid& g = m_volume.fields[m_selected];
        ImGui::TextDisabled("%d x %d x %d | %s",g.shape[0],g.shape[1],g.shape[2],g.unit.c_str());
        combo("Tool",&m_operation,
            "Total integral\0Add reference\0Subtract reference (difference density)\0Multiply reference\0Divide by reference\0Scale\0Gaussian smoothing\0Cartesian gradient\0Laplacian\0Energy-density conversion\0Line profile\0Planar average\0Macroscopic average\0Plane section\0Contour segments (z=0)\0Peak search\0Voronoi site integration\0Sphere integration\0Structure factors\0Fourier synthesis\0Patterson density\0Ewald site potentials\0Isosurface\0Resample onto reference\0Verify periodic endpoint planes\0");
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
        if (m_operation == 14 || m_operation == 22) inputFloat("Isovalue",&m_scalar);
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
        ImGui::Spacing();
        if (ImGui::Button("Calculate",ImVec2(-FLT_MIN,0)))
        {
            const Grid source = g, reference = referenceVolume.fields[m_reference];
            const auto sites = m_volume.sites;
            const int op = m_operation, count = m_count, axis = m_axis, window = m_window, radius = m_radius;
            const double scalar = m_scalar, sigma = m_sigma, alpha = m_alpha, realCutoff = m_realCutoff, reciprocalCutoff = m_reciprocalCutoff;
            const auto start = vector(m_start), end = vector(m_end), u = vector(m_u), v = vector(m_v);
            const std::string reflectionText(m_reflections), chargeText(m_charges);
            const bool colored = m_colorSurface;
            const bool appendSurface = m_appendSurface;
            const bool existingSurface = !m_surface.vertices.empty();
            const std::string existingColorUnit = m_colorUnit;
            m_error.clear();
            m_task.start([=]
            {
                Output out;
                out.volume.sites = sites;
                auto add = [&](Grid f) { out.volume.fields.push_back(std::move(f)); };
                auto profile = [&](const std::vector<glm::dvec2>& rows)
                { out.columns = 2; out.heading = "distance_A,value"; for (auto row : rows) out.table.insert(out.table.end(),{row.x,row.y}); };
                if (op == 0) { out.columns = 1; out.heading = "integral (field unit * A^3)"; out.table = {integrate(source)}; }
                else if (op >= 1 && op <= 4) add(arithmetic(source,reference,std::array<const char*,4>{"add","subtract","multiply","divide"}[op-1]));
                else if (op == 5) add(scale(source,scalar));
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
                return out;
            });
        }
        ImGui::Spacing();
        if (ImGui::CollapsingHeader("Export",ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (combo("Format",&m_exportFormat,"Grid: XSF\0Grid: Cube\0Grid: VASP\0Table: CSV\0Surface: OBJ\0Surface: PLY\0"))
            {
                auto output = std::filesystem::u8path(m_output);
                output.replace_extension(std::array<const char*,6>{".xsf",".cube",".vasp",".csv",".obj",".ply"}[m_exportFormat]);
                std::snprintf(m_output,sizeof(m_output),"%s",output.u8string().c_str());
            }
            if (ImGui::Button("Save as...",ImVec2(-FLT_MIN,0)))
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
    if (m_task.running()) ImGui::TextUnformatted("Calculating...");
    if (!m_error.empty()) ImGui::TextWrapped("%s",m_error.c_str());
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Appearance",ImGuiTreeNodeFlags_DefaultOpen))
    {
        combo("Colors",&m_palette,"Spectrum\0Blue-white-red\0Sequential blue\0");
        if (ImGui::Checkbox("Automatic range",&m_autoRange) && m_autoRange) { m_colorLow=m_autoLow; m_colorHigh=m_autoHigh; }
        if (!m_autoRange)
        {
            inputFloat("Minimum",&m_colorLow,0,0,"%.5g");
            inputFloat("Maximum",&m_colorHigh,0,0,"%.5g");
        }
        slider("Opacity",&m_opacity,.05f,1,"%.2f");
        if (ImGui::TreeNode("Lighting"))
        {
            slider("Specular",&m_specular,0,1,"%.2f");
            slider("Shininess",&m_shininess,4,128,"%.0f");
            ImGui::TreePop();
        }
    }
    ImGui::EndChild();
    ImGui::SameLine();
    ImGui::BeginChild("Electronic results",ImVec2(0,0),false,ImGuiWindowFlags_NoScrollWithMouse|ImGuiWindowFlags_NoScrollbar);
    if (!m_result.table.empty() && ImGui::CollapsingHeader("Results table",ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::TextWrapped("%s",m_result.heading.c_str());
        if (m_result.heading == "distance_A,value")
        {
            std::vector<float> values;
            for (std::size_t i = 1; i < m_result.table.size(); i += 2) values.push_back(static_cast<float>(m_result.table[i]));
            ImGui::PlotLines("Profile",values.data(),static_cast<int>(values.size()),0,nullptr,FLT_MAX,FLT_MAX,ImVec2(-1,140));
        }
        ImGui::BeginChild("Electronic table",ImVec2(0,100),true);
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
        if (m_pickerAction==2) save(*path);
        else load(*path,m_pickerAction==1);
    }
    ImGui::End();
}

void ElectronicPostProcessingDialog::load(const std::string& path, bool reference)
{
    const std::string q=std::array<const char*,5>{"auto","density","potential","elf","raw"}[m_quantity];
    const std::string units=m_cubeUnits ? "angstrom" : "bohr";
    m_error.clear();
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
    if (ImGui::Button("Fit view")) resetCamera();
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
    if (!m_surface.vertices.empty())
    {
        try
        {
            if (!std::isfinite(m_colorLow) || !std::isfinite(m_colorHigh) || m_colorHigh<m_colorLow)
                throw std::invalid_argument("Color range must be finite with maximum >= minimum.");
            const auto scale=ImGui::GetIO().DisplayFramebufferScale;
            const auto texture=m_viewport.render(static_cast<int>(size.x*scale.x),static_cast<int>(size.y*scale.y),m_yaw,m_pitch,m_zoom,m_pan,m_opacity,m_colorLow,m_colorHigh,m_palette,m_specular,m_shininess);
            draw->AddImage((ImTextureID)(intptr_t)texture,pos,ImVec2(pos.x+size.x,pos.y+size.y),ImVec2(0,1),ImVec2(1,0));
            const std::string caption=std::to_string(m_surface.vertices.size()/3)+" triangles";
            draw->AddText(ImVec2(pos.x+12,pos.y+12),IM_COL32(70,80,95,255),caption.c_str());
        }
        catch (const std::exception& e) { m_error=e.what(); }
    }
    else if (!m_volume.fields.empty())
    {
        const auto& g=m_volume.fields[m_selected];
        if (m_sliceField!=m_selected)
        {
            const auto range=std::minmax_element(g.values.begin(),g.values.end());
            m_sliceLow=static_cast<float>(*range.first); m_sliceHigh=static_cast<float>(*range.second); m_sliceField=m_selected;
        }
        m_autoLow=m_sliceLow; m_autoHigh=m_sliceHigh; m_colorUnit=g.unit;
        if(m_autoRange) { m_colorLow=m_autoLow; m_colorHigh=m_autoHigh; }
        const int z=g.name=="section" ? 0 : g.shape[2]/2;
        const int stride=std::max(1,(std::max(g.shape[0],g.shape[1])+255)/256);
        for(int y=0;y<g.shape[1];y+=stride) for(int x=0;x<g.shape[0];x+=stride)
        {
            const float t=m_colorHigh>m_colorLow ? static_cast<float>((g.values[g.index(x,y,z)]-m_colorLow)/(m_colorHigh-m_colorLow)) : .5f;
            const auto color=ElectronicViewport::color(t,m_palette);
            draw->AddRectFilled(ImVec2(pos.x+size.x*x/g.shape[0],pos.y+size.y*(1.0f-static_cast<float>(std::min(y+stride,g.shape[1]))/g.shape[1])),
                ImVec2(pos.x+size.x*std::min(x+stride,g.shape[0])/g.shape[0],pos.y+size.y*(1.0f-static_cast<float>(y)/g.shape[1])),
                ImGui::ColorConvertFloat4ToU32(ImVec4(color.x,color.y,color.z,1)));
        }
        const char* caption = "2D grid slice\nCalculate an isosurface to orbit in 3D";
        const float captionHeight = ImGui::CalcTextSize(caption,nullptr,false,size.x-24).y+24;
        draw->AddRectFilled(pos,ImVec2(pos.x+size.x,pos.y+captionHeight),IM_COL32(245,247,251,255));
        draw->AddText(nullptr,0,ImVec2(pos.x+12,pos.y+12),IM_COL32(30,40,55,255),caption,nullptr,size.x-24);
    }
    else draw->AddText(ImVec2(pos.x+20,pos.y+25),IM_COL32(90,100,115,255),"Open a charge-density, potential, Cube or XSF file.");
    draw->PopClipRect();
    const auto bar=ImGui::GetCursorScreenPos();
    const float width=std::max(100.0f,size.x-4);
    for(int i=0;i<128;++i)
    {
        const auto c=ElectronicViewport::color(static_cast<float>(i)/127,m_palette);
        draw->AddRectFilled(ImVec2(bar.x+width*i/128,bar.y),ImVec2(bar.x+width*(i+1)/128,bar.y+10),ImGui::ColorConvertFloat4ToU32(ImVec4(c.x,c.y,c.z,1)));
    }
    ImGui::Dummy(ImVec2(width,12));
    ImGui::Text("%.5g",m_colorLow); ImGui::SameLine(width*.4f); ImGui::TextUnformatted(m_colorUnit.c_str());
    ImGui::SameLine(std::max(width-85.0f,150.0f)); ImGui::Text("%.5g",m_colorHigh);
}
