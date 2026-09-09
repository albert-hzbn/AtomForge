#include "ui/ElectronicPostProcessing.h"
#include "imgui.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>

using namespace atomforge::electronic;

namespace
{
ImU32 scalarColor(double value, double low, double high, float alpha = 1)
{
    const float t = high > low ? static_cast<float>(std::clamp((value - low) / (high - low),0.0,1.0)) : 0.5f;
    return ImGui::ColorConvertFloat4ToU32(ImVec4(t, 0.25f + 0.5f * (1 - std::abs(2*t-1)), 1-t, alpha));
}
glm::dvec3 vector(const float* p) { return {p[0],p[1],p[2]}; }
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
            auto previousMesh = std::move(m_result.mesh);
            m_result = *m_task.result();
            if (m_result.appendMesh && previousMesh.vertices.size() + m_result.mesh.vertices.size() <= 12000000)
            {
                m_result.mesh.vertices.insert(m_result.mesh.vertices.end(),previousMesh.vertices.begin(),previousMesh.vertices.end());
                m_result.mesh.colors.insert(m_result.mesh.colors.end(),previousMesh.colors.begin(),previousMesh.colors.end());
            }
            if (m_result.referenceLoaded) { m_referenceVolume = m_result.volume; m_reference = 0; }
            else if (m_result.loaded) { m_volume = m_result.volume; m_selected = m_reference = 0; }
            else if (!m_result.volume.fields.empty())
            {
                m_selected = static_cast<int>(m_volume.fields.size());
                m_volume.fields.insert(m_volume.fields.end(),m_result.volume.fields.begin(),m_result.volume.fields.end());
            }
        }
        m_task.clearResult();
    }
    if (!m_open) return;
    ImGui::SetNextWindowSize(ImVec2(850,760),ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("Electronic Post-processing",&m_open)) { ImGui::End(); return; }
    ImGui::TextWrapped("VASP charge / potential / ELF, Gaussian Cube and XSF scalar grids. Coordinates: Angstrom. Imported and derived fields are kept separately from the structure editor.");
    ImGui::BeginDisabled(m_task.running());
    ImGui::InputText("Input file",m_input,sizeof(m_input));
    ImGui::Combo("Quantity",&m_quantity,"Auto (Cube/XSF: raw)\0Electron density\0Potential\0ELF / dimensionless\0Raw scalar\0");
    ImGui::Combo("Cube coordinates",&m_cubeUnits,"Bohr\0Angstrom (nonstandard files)\0");
    if (ImGui::Button("Load volume"))
    {
        const std::string path(m_input), q = std::array<const char*,5>{"auto","density","potential","elf","raw"}[m_quantity];
        const std::string units = m_cubeUnits ? "angstrom" : "bohr";
        m_error.clear();
        m_task.start([path,q,units] { Output out; out.volume = loadVolume(path,q,units); out.loaded = true; return out; });
    }
    ImGui::SameLine();
    if (ImGui::Button("Load as reference"))
    {
        const std::string path(m_input), q = std::array<const char*,5>{"auto","density","potential","elf","raw"}[m_quantity];
        const std::string units = m_cubeUnits ? "angstrom" : "bohr";
        m_error.clear();
        m_task.start([path,q,units] { Output out; out.volume = loadVolume(path,q,units); out.referenceLoaded = true; return out; });
    }
    ImGui::SameLine();
    if (ImGui::Button("Clear fields")) { m_volume = {}; m_referenceVolume = {}; m_result = {}; m_selected = m_reference = 0; }
    if (!m_volume.fields.empty())
    {
        auto selectField = [&](const char* label, int& selected, const Volume& volume)
        {
            selected = std::clamp(selected,0,static_cast<int>(volume.fields.size())-1);
            if (ImGui::BeginCombo(label,volume.fields[selected].name.c_str()))
            {
                for (int i = 0; i < static_cast<int>(volume.fields.size()); ++i)
                {
                    ImGui::PushID(i);
                    if (ImGui::Selectable(volume.fields[i].name.c_str(),i == selected)) selected = i;
                    ImGui::PopID();
                }
                ImGui::EndCombo();
            }
        };
        selectField("Field",m_selected,m_volume);
        const auto& referenceVolume = m_referenceVolume.fields.empty() ? m_volume : m_referenceVolume;
        m_reference = std::clamp(m_reference,0,static_cast<int>(referenceVolume.fields.size())-1);
        const Grid& g = m_volume.fields[m_selected];
        ImGui::Text("%d x %d x %d | %s | %s | %zu atoms",g.shape[0],g.shape[1],g.shape[2],g.unit.c_str(),g.periodic ? "periodic" : "finite, endpoints included",m_volume.sites.size());
        ImGui::Combo("Tool",&m_operation,
            "Total integral\0Add reference\0Subtract reference (difference density)\0Multiply reference\0Divide by reference\0Scale\0Gaussian smoothing\0Cartesian gradient\0Laplacian\0Energy-density conversion\0Line profile\0Planar average\0Macroscopic average\0Plane section\0Contour segments (z=0)\0Peak search\0Voronoi site integration\0Sphere integration\0Structure factors\0Fourier synthesis\0Patterson density\0Ewald site potentials\0Isosurface\0Resample onto reference\0Verify periodic endpoint planes\0");
        if ((m_operation >= 1 && m_operation <= 4) || m_operation == 23 || m_operation == 22) selectField("Reference field",m_reference,referenceVolume);
        if (m_operation == 5) ImGui::InputFloat("Scale factor",&m_scalar);
        if (m_operation == 6) { ImGui::InputFloat("Gaussian sigma (A)",&m_sigma); ImGui::InputInt("Kernel radius (grid steps)",&m_radius); }
        if (m_operation == 9) ImGui::TextWrapped("Gradient-expansion kinetic and local-virial energy densities (eV/A^3). Requires nonnegative e/A^3. Samples <= 1e-12 e/A^3 are masked to zero.");
        if (m_operation == 10 || m_operation == 13 || m_operation == 17) ImGui::InputFloat3("Start / origin / center (A)",m_start);
        if (m_operation == 10) ImGui::InputFloat3("End (A)",m_end);
        if (m_operation == 10 || m_operation == 13) ImGui::InputInt("Sample count",&m_count);
        if (m_operation == 11 || m_operation == 12) ImGui::Combo("Normal to lattice plane",&m_axis,"a*\0b*\0c*\0");
        if (m_operation == 12) ImGui::InputInt("Odd averaging window",&m_window);
        if (m_operation == 13) { ImGui::InputFloat3("Plane span u (A)",m_u); ImGui::InputFloat3("Plane span v (A)",m_v); }
        if (m_operation == 14 || m_operation == 22) ImGui::InputFloat("Isovalue",&m_scalar);
        if (m_operation == 17) ImGui::InputFloat("Sphere radius (A)",&m_scalar);
        if (m_operation == 15) ImGui::TextWrapped("26-neighbor grid maxima; adjacent equal values use a deterministic tie break. Refine the grid to check peak positions.");
        if (m_operation == 16) ImGui::TextWrapped("Integrates nearest-site Voronoi cells around imported atoms. Boundary samples share weight. These are geometric partitions, not Bader basins.");
        if (m_operation == 18 || m_operation == 19)
        {
            ImGui::TextUnformatted(m_operation == 18 ? "One h k l per line" : "One h k l real imag per line; include conjugate pairs");
            ImGui::InputTextMultiline("##reflections",m_reflections,sizeof(m_reflections),ImVec2(-1,100));
        }
        if (m_operation == 21)
        {
            ImGui::TextWrapped("Explicit charges in electron units, one per imported atom in file order. Requires a neutral cell; conducting boundary. Increase both cutoffs to check convergence.");
            ImGui::InputTextMultiline("Charges",m_charges,sizeof(m_charges),ImVec2(-1,70));
            ImGui::InputFloat("Ewald alpha (1/A)",&m_alpha);
            ImGui::InputFloat("Real cutoff (A)",&m_realCutoff);
            ImGui::InputFloat("Reciprocal cutoff (1/A, includes 2*pi)",&m_reciprocalCutoff);
        }
        if (m_operation == 22)
        {
            ImGui::Checkbox("Color surface using reference",&m_colorSurface);
            ImGui::Checkbox("Append to previous surface (multiple levels)",&m_appendSurface);
        }
        if (m_operation == 24) ImGui::TextWrapped("For Cube/XSF covering a complete periodic cell: checks all opposite endpoint values, then removes duplicate planes. Matching values alone do not establish physical periodicity.");
        if (ImGui::Button("Calculate"))
        {
            const Grid source = g, reference = referenceVolume.fields[m_reference];
            const auto sites = m_volume.sites;
            const int op = m_operation, count = m_count, axis = m_axis, window = m_window, radius = m_radius;
            const double scalar = m_scalar, sigma = m_sigma, alpha = m_alpha, realCutoff = m_realCutoff, reciprocalCutoff = m_reciprocalCutoff;
            const auto start = vector(m_start), end = vector(m_end), u = vector(m_u), v = vector(m_v);
            const std::string reflectionText(m_reflections), chargeText(m_charges);
            const bool colored = m_colorSurface;
            const bool appendSurface = m_appendSurface;
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
                else if (op == 22) { out.mesh = isosurface(source,scalar,colored ? &reference : nullptr); out.appendMesh = appendSurface; }
                else if (op == 23) add(resample(source,reference));
                else if (op == 24)
                {
                    add(periodicEndpoints(source));
                }
                return out;
            });
        }
        ImGui::Separator();
        ImGui::InputText("Output file",m_output,sizeof(m_output));
        ImGui::Combo("Export",&m_exportFormat,"Selected field: XSF\0Selected field: Cube\0Selected field: VASP\0Last table: CSV\0Last surface: OBJ / PLY\0");
        if (ImGui::Button("Save result"))
        {
            try
            {
                if (m_exportFormat < 3) saveVolume({m_volume.sites,{g}},m_output,std::array<const char*,3>{"xsf","cube","vasp"}[m_exportFormat]);
                else if (m_exportFormat == 4) saveMesh(m_result.mesh,m_output);
                else
                {
                    if (!m_result.columns) throw std::invalid_argument("No table to export");
                    std::ofstream stream{std::filesystem::u8path(m_output)};
                    if (!stream) throw std::runtime_error("Cannot open CSV output");
                    stream << "# " << m_result.heading << '\n' << std::setprecision(17);
                    for (std::size_t i = 0; i < m_result.table.size(); ++i) stream << m_result.table[i] << ((i+1)%m_result.columns ? ',' : '\n');
                    if (!stream) throw std::runtime_error("Failed writing CSV");
                }
                m_error = "Saved " + std::string(m_output);
            }
            catch (const std::exception& e) { m_error = e.what(); }
        }
    }
    ImGui::EndDisabled();
    if (m_task.running()) ImGui::TextUnformatted("Calculating...");
    if (!m_error.empty()) ImGui::TextWrapped("%s",m_error.c_str());
    if (!m_result.table.empty())
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
    ImGui::End();
}

void ElectronicPostProcessingDialog::drawPreview()
{
    if (!m_result.mesh.vertices.empty())
    {
        ImGui::Text("%zu triangles (preview displays up to 20000)",m_result.mesh.vertices.size()/3);
        ImGui::SliderFloat("Rotation",&m_yaw,-3.14f,3.14f);
        ImGui::SliderFloat("Tilt",&m_pitch,-1.57f,1.57f);
        ImGui::SliderFloat("Opacity",&m_opacity,0.05f,1);
        const auto pos = ImGui::GetCursorScreenPos();
        const ImVec2 size(std::max(100.0f,ImGui::GetContentRegionAvail().x),300);
        ImGui::InvisibleButton("Surface preview",size);
        auto* draw = ImGui::GetWindowDrawList();
        draw->AddRectFilled(pos,ImVec2(pos.x+size.x,pos.y+size.y),IM_COL32(25,28,35,255));
        const auto& vertices = m_result.mesh.vertices;
        glm::dvec3 low = vertices[0], high = low;
        for (auto p : vertices) { low = glm::min(low,p); high = glm::max(high,p); }
        const auto center = (low+high)*0.5;
        const double zoom = 0.8 * std::min(size.x,size.y) / std::max(1e-8,glm::length(high-low));
        auto rotate = [&](glm::dvec3 p)
        {
            p -= center;
            glm::dvec3 q(std::cos(m_yaw)*p.x+std::sin(m_yaw)*p.z,p.y,-std::sin(m_yaw)*p.x+std::cos(m_yaw)*p.z);
            return glm::dvec3(q.x,std::cos(m_pitch)*q.y-std::sin(m_pitch)*q.z,std::sin(m_pitch)*q.y+std::cos(m_pitch)*q.z);
        };
        struct Triangle { std::size_t index; double depth; };
        std::vector<Triangle> order;
        const std::size_t stride = std::max<std::size_t>(1,(vertices.size()/3+19999)/20000);
        for (std::size_t i = 0; i+2 < vertices.size(); i += 3*stride) order.push_back({i,rotate((vertices[i]+vertices[i+1]+vertices[i+2])/3.0).z});
        std::sort(order.begin(),order.end(),[](const Triangle& a,const Triangle& b) { return a.depth < b.depth; });
        const auto range = std::minmax_element(m_result.mesh.colors.begin(),m_result.mesh.colors.end());
        for (const auto& triangle : order)
        {
            ImVec2 points[3];
            for (int j = 0; j < 3; ++j) { const auto p = rotate(vertices[triangle.index+j]); points[j] = ImVec2(pos.x+size.x*0.5f+static_cast<float>(p.x*zoom),pos.y+size.y*0.5f-static_cast<float>(p.y*zoom)); }
            const auto a = rotate(vertices[triangle.index]), b = rotate(vertices[triangle.index+1]), c = rotate(vertices[triangle.index+2]);
            const auto normal = glm::normalize(glm::cross(b-a,c-a));
            const float light = static_cast<float>(0.35 + 0.65 * std::abs(glm::dot(normal,glm::normalize(glm::dvec3(0.3,0.5,1)))));
            auto color = ImGui::ColorConvertU32ToFloat4(scalarColor(m_result.mesh.colors[triangle.index],*range.first,*range.second,m_opacity));
            color.x *= light; color.y *= light; color.z *= light;
            draw->AddTriangleFilled(points[0],points[1],points[2],ImGui::ColorConvertFloat4ToU32(color));
        }
    }
    else if (!m_volume.fields.empty() && m_volume.fields[m_selected].name == "section")
    {
        const auto& g = m_volume.fields[m_selected];
        const auto range = std::minmax_element(g.values.begin(),g.values.end());
        ImGui::Text("Section range: %.6g to %.6g %s",*range.first,*range.second,g.unit.c_str());
        const auto pos = ImGui::GetCursorScreenPos();
        const ImVec2 size(std::max(100.0f,ImGui::GetContentRegionAvail().x),250);
        ImGui::InvisibleButton("Section preview",size);
        auto* draw = ImGui::GetWindowDrawList();
        const int stride = std::max(1,std::max(g.shape[0],g.shape[1])/200);
        for (int y = 0; y < g.shape[1]; y += stride)
            for (int x = 0; x < g.shape[0]; x += stride)
                draw->AddRectFilled(ImVec2(pos.x+size.x*x/g.shape[0],pos.y+size.y*y/g.shape[1]),ImVec2(pos.x+size.x*std::min(x+stride,g.shape[0])/g.shape[0],pos.y+size.y*std::min(y+stride,g.shape[1])/g.shape[1]),scalarColor(g.values[g.index(x,y,0)],*range.first,*range.second));
    }
}
