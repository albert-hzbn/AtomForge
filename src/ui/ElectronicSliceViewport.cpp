#include "ui/ResponsiveLayout.h"
#include "ui/ElectronicSliceViewport.h"
#include "electronic/DisplayRange.h"
#include "graphics/ElectronicViewport.h"
#include "imgui.h"

#include <algorithm>
#include <cmath>

using namespace atomforge::electronic;

SlicePlane ElectronicSliceViewport::plane(const Grid& grid) const
{
    if (m_axis==3) return slicePlane(grid,{m_point[0],m_point[1],m_point[2]},{m_normal[0],m_normal[1],m_normal[2]});
    const int a=m_axis==0 ? 1 : 0, b=m_axis==2 ? 1 : 2;
    const auto origin=grid.origin+grid.cell[m_axis]*static_cast<double>(m_position);
    return {origin,grid.cell[a],grid.cell[b],{origin,origin+grid.cell[a],origin+grid.cell[a]+grid.cell[b],origin+grid.cell[b]}};
}

void ElectronicSliceViewport::draw(const Grid& grid, float low, float high, int palette)
{
    ImGui::TextUnformatted("2D section");
    ImGui::SetNextItemWidth(responsive::dp(-1));
    m_dirty |= ImGui::Combo("##plane",&m_axis,"bc plane (normal a*)\0ac plane (normal b*)\0ab plane (normal c*)\0Arbitrary plane\0");
    if (m_axis==3)
    {
        ImGui::TextUnformatted("Point (fractional a, b, c)");
        ImGui::SetNextItemWidth(responsive::dp(-1));
        m_dirty |= ImGui::InputFloat3("##plane point",m_point,"%.4g");
        ImGui::TextUnformatted("Normal (Cartesian x, y, z)");
        ImGui::SetNextItemWidth(responsive::dp(-1));
        m_dirty |= ImGui::InputFloat3("##plane normal",m_normal,"%.4g");
    }
    else
    {
        ImGui::TextUnformatted("Slice position (fraction of cell)");
        ImGui::SetNextItemWidth(responsive::dp(-1));
        m_dirty |= ImGui::SliderFloat("##position",&m_position,0,1,"%.3f");
    }
    const auto geometry=plane(grid);
    if (geometry.boundary.size()<3)
    {
        ImGui::TextWrapped("Warning: this plane does not intersect the cell. Change the point or normal.");
        return;
    }
    if (m_dirty)
    {
        // Preview sampling is bounded; the native section/contour tools retain
        // their full user-selected resolution for numerical work and export.
        const int a=m_axis==0 ? 1 : 0, b=m_axis==2 ? 1 : 2;
        const int nu=m_axis==3 ? 129 : std::clamp(grid.shape[a]+(grid.periodic ? 1 : 0),2,129);
        const int nv=m_axis==3 ? 129 : std::clamp(grid.shape[b]+(grid.periodic ? 1 : 0),2,129);
        m_slice=sampleSlicePlane(grid,geometry,nu,nv);
        m_visibleValues.clear();
        const auto inverse=glm::inverse(grid.cell);
        for (int y=0;y<nv;++y) for (int x=0;x<nu;++x)
        {
            const auto f=inverse*(m_slice.position(x,y,0)-grid.origin);
            if (m_axis!=3 || (glm::all(glm::greaterThanEqual(f,glm::dvec3(-1e-9))) && glm::all(glm::lessThanEqual(f,glm::dvec3(1+1e-9)))))
                m_visibleValues.push_back(m_slice.values[m_slice.index(x,y,0)]);
        }
        // Include the intersection vertices even for a very small section.
        for (auto p : geometry.boundary) m_visibleValues.push_back(grid.sample(p));
    }
    if (m_estimateLevel) { m_level=static_cast<float>(displayRange(m_visibleValues).suggested); m_estimateLevel=false; }
    ImGui::Checkbox("Show contour line",&m_showContour);
    ImGui::TextUnformatted("2D isovalue (contour)");
    ImGui::SetNextItemWidth(responsive::dp(-1));
    bool levelChanged = ImGui::InputFloat("##level",&m_level,0,0,"%.5g");
    if (responsive::button("Estimate 2D level"))
    {
        m_level=static_cast<float>(displayRange(m_visibleValues).suggested);
        levelChanged=true;
    }
    if (responsive::button("Reset 2D")) reset();
    ImGui::SameLine();
    ImGui::SetNextItemWidth(responsive::dp(-1));
    ImGui::SliderFloat("##angle",&m_angle,-180,180,"%.0f deg");
    ImGui::TextWrapped("Drag: rotate in plane | Right drag: pan | Wheel: zoom");

    if (!std::isfinite(m_level))
    {
        ImGui::TextWrapped("Enter a finite 2D isovalue.");
        return;
    }
    if (m_dirty || levelChanged) m_contours=contours(m_slice,m_level);
    m_dirty=false;
    const auto bounds=std::minmax_element(m_visibleValues.begin(),m_visibleValues.end());
    if (m_level<*bounds.first || m_level>*bounds.second)
    {
        ImGui::PushStyleColor(ImGuiCol_Text,ImVec4(.85f,.32f,.08f,1));
        ImGui::TextWrapped("Warning: contour value is outside this slice's range (%.5g to %.5g %s).",*bounds.first,*bounds.second,grid.unit.c_str());
        ImGui::PopStyleColor();
    }
    else if (*bounds.first==*bounds.second) ImGui::TextWrapped("Constant slice: no distinct contour exists.");
    ImGui::TextWrapped("Range: %.5g to %.5g %s",low,high,grid.unit.c_str());

    const auto pos=ImGui::GetCursorScreenPos();
    const ImVec2 size(std::max(80.0f,ImGui::GetContentRegionAvail().x),std::max(80.0f,ImGui::GetContentRegionAvail().y));
    ImGui::InvisibleButton("2D canvas",size,ImGuiButtonFlags_MouseButtonLeft|ImGuiButtonFlags_MouseButtonRight);
    if (ImGui::IsItemHovered())
    {
        m_zoom=std::clamp(m_zoom*std::exp(ImGui::GetIO().MouseWheel*.12f),.1f,20.0f);
        if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) reset();
    }
    if (ImGui::IsItemActive())
    {
        const auto delta=ImGui::GetIO().MouseDelta;
        if (ImGui::IsMouseDown(ImGuiMouseButton_Left)) m_angle=std::remainder(m_angle+delta.x*.5f,360.0f);
        if (ImGui::IsMouseDown(ImGuiMouseButton_Right)) m_pan+=glm::vec2(delta.x/size.x,delta.y/size.y);
    }

    // Project onto an orthonormal basis of the actual lattice plane so skewed
    // cells retain their angles and lengths when rotated on screen.
    const auto ex=glm::normalize(m_slice.cell[0]);
    const auto ey=glm::normalize(m_slice.cell[1]-ex*glm::dot(ex,m_slice.cell[1]));
    const auto center=m_slice.origin+(m_slice.cell[0]+m_slice.cell[1])*.5;
    const double extent=std::max(glm::length(m_slice.cell[0]+m_slice.cell[1]),glm::length(m_slice.cell[0]-m_slice.cell[1]));
    const float scale=.85f*std::min(size.x,size.y)*m_zoom/static_cast<float>(extent);
    const float angle=m_angle*static_cast<float>(pi/180), cs=std::cos(angle), sn=std::sin(angle);
    auto screen=[&](glm::dvec3 point)
    {
        const auto p=point-center;
        const float x=static_cast<float>(glm::dot(p,ex)), y=static_cast<float>(glm::dot(p,ey));
        return ImVec2(pos.x+size.x*(.5f+m_pan.x)+scale*(cs*x-sn*y),pos.y+size.y*(.5f+m_pan.y)-scale*(sn*x+cs*y));
    };
    auto* draw=ImGui::GetWindowDrawList();
    draw->PushClipRect(pos,ImVec2(pos.x+size.x,pos.y+size.y),true);
    draw->AddRectFilled(pos,ImVec2(pos.x+size.x,pos.y+size.y),IM_COL32(245,247,251,255));
    const auto flags=draw->Flags;
    draw->Flags &= ~ImDrawListFlags_AntiAliasedFill;
    for (int y=0;y<m_slice.shape[1]-1;++y) for (int x=0;x<m_slice.shape[0]-1;++x)
    {
        const double value=(m_slice.values[m_slice.index(x,y,0)]+m_slice.values[m_slice.index(x+1,y,0)]
            +m_slice.values[m_slice.index(x,y+1,0)]+m_slice.values[m_slice.index(x+1,y+1,0)])*.25;
        const float t=high>low ? static_cast<float>((value-low)/(high-low)) : .5f;
        const auto c=ElectronicViewport::color(t,palette);
        auto polygon=std::vector<glm::dvec3>{m_slice.position(x,y,0),m_slice.position(x,y+1,0),m_slice.position(x+1,y+1,0),m_slice.position(x+1,y,0)};
        if (m_axis==3) polygon=clipToCell(grid,std::move(polygon));
        std::vector<ImVec2> points;
        for (auto p : polygon) points.push_back(screen(p));
        if (points.size()>=3) draw->AddConvexPolyFilled(points.data(),static_cast<int>(points.size()),ImGui::ColorConvertFloat4ToU32(ImVec4(c.x,c.y,c.z,1)));
    }
    draw->Flags=flags;
    for (std::size_t i=0;m_showContour && i+1<m_contours.size();i+=2)
    {
        auto start=m_contours[i], end=m_contours[i+1];
        if (m_axis==3 && !clipSegmentToCell(grid,start,end)) continue;
        const auto a=screen(start), b=screen(end);
        draw->AddLine(a,b,IM_COL32(30,40,55,255),1);
    }
    draw->PopClipRect();
}
