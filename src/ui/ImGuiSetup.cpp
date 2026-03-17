#include "ImGuiSetup.h"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

namespace
{
void applyAtomsEditorTheme()
{
    ImGuiStyle& style = ImGui::GetStyle();
    ImVec4* colors = style.Colors;

    style.WindowPadding = ImVec2(14.0f, 12.0f);
    style.FramePadding = ImVec2(10.0f, 7.0f);
    style.CellPadding = ImVec2(8.0f, 6.0f);
    style.ItemSpacing = ImVec2(10.0f, 8.0f);
    style.ItemInnerSpacing = ImVec2(8.0f, 6.0f);
    style.TouchExtraPadding = ImVec2(1.0f, 1.0f);
    style.IndentSpacing = 20.0f;
    style.ScrollbarSize = 14.0f;
    style.GrabMinSize = 10.0f;

    style.WindowRounding = 12.0f;
    style.ChildRounding = 10.0f;
    style.FrameRounding = 9.0f;
    style.PopupRounding = 12.0f;
    style.ScrollbarRounding = 12.0f;
    style.GrabRounding = 12.0f;
    style.TabRounding = 10.0f;

    style.WindowBorderSize = 1.0f;
    style.ChildBorderSize = 1.0f;
    style.PopupBorderSize = 1.0f;
    style.FrameBorderSize = 0.0f;
    style.TabBorderSize = 0.0f;

    colors[ImGuiCol_Text]                 = ImVec4(0.93f, 0.95f, 0.97f, 1.00f);
    colors[ImGuiCol_TextDisabled]         = ImVec4(0.53f, 0.58f, 0.63f, 1.00f);
    colors[ImGuiCol_WindowBg]             = ImVec4(0.08f, 0.10f, 0.13f, 1.00f);
    colors[ImGuiCol_ChildBg]              = ImVec4(0.10f, 0.12f, 0.16f, 1.00f);
    colors[ImGuiCol_PopupBg]              = ImVec4(0.09f, 0.11f, 0.15f, 1.00f);
    colors[ImGuiCol_Border]               = ImVec4(0.22f, 0.28f, 0.34f, 0.70f);
    colors[ImGuiCol_BorderShadow]         = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);

    colors[ImGuiCol_FrameBg]              = ImVec4(0.14f, 0.17f, 0.22f, 1.00f);
    colors[ImGuiCol_FrameBgHovered]       = ImVec4(0.18f, 0.23f, 0.29f, 1.00f);
    colors[ImGuiCol_FrameBgActive]        = ImVec4(0.21f, 0.28f, 0.35f, 1.00f);

    colors[ImGuiCol_TitleBg]              = ImVec4(0.07f, 0.09f, 0.12f, 1.00f);
    colors[ImGuiCol_TitleBgActive]        = ImVec4(0.10f, 0.13f, 0.17f, 1.00f);
    colors[ImGuiCol_TitleBgCollapsed]     = ImVec4(0.07f, 0.09f, 0.12f, 0.80f);
    colors[ImGuiCol_MenuBarBg]            = ImVec4(0.09f, 0.12f, 0.16f, 1.00f);

    colors[ImGuiCol_ScrollbarBg]          = ImVec4(0.07f, 0.08f, 0.11f, 0.80f);
    colors[ImGuiCol_ScrollbarGrab]        = ImVec4(0.27f, 0.33f, 0.40f, 0.95f);
    colors[ImGuiCol_ScrollbarGrabHovered] = ImVec4(0.34f, 0.42f, 0.50f, 0.95f);
    colors[ImGuiCol_ScrollbarGrabActive]  = ImVec4(0.43f, 0.53f, 0.63f, 0.95f);

    colors[ImGuiCol_CheckMark]            = ImVec4(0.89f, 0.73f, 0.31f, 1.00f);
    colors[ImGuiCol_SliderGrab]           = ImVec4(0.31f, 0.74f, 0.77f, 0.95f);
    colors[ImGuiCol_SliderGrabActive]     = ImVec4(0.42f, 0.82f, 0.84f, 1.00f);

    colors[ImGuiCol_Button]               = ImVec4(0.15f, 0.39f, 0.46f, 1.00f);
    colors[ImGuiCol_ButtonHovered]        = ImVec4(0.20f, 0.52f, 0.60f, 1.00f);
    colors[ImGuiCol_ButtonActive]         = ImVec4(0.24f, 0.60f, 0.68f, 1.00f);

    colors[ImGuiCol_Header]               = ImVec4(0.14f, 0.32f, 0.38f, 0.85f);
    colors[ImGuiCol_HeaderHovered]        = ImVec4(0.20f, 0.46f, 0.53f, 0.92f);
    colors[ImGuiCol_HeaderActive]         = ImVec4(0.25f, 0.56f, 0.63f, 1.00f);

    colors[ImGuiCol_Separator]            = ImVec4(0.22f, 0.28f, 0.34f, 0.75f);
    colors[ImGuiCol_SeparatorHovered]     = ImVec4(0.31f, 0.74f, 0.77f, 0.78f);
    colors[ImGuiCol_SeparatorActive]      = ImVec4(0.40f, 0.82f, 0.84f, 0.95f);

    colors[ImGuiCol_ResizeGrip]           = ImVec4(0.21f, 0.28f, 0.35f, 0.35f);
    colors[ImGuiCol_ResizeGripHovered]    = ImVec4(0.31f, 0.74f, 0.77f, 0.55f);
    colors[ImGuiCol_ResizeGripActive]     = ImVec4(0.40f, 0.82f, 0.84f, 0.85f);

    colors[ImGuiCol_Tab]                  = ImVec4(0.11f, 0.14f, 0.18f, 1.00f);
    colors[ImGuiCol_TabHovered]           = ImVec4(0.22f, 0.48f, 0.56f, 0.95f);
    colors[ImGuiCol_TabActive]            = ImVec4(0.17f, 0.36f, 0.43f, 1.00f);
    colors[ImGuiCol_TabUnfocused]         = ImVec4(0.09f, 0.11f, 0.14f, 1.00f);
    colors[ImGuiCol_TabUnfocusedActive]   = ImVec4(0.13f, 0.18f, 0.23f, 1.00f);

    colors[ImGuiCol_PlotLines]            = ImVec4(0.80f, 0.83f, 0.87f, 1.00f);
    colors[ImGuiCol_PlotLinesHovered]     = ImVec4(0.89f, 0.73f, 0.31f, 1.00f);
    colors[ImGuiCol_PlotHistogram]        = ImVec4(0.31f, 0.74f, 0.77f, 1.00f);
    colors[ImGuiCol_PlotHistogramHovered] = ImVec4(0.42f, 0.82f, 0.84f, 1.00f);

    colors[ImGuiCol_TableHeaderBg]        = ImVec4(0.11f, 0.15f, 0.19f, 1.00f);
    colors[ImGuiCol_TableBorderStrong]    = ImVec4(0.24f, 0.29f, 0.34f, 1.00f);
    colors[ImGuiCol_TableBorderLight]     = ImVec4(0.16f, 0.20f, 0.24f, 1.00f);
    colors[ImGuiCol_TableRowBg]           = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    colors[ImGuiCol_TableRowBgAlt]        = ImVec4(0.10f, 0.12f, 0.16f, 0.55f);

    colors[ImGuiCol_TextSelectedBg]       = ImVec4(0.31f, 0.74f, 0.77f, 0.28f);
    colors[ImGuiCol_DragDropTarget]       = ImVec4(0.89f, 0.73f, 0.31f, 0.95f);
    colors[ImGuiCol_NavHighlight]         = ImVec4(0.31f, 0.74f, 0.77f, 0.85f);
    colors[ImGuiCol_NavWindowingHighlight]= ImVec4(1.00f, 1.00f, 1.00f, 0.70f);
    colors[ImGuiCol_ModalWindowDimBg]     = ImVec4(0.03f, 0.04f, 0.06f, 0.66f);
}
} // namespace

void initImGui(GLFWwindow* window)
{
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    io.FontGlobalScale = 1.05f;
    io.ConfigWindowsMoveFromTitleBarOnly = true;

    applyAtomsEditorTheme();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");
}

void shutdownImGui()
{
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}
