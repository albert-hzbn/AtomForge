#pragma once

#include "app/EditorState.h"

#include <glm/glm.hpp>

struct Camera;
struct ImDrawList;

struct FrameActionRequests
{
    bool doDeleteSelected = false;
    bool requestMeasureDistance = false;
    bool requestMeasureAngle = false;
    bool requestAtomInfo = false;
    bool requestStructureInfo = false;
    bool requestUndo = false;
    bool requestRedo = false;
    bool requestViewAxisX = false;
    bool requestViewAxisY = false;
    bool requestViewAxisZ = false;
    bool requestViewLatticeA = false;
    bool requestRotateCrystalX = false;
    bool requestRotateCrystalY = false;
    bool requestRotateCrystalZ = false;
    bool requestViewLatticeB = false;
    bool requestViewLatticeC = false;
};

FrameActionRequests beginFrameActionRequests(EditorState& state);
void applyKeyboardShortcuts(EditorState& state, FrameActionRequests& requests);

void handlePendingAtomPick(
    Camera& camera,
    EditorState& state,
    const glm::vec3& cameraPosition,
    int windowWidth,
    int windowHeight,
    const glm::mat4& projection,
    const glm::mat4& view);

void handleBoxSelection(
    EditorState& state,
    int windowWidth,
    int windowHeight,
    const glm::mat4& projection,
    const glm::mat4& view,
    ImDrawList* drawList);

void handleRightClick(Camera& camera, EditorState& state);

// Blender-style grab mode: press G to grab selected atoms, move with mouse,
// optionally constrain to X/Y/Z axis, click to confirm, Escape/right-click to cancel.
void handleGrabMode(
    EditorState& state,
    Camera& camera,
    const glm::mat4& projection,
    const glm::mat4& view,
    int windowWidth,
    int windowHeight);

// Draw real-time position overlay for atoms being grabbed.
void drawGrabOverlay(
    const EditorState& state,
    ImDrawList* drawList,
    const glm::mat4& projection,
    const glm::mat4& view,
    int windowWidth,
    int windowHeight);