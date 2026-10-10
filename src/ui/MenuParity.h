#pragma once
// Where every structure-pipeline step appears in the menus, so each Build and
// Edit feature exists both as a menu option and as a pipeline step. Steps
// without a dedicated dialog are offered through the one-step operation dialog
// (OperationDialog) in the Build or Edit menu. Tests check the table covers
// every registered step.

#include <string>
#include <vector>

struct MenuLocation
{
    const char* step;     // pipeline step id
    const char* menu;     // "Build" or "Edit"
    const char* submenu;  // "" for top level, else the submenu
    const char* label;    // menu item text
    bool dialog;          // true: an existing dedicated dialog; false: the operation dialog
};

const std::vector<MenuLocation>& menuLocations();
// "Build > Vacancies", or empty when a step has no menu location.
std::string menuPathOf(const std::string& step);
