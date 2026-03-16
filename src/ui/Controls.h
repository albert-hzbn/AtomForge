#pragma once
#include "Camera.h"
#include <functional>

void drawControls(Camera& camera, std::function<void()> onDeleteSelected = nullptr);