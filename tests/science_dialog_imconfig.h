#pragma once
// ImGui configuration for the headless scientific-dialog test: any ImGui
// usage error (ID stack, Begin/End or disabled-scope mismatch) fails the test.
#include <stdexcept>
#define IM_ASSERT(_EXPR) do { if (!(_EXPR)) throw std::logic_error("ImGui assertion failed: " #_EXPR); } while (0)
