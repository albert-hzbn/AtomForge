#pragma once
// ImGui configuration for the headless scientific-dialog test: any ImGui
// usage error (ID stack, Begin/End or disabled-scope mismatch) fails the test.
#include <stdexcept>
// A function, so that an assertion in a destructor (noexcept) compiles
// without -Wterminate; one that fires there still stops the test.
[[noreturn]] inline void imguiTestAssertionFailed(const char* expression)
{
    throw std::logic_error(std::string("ImGui assertion failed: ") + expression);
}
#define IM_ASSERT(_EXPR) do { if (!(_EXPR)) imguiTestAssertionFailed(#_EXPR); } while (0)
