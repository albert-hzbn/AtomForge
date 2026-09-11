#include "util/ApplicationPaths.h"

#include <vector>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#else
#include <unistd.h>
#endif

std::filesystem::path applicationDirectory()
{
#ifdef _WIN32
    std::vector<wchar_t> buffer(32768);
    const DWORD length = GetModuleFileNameW(nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));
    if (length > 0 && length < buffer.size())
        return std::filesystem::path(std::wstring(buffer.data(), length)).parent_path();
#elif defined(__APPLE__)
    uint32_t length = 0;
    _NSGetExecutablePath(nullptr, &length);
    std::vector<char> buffer(length);
    if (_NSGetExecutablePath(buffer.data(), &length) == 0)
        return std::filesystem::path(buffer.data()).parent_path();
#else
    std::error_code error;
    const auto path = std::filesystem::read_symlink("/proc/self/exe", error);
    if (!error) return path.parent_path();
#endif
    return {};
}

std::filesystem::path findManualPdf(const std::filesystem::path& directory)
{
    if (directory.empty()) return {};
    const auto parent = directory.parent_path();
    for (const auto& candidate : {
        parent / "docs/manual/AtomForge-manual.pdf",
        parent.parent_path() / "docs/manual/AtomForge-manual.pdf",
        directory / "AtomForge-manual.pdf",
        parent / "AtomForge-manual.pdf",
        parent / "share/doc/atomforge/AtomForge-manual.pdf"})
    {
        std::error_code error;
        if (std::filesystem::is_regular_file(candidate, error)) return candidate;
    }
#ifdef ATOMFORGE_SOURCE_MANUAL
    // Out-of-source development builds still use the one published source PDF.
    // Packaged locations above take precedence and remain relocatable.
    if (directory == applicationDirectory())
    {
        const auto sourceManual = std::filesystem::path(ATOMFORGE_SOURCE_MANUAL);
        std::error_code error;
        if (std::filesystem::is_regular_file(sourceManual, error)) return sourceManual;
    }
#endif
    return {};
}
