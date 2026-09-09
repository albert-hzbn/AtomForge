#pragma once

#include <functional>
#include <string>
#include <utility>
#include <vector>

// Common path and directory utility functions shared across the application.

using DirectoryEntry = std::pair<std::string, bool>;

std::string normalizePathSeparators(const std::string& path);
std::string joinPath(const std::string& base, const std::string& name);
std::string parentPath(const std::string& path);
bool isDriveRootPath(const std::string& path);
std::string detectHomePath();
std::vector<std::string> getDriveRoots();
void appendUniquePath(std::vector<std::string>& paths, const std::string& value);

bool hasExtension(const std::string& name, const std::string& extension);
std::string replaceFileExtension(const std::string& filename,
                                 const std::string& extension,
                                 const std::string& fallbackBase);
void updateFilenameWithExtension(char* filenameBuffer,
                                 std::size_t bufferSize,
                                 const std::string& extension,
                                 const std::string& fallbackBase);

std::string toLowerStr(const std::string& s);

int atomicNumberFromSymbol(const std::string& symbol);

bool loadDirectoryEntries(const std::string& directory,
                          bool filterFiles,
                          const std::function<bool(const std::string&)>& includeFile,
                          std::vector<DirectoryEntry>& entries);

void drawDirectoryEntries(const std::vector<DirectoryEntry>& entries,
                          char* selectedFilename,
                          int idBase,
                          const std::function<void(const std::string&)>& onEnterDirectory,
                          bool* fileDoubleClicked = nullptr);

void pushDirectoryHistory(std::vector<std::string>& history,
                          int& historyIndex,
                          const std::string& dir);
