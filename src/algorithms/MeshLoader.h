#pragma once

#include <glm/glm.hpp>
#include <string>
#include <vector>

// Parse a Wavefront OBJ file into a flat triangle list.
// Returns true on success; sets error on failure.
bool parseObjMesh(const std::string& path,
                  std::vector<glm::vec3>& outVertices,
                  std::vector<unsigned int>& outIndices,
                  std::string& error);

// Parse a binary or ASCII STL file into a flat triangle list.
// Returns true on success; sets error on failure.
bool parseStlMesh(const std::string& path,
                  std::vector<glm::vec3>& outVertices,
                  std::vector<unsigned int>& outIndices,
                  std::string& error);
