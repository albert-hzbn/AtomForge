#pragma once
// Command-line argument helpers shared by the CLI modes.
#include <glm/glm.hpp>

#include <string>
#include <vector>

namespace cli
{
// Value following `flag`, or nullptr.
const char* findArg(int argc, char* argv[], const char* flag);
bool hasFlag(int argc, char* argv[], const char* flag);
// Values after every occurrence of a repeatable `flag`.
std::vector<std::string> findAllArgs(int argc, char* argv[], const char* flag);
// Numeric values; `def` when the flag is absent. Malformed or non-finite values throw
// std::invalid_argument naming the flag.
double argDouble(int argc, char* argv[], const char* flag, double def);
int argInt(int argc, char* argv[], const char* flag, int def);
// Open Babel format name from an output file extension (cif when there is none).
std::string detectFormat(const std::string& path);
// "x y z" triples (whitespace or comma separated).
bool parseIvec3(const char* text, glm::ivec3& out);
bool parseVec3(const char* text, glm::vec3& out);
}
