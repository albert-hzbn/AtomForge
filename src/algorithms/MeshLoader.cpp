#include "algorithms/MeshLoader.h"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

bool parseObjMesh(const std::string& path,
                  std::vector<glm::vec3>& outVertices,
                  std::vector<unsigned int>& outIndices,
                  std::string& error)
{
    std::ifstream in(path.c_str());
    if (!in)
    {
        error = "Cannot open OBJ file.";
        return false;
    }

    outVertices.clear();
    outIndices.clear();

    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;

        std::istringstream ls(line);
        std::string tag;
        ls >> tag;

        if (tag == "v")
        {
            float x = 0.0f, y = 0.0f, z = 0.0f;
            if (ls >> x >> y >> z)
                outVertices.push_back(glm::vec3(x, y, z));
        }
        else if (tag == "f")
        {
            std::vector<unsigned int> face;
            std::string token;
            while (ls >> token)
            {
                const size_t slashPos = token.find('/');
                const std::string idStr = (slashPos == std::string::npos) ? token : token.substr(0, slashPos);
                if (idStr.empty())
                    continue;

                const int idx = std::atoi(idStr.c_str());
                if (idx <= 0)
                    continue;

                face.push_back((unsigned int)(idx - 1));
            }

            if (face.size() >= 3)
            {
                for (size_t i = 1; i + 1 < face.size(); ++i)
                {
                    outIndices.push_back(face[0]);
                    outIndices.push_back(face[i]);
                    outIndices.push_back(face[i + 1]);
                }
            }
        }
    }

    if (outVertices.empty() || outIndices.empty())
    {
        error = "OBJ has no usable vertices/faces.";
        return false;
    }

    for (unsigned int idx : outIndices)
    {
        if (idx >= outVertices.size())
        {
            error = "OBJ face index out of range.";
            return false;
        }
    }

    return true;
}

bool parseStlMesh(const std::string& path,
                  std::vector<glm::vec3>& outVertices,
                  std::vector<unsigned int>& outIndices,
                  std::string& error)
{
    std::ifstream in(path.c_str(), std::ios::binary);
    if (!in)
    {
        error = "Cannot open STL file.";
        return false;
    }

    in.seekg(0, std::ios::end);
    const std::streamoff fileSize = in.tellg();
    in.seekg(0, std::ios::beg);

    char header[80] = {0};
    in.read(header, 80);
    unsigned int triCount = 0;
    in.read(reinterpret_cast<char*>(&triCount), sizeof(unsigned int));

    const std::streamoff expectedBinarySize = 84 + (std::streamoff)triCount * 50;
    const bool looksBinary = (fileSize == expectedBinarySize);

    outVertices.clear();
    outIndices.clear();

    if (looksBinary)
    {
        outVertices.reserve((size_t)triCount * 3);
        outIndices.reserve((size_t)triCount * 3);

        for (unsigned int t = 0; t < triCount; ++t)
        {
            float data[12] = {0.0f};
            unsigned short attr = 0;
            in.read(reinterpret_cast<char*>(data), sizeof(data));
            in.read(reinterpret_cast<char*>(&attr), sizeof(attr));
            if (!in)
            {
                error = "Binary STL parse failed.";
                return false;
            }

            const unsigned int base = (unsigned int)outVertices.size();
            outVertices.push_back(glm::vec3(data[3], data[4], data[5]));
            outVertices.push_back(glm::vec3(data[6], data[7], data[8]));
            outVertices.push_back(glm::vec3(data[9], data[10], data[11]));
            outIndices.push_back(base + 0);
            outIndices.push_back(base + 1);
            outIndices.push_back(base + 2);
        }

        return !outVertices.empty();
    }

    in.clear();
    in.seekg(0, std::ios::beg);

    std::string line;
    std::vector<glm::vec3> triVerts;
    triVerts.reserve(3);

    while (std::getline(in, line))
    {
        std::istringstream ls(line);
        std::string tag;
        ls >> tag;
        if (tag == "vertex")
        {
            float x = 0.0f, y = 0.0f, z = 0.0f;
            if (ls >> x >> y >> z)
            {
                triVerts.push_back(glm::vec3(x, y, z));
                if (triVerts.size() == 3)
                {
                    const unsigned int base = (unsigned int)outVertices.size();
                    outVertices.push_back(triVerts[0]);
                    outVertices.push_back(triVerts[1]);
                    outVertices.push_back(triVerts[2]);
                    outIndices.push_back(base + 0);
                    outIndices.push_back(base + 1);
                    outIndices.push_back(base + 2);
                    triVerts.clear();
                }
            }
        }
    }

    if (outVertices.empty() || outIndices.empty())
    {
        error = "ASCII STL has no usable triangles.";
        return false;
    }

    return true;
}
