#include "StructureLoader.h"

#include <openbabel3/openbabel/obconversion.h>
#include <openbabel3/openbabel/mol.h>
#include <openbabel3/openbabel/atom.h>
#include <openbabel3/openbabel/elements.h>
#include <openbabel3/openbabel/generic.h>
#include <openbabel3/openbabel/math/vector3.h>
#include <openbabel3/openbabel/plugin.h>
#include <openbabel3/openbabel/oberror.h>

#include <iostream>
#include <array>
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <sys/stat.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_inverse.hpp>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dirent.h>
#include <unistd.h>
#endif

namespace
{
constexpr const char* kIpfSidecarSuffix = ".atomforge-ipf";

struct IpfKey
{
    int atomicNumber = 0;
    long long qx = 0;
    long long qy = 0;
    long long qz = 0;

    bool operator==(const IpfKey& other) const
    {
        return atomicNumber == other.atomicNumber
            && qx == other.qx
            && qy == other.qy
            && qz == other.qz;
    }
};

struct IpfKeyHash
{
    size_t operator()(const IpfKey& key) const
    {
        size_t h = (size_t)key.atomicNumber;
        h ^= (size_t)(key.qx * 73856093ll);
        h ^= (size_t)(key.qy * 19349663ll);
        h ^= (size_t)(key.qz * 83492791ll);
        return h;
    }
};

struct IpfRecord
{
    IpfKey key;
    std::array<float, 3> color = {{0.0f, 0.0f, 0.0f}};
};

struct Vec3iHash
{
    size_t operator()(const glm::ivec3& value) const
    {
        return ((size_t)value.x * 73856093u)
             ^ ((size_t)value.y * 19349663u)
             ^ ((size_t)value.z * 83492791u);
    }
};

long long quantizeIpfCoord(double value)
{
    return (long long)std::llround(value * 10000.0);
}

IpfKey makeIpfKey(int atomicNumber, double x, double y, double z)
{
    IpfKey key;
    key.atomicNumber = atomicNumber;
    key.qx = quantizeIpfCoord(x);
    key.qy = quantizeIpfCoord(y);
    key.qz = quantizeIpfCoord(z);
    return key;
}

std::string ipfSidecarPath(const std::string& filename)
{
    const std::size_t slashPos = filename.find_last_of("/\\");
    const std::size_t dotPos = filename.find_last_of('.');
    const bool hasExtension = (dotPos != std::string::npos)
                           && (slashPos == std::string::npos || dotPos > slashPos);
    if (!hasExtension)
        return filename + kIpfSidecarSuffix;
    return filename.substr(0, dotPos) + kIpfSidecarSuffix;
}

std::string ipfLegacySidecarPath(const std::string& filename)
{
    return filename + kIpfSidecarSuffix;
}

void removeIpfSidecar(const std::string& filename)
{
    const std::string path = ipfSidecarPath(filename);
    const std::string legacyPath = ipfLegacySidecarPath(filename);
    std::remove(path.c_str());
    if (legacyPath != path)
        std::remove(legacyPath.c_str());
}

bool saveIpfSidecar(const Structure& structure, const std::string& filename)
{
    if (structure.grainColors.size() != structure.atoms.size() || structure.atoms.empty())
    {
        removeIpfSidecar(filename);
        return true;
    }

    const std::string path = ipfSidecarPath(filename);
    std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
    if (!out)
        return false;

    out << "ATOMFORGE_IPF_V1\n";
    out << structure.atoms.size() << "\n";
    out << std::setprecision(17);
    for (size_t i = 0; i < structure.atoms.size(); ++i)
    {
        const AtomSite& atom = structure.atoms[i];
        const std::array<float, 3>& color = structure.grainColors[i];
        out << atom.atomicNumber << ' '
            << atom.x << ' ' << atom.y << ' ' << atom.z << ' '
            << color[0] << ' ' << color[1] << ' ' << color[2] << '\n';
    }

    return out.good();
}

bool loadIpfSidecarRecords(const std::string& filename,
                          std::vector<IpfRecord>& records)
{
    records.clear();

    std::ifstream in;
    const std::string path = ipfSidecarPath(filename);
    in.open(path.c_str());
    if (!in)
    {
        in.clear();
        const std::string legacyPath = ipfLegacySidecarPath(filename);
        if (legacyPath == path)
            return false;
        in.open(legacyPath.c_str());
        if (!in)
            return false;
    }

    std::string header;
    std::getline(in, header);
    const bool isV1 = (header == "ATOMFORGE_IPF_V1");
    const bool isV2 = (header == "ATOMFORGE_IPF_V2");
    if (!isV1 && !isV2)
        return false;

    size_t count = 0;
    in >> count;
    if (!in)
        return false;

    records.reserve(count);
    for (size_t i = 0; i < count; ++i)
    {
        int atomicNumber = 0;
        double x = 0.0, y = 0.0, z = 0.0;
        float r = 0.0f, g = 0.0f, b = 0.0f;
        if (isV1)
        {
            in >> atomicNumber >> x >> y >> z >> r >> g >> b;
        }
        else
        {
            int ignoredRegionId = -1;
            in >> atomicNumber >> x >> y >> z >> r >> g >> b >> ignoredRegionId;
        }
        if (!in)
            return false;

        IpfRecord record;
        record.key = makeIpfKey(atomicNumber, x, y, z);
        record.color = {{r, g, b}};
        records.push_back(record);
    }

    return true;
}

bool restoreIpfSidecar(const std::string& filename, Structure& structure)
{
    if (structure.atoms.empty())
        return false;

    std::vector<IpfRecord> records;
    if (!loadIpfSidecarRecords(filename, records))
        return false;
    if (records.size() != structure.atoms.size())
        return false;

    std::unordered_multimap<IpfKey, size_t, IpfKeyHash> byKey;
    byKey.reserve(structure.atoms.size());
    for (size_t i = 0; i < structure.atoms.size(); ++i)
    {
        const AtomSite& atom = structure.atoms[i];
        byKey.insert(std::make_pair(makeIpfKey(atom.atomicNumber, atom.x, atom.y, atom.z), i));
    }

    std::vector<std::array<float, 3>> restored(structure.atoms.size(), {{0.0f, 0.0f, 0.0f}});
    std::vector<bool> assigned(structure.atoms.size(), false);
    size_t matched = 0;

    for (size_t i = 0; i < records.size(); ++i)
    {
        auto range = byKey.equal_range(records[i].key);
        size_t chosen = structure.atoms.size();
        for (auto it = range.first; it != range.second; ++it)
        {
            if (!assigned[it->second])
            {
                chosen = it->second;
                break;
            }
        }
        if (chosen == structure.atoms.size())
            break;

        restored[chosen] = records[i].color;
        assigned[chosen] = true;
        ++matched;
    }

    if (matched != structure.atoms.size())
    {
        // Fallback for formats that preserve atom order but perturb coordinates.
        for (size_t i = 0; i < structure.atoms.size(); ++i)
        {
            restored[i] = records[i].color;
        }
    }

    structure.grainColors.swap(restored);
    return true;
}

std::array<float, 3> ipfColorFromDirection(const glm::vec3& direction)
{
    float h = std::abs(direction.x);
    float k = std::abs(direction.y);
    float l = std::abs(direction.z);

    if (h < k) std::swap(h, k);
    if (h < l) std::swap(h, l);
    if (k < l) std::swap(k, l);

    const float len = std::sqrt(h * h + k * k + l * l);
    if (len < 1e-10f)
        return {{0.5f, 0.5f, 0.5f}};

    h /= len;
    k /= len;
    l /= len;

    float r = h - k;
    float g = k - l;
    float b = l * 1.7320508f;

    const float maxC = std::max({r, g, b, 1e-6f});
    r /= maxC;
    g /= maxC;
    b /= maxC;

    return {{r, g, b}};
}

bool recoverIpfFromGeometry(Structure& structure)
{
    if (structure.atoms.size() < 8)
        return false;

    constexpr float kCellSize = 4.0f;
    constexpr float kMaxNeighborDist = 4.2f;
    constexpr float kMinNeighborDist = 0.2f;
    constexpr int kTargetNeighbors = 16;

    std::vector<glm::vec3> positions(structure.atoms.size());
    positions.reserve(structure.atoms.size());
    for (size_t i = 0; i < structure.atoms.size(); ++i)
    {
        positions[i] = glm::vec3((float)structure.atoms[i].x,
                                 (float)structure.atoms[i].y,
                                 (float)structure.atoms[i].z);
    }

    auto gridCell = [&](const glm::vec3& pos) -> glm::ivec3 {
        return glm::ivec3((int)std::floor(pos.x / kCellSize),
                          (int)std::floor(pos.y / kCellSize),
                          (int)std::floor(pos.z / kCellSize));
    };

    std::unordered_map<glm::ivec3, std::vector<int>, Vec3iHash> grid;
    grid.reserve(structure.atoms.size());
    for (int i = 0; i < (int)positions.size(); ++i)
        grid[gridCell(positions[i])].push_back(i);

    std::vector<glm::vec3> crystalDirs(structure.atoms.size(), glm::vec3(0.0f, 0.0f, 1.0f));
    std::vector<std::vector<int>> neighborIds(structure.atoms.size());
    std::vector<bool> hasDir(structure.atoms.size(), false);
    int recoveredCount = 0;

    for (int i = 0; i < (int)positions.size(); ++i)
    {
        const glm::ivec3 c = gridCell(positions[i]);
        std::vector<std::pair<float, int>> neighbors;
        neighbors.reserve(64);

        for (int dx = -1; dx <= 1; ++dx)
        {
            for (int dy = -1; dy <= 1; ++dy)
            {
                for (int dz = -1; dz <= 1; ++dz)
                {
                    const glm::ivec3 cc(c.x + dx, c.y + dy, c.z + dz);
                    if (auto it = grid.find(cc); it != grid.end())
                    for (int idx : it->second)
                    {
                        if (idx == i)
                            continue;
                        const glm::vec3 d = positions[idx] - positions[i];
                        const float dist = glm::length(d);
                        if (dist < kMinNeighborDist || dist > kMaxNeighborDist)
                            continue;
                        neighbors.push_back(std::make_pair(dist, idx));
                    }
                }
            }
        }

        if (neighbors.size() < 3)
            continue;

        std::sort(neighbors.begin(), neighbors.end(),
                        [](const std::pair<float, int>& a,
                            const std::pair<float, int>& b) {
            return a.first < b.first;
        });

        const int useCount = std::min((int)neighbors.size(), kTargetNeighbors);
        neighborIds[i].reserve(useCount);
        std::vector<glm::vec3> unitVecs;
        unitVecs.reserve(useCount);
        for (int n = 0; n < useCount; ++n)
        {
            const int ni = neighbors[n].second;
            const glm::vec3 d = positions[ni] - positions[i];
            const float dn = glm::length(d);
            if (dn < 1e-8f)
                continue;
            const glm::vec3 u = d / dn;

            neighborIds[i].push_back(ni);
            unitVecs.push_back(u);
        }

        if (unitVecs.size() < 3)
            continue;

        // Build a local orthonormal frame from two non-collinear neighbor vectors.
        const glm::vec3 u1 = glm::normalize(unitVecs[0]);
        bool foundSecond = false;
        glm::vec3 u2raw(0.0f);
        for (size_t n = 1; n < unitVecs.size(); ++n)
        {
            const float cosine = std::abs(glm::dot(u1, unitVecs[n]));
            if (cosine < 0.85f)
            {
                u2raw = unitVecs[n];
                foundSecond = true;
                break;
            }
        }
        if (!foundSecond)
            continue;

        glm::vec3 u2 = u2raw - glm::dot(u2raw, u1) * u1;
        const float u2n = glm::length(u2);
        if (u2n < 1e-6f)
            continue;
        u2 /= u2n;

        glm::vec3 u3 = glm::cross(u1, u2);
        const float u3n = glm::length(u3);
        if (u3n < 1e-6f)
            continue;
        u3 /= u3n;

        const glm::mat3 R(u1, u2, u3);
        const glm::vec3 zDir(0.0f, 0.0f, 1.0f);
        const glm::vec3 crystalDir = glm::transpose(R) * zDir;
        if (glm::length(crystalDir) < 1e-6f)
            continue;

        crystalDirs[i] = glm::normalize(crystalDir);
        hasDir[i] = true;
        ++recoveredCount;
    }

    if (recoveredCount < (int)structure.atoms.size() / 4)
        return false;

    // Fill missing directions using local neighbors only.
    for (int pass = 0; pass < 3; ++pass)
    {
        bool changed = false;
        for (size_t i = 0; i < crystalDirs.size(); ++i)
        {
            if (hasDir[i])
                continue;

            glm::vec3 sum(0.0f, 0.0f, 0.0f);
            const std::vector<int>& nbrs = neighborIds[i];
            for (int ni : nbrs)
            {
                if (!hasDir[ni])
                    continue;
                glm::vec3 a = crystalDirs[ni];
                if (glm::dot(a, sum) < 0.0f)
                    a = -a;
                sum += a;
            }

            if (glm::length(sum) > 1e-8f)
            {
                crystalDirs[i] = glm::normalize(sum);
                hasDir[i] = true;
                changed = true;
            }
        }
        if (!changed)
            break;
    }

    for (size_t i = 0; i < crystalDirs.size(); ++i)
    {
        if (!hasDir[i])
            crystalDirs[i] = glm::vec3(0.0f, 0.0f, 1.0f);
    }

    // Smooth orientation field inside neighborhoods.
    for (int iter = 0; iter < 4; ++iter)
    {
        std::vector<glm::vec3> nextDirs = crystalDirs;
        for (size_t i = 0; i < crystalDirs.size(); ++i)
        {
            glm::vec3 sum = crystalDirs[i] * 2.0f;
            const std::vector<int>& nbrs = neighborIds[i];
            for (int ni : nbrs)
            {
                glm::vec3 a = crystalDirs[ni];
                if (glm::dot(a, sum) < 0.0f)
                    a = -a;
                sum += a;
            }

            if (glm::length(sum) > 1e-8f)
            {
                nextDirs[i] = glm::normalize(sum);
            }
        }
        crystalDirs.swap(nextDirs);
    }

    // Segment into connected orientation regions and color per region
    // to avoid speckled atom-wise noise.
    const float cosThreshold = std::cos(20.0f * 3.1415926535f / 180.0f);
    std::vector<int> regionId(crystalDirs.size(), -1);
    std::vector<glm::vec3> regionMeans;
    std::vector<std::vector<int>> regionMembers;

    for (size_t seed = 0; seed < crystalDirs.size(); ++seed)
    {
        if (regionId[seed] >= 0)
            continue;

        const int rid = (int)regionMeans.size();
        regionMeans.push_back(glm::normalize(crystalDirs[seed]));
        regionMembers.push_back(std::vector<int>());

        std::queue<int> q;
        q.push((int)seed);
        regionId[seed] = rid;
        regionMembers[rid].push_back((int)seed);

        glm::vec3 sum = crystalDirs[seed];

        while (!q.empty())
        {
            const int u = q.front();
            q.pop();

            glm::vec3 ref = (glm::length(sum) > 1e-8f)
                          ? glm::normalize(sum)
                          : glm::normalize(crystalDirs[u]);

            const std::vector<int>& nbrs = neighborIds[u];
            for (int v : nbrs)
            {
                if (regionId[v] >= 0)
                    continue;

                const glm::vec3 dv = glm::normalize(crystalDirs[v]);
                if (std::abs(glm::dot(dv, ref)) < cosThreshold)
                    continue;

                regionId[v] = rid;
                regionMembers[rid].push_back(v);
                q.push(v);

                glm::vec3 aligned = dv;
                if (glm::dot(aligned, sum) < 0.0f)
                    aligned = -aligned;
                sum += aligned;
            }
        }

        if (glm::length(sum) > 1e-8f)
            regionMeans[rid] = glm::normalize(sum);
    }

    const int nAtoms = (int)crystalDirs.size();
    const int minRegionSize = std::max(32, nAtoms / 180);
    std::vector<bool> isLarge(regionMeans.size(), false);
    for (size_t r = 0; r < regionMembers.size(); ++r)
        isLarge[r] = (int)regionMembers[r].size() >= minRegionSize;

    bool anyLarge = false;
    for (size_t r = 0; r < isLarge.size(); ++r)
        anyLarge = anyLarge || isLarge[r];

    if (!anyLarge && !regionMembers.empty())
    {
        size_t largest = 0;
        for (size_t r = 1; r < regionMembers.size(); ++r)
        {
            if (regionMembers[r].size() > regionMembers[largest].size())
                largest = r;
        }
        isLarge[largest] = true;
    }

    // Merge small noisy regions as whole components into large regions.
    std::vector<int> regionRemap(regionMeans.size(), -1);
    for (size_t r = 0; r < regionMeans.size(); ++r)
        regionRemap[r] = (int)r;

    std::vector<std::unordered_map<int, int>> regionAdjLarge(regionMeans.size());
    for (int i = 0; i < nAtoms; ++i)
    {
        const int ri = regionId[i];
        if (ri < 0)
            continue;

        const std::vector<int>& nbrs = neighborIds[i];
        for (int ni : nbrs)
        {
            if (ni <= i)
                continue;
            const int rj = regionId[ni];
            if (rj < 0 || rj == ri)
                continue;

            if (!isLarge[ri] && isLarge[rj])
                regionAdjLarge[ri][rj] += 1;
            if (!isLarge[rj] && isLarge[ri])
                regionAdjLarge[rj][ri] += 1;
        }
    }

    for (size_t r = 0; r < regionMeans.size(); ++r)
    {
        if (isLarge[r])
            continue;

        int bestRegion = -1;
        int bestCount = -1;
        float bestOrient = -1.0f;

        for (std::unordered_map<int, int>::const_iterator it = regionAdjLarge[r].begin();
             it != regionAdjLarge[r].end(); ++it)
        {
            const int nr = it->first;
            const int cnt = it->second;
            const float orient = std::abs(glm::dot(regionMeans[r], regionMeans[nr]));
            if (cnt > bestCount || (cnt == bestCount && orient > bestOrient))
            {
                bestRegion = nr;
                bestCount = cnt;
                bestOrient = orient;
            }
        }

        if (bestRegion < 0)
        {
            for (size_t nr = 0; nr < regionMeans.size(); ++nr)
            {
                if (!isLarge[nr])
                    continue;
                const float orient = std::abs(glm::dot(regionMeans[r], regionMeans[nr]));
                if (orient > bestOrient)
                {
                    bestOrient = orient;
                    bestRegion = (int)nr;
                }
            }
        }

        if (bestRegion >= 0)
            regionRemap[r] = bestRegion;
    }

    for (int i = 0; i < nAtoms; ++i)
    {
        const int rid = regionId[i];
        if (rid >= 0)
            regionId[i] = regionRemap[rid];
    }

    // Recompute stable region means after merging.
    std::vector<glm::vec3> finalMeans(regionMeans.size(), glm::vec3(0.0f));
    std::vector<int> finalCounts(regionMeans.size(), 0);
    for (int i = 0; i < nAtoms; ++i)
    {
        const int rid = regionId[i];
        if (rid < 0)
            continue;
        glm::vec3 d = glm::normalize(crystalDirs[i]);
        if (finalCounts[rid] > 0 && glm::dot(d, finalMeans[rid]) < 0.0f)
            d = -d;
        finalMeans[rid] += d;
        finalCounts[rid] += 1;
    }

    for (size_t r = 0; r < finalMeans.size(); ++r)
    {
        if (finalCounts[r] > 0 && glm::length(finalMeans[r]) > 1e-8f)
            finalMeans[r] = glm::normalize(finalMeans[r]);
        else
            finalMeans[r] = glm::vec3(0.0f, 0.0f, 1.0f);
    }

    std::vector<std::array<float, 3>> recovered(structure.atoms.size());
    for (size_t i = 0; i < crystalDirs.size(); ++i)
    {
        const int rid = (regionId[i] >= 0) ? regionId[i] : 0;
        const glm::vec3 dir = (rid >= 0 && rid < (int)finalMeans.size())
                            ? finalMeans[rid]
                            : glm::vec3(0.0f, 0.0f, 1.0f);
        recovered[i] = ipfColorFromDirection(dir);
    }

    structure.grainColors.swap(recovered);
    return true;
}

std::string normalizeSeparators(const std::string& path)
{
    std::string out = path;
    std::replace(out.begin(), out.end(), '\\', '/');
    return out;
}

std::string ensureTrailingSlash(const std::string& path)
{
    if (path.empty())
        return path;
    if (path.back() == '/' || path.back() == '\\')
        return path;
    return path + "/";
}

std::string parentDirectory(const std::string& path)
{
    if (path.empty())
        return std::string();

    std::string out = normalizeSeparators(path);
    while (out.size() > 1 && out.back() == '/')
        out.pop_back();

    std::size_t pos = out.find_last_of('/');
    if (pos == std::string::npos)
        return std::string();
    if (pos == 0)
        return out.substr(0, 1);
    return out.substr(0, pos);
}

std::string joinPath(const std::string& base, const std::string& name)
{
    if (base.empty())
        return name;
    if (base.back() == '/' || base.back() == '\\')
        return base + name;
    return base + "/" + name;
}

bool directoryExists(const std::string& path)
{
    struct stat st;
    if (stat(path.c_str(), &st) != 0)
        return false;
#ifdef _WIN32
    return (st.st_mode & _S_IFDIR) != 0;
#else
    return S_ISDIR(st.st_mode);
#endif
}

bool fileExists(const std::string& path)
{
    struct stat st;
    if (stat(path.c_str(), &st) != 0)
        return false;
#ifdef _WIN32
    return (st.st_mode & _S_IFREG) != 0;
#else
    return S_ISREG(st.st_mode);
#endif
}

bool hasSuffixCaseInsensitive(const std::string& value, const char* suffix)
{
    const std::size_t valueLen = value.size();
    const std::size_t suffixLen = std::strlen(suffix);
    if (valueLen < suffixLen)
        return false;

    const std::size_t offset = valueLen - suffixLen;
    for (std::size_t i = 0; i < suffixLen; ++i)
    {
        const char a = (char)std::tolower((unsigned char)value[offset + i]);
        const char b = (char)std::tolower((unsigned char)suffix[i]);
        if (a != b)
            return false;
    }
    return true;
}

bool directoryLooksLikeOpenBabelPluginDir(const std::string& path)
{
    // Fast path for older/known plugin names.
    if (fileExists(joinPath(path, "formats.obf")))
        return true;

#ifdef _WIN32
    WIN32_FIND_DATAA findData;
    HANDLE hFind = FindFirstFileA(joinPath(path, "*").c_str(), &findData);
    if (hFind == INVALID_HANDLE_VALUE)
        return false;

    bool foundPlugin = false;
    do
    {
        const char* name = findData.cFileName;
        if (!name || std::strcmp(name, ".") == 0 || std::strcmp(name, "..") == 0)
            continue;
        if ((findData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0)
            continue;

        const std::string filename(name);
        if (hasSuffixCaseInsensitive(filename, ".obf") || hasSuffixCaseInsensitive(filename, ".dll"))
        {
            foundPlugin = true;
            break;
        }
    } while (FindNextFileA(hFind, &findData) != 0);

    FindClose(hFind);
    return foundPlugin;
#else
    DIR* dir = opendir(path.c_str());
    if (!dir)
        return false;

    bool foundPlugin = false;
    struct dirent* entry = nullptr;
    while ((entry = readdir(dir)) != nullptr)
    {
        const char* name = entry->d_name;
        if (!name || std::strcmp(name, ".") == 0 || std::strcmp(name, "..") == 0)
            continue;

        const std::string filename(name);
        if (hasSuffixCaseInsensitive(filename, ".obf") ||
            hasSuffixCaseInsensitive(filename, ".so") ||
            hasSuffixCaseInsensitive(filename, ".dylib"))
        {
            foundPlugin = true;
            break;
        }
    }

    closedir(dir);
    return foundPlugin;
#endif
}

void appendUnique(std::vector<std::string>& values, const std::string& value)
{
    if (value.empty())
        return;
    const std::string normalized = normalizeSeparators(value);
    if (std::find(values.begin(), values.end(), normalized) == values.end())
        values.push_back(normalized);
}

void appendVersionedPluginCandidates(const std::string& openBabelRoot, std::vector<std::string>& candidates);

std::string executableDirectory()
{
#ifdef _WIN32
    char exePath[MAX_PATH] = {0};
    DWORD len = GetModuleFileNameA(nullptr, exePath, MAX_PATH);
    if (len > 0 && len < MAX_PATH)
        return parentDirectory(exePath);
#else
    char exePath[4096] = {0};
    const ssize_t len = readlink("/proc/self/exe", exePath, sizeof(exePath) - 1);
    if (len > 0)
    {
        exePath[len] = '\0';
        return parentDirectory(exePath);
    }
#endif
    return std::string();
}

#ifdef _WIN32
void appendWindowsPathVariants(const std::string& path, std::vector<std::string>& values)
{
    if (path.empty())
        return;

    const std::string normalized = normalizeSeparators(path);
    appendUnique(values, normalized);

    // Convert MSYS-style absolute paths to native Windows candidates.
    if (normalized.size() > 1 && normalized[0] == '/')
    {
        appendUnique(values, "C:/msys64" + normalized);
        appendUnique(values, "C:/msys2" + normalized);
    }
}

void appendCandidatesFromPathEnvironment(std::vector<std::string>& candidates)
{
    const char* pathEnv = std::getenv("PATH");
    if (!pathEnv || !*pathEnv)
        return;

    const std::string pathList(pathEnv);
    std::size_t start = 0;
    while (start <= pathList.size())
    {
        const std::size_t end = pathList.find(';', start);
        const std::string entry = normalizeSeparators(pathList.substr(start, end - start));

        if (!entry.empty())
        {
            // Typical Open Babel runtime path is .../<prefix>/bin.
            if (hasSuffixCaseInsensitive(entry, "/bin"))
            {
                const std::string prefix = parentDirectory(entry);
                if (!prefix.empty())
                    appendVersionedPluginCandidates(joinPath(joinPath(prefix, "lib"), "openbabel"), candidates);
            }

            // Sometimes plugin dir itself is directly on PATH.
            appendVersionedPluginCandidates(entry, candidates);
        }

        if (end == std::string::npos)
            break;
        start = end + 1;
    }
}
#endif

void appendVersionedPluginCandidates(const std::string& openBabelRoot, std::vector<std::string>& candidates)
{
    if (openBabelRoot.empty() || !directoryExists(openBabelRoot))
        return;

    // Allow direct plugin dir input as well.
    appendUnique(candidates, openBabelRoot);

#ifdef _WIN32
    WIN32_FIND_DATAA findData;
    HANDLE hFind = FindFirstFileA(joinPath(openBabelRoot, "*").c_str(), &findData);
    if (hFind == INVALID_HANDLE_VALUE)
        return;

    do
    {
        if ((findData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) == 0)
            continue;

        const char* name = findData.cFileName;
        if (!name || std::strcmp(name, ".") == 0 || std::strcmp(name, "..") == 0)
            continue;

        if (!std::isdigit((unsigned char)name[0]))
            continue;

        appendUnique(candidates, joinPath(openBabelRoot, name));
    } while (FindNextFileA(hFind, &findData) != 0);

    FindClose(hFind);
#else
    DIR* dir = opendir(openBabelRoot.c_str());
    if (!dir)
        return;

    struct dirent* entry = nullptr;
    while ((entry = readdir(dir)) != nullptr)
    {
        const char* name = entry->d_name;
        if (!name || std::strcmp(name, ".") == 0 || std::strcmp(name, "..") == 0)
            continue;

        if (!std::isdigit((unsigned char)name[0]))
            continue;

        appendUnique(candidates, joinPath(openBabelRoot, name));
    }

    closedir(dir);
#endif
}

void setBabelLibDir(const std::string& path)
{
    const std::string normalized = ensureTrailingSlash(normalizeSeparators(path));
#ifdef _WIN32
    _putenv_s("BABEL_LIBDIR", normalized.c_str());
#else
    setenv("BABEL_LIBDIR", normalized.c_str(), 1);
#endif
}

void ensureOpenBabelPlugins()
{
    static bool initialized = false;
    if (initialized)
        return;
    initialized = true;

    const char* currentLibDir = std::getenv("BABEL_LIBDIR");
    if (currentLibDir && *currentLibDir)
    {
        const std::string envPath = normalizeSeparators(currentLibDir);
        if (directoryExists(envPath) && directoryLooksLikeOpenBabelPluginDir(envPath))
        {
            // Respect valid user-provided plugin path first.
            setBabelLibDir(envPath);
            OpenBabel::OBPlugin::LoadAllPlugins();

            OpenBabel::OBConversion probe;
            if (probe.SetInFormat("cif") || probe.SetInFormat("vasp"))
                return;
        }
    }

    std::vector<std::string> candidates;

#ifdef ATOMFORGE_OPENBABEL_PLUGIN_DIR
#ifdef _WIN32
    appendWindowsPathVariants(ATOMFORGE_OPENBABEL_PLUGIN_DIR, candidates);
#else
    appendUnique(candidates, ATOMFORGE_OPENBABEL_PLUGIN_DIR);
#endif
#endif

    // Portable layout support: bundled plugins live next to the executable
    // under an "openbabel" directory.
    const std::string exeDir = executableDirectory();
    if (!exeDir.empty())
    {
        appendVersionedPluginCandidates(joinPath(exeDir, "openbabel"), candidates);

        // Typical system layout fallback when executable is in <prefix>/bin.
        const std::string prefixDir = parentDirectory(exeDir);
        if (!prefixDir.empty())
            appendVersionedPluginCandidates(joinPath(joinPath(prefixDir, "lib"), "openbabel"), candidates);
    }

#ifdef _WIN32
    HMODULE obModule = nullptr;
    if (GetModuleHandleExA(
            GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
            reinterpret_cast<LPCSTR>(&OpenBabel::OBPlugin::LoadAllPlugins),
            &obModule))
    {
        char modulePath[MAX_PATH] = {0};
        DWORD len = GetModuleFileNameA(obModule, modulePath, MAX_PATH);
        if (len > 0 && len < MAX_PATH)
        {
            const std::string binDir = parentDirectory(modulePath);
            const std::string prefixDir = parentDirectory(binDir);
            if (!prefixDir.empty())
                appendVersionedPluginCandidates(joinPath(joinPath(prefixDir, "lib"), "openbabel"), candidates);
        }
    }

    appendCandidatesFromPathEnvironment(candidates);

    appendVersionedPluginCandidates("C:/msys64/ucrt64/lib/openbabel", candidates);
    appendVersionedPluginCandidates("C:/msys64/mingw64/lib/openbabel", candidates);
    appendVersionedPluginCandidates("C:/msys2/ucrt64/lib/openbabel", candidates);
    appendVersionedPluginCandidates("C:/msys2/mingw64/lib/openbabel", candidates);

    // Also try MSYS-style paths in case the process can resolve them.
    appendVersionedPluginCandidates("/ucrt64/lib/openbabel", candidates);
    appendVersionedPluginCandidates("/mingw64/lib/openbabel", candidates);
#else
    appendVersionedPluginCandidates("/usr/lib/openbabel", candidates);
    appendVersionedPluginCandidates("/usr/lib64/openbabel", candidates);
    appendVersionedPluginCandidates("/usr/local/lib/openbabel", candidates);
    appendVersionedPluginCandidates("/usr/local/lib64/openbabel", candidates);
    appendVersionedPluginCandidates("/usr/lib/x86_64-linux-gnu/openbabel", candidates);
#endif

    bool pluginsReady = false;
    for (std::size_t i = 0; i < candidates.size(); ++i)
    {
        if (directoryExists(candidates[i]) && directoryLooksLikeOpenBabelPluginDir(candidates[i]))
        {
            setBabelLibDir(candidates[i]);
            OpenBabel::OBPlugin::LoadAllPlugins();

            OpenBabel::OBConversion probe;
            if (probe.SetInFormat("cif") || probe.SetInFormat("vasp"))
            {
                pluginsReady = true;
                break;
            }
        }
    }

    if (!pluginsReady)
        OpenBabel::OBPlugin::LoadAllPlugins();
}

class ScopedObWarningSilencer
{
public:
    ScopedObWarningSilencer()
        : m_previousLevel(OpenBabel::obErrorLog.GetOutputLevel())
    {
        // Keep Open Babel errors visible, hide warnings/info that are often
        // emitted for partially-defined space groups in valid CIF files.
        OpenBabel::obErrorLog.SetOutputLevel(OpenBabel::obError);
    }

    ~ScopedObWarningSilencer()
    {
        OpenBabel::obErrorLog.SetOutputLevel(m_previousLevel);
    }

private:
    OpenBabel::obMessageLevel m_previousLevel;
};

std::string toLowerCopy(const std::string& value)
{
    std::string out = value;
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = (char)std::tolower((unsigned char)out[i]);
    return out;
}

std::string extractExtension(const std::string& filename)
{
    const std::size_t slashPos = filename.find_last_of("/\\");
    const std::size_t dotPos = filename.find_last_of('.');
    if (dotPos == std::string::npos)
        return std::string();
    if (slashPos != std::string::npos && dotPos < slashPos)
        return std::string();
    return filename.substr(dotPos);
}

bool isSupportedExtension(const std::string& extLower)
{
    static const char* kSupportedExtensions[] = {
        ".cif",
        ".mol",
        ".pdb",
        ".xyz",
        ".sdf",
        ".vasp",
        ".mol2",
        ".pwi",
        ".gjf"
    };

    for (std::size_t i = 0; i < sizeof(kSupportedExtensions) / sizeof(kSupportedExtensions[0]); ++i)
    {
        if (extLower == kSupportedExtensions[i])
            return true;
    }
    return false;
}

std::string supportedExtensionsSummary()
{
    return ".cif, .mol, .pdb, .xyz, .sdf, .vasp, .mol2, .pwi, .gjf";
}

bool filenameIsPoscarOrContcar(const std::string& filename)
{
    // Extract basename (after last slash or backslash)
    const std::size_t slashPos = filename.find_last_of("/\\");
    const std::string basename = (slashPos != std::string::npos)
                               ? filename.substr(slashPos + 1)
                               : filename;
    return basename.compare(0, 6, "POSCAR") == 0
        || basename.compare(0, 7, "CONTCAR") == 0;
}

void wrapAtomsIntoPrimaryCell(Structure& structure)
{
    if (!structure.hasUnitCell || structure.atoms.empty())
        return;

    const glm::vec3 origin(0.0f, 0.0f, 0.0f);
    const glm::vec3 a((float)structure.cellVectors[0][0], (float)structure.cellVectors[0][1], (float)structure.cellVectors[0][2]);
    const glm::vec3 b((float)structure.cellVectors[1][0], (float)structure.cellVectors[1][1], (float)structure.cellVectors[1][2]);
    const glm::vec3 c((float)structure.cellVectors[2][0], (float)structure.cellVectors[2][1], (float)structure.cellVectors[2][2]);

    const glm::mat3 cellMat(a, b, c);
    const float det = glm::determinant(cellMat);
    if (std::abs(det) <= 1e-8f)
        return;

    const glm::mat3 invCellMat = glm::inverse(cellMat);
    constexpr float kWrapTol = 1e-5f;

    for (auto& atom : structure.atoms)
    {
        glm::vec3 pos((float)atom.x, (float)atom.y, (float)atom.z);
        glm::vec3 frac = invCellMat * (pos - origin);

        frac.x -= std::floor(frac.x);
        frac.y -= std::floor(frac.y);
        frac.z -= std::floor(frac.z);

        if (std::abs(frac.x) <= kWrapTol || std::abs(1.0f - frac.x) <= kWrapTol) frac.x = 0.0f;
        if (std::abs(frac.y) <= kWrapTol || std::abs(1.0f - frac.y) <= kWrapTol) frac.y = 0.0f;
        if (std::abs(frac.z) <= kWrapTol || std::abs(1.0f - frac.z) <= kWrapTol) frac.z = 0.0f;

        const glm::vec3 wrapped = origin + frac.x * a + frac.y * b + frac.z * c;
        atom.x = wrapped.x;
        atom.y = wrapped.y;
        atom.z = wrapped.z;
    }
}
}


void getDefaultElementColor(int Z,float& r,float& g,float& b)
{
    // CPK coloring scheme for all elements (1..118).
    // Values taken from common molecular visualization defaults.
    static const std::array<std::array<float,3>, 119> colors = {
        std::array<float,3>{0.0f, 0.0f, 0.0f}, // placeholder for index 0
        std::array<float,3>{1.00f, 1.00f, 1.00f}, // 1  H
        std::array<float,3>{0.85f, 1.00f, 1.00f}, // 2  He
        std::array<float,3>{0.80f, 0.50f, 1.00f}, // 3  Li
        std::array<float,3>{0.76f, 1.00f, 0.00f}, // 4  Be
        std::array<float,3>{1.00f, 0.71f, 0.71f}, // 5  B
        std::array<float,3>{0.20f, 0.20f, 0.20f}, // 6  C
        std::array<float,3>{0.00f, 0.00f, 1.00f}, // 7  N
        std::array<float,3>{1.00f, 0.00f, 0.00f}, // 8  O
        std::array<float,3>{0.00f, 1.00f, 0.00f}, // 9  F
        std::array<float,3>{0.70f, 0.89f, 0.96f}, // 10 Ne
        std::array<float,3>{0.67f, 0.36f, 0.95f}, // 11 Na
        std::array<float,3>{0.54f, 1.00f, 0.00f}, // 12 Mg
        std::array<float,3>{0.75f, 0.65f, 0.65f}, // 13 Al
        std::array<float,3>{0.94f, 0.78f, 0.62f}, // 14 Si
        std::array<float,3>{1.00f, 0.50f, 0.00f}, // 15 P
        std::array<float,3>{1.00f, 1.00f, 0.00f}, // 16 S
        std::array<float,3>{0.00f, 1.00f, 0.00f}, // 17 Cl
        std::array<float,3>{0.50f, 0.82f, 0.89f}, // 18 Ar
        std::array<float,3>{0.56f, 0.00f, 1.00f}, // 19 K
        std::array<float,3>{0.24f, 1.00f, 0.00f}, // 20 Ca
        std::array<float,3>{0.90f, 0.90f, 0.90f}, // 21 Sc
        std::array<float,3>{0.75f, 0.76f, 0.78f}, // 22 Ti
        std::array<float,3>{0.65f, 0.65f, 0.67f}, // 23 V
        std::array<float,3>{0.54f, 0.60f, 0.78f}, // 24 Cr
        std::array<float,3>{0.61f, 0.48f, 0.78f}, // 25 Mn
        std::array<float,3>{0.88f, 0.40f, 0.20f}, // 26 Fe
        std::array<float,3>{0.88f, 0.38f, 0.20f}, // 27 Co
        std::array<float,3>{0.78f, 0.79f, 0.78f}, // 28 Ni
        std::array<float,3>{0.78f, 0.50f, 0.20f}, // 29 Cu
        std::array<float,3>{0.49f, 0.50f, 0.69f}, // 30 Zn
        std::array<float,3>{0.76f, 0.56f, 0.56f}, // 31 Ga
        std::array<float,3>{0.40f, 0.56f, 0.56f}, // 32 Ge
        std::array<float,3>{0.74f, 0.50f, 0.89f}, // 33 As
        std::array<float,3>{1.00f, 0.63f, 0.00f}, // 34 Se
        std::array<float,3>{0.65f, 0.16f, 0.16f}, // 35 Br
        std::array<float,3>{0.36f, 0.72f, 0.82f}, // 36 Kr
        std::array<float,3>{0.44f, 0.18f, 0.69f}, // 37 Rb
        std::array<float,3>{0.00f, 1.00f, 0.00f}, // 38 Sr
        std::array<float,3>{0.58f, 1.00f, 1.00f}, // 39 Y
        std::array<float,3>{0.58f, 0.88f, 0.88f}, // 40 Zr
        std::array<float,3>{0.45f, 0.76f, 0.79f}, // 41 Nb
        std::array<float,3>{0.33f, 0.71f, 0.71f}, // 42 Mo
        std::array<float,3>{0.23f, 0.62f, 0.62f}, // 43 Tc
        std::array<float,3>{0.14f, 0.56f, 0.56f}, // 44 Ru
        std::array<float,3>{0.04f, 0.49f, 0.55f}, // 45 Rh
        std::array<float,3>{0.00f, 0.41f, 0.52f}, // 46 Pd
        std::array<float,3>{0.75f, 0.75f, 0.75f}, // 47 Ag
        std::array<float,3>{1.00f, 0.85f, 0.00f}, // 48 Cd
        std::array<float,3>{0.65f, 0.46f, 0.45f}, // 49 In
        std::array<float,3>{0.40f, 0.50f, 0.50f}, // 50 Sn
        std::array<float,3>{0.62f, 0.39f, 0.71f}, // 51 Sb
        std::array<float,3>{0.83f, 0.50f, 0.18f}, // 52 Te
        std::array<float,3>{0.58f, 0.00f, 0.58f}, // 53 I
        std::array<float,3>{0.26f, 0.62f, 0.69f}, // 54 Xe
        std::array<float,3>{0.34f, 0.09f, 0.56f}, // 55 Cs
        std::array<float,3>{0.00f, 0.79f, 0.00f}, // 56 Ba
        std::array<float,3>{0.44f, 0.83f, 1.00f}, // 57 La
        std::array<float,3>{1.00f, 1.00f, 0.78f}, // 58 Ce
        std::array<float,3>{0.85f, 1.00f, 0.78f}, // 59 Pr
        std::array<float,3>{0.78f, 1.00f, 0.78f}, // 60 Nd
        std::array<float,3>{0.64f, 1.00f, 0.78f}, // 61 Pm
        std::array<float,3>{0.56f, 1.00f, 0.78f}, // 62 Sm
        std::array<float,3>{0.38f, 1.00f, 0.78f}, // 63 Eu
        std::array<float,3>{0.27f, 1.00f, 0.78f}, // 64 Gd
        std::array<float,3>{0.19f, 1.00f, 0.78f}, // 65 Tb
        std::array<float,3>{0.12f, 1.00f, 0.78f}, // 66 Dy
        std::array<float,3>{0.10f, 1.00f, 0.61f}, // 67 Ho
        std::array<float,3>{0.10f, 0.86f, 0.56f}, // 68 Er
        std::array<float,3>{0.10f, 0.77f, 0.52f}, // 69 Tm
        std::array<float,3>{0.12f, 0.70f, 0.47f}, // 70 Yb
        std::array<float,3>{0.12f, 0.63f, 0.42f}, // 71 Lu
        std::array<float,3>{0.30f, 0.61f, 0.61f}, // 72 Hf
        std::array<float,3>{0.30f, 0.58f, 0.58f}, // 73 Ta
        std::array<float,3>{0.13f, 0.58f, 0.58f}, // 74 W
        std::array<float,3>{0.14f, 0.54f, 0.54f}, // 75 Re
        std::array<float,3>{0.14f, 0.51f, 0.51f}, // 76 Os
        std::array<float,3>{0.09f, 0.46f, 0.52f}, // 77 Ir
        std::array<float,3>{0.82f, 0.82f, 0.88f}, // 78 Pt
        std::array<float,3>{1.00f, 0.82f, 0.14f}, // 79 Au
        std::array<float,3>{0.72f, 0.72f, 0.82f}, // 80 Hg
        std::array<float,3>{0.65f, 0.33f, 0.30f}, // 81 Tl
        std::array<float,3>{0.34f, 0.35f, 0.38f}, // 82 Pb
        std::array<float,3>{0.62f, 0.31f, 0.71f}, // 83 Bi
        std::array<float,3>{0.67f, 0.36f, 0.95f}, // 84 Po
        std::array<float,3>{0.46f, 0.32f, 0.28f}, // 85 At
        std::array<float,3>{0.26f, 0.42f, 0.47f}, // 86 Rn
        std::array<float,3>{0.26f, 0.00f, 0.40f}, // 87 Fr
        std::array<float,3>{0.00f, 0.42f, 0.00f}, // 88 Ra
        std::array<float,3>{0.00f, 0.70f, 1.00f}, // 89 Ac
        std::array<float,3>{0.00f, 0.73f, 1.00f}, // 90 Th
        std::array<float,3>{0.00f, 0.73f, 1.00f}, // 91 Pa
        std::array<float,3>{0.00f, 0.50f, 1.00f}, // 92 U
        std::array<float,3>{0.00f, 0.47f, 1.00f}, // 93 Np
        std::array<float,3>{0.00f, 0.40f, 1.00f}, // 94 Pu
        std::array<float,3>{0.00f, 0.38f, 1.00f}, // 95 Am
        std::array<float,3>{0.00f, 0.35f, 1.00f}, // 96 Cm
        std::array<float,3>{0.00f, 0.32f, 1.00f}, // 97 Bk
        std::array<float,3>{0.00f, 0.30f, 1.00f}, // 98 Cf
        std::array<float,3>{0.18f, 0.22f, 0.88f}, // 99 Es
        std::array<float,3>{0.20f, 0.20f, 0.90f}, // 100 Fm
        std::array<float,3>{0.20f, 0.20f, 0.88f}, // 101 Md
        std::array<float,3>{0.20f, 0.20f, 0.87f}, // 102 No
        std::array<float,3>{0.20f, 0.20f, 0.86f}, // 103 Lr
        std::array<float,3>{0.39f, 0.13f, 0.67f}, // 104 Rf
        std::array<float,3>{0.46f, 0.13f, 0.64f}, // 105 Db
        std::array<float,3>{0.54f, 0.16f, 0.59f}, // 106 Sg
        std::array<float,3>{0.61f, 0.18f, 0.55f}, // 107 Bh
        std::array<float,3>{0.68f, 0.20f, 0.50f}, // 108 Hs
        std::array<float,3>{0.74f, 0.20f, 0.46f}, // 109 Mt
        std::array<float,3>{0.78f, 0.21f, 0.42f}, // 110 Ds
        std::array<float,3>{0.82f, 0.22f, 0.38f}, // 111 Rg
        std::array<float,3>{0.86f, 0.23f, 0.34f}, // 112 Cn
        std::array<float,3>{0.90f, 0.24f, 0.30f}, // 113 Nh
        std::array<float,3>{0.94f, 0.26f, 0.26f}, // 114 Fl
        std::array<float,3>{0.98f, 0.28f, 0.22f}, // 115 Mc
        std::array<float,3>{1.00f, 0.30f, 0.18f}, // 116 Lv
        std::array<float,3>{1.00f, 0.32f, 0.14f}, // 117 Ts
        std::array<float,3>{1.00f, 0.34f, 0.10f}, // 118 Og
    };

    if (Z < 1) Z = 1;
    if (Z >= (int)colors.size()) Z = (int)colors.size() - 1;

    r = colors[Z][0];
    g = colors[Z][1];
    b = colors[Z][2];
}

bool isSupportedStructureFile(const std::string& filename)
{
    const std::string extLower = toLowerCopy(extractExtension(filename));
    if (isSupportedExtension(extLower))
        return true;
    return filenameIsPoscarOrContcar(filename);
}

bool loadStructureFromFile(const std::string& filename, Structure& structure, std::string& errorMessage)
{
    structure = Structure();
    errorMessage.clear();

    if (filename.empty())
    {
        errorMessage = "No file selected.";
        return false;
    }

    const std::string ext = extractExtension(filename);
    const std::string extLower = toLowerCopy(ext);
    const bool detectedPoscar = !isSupportedExtension(extLower) && filenameIsPoscarOrContcar(filename);
    if (!isSupportedExtension(extLower) && !detectedPoscar)
    {
        if (extLower.empty())
            errorMessage = "Unsupported file format (missing extension). Supported formats: " + supportedExtensionsSummary();
        else
            errorMessage = "Unsupported file format '" + ext + "'. Supported formats: " + supportedExtensionsSummary();
        return false;
    }

    struct stat fileStat;
    if (stat(filename.c_str(), &fileStat) != 0)
    {
        errorMessage = "File not found: " + filename;
        return false;
    }

    ensureOpenBabelPlugins();

    ScopedObWarningSilencer silenceWarnings;

    OpenBabel::OBMol mol;
    OpenBabel::OBConversion conv;

    OpenBabel::OBFormat* inFmt = conv.FormatFromExt(filename.c_str());

    bool inFormatSet = false;
    if (inFmt)
        inFormatSet = conv.SetInFormat(inFmt);

    if (!inFormatSet)
    {
        const std::string extNoDot = (extLower.size() > 1) ? extLower.substr(1) : std::string();
        if (!extNoDot.empty())
            inFormatSet = conv.SetInFormat(extNoDot.c_str());

        if (!inFormatSet && (extLower == ".vasp" || detectedPoscar))
            inFormatSet = conv.SetInFormat("vasp");
    }

    if (!inFormatSet)
    {
        const char* babelLibDir = std::getenv("BABEL_LIBDIR");
        errorMessage = "Unsupported file format '" + ext + "'. Supported formats: " + supportedExtensionsSummary();
        if (!babelLibDir || !*babelLibDir)
            errorMessage += " (Open Babel plugins not found; set BABEL_LIBDIR to your openbabel/<version> plugin directory).";
        return false;
    }

    // Suppress automatic bond detection inside format readers (e.g. PDB calls
    // ConnectTheDots() internally).  Bond topology is not used by this
    // application so there is no value in the O(n²) work.
    conv.AddOption("b", OpenBabel::OBConversion::INOPTIONS);

    if (!conv.ReadFile(&mol, filename))
    {
        errorMessage = "Failed to load file. The file may be corrupted or unreadable.";
        return false;
    }

    // If the file provides unit cell / periodic information, make sure we
    // generate the full unit cell (symmetry-equivalent atoms) so that we
    // display the complete crystal structure.
    if (auto *data = mol.GetData(OpenBabel::OBGenericDataType::UnitCell))
    {
        if (auto *cell = dynamic_cast<OpenBabel::OBUnitCell*>(data))
        {
            cell->FillUnitCell(&mol);
            structure.hasUnitCell = true;

            auto vecs = cell->GetCellVectors();
            if (vecs.size() >= 3)
            {
                for (int i = 0; i < 3; ++i)
                {
                    structure.cellVectors[i][0] = vecs[i].GetX();
                    structure.cellVectors[i][1] = vecs[i].GetY();
                    structure.cellVectors[i][2] = vecs[i].GetZ();
                }
            }
        }
    }

    OpenBabel::OBAtomIterator ai;
    constexpr size_t kWarningAtomCount = 10000;
    constexpr size_t kMaxAtomCount = 500000;

    for(OpenBabel::OBAtom* atom = mol.BeginAtom(ai);
        atom;
        atom = mol.NextAtom(ai))
    {
        if (structure.atoms.size() >= kMaxAtomCount)
        {
            std::cerr << "Warning: Atom count capped at " << kMaxAtomCount 
                      << " to prevent memory exhaustion." << std::endl;
            break;
        }

        AtomSite site;

        site.atomicNumber = atom->GetAtomicNum();
        site.symbol = OpenBabel::OBElements::GetSymbol(site.atomicNumber);

        site.x = atom->GetX();
        site.y = atom->GetY();
        site.z = atom->GetZ();

        getDefaultElementColor(site.atomicNumber, site.r, site.g, site.b);

        structure.atoms.push_back(site);

        if (structure.atoms.size() == kWarningAtomCount)
        {
            std::cerr << "Loading large structure with " << kWarningAtomCount 
                      << " atoms (bond computation: O(n) with spatial hashing)..." << std::endl;
        }
    }

    if (structure.atoms.size() > kWarningAtomCount)
    {
        std::cerr << "Loaded structure with " << structure.atoms.size() 
                  << " atoms." << std::endl;
    }

    if (structure.hasUnitCell)
        wrapAtomsIntoPrimaryCell(structure);

    const bool restoredFromSidecar = restoreIpfSidecar(filename, structure);
    const bool recoveredFromGeometry = !restoredFromSidecar
                                   && recoverIpfFromGeometry(structure);

    if (restoredFromSidecar)
    {
        structure.ipfLoadStatus = "IPF loaded from saved metadata.";
    }
    else if (recoveredFromGeometry)
    {
        structure.ipfLoadStatus = "IPF reconstructed from geometry fallback.";
    }
    else
    {
        structure.ipfLoadStatus = "IPF not found; using non-orientation colors.";
    }

    if (structure.atoms.empty())
    {
        errorMessage = "File loaded but no atoms were found.";
        return false;
    }

    return true;
}

Structure loadStructure(const std::string& filename)
{
    Structure structure;
    std::string errorMessage;
    if (!filename.empty() && !loadStructureFromFile(filename, structure, errorMessage))
        std::cerr << "Warning: " << errorMessage << std::endl;

    return structure;
}

bool saveStructure(const Structure& structure, const std::string& filename, const std::string& format)
{
    ensureOpenBabelPlugins();

    // Certain formats (VASP POSCAR, LAMMPS data) require all atoms of the
    // same element to be contiguous.  If atoms are interleaved – as they are
    // after substitutional-solid-solution or other random builders – the
    // writer produces duplicate or malformed element headers.  Stable-sort
    // a working copy by atomic number before building the OBMol.
    const bool needsGrouping = (format == "vasp" || format == "VASP"
                             || format == "lmpdat" || format == "lammps");

    const Structure* src = &structure;
    Structure grouped;
    if (needsGrouping)
    {
        const std::size_t n = structure.atoms.size();
        std::vector<std::size_t> order(n);
        std::iota(order.begin(), order.end(), std::size_t(0));
        std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b)
        {
            return structure.atoms[a].atomicNumber < structure.atoms[b].atomicNumber;
        });

        grouped = structure;
        for (std::size_t i = 0; i < n; ++i)
            grouped.atoms[i] = structure.atoms[order[i]];

        if (structure.grainColors.size() == n)
        {
            grouped.grainColors.resize(n);
            for (std::size_t i = 0; i < n; ++i)
                grouped.grainColors[i] = structure.grainColors[order[i]];
        }
        if (structure.grainRegionIds.size() == n)
        {
            grouped.grainRegionIds.resize(n);
            for (std::size_t i = 0; i < n; ++i)
                grouped.grainRegionIds[i] = structure.grainRegionIds[order[i]];
        }
        src = &grouped;
    }

    OpenBabel::OBMol mol;
    mol.BeginModify();

    for (const auto& site : src->atoms)
    {
        OpenBabel::OBAtom* atom = mol.NewAtom();
        atom->SetAtomicNum(site.atomicNumber);
        atom->SetVector(site.x, site.y, site.z);
    }

    if (src->hasUnitCell)
    {
        OpenBabel::OBUnitCell* cell = new OpenBabel::OBUnitCell();
        OpenBabel::vector3 va(src->cellVectors[0][0], src->cellVectors[0][1], src->cellVectors[0][2]);
        OpenBabel::vector3 vb(src->cellVectors[1][0], src->cellVectors[1][1], src->cellVectors[1][2]);
        OpenBabel::vector3 vc(src->cellVectors[2][0], src->cellVectors[2][1], src->cellVectors[2][2]);
        cell->SetData(va, vb, vc);
        mol.SetData(cell);
    }

    // ConnectTheDots / PerceiveBondOrders are O(n²) and cause OOM crashes
    // for structures with more than a few hundred atoms, or for any periodic
    // structure where every atom is within bonding range of many neighbours.
    // For formats that actually store bonds (mol2, sdf) the data is
    // irrelevant for export; for all other formats (xyz, cif, vasp …) bonds
    // are not part of the file format at all.
    constexpr size_t kBondDetectionAtomLimit = 500;
    if (!src->hasUnitCell && src->atoms.size() <= kBondDetectionAtomLimit)
    {
        mol.ConnectTheDots();
        mol.PerceiveBondOrders();
    }
    mol.EndModify();

    OpenBabel::OBConversion conv;
    OpenBabel::OBFormat* outFmt = conv.FindFormat(format.c_str());
    if (!outFmt)
        return false;

    conv.SetOutFormat(outFmt);
    if (!conv.WriteFile(&mol, filename))
        return false;

    return saveIpfSidecar(structure, filename);
}
