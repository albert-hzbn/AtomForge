#include "StructureLoader.h"

#include <openbabel3/openbabel/obconversion.h>
#include <openbabel3/openbabel/mol.h>
#include <openbabel3/openbabel/atom.h>
#include <openbabel3/openbabel/elements.h>
#include <openbabel3/openbabel/generic.h>
#include <openbabel3/openbabel/math/vector3.h>
#include <openbabel3/openbabel/math/spacegroup.h>
#include <openbabel3/openbabel/plugin.h>
#include <openbabel3/openbabel/oberror.h>

#include <iostream>
#include <array>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
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
    std::array<double,3> position{};
    int region=-1;
    bool hasRegion=false;
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
    return (long long)std::floor(value * 500.0);
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

    const bool regions=structure.grainRegionIds.size()==structure.atoms.size();
    out << (regions ? "ATOMFORGE_IPF_V2\n" : "ATOMFORGE_IPF_V1\n");
    out << structure.atoms.size() << "\n";
    out << std::setprecision(17);
    for (size_t i = 0; i < structure.atoms.size(); ++i)
    {
        const AtomSite& atom = structure.atoms[i];
        const std::array<float, 3>& color = structure.grainColors[i];
        out << atom.atomicNumber << ' '
            << atom.x << ' ' << atom.y << ' ' << atom.z << ' '
            << color[0] << ' ' << color[1] << ' ' << color[2];
        if (regions) out << ' ' << structure.grainRegionIds[i];
        out << '\n';
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
    if (!in || count>10000000)
        return false;

    records.reserve(count);
    for (size_t i = 0; i < count; ++i)
    {
        int atomicNumber = 0;
        double x = 0.0, y = 0.0, z = 0.0;
        float r = 0.0f, g = 0.0f, b = 0.0f;
        int region=-1;
        if (isV1)
        {
            in >> atomicNumber >> x >> y >> z >> r >> g >> b;
        }
        else
        {
            in >> atomicNumber >> x >> y >> z >> r >> g >> b >> region;
        }
        if (!in)
            return false;

        if (atomicNumber<1 || atomicNumber>118 || !std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)
            || std::abs(x)>1e9 || std::abs(y)>1e9 || std::abs(z)>1e9
            || !std::isfinite(r) || !std::isfinite(g) || !std::isfinite(b)) return false;
        IpfRecord record;
        record.position={x,y,z};
        record.key = makeIpfKey(atomicNumber, x, y, z);
        record.color = {{r, g, b}};
        record.region=region; record.hasRegion=isV2;
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
    std::vector<int> regions(structure.atoms.size(),-1);
    std::vector<bool> assigned(structure.atoms.size(), false);
    size_t matched = 0;

    for (size_t i = 0; i < records.size(); ++i)
    {
        size_t chosen = structure.atoms.size();
        for (int dx=-1;dx<=1 && chosen==structure.atoms.size();++dx)
        for (int dy=-1;dy<=1 && chosen==structure.atoms.size();++dy)
        for (int dz=-1;dz<=1 && chosen==structure.atoms.size();++dz) {
            auto key=records[i].key; key.qx+=dx; key.qy+=dy; key.qz+=dz;
            const auto range=byKey.equal_range(key);
            for (auto it=range.first;it!=range.second;++it) {
                const auto& atom=structure.atoms[it->second];
                const auto& p=records[i].position;
                if (!assigned[it->second] && std::abs(atom.x-p[0])<.002 &&
                    std::abs(atom.y-p[1])<.002 && std::abs(atom.z-p[2])<.002) {
                    chosen=it->second; break;
                }
            }
        }
        if (chosen == structure.atoms.size())
            break;

        restored[chosen] = records[i].color;
        regions[chosen] = records[i].region;
        assigned[chosen] = true;
        ++matched;
    }

    if (matched != structure.atoms.size()) return false;

    structure.grainColors.swap(restored);
    if (!records.empty() && records.front().hasRegion) structure.grainRegionIds.swap(regions);
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

// Open Babel reads its data files (space-groups.txt among them) from
// BABEL_DATADIR. Without them a CIF that names its space group but lists no
// symmetry operations is not expanded, so point it at the data next to the
// plugins in use, or at a bundled copy.
void ensureOpenBabelData()
{
    const auto hasSpaceGroups = [](const std::string& dir) {
        std::ifstream probe(joinPath(dir, "space-groups.txt"));
        return probe.good();
    };
    const char* current = std::getenv("BABEL_DATADIR");
    if (current && *current && hasSpaceGroups(normalizeSeparators(current)))
        return;

    std::vector<std::string> candidates;
    const std::string exeDir = executableDirectory();
    if (!exeDir.empty())
    {
        appendVersionedPluginCandidates(joinPath(exeDir, "openbabel-data"), candidates);
        const std::string prefixDir = parentDirectory(exeDir);
        if (!prefixDir.empty())
            appendVersionedPluginCandidates(joinPath(joinPath(prefixDir, "share"), "openbabel"), candidates);
    }
    // <prefix>/lib/openbabel[/<version>] -> <prefix>/share/openbabel/<version>
    if (const char* libDir = std::getenv("BABEL_LIBDIR"); libDir && *libDir)
    {
        const std::string plugins = normalizeSeparators(libDir);
        const std::size_t lib = plugins.rfind("/lib/");
        if (lib != std::string::npos)
            appendVersionedPluginCandidates(joinPath(joinPath(plugins.substr(0, lib), "share"), "openbabel"), candidates);
    }
#ifdef _WIN32
    appendVersionedPluginCandidates("C:/msys64/ucrt64/share/openbabel", candidates);
    appendVersionedPluginCandidates("C:/msys64/mingw64/share/openbabel", candidates);
#else
    appendVersionedPluginCandidates("/usr/share/openbabel", candidates);
    appendVersionedPluginCandidates("/usr/local/share/openbabel", candidates);
#endif
    for (const std::string& candidate : candidates)
    {
        if (!hasSpaceGroups(candidate))
            continue;
        std::string dir = candidate;
        while (dir.size() > 1 && dir.back() == '/') dir.pop_back();
#ifdef _WIN32
        _putenv_s("BABEL_DATADIR", dir.c_str());
#else
        setenv("BABEL_DATADIR", dir.c_str(), 1);
#endif
        return;
    }
}

void ensureOpenBabelPlugins()
{
    static bool initialized = false;
    if (initialized)
        return;
    initialized = true;

    // Suppress OB diagnostic output during the entire plugin-discovery phase.
    // "Unable to find OpenBabel plugins" is emitted on every failed LoadAllPlugins()
    // attempt, which is noisy and unhelpful to the user.
    const OpenBabel::obMessageLevel savedLevel =
        OpenBabel::obErrorLog.GetOutputLevel();
    // Suppress all OB diagnostics during plugin discovery; obMessageLevel has
    // no "silent" constant so we cast -1 which is below obError (0).
    OpenBabel::obErrorLog.SetOutputLevel(
        static_cast<OpenBabel::obMessageLevel>(-1));

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
            {
                ensureOpenBabelData();
                OpenBabel::obErrorLog.SetOutputLevel(savedLevel);
                return;
            }
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

    ensureOpenBabelData();
    OpenBabel::obErrorLog.SetOutputLevel(savedLevel);
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
        ".gjf",
        ".extxyz",
        ".data",
        ".lmp"
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
    return ".cif, .xyz, .extxyz, .vasp (or POSCAR/CONTCAR), .data/.lmp (LAMMPS data), .pdb, .mol, .mol2, .sdf, .pwi, .gjf";
}

bool filenameIsPoscarOrContcar(const std::string& filename)
{
    // Extract basename (after last slash or backslash)
    const std::size_t slashPos = filename.find_last_of("/\\");
    const std::string basename = (slashPos != std::string::npos)
                               ? filename.substr(slashPos + 1)
                               : filename;
    // POSCAR, CONTCAR_relaxed, Cu_POSCAR, NaCl.poscar, POSCAR-001 (phonopy) ...
    std::string upper = basename;
    std::transform(upper.begin(), upper.end(), upper.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return upper.find("POSCAR") != std::string::npos
        || upper.find("CONTCAR") != std::string::npos;
}

void wrapAtomsIntoPrimaryCell(Structure& structure)
{
    if (!structure.hasUnitCell || structure.atoms.empty())
        return;

    // In double precision, relative to the cell's origin; an atom already
    // inside the cell keeps its coordinates exactly.
    const glm::dvec3 origin(structure.cellOffset[0], structure.cellOffset[1], structure.cellOffset[2]);
    const glm::dvec3 a(structure.cellVectors[0][0], structure.cellVectors[0][1], structure.cellVectors[0][2]);
    const glm::dvec3 b(structure.cellVectors[1][0], structure.cellVectors[1][1], structure.cellVectors[1][2]);
    const glm::dvec3 c(structure.cellVectors[2][0], structure.cellVectors[2][1], structure.cellVectors[2][2]);

    const glm::dmat3 cellMat(a, b, c);
    if (std::abs(glm::determinant(cellMat)) <= 1e-8)
        return;

    const glm::dmat3 invCellMat = glm::inverse(cellMat);
    constexpr double kWrapTol = 1e-5;

    for (auto& atom : structure.atoms)
    {
        const glm::dvec3 pos(atom.x, atom.y, atom.z);
        glm::dvec3 frac = invCellMat * (pos - origin);
        bool inside = true;
        for (int k = 0; k < 3; ++k)
            inside = inside && frac[k] >= -1e-12 && frac[k] < 1.0 - kWrapTol;
        if (inside)
            continue;

        for (int k = 0; k < 3; ++k)
        {
            frac[k] -= std::floor(frac[k]);
            // Just below 1 (or a rounding error below 0) is the 0 face.
            if (std::abs(frac[k]) <= kWrapTol || std::abs(1.0 - frac[k]) <= kWrapTol) frac[k] = 0.0;
        }

        const glm::dvec3 wrapped = origin + frac.x * a + frac.y * b + frac.z * c;
        atom.x = wrapped.x;
        atom.y = wrapped.y;
        atom.z = wrapped.z;
    }
}
}




// XYZ element symbols in standard case ("NI", "cu" -> "Ni", "Cu"): Open Babel
// only capitalises the first letter, so an all-capitals symbol would load as
// an unknown element. Returns false (text unchanged) if the frames do not
// follow the XYZ layout.
// The cell of an extended XYZ file: Lattice="ax ay az bx by bz cx cy cz" on
// the first frame's comment line (Open Babel does not read it).
static bool readExtendedXyzLattice(const std::string& text, std::array<std::array<double, 3>, 3>& cell)
{
    const std::size_t first = text.find('\n');
    if (first == std::string::npos) return false;
    const std::size_t second = text.find('\n', first + 1);
    const std::string comment = text.substr(first + 1, second == std::string::npos ? std::string::npos : second - first - 1);
    std::size_t key = std::string::npos;
    for (std::size_t at = comment.find("attice="); at != std::string::npos; at = comment.find("attice=", at + 1))
        if (at > 0 && (comment[at - 1] == 'L' || comment[at - 1] == 'l')) { key = at + 7; break; }
    if (key == std::string::npos || key >= comment.size()) return false;
    const char quote = comment[key];
    if (quote != '"' && quote != '\'') return false;
    const std::size_t close = comment.find(quote, key + 1);
    if (close == std::string::npos) return false;
    std::istringstream values(comment.substr(key + 1, close - key - 1));
    for (auto& row : cell)
        for (double& value : row)
            if (!(values >> value) || !std::isfinite(value)) return false;
    return true;
}

// POSCAR details Open Babel misreads: a negative scale factor is the cell
// volume, and a VASP 4 file has no element line (the species are named on
// the title line, as other tools write them). Returns false when the text
// is left unchanged.
static bool normalizePoscar(std::string& text)
{
    std::vector<std::string> lines;
    for (std::size_t pos = 0; pos < text.size();)
    {
        const std::size_t end = text.find('\n', pos);
        lines.push_back(text.substr(pos, end == std::string::npos ? std::string::npos : end - pos + 1));
        pos = end == std::string::npos ? text.size() : end + 1;
    }
    if (lines.size() < 7) return false;
    const auto tokens = [](const std::string& line) {
        std::istringstream in(line);
        std::vector<std::string> out;
        for (std::string token; in >> token;) out.push_back(token);
        return out;
    };
    bool changed = false;

    // Scale line: one negative value is the volume of the cell.
    const auto scale = tokens(lines[1]);
    if (scale.size() == 1)
    {
        char* end = nullptr;
        const double value = std::strtod(scale[0].c_str(), &end);
        if (end != scale[0].c_str() && value < 0.0)
        {
            double m[3][3];
            for (int r = 0; r < 3; ++r)
            {
                const auto row = tokens(lines[static_cast<std::size_t>(2 + r)]);
                if (row.size() < 3) return false;
                for (int c = 0; c < 3; ++c) m[r][c] = std::atof(row[static_cast<std::size_t>(c)].c_str());
            }
            const double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                             - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                             + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
            if (std::abs(det) < 1e-12) return false;
            char buffer[64];
            std::snprintf(buffer, sizeof(buffer), "%.12g\n", std::cbrt(-value / std::abs(det)));
            lines[1] = buffer;
            changed = true;
        }
    }

    // VASP 4: the line after the cell holds the counts; take the species from
    // the title ("Na Cl", "Fe_pv O") when it names one per count.
    const auto counts = tokens(lines[5]);
    bool allCounts = !counts.empty();
    for (const auto& count : counts)
        allCounts = allCounts && count.find_first_not_of("0123456789") == std::string::npos;
    if (!allCounts && !counts.empty())
    {
        // Element line with POTCAR labels ("Fe_pv", "O_s", "Ti/abc"): the element.
        std::string line;
        bool relabelled = false;
        for (const auto& word : counts)
        {
            std::size_t letters = 0;
            while (letters < word.size() && std::isalpha(static_cast<unsigned char>(word[letters]))) ++letters;
            const std::string symbol = letters >= 1 && letters <= 2 && letters < word.size() ? word.substr(0, letters) : word;
            relabelled = relabelled || symbol != word;
            line += (line.empty() ? "" : " ") + symbol;
        }
        if (relabelled)
        {
            lines[5] = line + "\n";
            changed = true;
        }
    }
    if (allCounts)
    {
        std::vector<std::string> species;
        for (const auto& word : tokens(lines[0]))
        {
            std::size_t letters = 0;
            while (letters < word.size() && std::isalpha(static_cast<unsigned char>(word[letters]))) ++letters;
            if (letters == 0 || letters > 2) { species.clear(); break; }
            std::string symbol = word.substr(0, letters);
            symbol[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(symbol[0])));
            if (symbol.size() == 2) symbol[1] = static_cast<char>(std::tolower(static_cast<unsigned char>(symbol[1])));
            species.push_back(symbol);
        }
        if (species.size() == counts.size())
        {
            std::string line;
            for (const auto& symbol : species) line += (line.empty() ? "" : " ") + symbol;
            lines.insert(lines.begin() + 5, line + "\n");
            changed = true;
        }
    }
    if (!changed) return false;
    text.clear();
    for (const auto& line : lines) text += line;
    return true;
}

// The columns of the species and the positions in an extended XYZ frame,
// from Properties=name:type:count:... on its comment line (for example
// id:I:1:species:S:1:pos:R:3). False without a usable Properties key.
static bool extendedXyzColumns(const std::string& comment, int& speciesColumn, int& positionColumn)
{
    std::size_t key = std::string::npos;
    for (std::size_t at = comment.find("roperties="); at != std::string::npos; at = comment.find("roperties=", at + 1))
        if (at > 0 && (comment[at - 1] == 'P' || comment[at - 1] == 'p')) { key = at + 10; break; }
    if (key == std::string::npos) return false;
    const std::size_t end = comment.find_first_of(" \t\r\n", key);
    std::istringstream spec(comment.substr(key, end == std::string::npos ? std::string::npos : end - key));
    speciesColumn = positionColumn = -1;
    int column = 0;
    for (std::string name, type, count; std::getline(spec, name, ':') && std::getline(spec, type, ':') && std::getline(spec, count, ':');)
    {
        int n = 0;
        try { n = std::stoi(count); } catch (const std::exception&) { return false; }
        if (n < 1) return false;
        std::string lower = name;
        for (char& c : lower) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (lower == "species" && n == 1) speciesColumn = column;
        if ((lower == "pos" || lower == "positions") && n == 3) positionColumn = column;
        column += n;
    }
    return speciesColumn >= 0 && positionColumn >= 0;
}

static bool normalizeXyzSymbols(std::string& text)
{
    std::string out;
    out.reserve(text.size());
    std::size_t pos = 0;
    const auto nextLine = [&](std::string& line) {
        if (pos >= text.size()) return false;
        const std::size_t end = text.find('\n', pos);
        line = text.substr(pos, end == std::string::npos ? std::string::npos : end - pos + 1);
        pos = end == std::string::npos ? text.size() : end + 1;
        return true;
    };
    std::string line;
    while (nextLine(line))
    {
        char* endPtr = nullptr;
        const long count = std::strtol(line.c_str(), &endPtr, 10);
        if (endPtr == line.c_str() || count < 0)
        {
            // Trailing blank lines are fine; anything else is not XYZ.
            if (line.find_first_not_of(" \t\r\n") != std::string::npos) return false;
            out += line;
            continue;
        }
        out += line;
        if (!nextLine(line)) return false;  // the comment line
        out += line;
        // Open Babel reads "symbol x y z": reorder the columns of an extended
        // XYZ frame that lists others first (an id, say).
        int speciesColumn = 0, positionColumn = 1;
        const bool reorder = extendedXyzColumns(line, speciesColumn, positionColumn) &&
                             (speciesColumn != 0 || positionColumn != 1);
        for (long n = 0; n < count; ++n)
        {
            if (!nextLine(line)) return false;
            if (reorder)
            {
                std::istringstream fields(line);
                std::vector<std::string> values;
                for (std::string value; fields >> value;) values.push_back(value);
                if ((int)values.size() <= std::max(speciesColumn, positionColumn + 2)) return false;
                line = values[(size_t)speciesColumn] + ' ' + values[(size_t)positionColumn] + ' ' +
                       values[(size_t)positionColumn + 1] + ' ' + values[(size_t)positionColumn + 2] + '\n';
            }
            const std::size_t start = line.find_first_not_of(" \t");
            std::size_t stop = start;
            while (stop < line.size() && std::isalpha(static_cast<unsigned char>(line[stop]))) ++stop;
            if (start != std::string::npos && stop > start)
            {
                line[start] = static_cast<char>(std::toupper(static_cast<unsigned char>(line[start])));
                for (std::size_t k = start + 1; k < stop; ++k)
                    line[k] = static_cast<char>(std::tolower(static_cast<unsigned char>(line[k])));
            }
            out += line;
        }
    }
    text.swap(out);
    return true;
}

bool isSupportedStructureFile(const std::string& filename)
{
    const std::string extLower = toLowerCopy(extractExtension(filename));
    if (isSupportedExtension(extLower))
        return true;
    return filenameIsPoscarOrContcar(filename);
}

bool isLargeDataFile(const std::string& filename)
{
    const std::string extLower = toLowerCopy(extractExtension(filename));
    return extLower == ".afcloud" || extLower == ".dump" || extLower == ".lammpstrj";
}

static bool finishLoadedStructure(const std::string& filename, Structure& structure, std::string& errorMessage);

// A disordered CIF lists several species on one site with occupancies below
// 1 (Cu 0.5 and Ni 0.5 at the origin). An editable structure holds one atom
// per site: keep the species with the largest occupancy (the first listed on
// a tie) and say so in the load note. Sites with one species and an
// occupancy below 1 (partial vacancies) are kept and counted.
static void resolvePartialOccupancy(Structure& structure, const std::vector<double>& occupancies)
{
    if (occupancies.size() != structure.atoms.size()) return;
    if (std::none_of(occupancies.begin(), occupancies.end(), [](double o) { return o < 1.0 - 1e-6; })) return;
    if (structure.hasUnitCell) wrapAtomsIntoPrimaryCell(structure);

    // Group atoms by position, to 0.01 A; atoms at fractional 0 and 1 are
    // the same site after wrapping.
    std::map<std::array<long long, 3>, std::vector<std::size_t>> sites;
    for (std::size_t i = 0; i < structure.atoms.size(); ++i)
    {
        const auto& a = structure.atoms[i];
        sites[{std::llround(a.x * 100.0), std::llround(a.y * 100.0), std::llround(a.z * 100.0)}].push_back(i);
    }
    std::vector<char> keep(structure.atoms.size(), 1);
    std::size_t shared = 0, partial = 0;
    for (const auto& [key, members] : sites)
    {
        if (members.size() > 1)
        {
            ++shared;
            std::size_t best = members[0];
            for (std::size_t m : members)
                if (occupancies[m] > occupancies[best] + 1e-9) best = m;
            for (std::size_t m : members) keep[m] = m == best;
        }
        else if (occupancies[members[0]] < 1.0 - 1e-6)
            ++partial;
    }
    if (shared)
    {
        std::vector<AtomSite> kept;
        kept.reserve(structure.atoms.size());
        for (std::size_t i = 0; i < structure.atoms.size(); ++i)
            if (keep[i]) kept.push_back(structure.atoms[i]);
        structure.atoms.swap(kept);
    }
    std::string note;
    if (shared)
        note = std::to_string(shared) + (shared == 1 ? " site holds" : " sites hold") +
               " more than one species (partial occupancy); each keeps its majority species "
               "(Build > Substitutional Solid Solution makes a random mixture).";
    if (partial)
        note += std::string(note.empty() ? "" : " ") + std::to_string(partial) +
                (partial == 1 ? " site is" : " sites are") + " partially occupied; its atom is kept.";
    structure.loadNote = note;
}

// The element of a LAMMPS atom type: a label after its mass ("1 63.546 # Cu",
// as AtomForge, ASE and others write), a label from an "Atom Type Labels"
// section, or else the element whose standard mass is nearest (within 0.5 u).
static int lammpsTypeElement(const std::string& label, double mass)
{
    if (!label.empty())
    {
        std::string symbol = label;
        symbol[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(symbol[0])));
        for (std::size_t k = 1; k < symbol.size(); ++k) symbol[k] = static_cast<char>(std::tolower(static_cast<unsigned char>(symbol[k])));
        const int z = static_cast<int>(OpenBabel::OBElements::GetAtomicNum(symbol.c_str()));
        if (z > 0) return z;
    }
    int best = 0;
    double bestDifference = 0.5;
    for (int z = 1; z <= 118 && mass > 0.0; ++z)
        if (const double difference = std::abs(elementAtomicMass(z) - mass); difference < bestDifference)
        {
            best = z;
            bestDifference = difference;
        }
    return best;
}

// A LAMMPS data file (read_data): the box, the Masses (or type labels) and
// the Atoms section in the atomic, charge, full or molecular style (from the
// "Atoms # style" comment, else from the number of columns).
static bool readLammpsData(const std::string& text, Structure& structure, std::string& errorMessage)
{
    std::istringstream in(text);
    std::string line;
    std::getline(in, line);  // title
    long long atomCount = -1;
    int typeCount = 0;
    double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
    double xy = 0, xz = 0, yz = 0;
    bool box[3] = {false, false, false};
    std::map<int, double> masses;
    std::map<int, std::string> labels;
    std::vector<std::vector<std::string>> atomRows;
    std::string style;
    std::string section;
    const auto stripComment = [](const std::string& value, std::string* comment) {
        const std::size_t hash = value.find('#');
        if (comment) *comment = hash == std::string::npos ? std::string() : value.substr(hash + 1);
        return hash == std::string::npos ? value : value.substr(0, hash);
    };
    const auto words = [](const std::string& value) {
        std::istringstream fields(value);
        std::vector<std::string> out;
        for (std::string word; fields >> word;) out.push_back(word);
        return out;
    };
    const auto trim = [](std::string value) {
        const std::size_t a = value.find_first_not_of(" \t\r\n");
        const std::size_t b = value.find_last_not_of(" \t\r\n");
        return a == std::string::npos ? std::string() : value.substr(a, b - a + 1);
    };
    while (std::getline(in, line))
    {
        std::string comment;
        const std::string body = trim(stripComment(line, &comment));
        if (body.empty()) continue;
        const auto w = words(body);
        if (!std::isdigit(static_cast<unsigned char>(w[0][0])) && w[0][0] != '-' && w[0][0] != '+' && w[0][0] != '.')
        {
            // A section header (Masses, Atoms, Velocities, Bonds, ...).
            section = body;
            if (section == "Atoms") style = trim(comment);
            continue;
        }
        try
        {
            if (section.empty())
            {
                if (w.size() >= 2 && w[1] == "atoms") atomCount = std::stoll(w[0]);
                else if (w.size() >= 3 && w[1] == "atom" && w[2] == "types") typeCount = std::stoi(w[0]);
                else if (w.size() >= 4 && w[2] == "xlo" && w[3] == "xhi") { lo[0] = std::stod(w[0]); hi[0] = std::stod(w[1]); box[0] = true; }
                else if (w.size() >= 4 && w[2] == "ylo" && w[3] == "yhi") { lo[1] = std::stod(w[0]); hi[1] = std::stod(w[1]); box[1] = true; }
                else if (w.size() >= 4 && w[2] == "zlo" && w[3] == "zhi") { lo[2] = std::stod(w[0]); hi[2] = std::stod(w[1]); box[2] = true; }
                else if (w.size() >= 6 && w[3] == "xy" && w[4] == "xz" && w[5] == "yz") { xy = std::stod(w[0]); xz = std::stod(w[1]); yz = std::stod(w[2]); }
            }
            else if (section == "Masses" && w.size() >= 2)
            {
                masses[std::stoi(w[0])] = std::stod(w[1]);
                const auto label = words(comment);
                if (!label.empty()) labels.emplace(std::stoi(w[0]), label[0]);
            }
            else if (section == "Atom Type Labels" && w.size() >= 2)
                labels[std::stoi(w[0])] = w[1];
            else if (section == "Atoms")
                atomRows.push_back(w);
        }
        catch (const std::exception&)
        {
            errorMessage = "Cannot read this line of the LAMMPS data file: " + trim(line);
            return false;
        }
    }
    if (atomCount < 0 || !box[0] || !box[1] || !box[2])
    {
        errorMessage = "Not a LAMMPS data file: it needs the atom count and the xlo xhi, ylo yhi and zlo zhi box lines.";
        return false;
    }
    if (static_cast<long long>(atomRows.size()) != atomCount)
    {
        errorMessage = "The LAMMPS data file declares " + std::to_string(atomCount) + " atoms but its Atoms section lists " +
                       std::to_string(atomRows.size()) + ".";
        return false;
    }
    // Columns of the type and x in each style; image flags may follow z.
    int typeColumn = 1, xColumn = 2;
    const std::string kind = words(style).empty() ? std::string() : words(style)[0];
    const std::size_t columns = atomRows.empty() ? 0 : atomRows[0].size();
    if (kind == "atomic") { typeColumn = 1; xColumn = 2; }
    else if (kind == "charge") { typeColumn = 1; xColumn = 3; }
    else if (kind == "full") { typeColumn = 2; xColumn = 4; }
    else if (kind == "molecular" || kind == "bond" || kind == "angle") { typeColumn = 2; xColumn = 3; }
    else if (!kind.empty())
    {
        errorMessage = "LAMMPS atom style '" + kind + "' is not supported; use atomic, charge, molecular or full.";
        return false;
    }
    else if (columns == 5 || columns == 8) { typeColumn = 1; xColumn = 2; }
    else if (columns == 6 || columns == 9) { typeColumn = 1; xColumn = 3; }
    else if (columns == 7 || columns == 10) { typeColumn = 2; xColumn = 4; }
    else
    {
        errorMessage = "Cannot tell the LAMMPS atom style; add it after the section name, e.g. \"Atoms # atomic\".";
        return false;
    }

    std::map<int, int> elementOfType;
    std::string unknownTypes;
    for (int type = 1; type <= std::max(typeCount, 1); ++type)
    {
        const auto mass = masses.find(type);
        const auto label = labels.find(type);
        const int z = lammpsTypeElement(label == labels.end() ? std::string() : label->second, mass == masses.end() ? 0.0 : mass->second);
        if (z > 0) elementOfType[type] = z;
        else unknownTypes += (unknownTypes.empty() ? "" : ", ") + std::to_string(type);
    }
    if (!unknownTypes.empty())
    {
        errorMessage = "The element of LAMMPS atom type " + unknownTypes +
                       " is unknown: give its mass in a Masses section, or its symbol after the mass (\"1 63.546 # Cu\").";
        return false;
    }

    // Sorted by atom id, as LAMMPS numbers them.
    std::vector<std::pair<long long, AtomSite>> sites;
    sites.reserve(atomRows.size());
    for (const auto& row : atomRows)
    {
        if ((int)row.size() < xColumn + 3)
        {
            errorMessage = "A LAMMPS atom line has too few columns for the " + (kind.empty() ? std::string("detected") : kind) + " style.";
            return false;
        }
        try
        {
            const int type = std::stoi(row[(size_t)typeColumn]);
            const auto element = elementOfType.find(type);
            if (element == elementOfType.end())
            {
                errorMessage = "LAMMPS atom " + row[0] + " has type " + row[(size_t)typeColumn] + ", which is not among the declared types.";
                return false;
            }
            AtomSite site;
            site.atomicNumber = element->second;
            site.symbol = elementSymbol(site.atomicNumber);
            site.x = std::stod(row[(size_t)xColumn]);
            site.y = std::stod(row[(size_t)xColumn + 1]);
            site.z = std::stod(row[(size_t)xColumn + 2]);
            getDefaultElementColor(site.atomicNumber, site.r, site.g, site.b);
            sites.push_back({std::stoll(row[0]), site});
        }
        catch (const std::exception&)
        {
            errorMessage = "Cannot read LAMMPS atom line starting " + row[0] + ".";
            return false;
        }
    }
    std::stable_sort(sites.begin(), sites.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    structure.atoms.reserve(sites.size());
    for (auto& entry : sites) structure.atoms.push_back(std::move(entry.second));

    structure.hasUnitCell = true;
    structure.cellVectors = {{{hi[0] - lo[0], 0.0, 0.0}, {xy, hi[1] - lo[1], 0.0}, {xz, yz, hi[2] - lo[2]}}};
    structure.cellOffset = {lo[0], lo[1], lo[2]};
    return true;
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

    if (extLower == ".data" || extLower == ".lmp")
    {
        std::ifstream in(filename, std::ios::binary);
        const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        if (!readLammpsData(text, structure, errorMessage)) return false;
        return finishLoadedStructure(filename, structure, errorMessage);
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

        if (!inFormatSet && extLower == ".extxyz")
            inFormatSet = conv.SetInFormat("xyz");
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

    bool read = false;
    bool selectiveDynamics = false;
    bool xyzLattice = false;
    std::array<std::array<double, 3>, 3> xyzCell{};
    if (extLower == ".xyz" || extLower == ".extxyz")
    {
        std::ifstream in(filename, std::ios::binary);
        std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        xyzLattice = readExtendedXyzLattice(text, xyzCell);
        read = normalizeXyzSymbols(text) ? conv.ReadString(&mol, text) : conv.ReadFile(&mol, filename);
    }
    else if (extLower == ".vasp" || extLower == ".poscar" || detectedPoscar)
    {
        std::ifstream in(filename, std::ios::binary);
        std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        // Selective dynamics: the line after the counts starts with S. The
        // flags are not kept, which a relaxation of the saved file would miss.
        {
            std::istringstream lines(text);
            std::string line;
            for (int n = 0; n < 9 && std::getline(lines, line); ++n)
            {
                const std::size_t first = line.find_first_not_of(" \t");
                if (n >= 6 && first != std::string::npos && (line[first] == 'S' || line[first] == 's'))
                {
                    selectiveDynamics = true;
                    break;
                }
            }
        }
        if (normalizePoscar(text))
        {
            // Open Babel's VASP reader needs a file, not a string.
            std::error_code error;
            const auto temporary = std::filesystem::temp_directory_path(error) /
                                   ("atomforge-" + std::to_string(std::hash<std::string>{}(filename)) + "-POSCAR");
            {
                std::ofstream out(temporary, std::ios::binary);
                out << text;
            }
            read = conv.ReadFile(&mol, temporary.string());
            std::filesystem::remove(temporary, error);
        }
        else
            read = conv.ReadFile(&mol, filename);
    }
    else
        read = conv.ReadFile(&mol, filename);
    if (!read)
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
            // FillUnitCell also deletes atoms it judges duplicates, which
            // corrupts P1 data such as dislocation cores (Open Babel versions
            // differ in tolerance). Only expand a declared nontrivial group.
            const OpenBabel::SpaceGroup* group = cell->GetSpaceGroup();
            if (group && group->GetId() > 1)
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

    if (!structure.hasUnitCell && xyzLattice)
    {
        structure.hasUnitCell = true;
        structure.cellVectors = xyzCell;
    }

    OpenBabel::OBAtomIterator ai;
    // Larger XYZ files open as atom clouds (app/LargeData.h, same limit);
    // never keep only part of a structure, which a save would then write.
    constexpr size_t kMaxAtomCount = 5000000;
    if (mol.NumAtoms() > kMaxAtomCount)
    {
        errorMessage = "The file has " + std::to_string(mol.NumAtoms()) + " atoms, more than the " +
                       std::to_string(kMaxAtomCount) + " that can be edited; save it as extended XYZ to view it as an atom cloud.";
        return false;
    }
    structure.atoms.reserve(mol.NumAtoms());
    std::vector<double> occupancies;
    occupancies.reserve(mol.NumAtoms());

    for(OpenBabel::OBAtom* atom = mol.BeginAtom(ai);
        atom;
        atom = mol.NextAtom(ai))
    {

        AtomSite site;

        site.atomicNumber = atom->GetAtomicNum();
        site.symbol = OpenBabel::OBElements::GetSymbol(site.atomicNumber);

        site.x = atom->GetX();
        site.y = atom->GetY();
        site.z = atom->GetZ();

        getDefaultElementColor(site.atomicNumber, site.r, site.g, site.b);

        double occupancy = 1.0;
        if (auto* data = dynamic_cast<OpenBabel::OBPairFloatingPoint*>(atom->GetData("_atom_site_occupancy")))
            occupancy = data->GetGenericValue();
        occupancies.push_back(occupancy);

        structure.atoms.push_back(site);
    }

    resolvePartialOccupancy(structure, occupancies);
    if (selectiveDynamics)
        structure.loadNote += std::string(structure.loadNote.empty() ? "" : " ") +
                              "The file's selective-dynamics flags (fixed atoms) are not kept; saving writes every atom free to move.";
    return finishLoadedStructure(filename, structure, errorMessage);
}

// Wrapping, orientation colours and the empty check, for every reader.
static bool finishLoadedStructure(const std::string& filename, Structure& structure, std::string& errorMessage)
{
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

namespace
{
thread_local std::string g_lastSaveError;
}

const std::string& lastSaveError()
{
    return g_lastSaveError;
}

bool saveStructure(const Structure& structure, const std::string& filename, const std::string& format)
{
    g_lastSaveError.clear();
    ensureOpenBabelPlugins();

    // CIF and POSCAR hold positions relative to a cell: without one Open
    // Babel writes Cartesian values as fractional (CIF) or a zero lattice.
    std::string formatLower = format;
    std::transform(formatLower.begin(), formatLower.end(), formatLower.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (!structure.hasUnitCell && (formatLower == "cif" || formatLower == "vasp" || formatLower == "poscar"))
    {
        g_lastSaveError = (formatLower == "cif" ? std::string("CIF") : std::string("VASP")) +
            " files need a unit cell, and this structure has none. Save it as XYZ, PDB or MOL, "
            "or give it a cell with Edit > Structure Operations > Add Vacuum.";
        return false;
    }

    if (formatLower == "extxyz") formatLower = "xyz";

    // CIF is written here, in P1 with the cell as lengths and angles and
    // fractional coordinates to twelve decimals (Open Babel keeps five, which
    // moves atoms by up to 5e-6 of the cell edge).
    if (formatLower == "cif")
    {
        std::ofstream out(filename, std::ios::binary);
        if (!out)
        {
            g_lastSaveError = "Could not write " + filename + "; check that the folder exists and is writable.";
            return false;
        }
        glm::dmat3 cell(1.0);  // columns a, b, c
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) cell[r][c] = structure.cellVectors[(size_t)r][(size_t)c];
        if (std::abs(glm::determinant(cell)) < 1e-9)
        {
            g_lastSaveError = "The unit cell has no volume, so it cannot be written as CIF.";
            return false;
        }
        const glm::dmat3 toFractional = glm::inverse(cell);
        const double la = glm::length(cell[0]), lb = glm::length(cell[1]), lc = glm::length(cell[2]);
        const auto angle = [](const glm::dvec3& u, const glm::dvec3& v) {
            const double cosine = std::max(-1.0, std::min(1.0, glm::dot(u, v) / (glm::length(u) * glm::length(v))));
            return std::acos(cosine) * 180.0 / 3.14159265358979323846;
        };
        const auto clean = [](double v) { return std::abs(v) < 5e-11 ? 0.0 : v; };
        char buffer[512];  // six cell parameters of up to 24 characters each
        out << "# CIF file written by AtomForge\ndata_AtomForge\n";
        std::snprintf(buffer, sizeof(buffer),
                      "_cell_length_a %.10g\n_cell_length_b %.10g\n_cell_length_c %.10g\n"
                      "_cell_angle_alpha %.10g\n_cell_angle_beta %.10g\n_cell_angle_gamma %.10g\n",
                      la, lb, lc, angle(cell[1], cell[2]), angle(cell[0], cell[2]), angle(cell[0], cell[1]));
        out << buffer;
        out << "_symmetry_space_group_name_H-M 'P 1'\n_symmetry_Int_Tables_number 1\n"
               "loop_\n_symmetry_equiv_pos_as_xyz\n'x, y, z'\n"
               "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n";
        std::map<std::string, int> serial;
        for (const auto& site : structure.atoms)
        {
            const std::string symbol = site.symbol.empty() ? std::string(elementSymbol(site.atomicNumber)) : site.symbol;
            const glm::dvec3 f = toFractional * glm::dvec3(site.x, site.y, site.z);
            const std::string label = symbol + std::to_string(++serial[symbol]);
            std::snprintf(buffer, sizeof(buffer), "%-7s %-3s %16.12f %16.12f %16.12f 1\n", label.c_str(), symbol.c_str(),
                          clean(f.x), clean(f.y), clean(f.z));
            out << buffer;
        }
        out.close();
        if (!out)
        {
            g_lastSaveError = "Could not write " + filename + "; check that the folder exists and is writable.";
            return false;
        }
        return saveIpfSidecar(structure, filename);
    }

    // LAMMPS data (read_data, atomic style) is written here. LAMMPS boxes
    // have a along x and b in the xy plane, so the cell and positions are
    // rotated into that frame; each mass names its element ("1 63.546 # Cu").
    // A structure without a cell gets a box around it with 10 A of space.
    if (formatLower == "data" || formatLower == "lmp" || formatLower == "lammps" || formatLower == "lmpdat")
    {
        std::ofstream out(filename, std::ios::binary);
        if (!out)
        {
            g_lastSaveError = "Could not write " + filename + "; check that the folder exists and is writable.";
            return false;
        }
        glm::dmat3 cell(1.0);  // columns a, b, c
        glm::dvec3 origin(0.0);
        if (structure.hasUnitCell)
        {
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c) cell[r][c] = structure.cellVectors[(size_t)r][(size_t)c];
            origin = glm::dvec3(structure.cellOffset[0], structure.cellOffset[1], structure.cellOffset[2]);
        }
        else
        {
            glm::dvec3 low(0.0), high(0.0);
            for (std::size_t i = 0; i < structure.atoms.size(); ++i)
            {
                const glm::dvec3 q(structure.atoms[i].x, structure.atoms[i].y, structure.atoms[i].z);
                low = i ? glm::min(low, q) : q;
                high = i ? glm::max(high, q) : q;
            }
            origin = low - glm::dvec3(10.0);
            cell = glm::dmat3(1.0);
            for (int k = 0; k < 3; ++k) cell[k][k] = high[k] - low[k] + 20.0;
        }
        if (std::abs(glm::determinant(cell)) < 1e-9)
        {
            g_lastSaveError = "The unit cell has no volume, so it cannot be written as a LAMMPS box.";
            return false;
        }
        const glm::dvec3 a = cell[0], b = cell[1], c = cell[2];
        const double ax = glm::length(a);
        const glm::dvec3 ahat = a / ax;
        const double bx = glm::dot(b, ahat);
        const double by = glm::length(glm::cross(ahat, b));
        const double cx = glm::dot(c, ahat);
        const double cy = (glm::dot(b, c) - bx * cx) / by;
        const double cz2 = glm::dot(c, c) - cx * cx - cy * cy;
        const double cz = std::sqrt(std::max(0.0, cz2)) * (glm::determinant(cell) < 0.0 ? -1.0 : 1.0);
        glm::dmat3 lammps(1.0);
        lammps[0] = glm::dvec3(ax, 0.0, 0.0);
        lammps[1] = glm::dvec3(bx, by, 0.0);
        lammps[2] = glm::dvec3(cx, cy, cz);
        const glm::dmat3 toFractional = glm::inverse(cell);
        if (cz < 0.0)
        {
            g_lastSaveError = "The cell is left-handed; LAMMPS needs a right-handed box (swap two cell vectors first).";
            return false;
        }

        std::vector<int> elements;
        for (const auto& site : structure.atoms)
            if (std::find(elements.begin(), elements.end(), site.atomicNumber) == elements.end()) elements.push_back(site.atomicNumber);
        std::sort(elements.begin(), elements.end());

        const auto clean = [](double v) { return std::abs(v) < 5e-11 ? 0.0 : v; };
        char buffer[200];
        out << "LAMMPS data file written by AtomForge\n\n" << structure.atoms.size() << " atoms\n" << elements.size() << " atom types\n\n";
        std::snprintf(buffer, sizeof(buffer), "0.0 %.10g xlo xhi\n0.0 %.10g ylo yhi\n0.0 %.10g zlo zhi\n", clean(ax), clean(by), clean(cz));
        out << buffer;
        if (std::abs(bx) > 1e-10 || std::abs(cx) > 1e-10 || std::abs(cy) > 1e-10)
        {
            std::snprintf(buffer, sizeof(buffer), "%.10g %.10g %.10g xy xz yz\n", clean(bx), clean(cx), clean(cy));
            out << buffer;
        }
        out << "\nMasses\n\n";
        for (std::size_t t = 0; t < elements.size(); ++t)
        {
            std::snprintf(buffer, sizeof(buffer), "%zu %.6g  # %s\n", t + 1, elementAtomicMass(elements[t]), elementSymbol(elements[t]));
            out << buffer;
        }
        out << "\nAtoms  # atomic\n\n";
        for (std::size_t i = 0; i < structure.atoms.size(); ++i)
        {
            const auto& site = structure.atoms[i];
            const glm::dvec3 f = toFractional * (glm::dvec3(site.x, site.y, site.z) - origin);
            const glm::dvec3 q = lammps * f;
            const std::size_t type = static_cast<std::size_t>(std::find(elements.begin(), elements.end(), site.atomicNumber) - elements.begin()) + 1;
            std::snprintf(buffer, sizeof(buffer), "%zu %zu %.10f %.10f %.10f\n", i + 1, type, clean(q.x), clean(q.y), clean(q.z));
            out << buffer;
        }
        out.close();
        if (!out)
        {
            g_lastSaveError = "Could not write " + filename + "; check that the folder exists and is writable.";
            return false;
        }
        return saveIpfSidecar(structure, filename);
    }

    // XYZ is written here: full precision (Open Babel keeps five decimals)
    // and, for a periodic structure, the cell as extended XYZ (ASE).
    if (formatLower == "xyz")
    {
        std::ofstream out(filename, std::ios::binary);
        if (!out)
        {
            g_lastSaveError = "Could not write " + filename + "; check that the folder exists and is writable.";
            return false;
        }
        const auto clean = [](double v) { return std::abs(v) < 5e-11 ? 0.0 : v; };
        char buffer[160];
        out << structure.atoms.size() << "\n";
        if (structure.hasUnitCell)
        {
            out << "Lattice=\"";
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                {
                    std::snprintf(buffer, sizeof(buffer), "%s%.10g", r || c ? " " : "", clean(structure.cellVectors[r][c]));
                    out << buffer;
                }
            out << "\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n";
        }
        else
            out << "AtomForge structure\n";
        for (const auto& site : structure.atoms)
        {
            const std::string symbol = site.symbol.empty() ? std::string(elementSymbol(site.atomicNumber)) : site.symbol;
            std::snprintf(buffer, sizeof(buffer), "%-3s %16.10f %16.10f %16.10f\n", symbol.c_str(), clean(site.x), clean(site.y), clean(site.z));
            out << buffer;
        }
        out.close();
        if (!out)
        {
            g_lastSaveError = "Could not write " + filename + "; check that the folder exists and is writable.";
            return false;
        }
        return saveIpfSidecar(structure, filename);
    }

    // Certain formats (VASP POSCAR, LAMMPS data) require all atoms of the
    // same element to be contiguous.  If atoms are interleaved – as they are
    // after substitutional-solid-solution or other random builders – the
    // writer produces duplicate or malformed element headers.  Stable-sort
    // a working copy by atomic number before building the OBMol.
    const bool needsGrouping = (formatLower == "vasp" || formatLower == "poscar"
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

    // POSCAR is written here: the species as the title (VASP 4 readers take
    // them from there), the lattice as given and Direct coordinates to
    // twelve decimals, atoms grouped by element (stable within each).
    if (formatLower == "vasp" || formatLower == "poscar")
    {
        const Structure& written = *src;
        std::ofstream out(filename, std::ios::binary);
        if (!out)
        {
            g_lastSaveError = "Could not write " + filename + "; check that the folder exists and is writable.";
            return false;
        }
        glm::dmat3 cell(1.0);  // columns a, b, c
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) cell[r][c] = written.cellVectors[(size_t)r][(size_t)c];
        if (std::abs(glm::determinant(cell)) < 1e-9)
        {
            g_lastSaveError = "The unit cell has no volume, so it cannot be written as POSCAR.";
            return false;
        }
        const glm::dmat3 toFractional = glm::inverse(cell);
        std::vector<int> elements;
        for (const auto& site : written.atoms)
            if (std::find(elements.begin(), elements.end(), site.atomicNumber) == elements.end()) elements.push_back(site.atomicNumber);
        std::sort(elements.begin(), elements.end());
        const auto clean = [](double v) { return std::abs(v) < 5e-11 ? 0.0 : v; };
        char buffer[200];
        std::string symbols, counts;
        for (int z : elements)
        {
            const auto count = std::count_if(written.atoms.begin(), written.atoms.end(), [z](const AtomSite& site) { return site.atomicNumber == z; });
            symbols += std::string(symbols.empty() ? "" : " ") + elementSymbol(z);
            counts += std::string(counts.empty() ? "" : " ") + std::to_string(count);
        }
        out << symbols << "\n1.0\n";
        for (int r = 0; r < 3; ++r)
        {
            std::snprintf(buffer, sizeof(buffer), "  %20.12f %20.12f %20.12f\n", clean(written.cellVectors[(size_t)r][0]),
                          clean(written.cellVectors[(size_t)r][1]), clean(written.cellVectors[(size_t)r][2]));
            out << buffer;
        }
        out << symbols << "\n" << counts << "\nDirect\n";
        for (int z : elements)
            for (const auto& site : written.atoms)
            {
                if (site.atomicNumber != z) continue;
                const glm::dvec3 f = toFractional * glm::dvec3(site.x, site.y, site.z);
                std::snprintf(buffer, sizeof(buffer), "  %16.12f %16.12f %16.12f\n", clean(f.x), clean(f.y), clean(f.z));
                out << buffer;
            }
        out.close();
        if (!out)
        {
            g_lastSaveError = "Could not write " + filename + "; check that the folder exists and is writable.";
            return false;
        }
        return saveIpfSidecar(written, filename);
    }

    OpenBabel::OBMol mol;
    mol.BeginModify();

    for (const auto& site : src->atoms)
    {
        OpenBabel::OBAtom* atom = mol.NewAtom();
        atom->SetAtomicNum(site.atomicNumber);
        // Rounding noise would print as "-0.00000".
        const auto clean = [](double v) { return std::abs(v) < 5e-6 ? 0.0 : v; };
        atom->SetVector(clean(site.x), clean(site.y), clean(site.z));
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
    {
        g_lastSaveError = "The " + format + " format is not available (Open Babel plugins missing?).";
        return false;
    }

    conv.SetOutFormat(outFmt);
    if (!conv.WriteFile(&mol, filename))
    {
        g_lastSaveError = "Could not write " + filename + "; check that the folder exists and is writable.";
        return false;
    }

    return saveIpfSidecar(structure, filename);
}
