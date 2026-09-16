#include "model/Workspace.h"
#include "util/TaskControl.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace atomforge
{
namespace
{
class Archive
{
public:
    explicit Archive(std::istream& input, std::uint64_t bytes) : in(&input), remaining(bytes) {}
    explicit Archive(std::ostream& output) : out(&output) {}
    bool reading() const { return in != nullptr; }
    void bytes(char* data, std::size_t count)
    {
        taskCheckpoint();
        if (in) {
            if (count>remaining) throw std::runtime_error("Truncated project");
            in->read(data,static_cast<std::streamsize>(count)); remaining-=count;
            if (!*in) throw std::runtime_error("Cannot read project");
        } else {
            out->write(data,static_cast<std::streamsize>(count));
            if (!*out) throw std::runtime_error("Cannot write project");
        }
    }
    void value(std::uint64_t& number)
    {
        char buffer[8];
        for (int i=0;i<8;++i) buffer[i]=static_cast<char>((number>>(8*i))&255);
        bytes(buffer,8);
        if (in) { number=0; for (int i=0;i<8;++i) number|=std::uint64_t(static_cast<unsigned char>(buffer[i]))<<(8*i); }
    }
    void value(double& number)
    {
        static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559);
        std::uint64_t bits=0;
        std::memcpy(&bits,&number,8); value(bits); if (in) std::memcpy(&number,&bits,8);
        if (!std::isfinite(number)) throw std::runtime_error("Non-finite project value");
    }
    void value(float& number) { double copy=number; value(copy); if (in) number=static_cast<float>(copy); if (!std::isfinite(number)) throw std::runtime_error("Float out of range"); }
    void value(int& number) { double copy=number; value(copy); if (copy!=std::floor(copy) || copy<std::numeric_limits<int>::min() || copy>std::numeric_limits<int>::max()) throw std::runtime_error("Invalid project integer"); if (in) number=static_cast<int>(copy); }
    void value(bool& number) { int copy=number; value(copy); if (copy!=0 && copy!=1) throw std::runtime_error("Invalid project boolean"); if (in) number=copy!=0; }
    void value(std::string& text)
    {
        std::uint64_t count=text.size(); value(count);
        if (count>1048576 || (in && count>remaining)) throw std::runtime_error("Project string too large");
        if (in) text.resize(static_cast<std::size_t>(count));
        bytes(text.data(),text.size());
    }
    template<class T,std::size_t N> void value(std::array<T,N>& values) { for (auto& item:values) value(item); }
    void value(std::map<std::string,double>& values)
    {
        std::uint64_t count=values.size(); value(count);
        if (count>10000) throw std::runtime_error("Too many project settings");
        if (in) {
            for (std::uint64_t i=0;i<count;++i) {
                std::string key; double number=0; value(key); value(number);
                if (!values.emplace(key,number).second) throw std::runtime_error("Duplicate project setting");
            }
        } else {
            for (const auto& pair:values) {
                std::string key=pair.first; double number=pair.second; value(key); value(number);
            }
        }
    }
    template<glm::length_t N,class T,glm::qualifier Q> void value(glm::vec<N,T,Q>& values) { for (glm::length_t i=0;i<N;++i) value(values[i]); }
    void value(glm::dmat3& matrix) { for (int i=0;i<3;++i) value(matrix[i]); }
    template<class T> void value(std::vector<T>& values)
    {
        std::uint64_t count=values.size(); value(count);
        if (count>100000000 || count>2147483648ull/sizeof(T) || (in && count>remaining/8)) throw std::runtime_error("Project array too large");
        if (in) values.resize(static_cast<std::size_t>(count));
        for (auto& item:values) value(item);
    }
    void value(AtomSite& a) { value(a.symbol); value(a.atomicNumber); value(a.x); value(a.y); value(a.z); value(a.r); value(a.g); value(a.b); }
    void value(Structure& s) {
        value(s.atoms); value(s.hasUnitCell); value(s.cellVectors); value(s.cellOffset);
        value(s.pbcBoundaryTol); value(s.grainColors); value(s.grainRegionIds); value(s.ipfLoadStatus);
        value(s.dislocationLoopPoints); value(s.dislocationDetectionDone);
        if ((!s.grainColors.empty() && s.grainColors.size()!=s.atoms.size()) ||
            (!s.grainRegionIds.empty() && s.grainRegionIds.size()!=s.atoms.size()))
            throw std::runtime_error("Misaligned project atom metadata");
    }
    void value(electronic::Site& s) { value(s.number); value(s.position); }
    void value(electronic::Grid& g) { value(g.shape); value(g.cell); value(g.origin); value(g.periodic); value(g.name); value(g.unit); value(g.values); g.validate(); }
    void value(electronic::Volume& v) { value(v.sites); value(v.fields); }
    void value(electronic::Mesh& m) {
        value(m.vertices); value(m.colors);
        if (m.vertices.size()%3 || (!m.colors.empty() && m.colors.size()!=m.vertices.size())) throw std::runtime_error("Invalid project mesh");
    }
    void value(Workspace& w) {
        value(w.structure); value(w.volume); value(w.reference); value(w.surface);
        value(w.table); value(w.columns); value(w.heading); value(w.settings); value(w.sliceSettings); value(w.camera);
        value(w.sourcePath); value(w.referencePath); value(w.title); value(w.history);
        if (!w.camera.empty() && w.camera.size()!=7) throw std::runtime_error("Invalid project camera");
        for (double coordinate:w.camera)
            if (std::abs(coordinate)>1e9) throw std::runtime_error("Project camera outside supported range");
        if (w.columns<0 || w.columns>10000 || (!w.table.empty() && (w.columns==0 || w.table.size()%w.columns))) throw std::runtime_error("Invalid project table");
    }
    std::uint64_t left() const { return remaining; }
private:
    std::istream* in=nullptr;
    std::ostream* out=nullptr;
    std::uint64_t remaining=0;
};
}

void saveWorkspace(const std::vector<Workspace>& tabs,const std::string& filename)
{
    if (tabs.empty() || tabs.size()>1000) throw std::runtime_error("A project requires 1-1000 tabs");
    const auto path=std::filesystem::u8path(filename);
    auto temporary=path;
    temporary += ".tmp-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    try {
        std::ofstream stream(temporary,std::ios::binary);
        Archive archive(static_cast<std::ostream&>(stream));
        std::string magic="ATOMFORGE_PROJECT_1"; archive.value(magic);
        // The serializer does not mutate values while writing.
        auto& data=const_cast<std::vector<Workspace>&>(tabs);
        archive.value(data);
        stream.close();
        if (!stream) throw std::runtime_error("Failed to flush project");
#ifdef _WIN32
        if (!MoveFileExW(temporary.c_str(),path.c_str(),MOVEFILE_REPLACE_EXISTING|MOVEFILE_WRITE_THROUGH))
            throw std::runtime_error("Cannot replace project file");
#else
        std::filesystem::rename(temporary,path);
#endif
    } catch (...) {
        std::error_code error; std::filesystem::remove(temporary,error); throw;
    }
}

std::vector<Workspace> loadWorkspace(const std::string& filename)
{
    const auto path=std::filesystem::u8path(filename);
    const auto size=std::filesystem::file_size(path);
    if (size>8ull*1024*1024*1024) throw std::runtime_error("Project exceeds 8 GiB limit");
    std::ifstream stream(path,std::ios::binary);
    Archive archive(stream,size);
    std::string magic; archive.value(magic);
    if (magic!="ATOMFORGE_PROJECT_1") throw std::runtime_error("Unsupported project version");
    std::uint64_t count=0; archive.value(count);
    if (count==0 || count>1000) throw std::runtime_error("Invalid project tab count");
    std::vector<Workspace> result(static_cast<std::size_t>(count));
    for (auto& tab:result) archive.value(tab);
    if (archive.left()) throw std::runtime_error("Unexpected trailing project data");
    return result;
}
}
