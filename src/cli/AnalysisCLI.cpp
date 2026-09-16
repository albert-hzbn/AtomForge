#include "cli/AnalysisCLI.h"
#include "algorithms/CommonNeighbourAnalysis.h"
#include "algorithms/AngularDistributionAnalysis.h"
#include "algorithms/RadialDistributionAnalysis.h"
#include "algorithms/ShortRangeOrderAnalysis.h"
#include "algorithms/InterstitialVoidAnalysis.h"
#include "algorithms/CellSculptorAlgo.h"
#include "io/StructureLoader.h"

#include <cmath>
#include <cstdint>
#include <set>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

int runAnalysisCLI(int argc, char* argv[])
{
    for (int i=1;i<argc;++i) if (std::string(argv[i])=="--help") {
        std::cout << "Analysis: --analyze MODE --input FILE --output FILE [options]\n"
            "cna: --cutoff R (0 uses --scale 1.18), --no-pbc\n"
            "rdf: --minimum 0 --cutoff 8 --bins 200 --smooth 0 --centre-z 0 --neighbor-z 0 --raw --no-pbc\n"
            "adf: --cutoff 3.5 --bins 180 --smooth 0 --centre SYMBOL --neighbor-a SYMBOL --neighbor-b SYMBOL --raw --no-pbc\n"
            "sro: --shells 3 --tolerance 0.1\n"
            "interstitial: --resolution 14 --limit 0 --clearance 0.6 --separation 0.8\n"
            "sculpt: --slabs \"h k l lower upper periodic;...\" --nx 1 --ny 1 --nz 1\n"
            "Lengths are Angstrom. Periodic sculpt uses integer start/count bounds.\n"
            "Conversion: --convert --input FILE --output FILE [--format FORMAT]\n"
            "CSV is used for analysis results; sculpt writes VASP.\n";
        return 0;
    }
    try {
        std::map<std::string,std::string> args;
        for (int i=1;i<argc;++i) {
            const std::string key=argv[i];
            if (key=="--convert" || key=="--no-pbc" || key=="--raw") {
                if (!args.emplace(key,"1").second) throw std::invalid_argument("Repeated option "+key);
            }
            else {
                if (i+1==argc) throw std::invalid_argument("Missing value for "+key);
                if (!args.emplace(key,argv[++i]).second) throw std::invalid_argument("Repeated option "+key);
            }
        }
        const auto number=[&](const char* key,double fallback) {
            if (!args.count(key)) return fallback;
            std::size_t used=0;
            const double value=std::stod(args.at(key),&used);
            if (used!=args.at(key).size() || !std::isfinite(value)) throw std::invalid_argument(key);
            return value;
        };
        const auto integer=[&](const char* key,int fallback,int low,int high) {
            const double value=number(key,fallback);
            if (value!=std::floor(value) || value<low || value>high) throw std::invalid_argument(key);
            return static_cast<int>(value);
        };
        if (!args.count("--input") || !args.count("--output"))
            throw std::invalid_argument("--input FILE and --output FILE are required");
        const std::map<std::string,std::set<std::string>> options={
            {"convert",{"--convert","--format"}},
            {"cna",{"--cutoff","--scale","--no-pbc"}},
            {"rdf",{"--cutoff","--minimum","--bins","--smooth","--centre-z","--neighbor-z","--raw","--no-pbc"}},
            {"adf",{"--cutoff","--bins","--smooth","--centre","--neighbor-a","--neighbor-b","--raw","--no-pbc"}},
            {"sro",{"--shells","--tolerance"}},
            {"interstitial",{"--resolution","--limit","--clearance","--separation"}},
            {"sculpt",{"--slabs","--nx","--ny","--nz"}}};
        const auto selected=args.count("--convert") ? "convert" : args.at("--analyze");
        if (!options.count(selected)) throw std::invalid_argument("Unknown analysis: "+selected);
        for (const auto& arg:args)
            if (arg.first!="--input" && arg.first!="--output" &&
                !(arg.first=="--analyze" && selected!="convert") && !options.at(selected).count(arg.first))
                throw std::invalid_argument("Unsupported option for "+selected+": "+arg.first);
        Structure structure; std::string error;
        if (!loadStructureFromFile(args.at("--input"),structure,error)) throw std::runtime_error(error);
        if (args.count("--convert")) {
            std::string format=args.count("--format") ? args.at("--format")
                : std::filesystem::path(args.at("--output")).extension().string().substr(1);
            if (format=="pwi") format="pwscf";
            if (format=="extxyz") {
                std::ofstream file(std::filesystem::u8path(args.at("--output")));
                file << std::setprecision(17) << structure.atoms.size() << '\n';
                if (structure.hasUnitCell) {
                    file << "Lattice=\"";
                    for (const auto& vector:structure.cellVectors) for (double value:vector) file << value << ' ';
                    file << "\" ";
                }
                file << "Properties=species:S:1:pos:R:3\n";
                for (const auto& atom:structure.atoms)
                    file << atom.symbol << ' ' << atom.x << ' ' << atom.y << ' ' << atom.z << '\n';
                file.close();
                if (!file) throw std::runtime_error("Cannot write extended XYZ");
                return 0;
            }
            if (!saveStructure(structure,args.at("--output"),format)) throw std::runtime_error("Conversion failed");
            return 0;
        }
        std::ostringstream out;
        out << std::setprecision(17);
        const auto mode=args.at("--analyze");
        if (mode=="sculpt") {
            std::vector<CellSlabPlane> slabs;
            std::istringstream input(args.at("--slabs"));
            std::string row;
            while (std::getline(input,row,';')) {
                CellSlabPlane slab;
                int periodic=0;
                std::istringstream values(row);
                if (!(values>>slab.h>>slab.k>>slab.l>>slab.d1>>slab.d2>>periodic) ||
                    (slab.h==0 && slab.k==0 && slab.l==0) || !std::isfinite(slab.d1) || !std::isfinite(slab.d2) ||
                    slab.h < -1000 || slab.h > 1000 || slab.k < -1000 || slab.k > 1000 ||
                    slab.l < -1000 || slab.l > 1000 || (periodic!=0 && periodic!=1))
                    throw std::invalid_argument("Slabs use h k l lower upper periodic, separated by semicolons");
                slab.usePeriodic=periodic!=0;
                if (slab.usePeriodic) {
                    if (!structure.hasUnitCell || slab.d1!=std::floor(slab.d1) || slab.d2!=std::floor(slab.d2) ||
                        std::abs(slab.d1)>10000 || slab.d2<1 || slab.d2>10000)
                        throw std::invalid_argument("Periodic slabs need a cell, integer start and positive period count");
                    slab.startPlane=static_cast<int>(slab.d1); slab.nPeriods=static_cast<int>(slab.d2);
                } else if (slab.d2<slab.d1) throw std::invalid_argument("Slab upper bound is below its lower bound");
                slabs.push_back(slab);
            }
            if (slabs.empty()) throw std::invalid_argument("At least one slab is required");
            const int nx=integer("--nx",1,1,100),ny=integer("--ny",1,1,100),nz=integer("--nz",1,1,100);
            if (static_cast<std::uint64_t>(nx)*ny*nz*structure.atoms.size()>10000000)
                throw std::invalid_argument("Sculpted supercell exceeds 10 million atoms");
            const auto supercell=cscBuildSupercell(structure,nx,ny,nz);
            const auto result=cscApplySlabs(supercell,slabs,structure);
            if (!saveStructure(result,args.at("--output"),"vasp")) throw std::runtime_error("Cannot save sculpted structure");
            return 0;
        }
        const bool periodic=!args.count("--no-pbc");
        if (mode=="cna") {
            atomforge::analysis::CnaParams p;
            p.usePbc=periodic; p.cutoffRadius=number("--cutoff",0); p.cutoffScale=number("--scale",1.18);
            const auto result=atomforge::analysis::computeCna(structure,p);
            if (!result.valid) throw std::runtime_error(result.message);
            out << "index,symbol,coordination,common,bonds,chain,environment\n";
            for (const auto& row:result.atomRows)
                out << row.index << ',' << row.symbol << ',' << row.coordination << ','
                    << row.dominantSignature.common << ',' << row.dominantSignature.bonds << ','
                    << row.dominantSignature.chain << ',' << row.environment << '\n';
        } else if (mode=="rdf") {
            atomforge::analysis::RdfParams p;
            p.usePbc=periodic; p.normalize=!args.count("--raw");
            p.rMin=number("--minimum",0); p.rMax=number("--cutoff",8);
            p.binCount=integer("--bins",200,1,100000);
            p.refAtomicNumberFilter=integer("--centre-z",0,0,118);
            p.targetAtomicNumberFilter=integer("--neighbor-z",0,0,118);
            p.smoothingPasses=integer("--smooth",0,0,1000);
            const auto result=atomforge::analysis::computeRdf(structure,p);
            if (!result.valid) throw std::runtime_error(result.message);
            out << "radius_A,g,raw_count,coordination\n";
            for (const auto& row:result.bins) out << row.rCenter << ',' << row.g << ',' << row.rawCount << ',' << row.cumulative << '\n';
        } else if (mode=="adf") {
            AdfParams p;
            p.usePbc=periodic; p.normalize=!args.count("--raw"); p.rCutoff=number("--cutoff",3.5);
            p.binCount=integer("--bins",180,1,100000); p.smoothPasses=integer("--smooth",0,0,1000);
            if (args.count("--centre")) { p.centreMode=AdfCentreMode::ByElement; p.centreSymbol=args.at("--centre"); }
            if (args.count("--neighbor-a") || args.count("--neighbor-b")) {
                if (p.centreSymbol.empty()) throw std::invalid_argument("Pair selection requires --centre");
                p.centreMode=AdfCentreMode::ByPair;
                p.neighSymbol1=args["--neighbor-a"]; p.neighSymbol2=args["--neighbor-b"];
            }
            const auto result=computeADF(structure,p);
            if (!result.valid) throw std::runtime_error(result.message);
            out << "angle_deg,raw_count,value\n";
            for (const auto& row:result.bins) out << row.angleDeg << ',' << row.count << ',' << row.value << '\n';
        } else if (mode=="sro") {
            const double tolerance=number("--tolerance",.1);
            if (tolerance<=0) throw std::invalid_argument("Tolerance must be positive");
            const auto result=analyzeSRO(structure,integer("--shells",3,1,100),tolerance);
            if (!result.success) throw std::runtime_error(result.message);
            out << "method,shell,distance_A,centre,neighbor,alpha\n";
            for (const auto& row:result.warrenCowley.entries)
                out << "Warren-Cowley," << row.shell << ',' << row.shellDistance << ',' << row.centerElement << ',' << row.neighborElement << ',' << row.alpha << '\n';
            for (const auto& row:result.raoCurtin.entries)
                out << "Rao-Curtin," << row.shell << ',' << row.shellDistance << ',' << row.elem_i << ',' << row.elem_j << ',' << row.alpha << '\n';
        } else if (mode=="interstitial") {
            InterstitialVoidDetectionParams p;
            p.gridResolution=integer("--resolution",14,2,100);
            p.maxVoids=integer("--limit",0,0,100000);
            p.minClearance=number("--clearance",.6); p.minSeparation=number("--separation",.8);
            if (p.minClearance<0 || p.minSeparation<=0) throw std::invalid_argument("Invalid clearance/separation");
            const auto result=detectInterstitialVoidRegions(structure,p);
            if (!result.success) throw std::runtime_error(result.message);
            out << "x_A,y_A,z_A,clearance_A,coordination,volume_A3,kind\n";
            for (const auto& row:result.regions)
                out << row.position.x << ',' << row.position.y << ',' << row.position.z << ',' << row.clearance << ',' << row.coordination << ',' << row.volume << ',' << static_cast<int>(row.kind) << '\n';
        } else throw std::invalid_argument("Analysis must be cna, rdf, adf, sro or interstitial");
        std::ofstream file(args.at("--output"));
        file << out.str(); file.close();
        if (!file) throw std::runtime_error("Cannot write analysis output");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
