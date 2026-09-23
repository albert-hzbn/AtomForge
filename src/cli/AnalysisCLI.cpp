#include "cli/AnalysisCLI.h"
#include "algorithms/CommonNeighbourAnalysis.h"
#include "algorithms/AngularDistributionAnalysis.h"
#include "algorithms/RadialDistributionAnalysis.h"
#include "algorithms/ShortRangeOrderAnalysis.h"
#include "algorithms/InterstitialVoidAnalysis.h"
#include "algorithms/CellSculptorAlgo.h"
#include "algorithms/DislocationFit.h"
#include "algorithms/DragPrep.h"
#include "algorithms/NyeTensor.h"
#include "algorithms/PatternMatch.h"
#include "algorithms/VitekMap.h"
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
            "nye: --reference FILE --cutoff 3.0 --no-pbc  (dislocation density tensor;\n"
            "     --input is the deformed structure, --reference the perfect one, same\n"
            "     atom count/order in both; see algorithms/NyeTensor.h)\n"
            "vitek: --reference FILE --line \"u v w\" --burgers M --cutoff 0 --no-pbc\n"
            "     (differential-displacement map; --line is the dislocation line/Burgers\n"
            "     direction, --burgers its magnitude in Angstrom -- 0 disables screw-\n"
            "     component wrapping; --cutoff 0 auto-detects the neighbor shell;\n"
            "     see algorithms/VitekMap.h)\n"
            "pattern-init: --cutoff C --no-pbc  (--input is the perfect reference\n"
            "     structure; writes a dx,dy,dz pattern file to --output; see\n"
            "     algorithms/PatternMatch.h)\n"
            "pattern-detect: --pattern FILE --cutoff C --angle-threshold 10 --no-pbc\n"
            "     (--input is the possibly-defective structure, --pattern a file\n"
            "     written by pattern-init)\n"
            "drag: --final FILE --zeta Z --clip-displacement --no-pbc  (prepares a\n"
            "     constrained-minimization/migration-barrier calculation; --input is\n"
            "     the initial configuration, --final the final one, same atom count/\n"
            "     order; writes the interpolated structure to --output and the\n"
            "     per-atom reaction-coordinate direction to --output + '.constraint.csv';\n"
            "     see algorithms/DragPrep.h)\n"
            "fit-dislocation: --reference FILE --line \"u v w\" --cutoff 3.0 --area 1.0\n"
            "     --no-pbc  (recovers a dislocation's position/Burgers vector from a\n"
            "     measured field, by Nye-tensor centroid/integral; --input is the\n"
            "     deformed structure, --reference the perfect one; --area is the\n"
            "     cross-sectional area per atom (A^2) used to scale the Burgers vector\n"
            "     estimate; see algorithms/DislocationFit.h)\n"
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
            if (key=="--convert" || key=="--no-pbc" || key=="--raw" || key=="--clip-displacement") {
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
            {"nye",{"--reference","--cutoff","--no-pbc"}},
            {"vitek",{"--reference","--line","--burgers","--cutoff","--no-pbc"}},
            {"pattern-init",{"--cutoff","--no-pbc"}},
            {"pattern-detect",{"--pattern","--cutoff","--angle-threshold","--no-pbc"}},
            {"drag",{"--final","--zeta","--clip-displacement","--no-pbc"}},
            {"fit-dislocation",{"--reference","--line","--cutoff","--area","--no-pbc"}},
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
        if (mode=="drag") {
            if (!args.count("--final")) throw std::invalid_argument("drag requires --final FILE (the final configuration)");
            Structure finalStructure; std::string finalError;
            if (!loadStructureFromFile(args.at("--final"),finalStructure,finalError)) throw std::runtime_error(finalError);
            const double zeta=number("--zeta",0.5);
            const auto result=atomforge::prepareDrag(structure,finalStructure,zeta,args.count("--clip-displacement")>0,periodic);
            if (!result.success) throw std::runtime_error(result.message);
            std::string fmt=std::filesystem::path(args.at("--output")).extension().string();
            if (!fmt.empty() && fmt[0]=='.') fmt=fmt.substr(1);
            if (fmt.empty()) fmt="vasp";
            if (!saveStructure(result.interpolated,args.at("--output"),fmt)) throw std::runtime_error("Cannot save interpolated structure");
            std::ofstream constraintFile(args.at("--output")+".constraint.csv");
            constraintFile << std::setprecision(17) << "index,dx,dy,dz\n";
            for (std::size_t i=0;i<result.constraintDirection.size();++i)
                constraintFile << i << ',' << result.constraintDirection[i].x << ',' << result.constraintDirection[i].y << ',' << result.constraintDirection[i].z << '\n';
            constraintFile.close();
            if (!constraintFile) throw std::runtime_error("Cannot write constraint file");
            std::cout << result.message << "\n";
            return 0;
        }
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
        } else if (mode=="nye") {
            if (!args.count("--reference")) throw std::invalid_argument("nye requires --reference FILE (the undeformed structure)");
            Structure reference; std::string refError;
            if (!loadStructureFromFile(args.at("--reference"),reference,refError)) throw std::runtime_error(refError);
            const float cutoff=static_cast<float>(number("--cutoff",3.0));
            if (!(cutoff>0)) throw std::invalid_argument("--cutoff must be positive");
            const auto result=atomforge::computeNyeTensor(reference,structure,cutoff,periodic);
            if (!result.success) throw std::runtime_error(result.message);
            out << "index,symbol,alpha_xx,alpha_xy,alpha_xz,alpha_yx,alpha_yy,alpha_yz,alpha_zx,alpha_zy,alpha_zz,norm\n";
            for (std::size_t i=0;i<result.alpha.size();++i) {
                out << i << ',' << (i<structure.atoms.size() ? structure.atoms[i].symbol : std::string());
                for (int a=0;a<3;++a) for (int b=0;b<3;++b) out << ',' << result.alpha[i][a][b];
                out << ',' << result.alphaNorm[i] << '\n';
            }
        } else if (mode=="fit-dislocation") {
            if (!args.count("--reference")) throw std::invalid_argument("fit-dislocation requires --reference FILE (the undeformed structure)");
            if (!args.count("--line")) throw std::invalid_argument("fit-dislocation requires --line \"u v w\" (Cartesian assumed line direction)");
            Structure reference; std::string refError;
            if (!loadStructureFromFile(args.at("--reference"),reference,refError)) throw std::runtime_error(refError);
            std::istringstream lineStream(args.at("--line"));
            double lx=0,ly=0,lz=0;
            if (!(lineStream>>lx>>ly>>lz)) throw std::invalid_argument("Cannot parse --line \"u v w\"");
            const float cutoff=static_cast<float>(number("--cutoff",3.0));
            if (!(cutoff>0)) throw std::invalid_argument("--cutoff must be positive");
            const double area=number("--area",1.0);
            const auto nyeResult=atomforge::computeNyeTensor(reference,structure,cutoff,periodic);
            if (!nyeResult.success) throw std::runtime_error(nyeResult.message);
            std::vector<glm::dvec3> positions(reference.atoms.size());
            for (std::size_t i=0;i<positions.size();++i)
                positions[i]=glm::dvec3(reference.atoms[i].x,reference.atoms[i].y,reference.atoms[i].z);
            const auto fit=atomforge::fitDislocationFromNye(positions,nyeResult,glm::dvec3(lx,ly,lz),area);
            if (!fit.success) throw std::runtime_error(fit.message);
            out << "quantity,x,y,z\n";
            out << "line_position," << fit.linePosition.x << ',' << fit.linePosition.y << ',' << fit.linePosition.z << '\n';
            out << "burgers_vector," << fit.burgersVector.x << ',' << fit.burgersVector.y << ',' << fit.burgersVector.z << '\n';
        } else if (mode=="vitek") {
            if (!args.count("--reference")) throw std::invalid_argument("vitek requires --reference FILE (the undeformed structure)");
            if (!args.count("--line")) throw std::invalid_argument("vitek requires --line \"u v w\" (Cartesian dislocation line/Burgers direction)");
            Structure reference; std::string refError;
            if (!loadStructureFromFile(args.at("--reference"),reference,refError)) throw std::runtime_error(refError);
            std::istringstream lineStream(args.at("--line"));
            double lx=0,ly=0,lz=0;
            if (!(lineStream>>lx>>ly>>lz)) throw std::invalid_argument("Cannot parse --line \"u v w\"");
            const double burgers=number("--burgers",0);
            const float cutoff=static_cast<float>(number("--cutoff",0));
            const auto result=atomforge::computeVitekMap(reference,structure,glm::dvec3(lx,ly,lz),burgers,cutoff,periodic);
            if (!result.success) throw std::runtime_error(result.message);
            out << "index_i,index_j,xi_A,yi_A,zi_A,xj_A,yj_A,zj_A,screw_A,edge_1_A,edge_2_A\n";
            for (const auto& p:result.pairs)
                out << p.indexI << ',' << p.indexJ << ','
                    << p.positionI.x << ',' << p.positionI.y << ',' << p.positionI.z << ','
                    << p.positionJ.x << ',' << p.positionJ.y << ',' << p.positionJ.z << ','
                    << p.screwComponent << ',' << p.edgeComponent.x << ',' << p.edgeComponent.y << '\n';
        } else if (mode=="pattern-init") {
            const float cutoff=static_cast<float>(number("--cutoff",0));
            if (!(cutoff>0)) throw std::invalid_argument("--cutoff must be positive");
            const auto result=atomforge::buildPattern(structure,cutoff,periodic);
            if (!result.success) throw std::runtime_error(result.message);
            out << "dx,dy,dz\n";
            for (const auto& d:result.directions) out << d.x << ',' << d.y << ',' << d.z << '\n';
        } else if (mode=="pattern-detect") {
            if (!args.count("--pattern")) throw std::invalid_argument("pattern-detect requires --pattern FILE (written by pattern-init)");
            const float cutoff=static_cast<float>(number("--cutoff",0));
            if (!(cutoff>0)) throw std::invalid_argument("--cutoff must be positive");
            const double angleThreshold=number("--angle-threshold",10.0);
            atomforge::CrystalPattern pattern;
            std::ifstream patternFile(args.at("--pattern"));
            if (!patternFile) throw std::runtime_error("Cannot read pattern file '"+args.at("--pattern")+"'");
            std::string header; std::getline(patternFile,header);
            std::string row;
            while (std::getline(patternFile,row)) {
                if (row.empty()) continue;
                std::istringstream rs(row);
                double dx,dy,dz; char comma;
                if (!(rs>>dx>>comma>>dy>>comma>>dz)) throw std::runtime_error("Malformed pattern file row: "+row);
                pattern.directions.push_back(glm::dvec3(dx,dy,dz));
            }
            pattern.success=!pattern.directions.empty();
            if (!pattern.success) throw std::runtime_error("Pattern file contains no directions");
            const auto result=atomforge::detectPattern(structure,pattern,cutoff,angleThreshold,periodic);
            if (!result.success) throw std::runtime_error(result.message);
            out << "index,matched,neighbor_count,max_angle_deviation_deg\n";
            for (const auto& matchRow:result.rows)
                out << matchRow.index << ',' << (matchRow.matched?1:0) << ',' << matchRow.neighborCount << ',' << matchRow.maxAngleDeviationDeg << '\n';
        } else throw std::invalid_argument("Analysis must be cna, rdf, adf, sro, interstitial, nye, vitek, pattern-init or pattern-detect");
        std::ofstream file(args.at("--output"));
        file << out.str(); file.close();
        if (!file) throw std::runtime_error("Cannot write analysis output");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
