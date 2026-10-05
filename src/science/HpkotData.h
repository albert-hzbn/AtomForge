// HPKOT band-path tables from SeeK-path (band_path_data, MIT licence; point
// coordinates provided by Y. Hinuma). Hinuma, Pizzi, Kumagai, Oba, Tanaka,
// Comp. Mat. Sci. 128, 140 (2017), doi:10.1016/j.commatsci.2016.10.015.
#pragma once

namespace atomforge::science::hpkot
{
struct Parameter { const char* name; const char* expression; };
struct Point { const char* label; const char* x; const char* y; const char* z; };
struct Segment { const char* from; const char* to; };
struct Lattice { const char* name; const Parameter* parameters; int parameterCount; const Point* points; int pointCount; const Segment* path; int segmentCount; };

inline const Parameter aP2Parameters[] = {{nullptr, nullptr}};
inline const Point aP2Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"Z", "0", "0", "1/2"},
    {"Y", "0", "1/2", "0"},
    {"X", "1/2", "0", "0"},
    {"V", "1/2", "1/2", "0"},
    {"U", "1/2", "0", "1/2"},
    {"T", "0", "1/2", "1/2"},
    {"R", "1/2", "1/2", "1/2"},
};
inline const Segment aP2Path[] = {{"GAMMA", "X"}, {"Y", "GAMMA"}, {"GAMMA", "Z"}, {"R", "GAMMA"}, {"GAMMA", "T"}, {"U", "GAMMA"}, {"GAMMA", "V"}};

inline const Parameter aP3Parameters[] = {{nullptr, nullptr}};
inline const Point aP3Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"Z", "0", "0", "1/2"},
    {"Y", "0", "1/2", "0"},
    {"Y_2", "0", "-1/2", "0"},
    {"X", "1/2", "0", "0"},
    {"V_2", "1/2", "-1/2", "0"},
    {"U_2", "-1/2", "0", "1/2"},
    {"T_2", "0", "-1/2", "1/2"},
    {"R_2", "-1/2", "-1/2", "1/2"},
};
inline const Segment aP3Path[] = {{"GAMMA", "X"}, {"Y", "GAMMA"}, {"GAMMA", "Z"}, {"R_2", "GAMMA"}, {"GAMMA", "T_2"}, {"U_2", "GAMMA"}, {"GAMMA", "V_2"}};

inline const Parameter cF1Parameters[] = {{nullptr, nullptr}};
inline const Point cF1Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"X", "1/2", "0", "1/2"},
    {"L", "1/2", "1/2", "1/2"},
    {"W", "1/2", "1/4", "3/4"},
    {"W_2", "3/4", "1/4", "1/2"},
    {"K", "3/8", "3/8", "3/4"},
    {"U", "5/8", "1/4", "5/8"},
};
inline const Segment cF1Path[] = {{"GAMMA", "X"}, {"X", "U"}, {"K", "GAMMA"}, {"GAMMA", "L"}, {"L", "W"}, {"W", "X"}, {"X", "W_2"}};

inline const Parameter cF2Parameters[] = {{nullptr, nullptr}};
inline const Point cF2Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"X", "1/2", "0", "1/2"},
    {"L", "1/2", "1/2", "1/2"},
    {"W", "1/2", "1/4", "3/4"},
    {"W_2", "3/4", "1/4", "1/2"},
    {"K", "3/8", "3/8", "3/4"},
    {"U", "5/8", "1/4", "5/8"},
};
inline const Segment cF2Path[] = {{"GAMMA", "X"}, {"X", "U"}, {"K", "GAMMA"}, {"GAMMA", "L"}, {"L", "W"}, {"W", "X"}};

inline const Parameter cI1Parameters[] = {{nullptr, nullptr}};
inline const Point cI1Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"H", "1/2", "-1/2", "1/2"},
    {"P", "1/4", "1/4", "1/4"},
    {"N", "0", "0", "1/2"},
};
inline const Segment cI1Path[] = {{"GAMMA", "H"}, {"H", "N"}, {"N", "GAMMA"}, {"GAMMA", "P"}, {"P", "H"}, {"P", "N"}};

inline const Parameter cP1Parameters[] = {{nullptr, nullptr}};
inline const Point cP1Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"R", "1/2", "1/2", "1/2"},
    {"M", "1/2", "1/2", "0"},
    {"X", "0", "1/2", "0"},
    {"X_1", "1/2", "0", "0"},
};
inline const Segment cP1Path[] = {{"GAMMA", "X"}, {"X", "M"}, {"M", "GAMMA"}, {"GAMMA", "R"}, {"R", "X"}, {"R", "M"}, {"M", "X_1"}};

inline const Parameter cP2Parameters[] = {{nullptr, nullptr}};
inline const Point cP2Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"R", "1/2", "1/2", "1/2"},
    {"M", "1/2", "1/2", "0"},
    {"X", "0", "1/2", "0"},
    {"X_1", "1/2", "0", "0"},
};
inline const Segment cP2Path[] = {{"GAMMA", "X"}, {"X", "M"}, {"M", "GAMMA"}, {"GAMMA", "R"}, {"R", "X"}, {"R", "M"}};

inline const Parameter hP1Parameters[] = {{nullptr, nullptr}};
inline const Point hP1Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"A", "0", "0", "1/2"},
    {"K", "1/3", "1/3", "0"},
    {"H", "1/3", "1/3", "1/2"},
    {"H_2", "1/3", "1/3", "-1/2"},
    {"M", "1/2", "0", "0"},
    {"L", "1/2", "0", "1/2"},
};
inline const Segment hP1Path[] = {{"GAMMA", "M"}, {"M", "K"}, {"K", "GAMMA"}, {"GAMMA", "A"}, {"A", "L"}, {"L", "H"}, {"H", "A"}, {"L", "M"}, {"H", "K"}, {"K", "H_2"}};

inline const Parameter hP2Parameters[] = {{nullptr, nullptr}};
inline const Point hP2Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"A", "0", "0", "1/2"},
    {"K", "1/3", "1/3", "0"},
    {"H", "1/3", "1/3", "1/2"},
    {"H_2", "1/3", "1/3", "-1/2"},
    {"M", "1/2", "0", "0"},
    {"L", "1/2", "0", "1/2"},
};
inline const Segment hP2Path[] = {{"GAMMA", "M"}, {"M", "K"}, {"K", "GAMMA"}, {"GAMMA", "A"}, {"A", "L"}, {"L", "H"}, {"H", "A"}, {"L", "M"}, {"H", "K"}};

inline const Parameter hR1Parameters[] = {{"D", "a*a/4/c/c"}, {"Y", "5/6-2*D"}, {"N", "1/3+D"} };
inline const Point hR1Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"T", "1/2", "1/2", "1/2"},
    {"L", "1/2", "0", "0"},
    {"L_2", "0", "-1/2", "0"},
    {"L_4", "0", "0", "-1/2"},
    {"F", "1/2", "0", "1/2"},
    {"F_2", "1/2", "1/2", "0"},
    {"S_0", "N", "-N", "0"},
    {"S_2", "1-N", "0", "N"},
    {"S_4", "N", "0", "-N"},
    {"S_6", "1-N", "N", "0"},
    {"H_0", "1/2", "-1+Y", "1-Y"},
    {"H_2", "Y", "1-Y", "1/2"},
    {"H_4", "Y", "1/2", "1-Y"},
    {"H_6", "1/2", "1-Y", "-1+Y"},
    {"M_0", "N", "-1+Y", "N"},
    {"M_2", "1-N", "1-Y", "1-N"},
    {"M_4", "Y", "N", "N"},
    {"M_6", "1-N", "1-N", "1-Y"},
    {"M_8", "N", "N", "-1+Y"},
};
inline const Segment hR1Path[] = {{"GAMMA", "T"}, {"T", "H_2"}, {"H_0", "L"}, {"L", "GAMMA"}, {"GAMMA", "S_0"}, {"S_2", "F"}, {"F", "GAMMA"}};

inline const Parameter hR2Parameters[] = {{"Z", "1/6-c*c/9/a/a"}, {"H", "1/2-2*Z"}, {"N", "1/2+Z"} };
inline const Point hR2Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"T", "1/2", "-1/2", "1/2"},
    {"P_0", "H", "-1+H", "H"},
    {"P_2", "H", "H", "H"},
    {"R_0", "1-H", "-H", "-H"},
    {"M", "1-N", "-N", "1-N"},
    {"M_2", "N", "-1+N", "-1+N"},
    {"L", "1/2", "0", "0"},
    {"F", "1/2", "-1/2", "0"},
};
inline const Segment hR2Path[] = {{"GAMMA", "L"}, {"L", "T"}, {"T", "P_0"}, {"P_2", "GAMMA"}, {"GAMMA", "F"}};

inline const Parameter mC1Parameters[] = {{"Z", "(2+a/c*cosbeta)/4/sinbeta/sinbeta"}, {"H", "1/2-2*Z*c*cosbeta/a"}, {"S", "3/4-b*b/4/a/a/sinbeta/sinbeta"}, {"P", "S-(3/4-S)*a*cosbeta/c"} };
inline const Point mC1Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"Y_2", "-1/2", "1/2", "0"},
    {"Y_4", "1/2", "-1/2", "0"},
    {"A", "0", "0", "1/2"},
    {"M_2", "-1/2", "1/2", "1/2"},
    {"V", "1/2", "0", "0"},
    {"V_2", "0", "1/2", "0"},
    {"L_2", "0", "1/2", "1/2"},
    {"C", "1-S", "1-S", "0"},
    {"C_2", "-1+S", "S", "0"},
    {"C_4", "S", "-1+S", "0"},
    {"D", "-1+P", "P", "1/2"},
    {"D_2", "1-P", "1-P", "1/2"},
    {"E", "-1+Z", "1-Z", "1-H"},
    {"E_2", "-Z", "Z", "H"},
    {"E_4", "Z", "-Z", "1-H"},
};
inline const Segment mC1Path[] = {{"GAMMA", "C"}, {"C_2", "Y_2"}, {"Y_2", "GAMMA"}, {"GAMMA", "M_2"}, {"M_2", "D"}, {"D_2", "A"}, {"A", "GAMMA"}, {"L_2", "GAMMA"}, {"GAMMA", "V_2"}};

inline const Parameter mC2Parameters[] = {{"Z", "(a*a/b/b+(1+a/c*cosbeta)/sinbeta/sinbeta)/4"}, {"M", "(1+a*a/b/b)/4"}, {"D", "-a*c*cosbeta/2/b/b"}, {"X", "1/2-2*Z*c*cosbeta/a"}, {"P", "1+Z-2*M"}, {"S", "X-2*D"} };
inline const Point mC2Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"Y", "1/2", "1/2", "0"},
    {"A", "0", "0", "1/2"},
    {"M", "1/2", "1/2", "1/2"},
    {"V_2", "0", "1/2", "0"},
    {"L_2", "0", "1/2", "1/2"},
    {"F", "-1+P", "1-P", "1-S"},
    {"F_2", "1-P", "P", "S"},
    {"F_4", "P", "1-P", "1-S"},
    {"H", "-Z", "Z", "X"},
    {"H_2", "Z", "1-Z", "1-X"},
    {"H_4", "Z", "-Z", "1-X"},
    {"G", "-M", "M", "D"},
    {"G_2", "M", "1-M", "-D"},
    {"G_4", "M", "-M", "-D"},
    {"G_6", "1-M", "M", "D"},
};
inline const Segment mC2Path[] = {{"GAMMA", "Y"}, {"Y", "M"}, {"M", "A"}, {"A", "GAMMA"}, {"L_2", "GAMMA"}, {"GAMMA", "V_2"}};

inline const Parameter mC3Parameters[] = {{"Z", "(a*a/b/b+(1+a/c*cosbeta)/sinbeta/sinbeta)/4"}, {"R", "1-Z*b*b/a/a"}, {"E", "1/2-2*Z*c*cosbeta/a"}, {"F", "E/2+a*a/4/b/b+a*c*cosbeta/2/b/b"}, {"U", "2*F-Z"}, {"W", "c/2/a/cosbeta*(1-4*U+a*a*sinbeta*sinbeta/b/b)"}, {"D", "-1/4+W/2-Z*c*cosbeta/a"} };
inline const Point mC3Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"Y", "1/2", "1/2", "0"},
    {"A", "0", "0", "1/2"},
    {"M_2", "-1/2", "1/2", "1/2"},
    {"V", "1/2", "0", "0"},
    {"V_2", "0", "1/2", "0"},
    {"L_2", "0", "1/2", "1/2"},
    {"I", "-1+R", "R", "1/2"},
    {"I_2", "1-R", "1-R", "1/2"},
    {"K", "-U", "U", "W"},
    {"K_2", "-1+U", "1-U", "1-W"},
    {"K_4", "1-U", "U", "W"},
    {"H", "-Z", "Z", "E"},
    {"H_2", "Z", "1-Z", "1-E"},
    {"H_4", "Z", "-Z", "1-E"},
    {"N", "-F", "F", "D"},
    {"N_2", "F", "1-F", "-D"},
    {"N_4", "F", "-F", "-D"},
    {"N_6", "1-F", "F", "D"},
};
inline const Segment mC3Path[] = {{"GAMMA", "A"}, {"A", "I_2"}, {"I", "M_2"}, {"M_2", "GAMMA"}, {"GAMMA", "Y"}, {"L_2", "GAMMA"}, {"GAMMA", "V_2"}};

inline const Parameter mP1Parameters[] = {{"Y", "(1+a/c*cosbeta)/2/sinbeta/sinbeta"}, {"N", "1/2+Y*c*cosbeta/a"} };
inline const Point mP1Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"Z", "0", "1/2", "0"},
    {"B", "0", "0", "1/2"},
    {"B_2", "0", "0", "-1/2"},
    {"Y", "1/2", "0", "0"},
    {"Y_2", "-1/2", "0", "0"},
    {"C", "1/2", "1/2", "0"},
    {"C_2", "-1/2", "1/2", "0"},
    {"D", "0", "1/2", "1/2"},
    {"D_2", "0", "1/2", "-1/2"},
    {"A", "-1/2", "0", "1/2"},
    {"E", "-1/2", "1/2", "1/2"},
    {"H", "-Y", "0", "1-N"},
    {"H_2", "-1+Y", "0", "N"},
    {"H_4", "-Y", "0", "-N"},
    {"M", "-Y", "1/2", "1-N"},
    {"M_2", "-1+Y", "1/2", "N"},
    {"M_4", "-Y", "1/2", "-N"},
};
inline const Segment mP1Path[] = {{"GAMMA", "Z"}, {"Z", "D"}, {"D", "B"}, {"B", "GAMMA"}, {"GAMMA", "A"}, {"A", "E"}, {"E", "Z"}, {"Z", "C_2"}, {"C_2", "Y_2"}, {"Y_2", "GAMMA"}};

inline const Parameter oA1Parameters[] = {{"X", "(1+b*b/c/c)/4"} };
inline const Point oA1Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"Y", "-1/2", "1/2", "0"},
    {"T", "-1/2", "1/2", "1/2"},
    {"Z", "0", "0", "1/2"},
    {"S", "0", "1/2", "0"},
    {"R", "0", "1/2", "1/2"},
    {"SIGMA_0", "X", "X", "0"},
    {"C_0", "-X", "1-X", "0"},
    {"A_0", "X", "X", "1/2"},
    {"E_0", "-X", "1-X", "1/2"},
};
inline const Segment oA1Path[] = {{"GAMMA", "Y"}, {"Y", "C_0"}, {"SIGMA_0", "GAMMA"}, {"GAMMA", "Z"}, {"Z", "A_0"}, {"E_0", "T"}, {"T", "Y"}, {"GAMMA", "S"}, {"S", "R"}, {"R", "Z"}, {"Z", "T"}};

inline const Parameter oA2Parameters[] = {{"X", "(1+c*c/b/b)/4"} };
inline const Point oA2Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"Y", "1/2", "1/2", "0"},
    {"T", "1/2", "1/2", "1/2"},
    {"T_2", "1/2", "1/2", "-1/2"},
    {"Z", "0", "0", "1/2"},
    {"Z_2", "0", "0", "-1/2"},
    {"S", "0", "1/2", "0"},
    {"R", "0", "1/2", "1/2"},
    {"R_2", "0", "1/2", "-1/2"},
    {"DELTA_0", "-X", "X", "0"},
    {"F_0", "X", "1-X", "0"},
    {"B_0", "-X", "X", "1/2"},
    {"B_2", "-X", "X", "-1/2"},
    {"G_0", "X", "1-X", "1/2"},
    {"G_2", "X", "1-X", "-1/2"},
};
inline const Segment oA2Path[] = {{"GAMMA", "Y"}, {"Y", "F_0"}, {"DELTA_0", "GAMMA"}, {"GAMMA", "Z"}, {"Z", "B_0"}, {"G_0", "T"}, {"T", "Y"}, {"GAMMA", "S"}, {"S", "R"}, {"R", "Z"}, {"Z", "T"}};

inline const Parameter oC1Parameters[] = {{"X", "(1+a*a/b/b)/4"} };
inline const Point oC1Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"Y", "-1/2", "1/2", "0"},
    {"T", "-1/2", "1/2", "1/2"},
    {"Z", "0", "0", "1/2"},
    {"S", "0", "1/2", "0"},
    {"R", "0", "1/2", "1/2"},
    {"SIGMA_0", "X", "X", "0"},
    {"C_0", "-X", "1-X", "0"},
    {"A_0", "X", "X", "1/2"},
    {"E_0", "-X", "1-X", "1/2"},
};
inline const Segment oC1Path[] = {{"GAMMA", "Y"}, {"Y", "C_0"}, {"SIGMA_0", "GAMMA"}, {"GAMMA", "Z"}, {"Z", "A_0"}, {"E_0", "T"}, {"T", "Y"}, {"GAMMA", "S"}, {"S", "R"}, {"R", "Z"}, {"Z", "T"}};

inline const Parameter oC2Parameters[] = {{"X", "(1+b*b/a/a)/4"} };
inline const Point oC2Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"Y", "1/2", "1/2", "0"},
    {"T", "1/2", "1/2", "1/2"},
    {"T_2", "1/2", "1/2", "-1/2"},
    {"Z", "0", "0", "1/2"},
    {"Z_2", "0", "0", "-1/2"},
    {"S", "0", "1/2", "0"},
    {"R", "0", "1/2", "1/2"},
    {"R_2", "0", "1/2", "-1/2"},
    {"DELTA_0", "-X", "X", "0"},
    {"F_0", "X", "1-X", "0"},
    {"B_0", "-X", "X", "1/2"},
    {"B_2", "-X", "X", "-1/2"},
    {"G_0", "X", "1-X", "1/2"},
    {"G_2", "X", "1-X", "-1/2"},
};
inline const Segment oC2Path[] = {{"GAMMA", "Y"}, {"Y", "F_0"}, {"DELTA_0", "GAMMA"}, {"GAMMA", "Z"}, {"Z", "B_0"}, {"G_0", "T"}, {"T", "Y"}, {"GAMMA", "S"}, {"S", "R"}, {"R", "Z"}, {"Z", "T"}};

inline const Parameter oF1Parameters[] = {{"J", "(1+a*a/b/b-a*a/c/c)/4"}, {"H", "(1+a*a/b/b+a*a/c/c)/4"} };
inline const Point oF1Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"T", "1", "1/2", "1/2"},
    {"Z", "1/2", "1/2", "0"},
    {"Y", "1/2", "0", "1/2"},
    {"SIGMA_0", "0", "H", "H"},
    {"U_0", "1", "1-H", "1-H"},
    {"A_0", "1/2", "1/2+J", "J"},
    {"C_0", "1/2", "1/2-J", "1-J"},
    {"L", "1/2", "1/2", "1/2"},
};
inline const Segment oF1Path[] = {{"GAMMA", "Y"}, {"Y", "T"}, {"T", "Z"}, {"Z", "GAMMA"}, {"GAMMA", "SIGMA_0"}, {"U_0", "T"}, {"Y", "C_0"}, {"A_0", "Z"}, {"GAMMA", "L"}};

inline const Parameter oF2Parameters[] = {{"J", "(1+c*c/a/a-c*c/b/b)/4"}, {"K", "(1+c*c/a/a+c*c/b/b)/4"} };
inline const Point oF2Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"T", "0", "1/2", "1/2"},
    {"Z", "1/2", "1/2", "1"},
    {"Y", "1/2", "0", "1/2"},
    {"LAMBDA_0", "K", "K", "0"},
    {"Q_0", "1-K", "1-K", "1"},
    {"G_0", "1/2-J", "1-J", "1/2"},
    {"H_0", "1/2+J", "J", "1/2"},
    {"L", "1/2", "1/2", "1/2"},
};
inline const Segment oF2Path[] = {{"GAMMA", "T"}, {"T", "Z"}, {"Z", "Y"}, {"Y", "GAMMA"}, {"GAMMA", "LAMBDA_0"}, {"Q_0", "Z"}, {"T", "G_0"}, {"H_0", "Y"}, {"GAMMA", "L"}};

inline const Parameter oF3Parameters[] = {{"H", "(1+a*a/b/b-a*a/c/c)/4"}, {"K", "(1+b*b/a/a-b*b/c/c)/4"}, {"P", "(1+c*c/b/b-c*c/a/a)/4"} };
inline const Point oF3Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"T", "0", "1/2", "1/2"},
    {"Z", "1/2", "1/2", "0"},
    {"Y", "1/2", "0", "1/2"},
    {"A_0", "1/2", "1/2+H", "H"},
    {"C_0", "1/2", "1/2-H", "1-H"},
    {"B_0", "1/2+K", "1/2", "K"},
    {"D_0", "1/2-K", "1/2", "1-K"},
    {"G_0", "P", "1/2+P", "1/2"},
    {"H_0", "1-P", "1/2-P", "1/2"},
    {"L", "1/2", "1/2", "1/2"},
};
inline const Segment oF3Path[] = {{"GAMMA", "Y"}, {"Y", "C_0"}, {"A_0", "Z"}, {"Z", "B_0"}, {"D_0", "T"}, {"T", "G_0"}, {"H_0", "Y"}, {"T", "GAMMA"}, {"GAMMA", "Z"}, {"GAMMA", "L"}};

inline const Parameter oI1Parameters[] = {{"Z", "(1+a*a/c/c)/4"}, {"H", "(1+b*b/c/c)/4"}, {"D", "(b*b-a*a)/4/c/c"}, {"N", "(a*a+b*b)/4/c/c"} };
inline const Point oI1Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"X", "1/2", "1/2", "-1/2"},
    {"S", "1/2", "0", "0"},
    {"R", "0", "1/2", "0"},
    {"T", "0", "0", "1/2"},
    {"W", "1/4", "1/4", "1/4"},
    {"SIGMA_0", "-Z", "Z", "Z"},
    {"F_2", "Z", "1-Z", "-Z"},
    {"Y_0", "H", "-H", "H"},
    {"U_0", "1-H", "H", "-H"},
    {"L_0", "-N", "N", "1/2-D"},
    {"M_0", "N", "-N", "1/2+D"},
    {"J_0", "1/2-D", "1/2+D", "-N"},
};
inline const Segment oI1Path[] = {{"GAMMA", "X"}, {"X", "F_2"}, {"SIGMA_0", "GAMMA"}, {"GAMMA", "Y_0"}, {"U_0", "X"}, {"GAMMA", "R"}, {"R", "W"}, {"W", "S"}, {"S", "GAMMA"}, {"GAMMA", "T"}, {"T", "W"}};

inline const Parameter oI2Parameters[] = {{"Z", "(1+b*b/a/a)/4"}, {"H", "(1+c*c/a/a)/4"}, {"D", "(c*c-b*b)/4/a/a"}, {"N", "(b*b+c*c)/4/a/a"} };
inline const Point oI2Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"X", "-1/2", "1/2", "1/2"},
    {"S", "1/2", "0", "0"},
    {"R", "0", "1/2", "0"},
    {"T", "0", "0", "1/2"},
    {"W", "1/4", "1/4", "1/4"},
    {"Y_0", "Z", "-Z", "Z"},
    {"U_2", "-Z", "Z", "1-Z"},
    {"LAMBDA_0", "H", "H", "-H"},
    {"G_2", "-H", "1-H", "H"},
    {"K", "1/2-D", "-N", "N"},
    {"K_2", "1/2+D", "N", "-N"},
    {"K_4", "-N", "1/2-D", "1/2+D"},
};
inline const Segment oI2Path[] = {{"GAMMA", "X"}, {"X", "U_2"}, {"Y_0", "GAMMA"}, {"GAMMA", "LAMBDA_0"}, {"G_2", "X"}, {"GAMMA", "R"}, {"R", "W"}, {"W", "S"}, {"S", "GAMMA"}, {"GAMMA", "T"}, {"T", "W"}};

inline const Parameter oI3Parameters[] = {{"Z", "(1+c*c/b/b)/4"}, {"Y", "(1+a*a/b/b)/4"}, {"D", "(a*a-c*c)/4/b/b"}, {"M", "(c*c+a*a)/4/b/b"} };
inline const Point oI3Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"X", "1/2", "-1/2", "1/2"},
    {"S", "1/2", "0", "0"},
    {"R", "0", "1/2", "0"},
    {"T", "0", "0", "1/2"},
    {"W", "1/4", "1/4", "1/4"},
    {"SIGMA_0", "-Y", "Y", "Y"},
    {"F_0", "Y", "-Y", "1-Y"},
    {"LAMBDA_0", "Z", "Z", "-Z"},
    {"G_0", "1-Z", "-Z", "Z"},
    {"V_0", "M", "1/2-D", "-M"},
    {"H_0", "-M", "1/2+D", "M"},
    {"H_2", "1/2+D", "-M", "1/2-D"},
};
inline const Segment oI3Path[] = {{"GAMMA", "X"}, {"X", "F_0"}, {"SIGMA_0", "GAMMA"}, {"GAMMA", "LAMBDA_0"}, {"G_0", "X"}, {"GAMMA", "R"}, {"R", "W"}, {"W", "S"}, {"S", "GAMMA"}, {"GAMMA", "T"}, {"T", "W"}};

inline const Parameter oP1Parameters[] = {{nullptr, nullptr}};
inline const Point oP1Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"X", "1/2", "0", "0"},
    {"Z", "0", "0", "1/2"},
    {"U", "1/2", "0", "1/2"},
    {"Y", "0", "1/2", "0"},
    {"S", "1/2", "1/2", "0"},
    {"T", "0", "1/2", "1/2"},
    {"R", "1/2", "1/2", "1/2"},
};
inline const Segment oP1Path[] = {{"GAMMA", "X"}, {"X", "S"}, {"S", "Y"}, {"Y", "GAMMA"}, {"GAMMA", "Z"}, {"Z", "U"}, {"U", "R"}, {"R", "T"}, {"T", "Z"}, {"X", "U"}, {"Y", "T"}, {"S", "R"}};

inline const Parameter tI1Parameters[] = {{"H", "(1+c*c/a/a)/4"} };
inline const Point tI1Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"M", "-1/2", "1/2", "1/2"},
    {"X", "0", "0", "1/2"},
    {"P", "1/4", "1/4", "1/4"},
    {"Z", "H", "H", "-H"},
    {"Z_0", "-H", "1-H", "H"},
    {"N", "0", "1/2", "0"},
};
inline const Segment tI1Path[] = {{"GAMMA", "X"}, {"X", "M"}, {"M", "GAMMA"}, {"GAMMA", "Z"}, {"Z_0", "M"}, {"X", "P"}, {"P", "N"}, {"N", "GAMMA"}};

inline const Parameter tI2Parameters[] = {{"H", "(1+a*a/c/c)/4"}, {"Z", "a*a/2/c/c"} };
inline const Point tI2Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"M", "1/2", "1/2", "-1/2"},
    {"X", "0", "0", "1/2"},
    {"P", "1/4", "1/4", "1/4"},
    {"N", "0", "1/2", "0"},
    {"S_0", "-H", "H", "H"},
    {"S", "H", "1-H", "-H"},
    {"R", "-Z", "Z", "1/2"},
    {"G", "1/2", "1/2", "-Z"},
};
inline const Segment tI2Path[] = {{"GAMMA", "X"}, {"X", "P"}, {"P", "N"}, {"N", "GAMMA"}, {"GAMMA", "M"}, {"M", "S"}, {"S_0", "GAMMA"}, {"X", "R"}, {"G", "M"}};

inline const Parameter tP1Parameters[] = {{nullptr, nullptr}};
inline const Point tP1Points[] = {
    {"GAMMA", "0", "0", "0"},
    {"Z", "0", "0", "1/2"},
    {"M", "1/2", "1/2", "0"},
    {"A", "1/2", "1/2", "1/2"},
    {"R", "0", "1/2", "1/2"},
    {"X", "0", "1/2", "0"},
};
inline const Segment tP1Path[] = {{"GAMMA", "X"}, {"X", "M"}, {"M", "GAMMA"}, {"GAMMA", "Z"}, {"Z", "R"}, {"R", "A"}, {"A", "Z"}, {"X", "R"}, {"M", "A"}};

inline const Lattice lattices[] = {
    {"aP2", aP2Parameters, 0, aP2Points, 8, aP2Path, 7},
    {"aP3", aP3Parameters, 0, aP3Points, 9, aP3Path, 7},
    {"cF1", cF1Parameters, 0, cF1Points, 7, cF1Path, 7},
    {"cF2", cF2Parameters, 0, cF2Points, 7, cF2Path, 6},
    {"cI1", cI1Parameters, 0, cI1Points, 4, cI1Path, 6},
    {"cP1", cP1Parameters, 0, cP1Points, 5, cP1Path, 7},
    {"cP2", cP2Parameters, 0, cP2Points, 5, cP2Path, 6},
    {"hP1", hP1Parameters, 0, hP1Points, 7, hP1Path, 10},
    {"hP2", hP2Parameters, 0, hP2Points, 7, hP2Path, 9},
    {"hR1", hR1Parameters, 3, hR1Points, 20, hR1Path, 7},
    {"hR2", hR2Parameters, 3, hR2Points, 9, hR2Path, 5},
    {"mC1", mC1Parameters, 4, mC1Points, 16, mC1Path, 9},
    {"mC2", mC2Parameters, 6, mC2Points, 16, mC2Path, 6},
    {"mC3", mC3Parameters, 7, mC3Points, 19, mC3Path, 7},
    {"mP1", mP1Parameters, 2, mP1Points, 18, mP1Path, 10},
    {"oA1", oA1Parameters, 1, oA1Points, 10, oA1Path, 11},
    {"oA2", oA2Parameters, 1, oA2Points, 15, oA2Path, 11},
    {"oC1", oC1Parameters, 1, oC1Points, 10, oC1Path, 11},
    {"oC2", oC2Parameters, 1, oC2Points, 15, oC2Path, 11},
    {"oF1", oF1Parameters, 2, oF1Points, 9, oF1Path, 9},
    {"oF2", oF2Parameters, 2, oF2Points, 9, oF2Path, 9},
    {"oF3", oF3Parameters, 3, oF3Points, 11, oF3Path, 10},
    {"oI1", oI1Parameters, 4, oI1Points, 13, oI1Path, 11},
    {"oI2", oI2Parameters, 4, oI2Points, 13, oI2Path, 11},
    {"oI3", oI3Parameters, 4, oI3Points, 13, oI3Path, 11},
    {"oP1", oP1Parameters, 0, oP1Points, 8, oP1Path, 12},
    {"tI1", tI1Parameters, 1, tI1Points, 7, tI1Path, 8},
    {"tI2", tI2Parameters, 2, tI2Points, 9, tI2Path, 9},
    {"tP1", tP1Parameters, 0, tP1Points, 6, tP1Path, 9},
};
}
