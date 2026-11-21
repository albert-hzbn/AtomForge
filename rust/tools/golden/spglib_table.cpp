// Exports spglib's space-group database (the 530 Hall settings: names, setting
// choice and symmetry operations in spglib's order) for af-io's CIF reader:
// one line per Hall number,
//   hall|number|hall_symbol|international_short|international_full|international|choice|r00 ... r22 t0 t1 t2;...
#include <spglib.h>

#include <cstdio>

int main(int argc, char** argv)
{
    FILE* out = argc > 1 ? std::fopen(argv[1], "w") : stdout;
    std::fprintf(out, "# spglib %d.%d.%d space-group database (Hall settings 1-530), exported by tools/golden/spglib_table.cpp\n",
                 spg_get_major_version(), spg_get_minor_version(), spg_get_micro_version());
    static int rotations[192][3][3];
    static double translations[192][3];
    for (int hall = 1; hall <= 530; ++hall) {
        const SpglibSpacegroupType t = spg_get_spacegroup_type(hall);
        std::fprintf(out, "%d|%d|%s|%s|%s|%s|%s|", hall, t.number, t.hall_symbol, t.international_short, t.international_full, t.international, t.choice);
        const int count = spg_get_symmetry_from_database(rotations, translations, hall);
        for (int n = 0; n < count; ++n) {
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c) std::fprintf(out, "%d ", rotations[n][r][c]);
            std::fprintf(out, "%.17g %.17g %.17g;", translations[n][0], translations[n][1], translations[n][2]);
        }
        std::fprintf(out, "\n");
    }
    if (out != stdout) std::fclose(out);
    return 0;
}
