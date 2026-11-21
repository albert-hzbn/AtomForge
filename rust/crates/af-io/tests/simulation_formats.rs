//! Port of tests/simulation_formats.cpp.

use af_core::{AtomSite, Structure};
use af_io::simulation::{read_simulation_format, write_simulation_format};

fn site(symbol: &str, z: i32, p: [f64; 3], fixed: u8) -> AtomSite {
    AtomSite { symbol: symbol.into(), atomic_number: z, x: p[0], y: p[1], z: p[2], fixed, ..Default::default() }
}

fn sample() -> Structure {
    let mut s = Structure::default();
    s.has_unit_cell = true;
    s.cell_vectors = [[4.1, 0.0, 0.0], [1.3, 3.9, 0.0], [0.7, -0.4, 6.2]];
    s.atoms = vec![
        site("Fe", 26, [0.1, 0.2, 0.3], 0b101),
        site("O", 8, [1.234567891, 2.5, 3.75], 0),
        site("Fe", 26, [3.0, 1.0, 5.0], 0b010),
        site("Mo", 42, [-0.5, 0.25, 1.125], 0),
    ];
    s
}

fn gap(a: &Structure, b: &Structure) -> f64 {
    let mut worst: f64 = 0.0;
    for (x, y) in a.atoms.iter().zip(&b.atoms) {
        worst = worst.max((x.x - y.x).abs()).max((x.y - y.y).abs()).max((x.z - y.z).abs());
    }
    for r in 0..3 {
        for c in 0..3 {
            worst = worst.max((a.cell_vectors[r][c] - b.cell_vectors[r][c]).abs());
        }
    }
    worst
}

fn round_trip(s: &Structure, format: &str) -> Structure {
    let text = write_simulation_format(s, format).unwrap_or_else(|e| panic!("{format} write: {e}"));
    read_simulation_format(&text, format).unwrap_or_else(|e| panic!("{format} read: {e}\n{text}"))
}

#[test]
fn round_trips() {
    let s = sample();
    for format in ["pwi", "cell", "gen", "xsf"] {
        let back = round_trip(&s, format);
        assert_eq!(back.atoms.len(), s.atoms.len(), "{format}: atom count");
        assert!(back.has_unit_cell, "{format}: cell kept");
        assert!(gap(&s, &back) < 1e-9, "{format}: positions and cell within 1e-9 A");
        for (a, b) in s.atoms.iter().zip(&back.atoms) {
            assert!(a.symbol == b.symbol && a.atomic_number == b.atomic_number, "{format}: species");
            if format == "pwi" || format == "cell" {
                assert_eq!(a.fixed, b.fixed, "{format}: fixed coordinates");
            }
        }
    }
}

#[test]
fn magnetic_moments() {
    let mut magnetic = sample();
    magnetic.atoms[0].magmom = 2.5;
    magnetic.atoms[2].magmom = -2.5;
    let back = round_trip(&magnetic, "cell");
    assert!(back.atoms[0].magmom == 2.5 && back.atoms[2].magmom == -2.5 && back.atoms[1].magmom == 0.0, "cell SPIN= round trip");
    let spaced = read_simulation_format(
        "%BLOCK LATTICE_CART\n3 0 0\n0 3 0\n0 0 3\n%ENDBLOCK LATTICE_CART\n%BLOCK POSITIONS_ABS\nNi 0 0 0 SPIN = -1.5\n%ENDBLOCK POSITIONS_ABS\n",
        "cell",
    )
    .unwrap();
    assert_eq!(spaced.atoms[0].magmom, -1.5, "cell SPIN = x");
    let text = write_simulation_format(&magnetic, "pwi").unwrap();
    assert!(text.contains("nspin = 2") && text.contains("ntyp = 4"), "pwi species per moment: {text}");
    assert!(text.contains("starting_magnetization(1) = 0.5") && text.contains("starting_magnetization(3) = -0.5"), "pwi starting magnetization");
    assert!(text.contains("Fe1 ") && text.contains("Fe2 "), "pwi numbered species labels");
    let parsed = read_simulation_format(&text, "pwi").unwrap();
    assert_eq!(parsed.atoms[2].symbol, "Fe", "pwi labels read as elements");
}

#[test]
fn periodicity_molecules_and_cp2k() {
    let mut slab = sample();
    slab.pbc = [true, true, false];
    assert_eq!(round_trip(&slab, "xsf").pbc, slab.pbc, "xsf slab periodicity");
    let mut wire = sample();
    wire.pbc = [true, false, false];
    assert_eq!(round_trip(&wire, "xsf").pbc, wire.pbc, "xsf polymer periodicity");
    let mut molecule = Structure::default();
    molecule.atoms = vec![site("O", 8, [0.0, 0.0, 0.1197], 0), site("H", 1, [0.0, 0.7572, -0.4786], 0), site("H", 1, [0.0, -0.7572, -0.4786], 0)];
    for format in ["gen", "xsf"] {
        let back = round_trip(&molecule, format);
        assert!(!back.has_unit_cell && back.atoms.len() == 3 && (back.atoms[1].y - 0.7572).abs() < 1e-9, "{format}: molecule");
    }
    for format in ["pwi", "cell", "cp2k"] {
        let error = write_simulation_format(&molecule, format).unwrap_err();
        assert!(error.contains("unit cell"), "{format}: needs a cell");
    }
    let text = write_simulation_format(&slab, "cp2k").unwrap();
    assert!(text.contains("PERIODIC XY\n") && text.contains("&KIND Mo") && text.contains("  &COORD"), "cp2k sections");
}

#[test]
fn pwi_cell_gen_xsf_details() {
    let text = "&SYSTEM\n ibrav = 0, celldm(1) = 6.8306, nat = 2, ntyp = 1\n/\n\
                CELL_PARAMETERS alat\n 0.5 0.5 0.0\n 0.0 0.5 0.5\n 0.5 0.0 0.5\n\
                ATOMIC_POSITIONS {crystal}\n Cu 0.0 0.0 0.0\n Cu 0.25d0 0.25d0 0.25d0 0 1 1\n";
    let cu = read_simulation_format(text, "pwi").unwrap();
    let alat = 6.8306 * 0.529177210903;
    assert!((cu.cell_vectors[0][0] - 0.5 * alat).abs() < 1e-9, "pwi alat cell");
    assert!((cu.atoms[1].x - 0.25 * (0.5 + 0.0 + 0.5) * alat).abs() < 1e-9, "pwi crystal coordinates");
    assert_eq!(cu.atoms[1].fixed, 1, "pwi if_pos 0 fixes x");
    let error = read_simulation_format("&SYSTEM\n ibrav = 2\n/\n", "pwi").unwrap_err();
    assert!(error.contains("ibrav"), "pwi refuses other ibrav");

    let text = "%block lattice_abc\nang\n3 4 5\n90 90 120\n%endblock lattice_abc\n%BLOCK POSITIONS_FRAC\nSi 0 0 0\nSi:tag 0.5 0.5 0.5\n%ENDBLOCK POSITIONS_FRAC\n";
    let si = read_simulation_format(text, "cell").unwrap();
    assert!((si.cell_vectors[1][0] + 2.0).abs() < 1e-9 && (si.cell_vectors[1][1] - 4.0 * 3f64.sqrt() / 2.0).abs() < 1e-9, "cell gamma 120");
    assert!((si.atoms[1].z - 2.5).abs() < 1e-9 && si.atoms[1].symbol == "Si", "cell fractional and labels");
    let bohr = "%BLOCK LATTICE_CART\nbohr\n10 0 0\n0 10 0\n0 0 10\n%ENDBLOCK LATTICE_CART\n%BLOCK POSITIONS_ABS\nbohr\nC 1 0 0\n%ENDBLOCK POSITIONS_ABS\n";
    let si = read_simulation_format(bohr, "cell").unwrap();
    assert!((si.cell_vectors[0][0] - 5.29177210903).abs() < 1e-9 && (si.atoms[0].x - 0.529177210903).abs() < 1e-12, "cell bohr units");

    let text = "2 F\n Ga As\n 1 1 0.0 0.0 0.0\n 2 2 0.25 0.25 0.25\n 0.0 0.0 0.0\n 0.0 2.8 2.8\n 2.8 0.0 2.8\n 2.8 2.8 0.0\n";
    let gaas = read_simulation_format(text, "gen").unwrap();
    assert!((gaas.atoms[1].x - 1.4).abs() < 1e-12 && gaas.atoms[1].symbol == "As", "gen fractional");
    assert!(read_simulation_format("1 C\n Cu\n 1 3 0 0 0\n", "gen").is_err(), "gen species index checked");

    let x = read_simulation_format("ATOMS\n 6 0 0 0\n O 0 0 1.128\n", "xsf").unwrap();
    assert!(x.atoms.len() == 2 && x.atoms[0].symbol == "C" && x.atoms[1].symbol == "O", "xsf atomic numbers and symbols");
    assert!(read_simulation_format("CRYSTAL\nPRIMVEC\n1 0 0\n", "xsf").is_err(), "truncated xsf reports");
}
