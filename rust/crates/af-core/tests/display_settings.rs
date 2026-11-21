//! Port of tests/display_settings.cpp.

use af_core::display::*;
use af_core::elements::{colour_schemes, default_colours, matching_scheme, scheme_colours};
use af_core::occlusion::{ambient_occlusion, darken_occluded, occlusion_directions};
use glam::Vec3;

fn close(a: [f32; 3], b: [f32; 3], tolerance: f32) -> bool {
    (Vec3::from(a) - Vec3::from(b)).length() <= tolerance
}

#[test]
fn light_direction_and_presets() {
    let defaults = DisplaySettings::default();
    assert!(close(light_direction(&defaults), Vec3::ONE.normalize().to_array(), 1e-5), "default light along (1, 1, 1)");
    let mut s = DisplaySettings::default();
    s.light_azimuth = 0.0;
    s.light_elevation = 0.0;
    assert!(close(light_direction(&s), [1.0, 0.0, 0.0], 1e-5), "azimuth 0 lights along +x");
    s.light_azimuth = 90.0;
    assert!(close(light_direction(&s), [0.0, 0.0, 1.0], 1e-5), "azimuth 90 lights along +z");
    s.light_elevation = 89.0;
    assert!((light_direction(&s)[1] - (85.0f32 * 3.14159265 / 180.0).sin()).abs() < 1e-5, "elevation clamped to 85");
    s.light_follows_camera = true;
    s.light_azimuth = 0.0;
    s.light_elevation = 0.0;
    assert!(close(light_direction(&s), [0.0, 0.0, 1.0], 1e-5), "a headlight at azimuth 0 comes from the viewer");
    s.light_azimuth = 90.0;
    assert!(close(light_direction(&s), [1.0, 0.0, 0.0], 1e-5), "a headlight at azimuth 90 comes from the right");

    assert_eq!(matching_preset(&defaults), 0, "the defaults match the Default preset");
    for (i, p) in PRESETS.iter().enumerate() {
        let mut look = DisplaySettings::default();
        look.atom_scale = 1.4;
        apply_preset(&mut look, p);
        assert_eq!(matching_preset(&look), i as i32, "preset recognised: {}", p.name);
        assert_eq!(look.atom_scale, 1.4, "a preset keeps the atom size");
    }
    let mut publication = DisplaySettings::default();
    apply_preset(&mut publication, PRESETS.last().unwrap());
    assert!(publication.custom_background && publication.background == [1.0; 3] && publication.custom_cell_colour && publication.cell_colour == [0.0; 3]);
    apply_preset(&mut publication, &PRESETS[1]);
    assert!(!publication.custom_background && !publication.custom_cell_colour, "other presets return to the theme's colours");
    let mut customised = DisplaySettings::default();
    customised.ambient = 0.41;
    assert_eq!(matching_preset(&customised), -1, "a changed slider is a custom look");
    customised = DisplaySettings::default();
    customised.custom_background = true;
    assert_eq!(matching_preset(&customised), -1, "a custom background is not the Default look");
}

fn all_settings() -> DisplaySettings {
    DisplaySettings {
        ambient: 0.123456,
        shadow_strength: 0.9,
        contrast: 2.1,
        saturation: 0.3,
        light_azimuth: 271.5,
        light_elevation: -20.25,
        light_follows_camera: true,
        specular_intensity: 1.7,
        shininess_scale: 0.2,
        shininess_floor: 200.0,
        atom_scale: 0.35,
        bond_scale: 2.5,
        uniform_bond_colour: true,
        bond_colour: [0.1, 0.2, 0.3],
        custom_background: true,
        background: [0.0, 0.0, 0.0],
        custom_cell_colour: true,
        cell_colour: [0.9, 0.5, 0.25],
        cell_line_width: 6.0,
        depth_cueing: 0.7,
        ambient_occlusion: 0.4,
        outline_width: 0.25,
        outline_colour: [0.2, 0.1, 0.0],
        selection_colour: [0.0, 0.8, 1.0],
        bond_tolerance: 1.3,
        metal_bonds: true,
        bond_atom_limit: 1200000.0,
    }
}

#[test]
fn text_round_trips_and_reading() {
    let defaults = DisplaySettings::default();
    let all = all_settings();
    assert_eq!(from_text(&to_text(&all), &defaults), all, "every setting round-trips through text");
    assert_eq!(from_text(&to_text(&defaults), &defaults), defaults, "the defaults round-trip");
    assert!(to_text(&defaults).contains("background=theme"));
    let read = from_text(
        "# comment\nambient=7\natom_scale=-1\nlight_azimuth=-30\nlight_elevation=120\ncontrast=nonsense\nunknown_key=3\n\
         bond_colour=1 0 2\nbackground = 0.5 0.5\n cell_colour = theme \n",
        &defaults,
    );
    assert_eq!(read.ambient, AMBIENT.high, "ambient clamped");
    assert_eq!(read.atom_scale, ATOM_SCALE.low, "atom size clamped");
    assert!((read.light_azimuth - 330.0).abs() < 1e-4, "azimuth wrapped into 0-360");
    assert_eq!(read.light_elevation, ELEVATION.high, "elevation clamped");
    assert_eq!(read.contrast, defaults.contrast, "an unreadable number keeps the default");
    assert!(read.uniform_bond_colour && read.bond_colour == [1.0, 0.0, 1.0], "colour channels clamped to 0-1");
    assert!(!read.custom_background, "an incomplete colour is ignored");
    assert!(!read.custom_cell_colour, "'theme' keeps the theme's cell colour");
    let mut base = DisplaySettings::default();
    base.bond_scale = 2.0;
    assert_eq!(from_text("ambient=0.5\n", &base).bond_scale, 2.0, "missing keys keep the base value");
    assert_eq!(from_text("ambient=0.5\n", &base).ambient, 0.5, "present keys are read");

    let mut look = default_look();
    assert!(element_text(&look.colours, &look.radii).is_empty(), "the built-in tables write nothing");
    look.display = all;
    look.colours[29] = [0.1, 0.2, 0.3];
    look.radii[8] = 0.9;
    look.radii[118] = 2.25;
    let elements = element_text(&look.colours, &look.radii);
    assert!(elements.contains("colour.Cu=") && elements.contains("radius.O=") && elements.contains("radius.Og="));
    assert_eq!(elements.matches('\n').count(), 3, "unchanged elements not written");
    let back = read_look(&look_text(&look));
    assert!(back.display == all && back.colours == look.colours && back.radii == look.radii, "a look round-trips through text");
    assert_eq!(from_text(&look_text(&look), &defaults), all, "element lines do not disturb the settings");
    let mut colours = default_look().colours;
    let mut radii = default_look().radii;
    read_element_text("colour.Xx=1 1 1\nradius.Fe=nonsense\nradius.Fe=99\ncolour.H=0.5 0.5\n", &mut colours, &mut radii);
    assert_eq!(colours, default_look().colours, "unknown elements and incomplete colours are skipped");
    assert_eq!(radii[26], 5.0, "a radius is clamped");

    let extremes = from_text("cell_line_width=40\noutline_width=-1\nbond_tolerance=9\nbond_atom_limit=12.6\nmetal_bonds=yes\n", &defaults);
    assert!(
        extremes.cell_line_width == CELL_LINE_WIDTH.high
            && extremes.outline_width == 0.0
            && extremes.bond_tolerance == BOND_TOLERANCE.high
            && extremes.bond_atom_limit == 13.0
            && extremes.metal_bonds
    );
    assert!(
        defaults.depth_cueing == 0.0
            && defaults.ambient_occlusion == 0.0
            && defaults.outline_width == 0.0
            && defaults.bond_tolerance == 1.18
            && !defaults.metal_bonds
            && defaults.bond_atom_limit == 500000.0
    );

    let path = std::env::temp_dir().join("atomforge_rs_display_settings_test.txt");
    assert!(save_defaults(&look, &path), "settings saved");
    let loaded = load_defaults(&path);
    assert!(loaded.display == all && loaded.colours == look.colours && loaded.radii == look.radii, "saved look loads back");
    std::fs::remove_file(&path).unwrap();
    let missing = load_defaults(&path);
    assert!(missing.display == DisplaySettings::default() && missing.colours == default_look().colours, "a missing file gives the defaults");
}

#[test]
fn ambient_occlusion_behaviour() {
    assert_eq!(occlusion_directions().len(), 48);
    let mut sum = Vec3::ZERO;
    for d in occlusion_directions() {
        assert!((d.length() - 1.0).abs() < 1e-5, "unit directions");
        sum += *d;
    }
    assert!(sum.length() < 0.05 * 48.0, "directions spread evenly");
    assert_eq!(ambient_occlusion(&[Vec3::ZERO], &[1.0], 2.0, 2_000_000), vec![0.0], "a lone atom is open");
    assert_eq!(ambient_occlusion(&[Vec3::ZERO, Vec3::new(10.0, 0.0, 0.0)], &[1.0, 1.0], 2.0, 2_000_000), vec![0.0, 0.0], "out of reach");
    let mut centres = vec![Vec3::ZERO];
    for x in -1..=1 {
        for y in -1..=1 {
            for z in -1..=1 {
                if x != 0 || y != 0 || z != 0 {
                    centres.push(Vec3::new(1.6 * x as f32, 1.6 * y as f32, 1.6 * z as f32));
                }
            }
        }
    }
    let shell = ambient_occlusion(&centres, &vec![1.0; centres.len()], 2.0, 2_000_000);
    assert!(shell[0] > 0.9, "an atom in a closed shell is enclosed: {}", shell[0]);
    assert!(shell[1] < shell[0], "a shell atom sees more sky than the centre");

    let mut slab = Vec::new();
    for x in 0..8 {
        for y in 0..8 {
            for z in 0..5 {
                slab.push(Vec3::new(2.0 * x as f32, 2.0 * y as f32, 2.0 * z as f32));
            }
        }
    }
    let ones = vec![1.0; slab.len()];
    let occ = ambient_occlusion(&slab, &ones, 2.0, 2_000_000);
    let at = |x: usize, y: usize, z: usize| (x * 8 + y) * 5 + z;
    assert!(occ[at(4, 4, 0)] < occ[at(4, 4, 2)], "the surface is lighter than the middle");
    assert!((occ[at(4, 4, 0)] - occ[at(4, 4, 4)]).abs() < 0.1, "both surfaces alike");
    assert!(ambient_occlusion(&slab, &ones, 2.0, 100).is_empty(), "max_atoms skips large sets");
    let mut big = Vec::new();
    for x in 0..16 {
        for y in 0..16 {
            for z in 0..12 {
                big.push(Vec3::new(2.0 * x as f32, 2.0 * y as f32, 2.0 * z as f32));
            }
        }
    }
    let big_occ = ambient_occlusion(&big, &vec![1.0; big.len()], 2.0, 2_000_000);
    let one = ambient_occlusion(&[big[0], big[1], big[12], big[192]], &[1.0; 4], 2.0, 2_000_000);
    assert!(big_occ.len() == big.len() && big_occ[0] > 0.0, "threaded occlusion computed");
    assert!((one[0] - big_occ[0]).abs() < 0.25, "a corner atom's occlusion comes from its own neighbours");
    let mut colours = vec![[1.0f32; 3]; slab.len()];
    darken_occluded(&slab, &ones, &mut colours, 0.0);
    assert!(colours.iter().all(|c| *c == [1.0; 3]), "strength 0 changes nothing");
    darken_occluded(&slab, &ones, &mut colours, 1.0);
    assert!(colours[at(4, 4, 2)][0] < colours[at(4, 4, 0)][0] && colours[at(4, 4, 0)][0] < 1.0, "buried atoms darken most");
}

#[test]
fn colour_schemes_behaviour() {
    assert_eq!(colour_schemes().len(), 4);
    for i in 0..4 {
        let colours = scheme_colours(i);
        assert_eq!(colours.len(), 119);
        assert_eq!(matching_scheme(&colours), i);
        assert!(colours.iter().all(|c| c.iter().all(|&v| (0.0..=1.0).contains(&v))));
    }
    assert_eq!(scheme_colours(0), default_colours(), "the AtomForge scheme is the built-in table");
    assert!(close(scheme_colours(1)[6], [0.565, 0.565, 0.565], 1e-3), "Jmol carbon");
    assert!(close(scheme_colours(1)[29], [0.784, 0.502, 0.2], 1e-3), "Jmol copper");
    assert_eq!(scheme_colours(1)[118], default_colours()[118], "elements past the Jmol table keep the AtomForge colour");
    assert!(close(scheme_colours(2)[8], [0.941, 0.0, 0.0], 1e-3), "RasMol oxygen");
    let grey = scheme_colours(3);
    assert!(grey[29][0] == grey[29][1] && grey[29][1] == grey[29][2] && grey[1][0] <= 0.9501 && grey[6][0] >= 0.1499);
    let mut edited = scheme_colours(1);
    edited[26] = [0.0, 0.0, 0.0];
    assert_eq!(matching_scheme(&edited), -1, "an edited scheme is no scheme");
}
