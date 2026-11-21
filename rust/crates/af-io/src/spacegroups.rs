//! spglib's space-group database (530 Hall settings) as exported by
//! tools/golden/spglib_table.cpp: names, setting choices and symmetry
//! operations in spglib's order, for reading CIF files that give a space
//! group without listing its operations.

use std::sync::OnceLock;

use glam::{DMat3, DVec3};

pub struct SpaceGroupType {
    pub hall_number: i32,
    pub number: i32,
    pub hall_symbol: String,
    pub international_short: String,
    pub international_full: String,
    pub international: String,
    pub choice: String,
    /// Rotation (rows as read, applied to column vectors) and translation.
    pub operations: Vec<(DMat3, DVec3)>,
}

pub fn types() -> &'static [SpaceGroupType] {
    static TABLE: OnceLock<Vec<SpaceGroupType>> = OnceLock::new();
    TABLE.get_or_init(|| {
        include_str!("../data/spglib_hall.txt")
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
            .map(|line| {
                let f: Vec<&str> = line.splitn(8, '|').collect();
                let operations = f[7]
                    .split(';')
                    .filter(|op| !op.trim().is_empty())
                    .map(|op| {
                        let v: Vec<f64> = op.split_whitespace().map(|x| x.parse().unwrap()).collect();
                        // r[row][col]: column-major from_cols takes columns.
                        let rotation = DMat3::from_cols(
                            DVec3::new(v[0], v[3], v[6]),
                            DVec3::new(v[1], v[4], v[7]),
                            DVec3::new(v[2], v[5], v[8]),
                        );
                        (rotation, DVec3::new(v[9], v[10], v[11]))
                    })
                    .collect();
                SpaceGroupType {
                    hall_number: f[0].parse().unwrap(),
                    number: f[1].parse().unwrap(),
                    hall_symbol: f[2].to_string(),
                    international_short: f[3].to_string(),
                    international_full: f[4].to_string(),
                    international: f[5].to_string(),
                    choice: f[6].to_string(),
                    operations,
                }
            })
            .collect()
    })
}

/// The type of Hall number 1-530.
pub fn of_hall(hall: i32) -> &'static SpaceGroupType {
    &types()[(hall - 1) as usize]
}
