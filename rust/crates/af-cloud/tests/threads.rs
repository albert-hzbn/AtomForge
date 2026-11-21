//! A multi-threaded build lays chunks out in completion order but holds the
//! same chunks (header order, counts, boxes) and the same atoms as a
//! single-threaded one.

use af_cloud::{generate_crystal, BuildOptions, CloudFile};

#[test]
fn threaded_build_equals_single_threaded() {
    let folder = std::env::temp_dir().join(format!("af-cloud-threads-{}", std::process::id()));
    std::fs::create_dir_all(&folder).unwrap();
    let build = |threads: u32| {
        let path = folder.join(format!("fcc-{threads}.afcloud"));
        let options = BuildOptions { chunk_atoms: 1024, threads, seed: 5, temporary: None };
        generate_crystal(&path, "fcc", 3.6, 0.0, [12, 11, 10], "Cu", "Ni", None, options).unwrap();
        CloudFile::open(&path).unwrap()
    };
    let (one, many) = (build(1), build(4));
    let (a, b) = (one.info(), many.info());
    assert_eq!(a.atoms, 5280);
    assert_eq!((a.atoms, a.lower, a.upper, &a.species, &a.source), (b.atoms, b.lower, b.upper, &b.species, &b.source));
    assert_eq!(a.chunks.len(), b.chunks.len());
    for (c, (x, y)) in a.chunks.iter().zip(&b.chunks).enumerate() {
        assert_eq!((x.count, x.lower, x.upper), (y.count, y.lower, y.upper), "chunk {c}");
        assert_eq!(one.read(c, 0, x.count).unwrap(), many.read(c, 0, y.count).unwrap(), "atoms of chunk {c}");
    }
    let _ = std::fs::remove_dir_all(&folder);
}
