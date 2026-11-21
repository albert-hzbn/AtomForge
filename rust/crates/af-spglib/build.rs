//! Links the vendored spglib 2.7.0 (vendor/lib, MSVC build from the spglib
//! Python wheel) and copies symspg.dll next to the build outputs so tests,
//! examples and binaries find it at run time.

use std::path::PathBuf;

fn main() {
    let manifest = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let vendor = manifest.join("vendor").join("lib");
    println!("cargo:rerun-if-changed=vendor/lib/symspg.lib");
    println!("cargo:rerun-if-changed=vendor/lib/symspg.dll");
    let target = std::env::var("TARGET").unwrap();
    if !target.contains("windows-msvc") {
        // Elsewhere link the system spglib (libsymspg).
        println!("cargo:rustc-link-lib=dylib=symspg");
        return;
    }
    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    std::fs::copy(vendor.join("symspg.lib"), out.join("symspg.lib")).expect("copy symspg.lib");
    std::fs::copy(vendor.join("symspg.dll"), out.join("symspg.dll")).expect("copy symspg.dll");
    // OUT_DIR is inside target/, so cargo also puts it on PATH for `cargo test` / `cargo run`.
    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=dylib=symspg");
    // Also next to the final executables (target/<profile>/ and its deps/).
    if let Some(profile_dir) = out.ancestors().nth(3) {
        let _ = std::fs::copy(vendor.join("symspg.dll"), profile_dir.join("symspg.dll"));
        let _ = std::fs::copy(vendor.join("symspg.dll"), profile_dir.join("deps").join("symspg.dll"));
    }
}
