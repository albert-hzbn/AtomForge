//! Windows: embeds the icon, the application manifest (no elevation prompt,
//! per-monitor DPI awareness) and the version information in AtomForge.exe,
//! as the C++ src/app/AtomForge.rc did, with the Windows SDK's rc.exe (MSVC)
//! or windres (GNU). Linux: the program finds its libraries (spglib) beside it.

use std::path::{Path, PathBuf};

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=res/AtomForge.manifest");
    let target = std::env::var("TARGET").unwrap_or_default();
    if target.contains("linux") {
        println!("cargo:rustc-link-arg-bins=-Wl,-rpath,$ORIGIN");
    }
    if !target.contains("windows") {
        return;
    }
    let manifest_dir = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let icon = manifest_dir.join("../../../assets/icon/atomforge.ico");
    println!("cargo:rerun-if-changed={}", icon.display());
    let version = std::env::var("CARGO_PKG_VERSION").unwrap_or_else(|_| "0.0.0".into());
    let numbers: Vec<u32> = version.split('.').map(|p| p.parse().unwrap_or(0)).chain(std::iter::repeat(0)).take(4).collect();
    let comma = numbers.iter().map(|n| n.to_string()).collect::<Vec<_>>().join(",");
    let dotted = numbers.iter().map(|n| n.to_string()).collect::<Vec<_>>().join(".");
    // No #include: constants spelt out, so rc.exe needs no SDK include path.
    let rc = format!(
        r#"1 ICON "{icon}"
1 24 "{manifest}"
1 VERSIONINFO
FILEVERSION {comma}
PRODUCTVERSION {comma}
FILEFLAGSMASK 0x3fL
FILEFLAGS 0x0L
FILEOS 0x4L
FILETYPE 0x1L
FILESUBTYPE 0x0L
BEGIN
    BLOCK "StringFileInfo"
    BEGIN
        BLOCK "040904B0"
        BEGIN
            VALUE "CompanyName", "Albert Linda"
            VALUE "FileDescription", "AtomForge - Atomistic structure builder and viewer"
            VALUE "FileVersion", "{dotted}"
            VALUE "InternalName", "AtomForge"
            VALUE "LegalCopyright", "Copyright (C) 2024-2026 Albert Linda"
            VALUE "OriginalFilename", "AtomForge.exe"
            VALUE "ProductName", "AtomForge"
            VALUE "ProductVersion", "{version}"
        END
    END
    BLOCK "VarFileInfo"
    BEGIN
        VALUE "Translation", 0x0409, 1200
    END
END
"#,
        icon = escape(&icon),
        manifest = escape(&manifest_dir.join("res/AtomForge.manifest")),
    );
    let out = PathBuf::from(std::env::var("OUT_DIR").unwrap());
    let rc_path = out.join("AtomForge.rc");
    std::fs::write(&rc_path, rc).expect("write AtomForge.rc");
    if target.contains("msvc") {
        let res = out.join("AtomForge.res");
        let Some(rc_exe) = find_rc() else {
            println!("cargo:warning=rc.exe not found (set RC to it): AtomForge.exe is built without its icon and manifest");
            return;
        };
        let status = std::process::Command::new(rc_exe).arg("/nologo").arg("/fo").arg(&res).arg(&rc_path).status();
        if !matches!(status, Ok(s) if s.success()) {
            println!("cargo:warning=rc.exe failed: AtomForge.exe is built without its icon and manifest");
            return;
        }
        println!("cargo:rustc-link-arg-bins={}", res.display());
    } else {
        let object = out.join("AtomForge-res.o");
        let windres = std::env::var("WINDRES").unwrap_or_else(|_| "windres".into());
        let status = std::process::Command::new(windres).arg(&rc_path).arg("-O").arg("coff").arg("-o").arg(&object).status();
        if !matches!(status, Ok(s) if s.success()) {
            println!("cargo:warning=windres failed: AtomForge.exe is built without its icon and manifest");
            return;
        }
        println!("cargo:rustc-link-arg-bins={}", object.display());
    }
}

fn escape(p: &Path) -> String {
    p.canonicalize().unwrap_or_else(|_| p.to_path_buf()).to_string_lossy().trim_start_matches(r"\\?\").replace('\\', "\\\\")
}

/// RC from the environment, else the newest Windows 10/11 SDK's x64 rc.exe.
fn find_rc() -> Option<PathBuf> {
    if let Ok(rc) = std::env::var("RC") {
        return Some(PathBuf::from(rc));
    }
    let base = PathBuf::from(std::env::var("ProgramFiles(x86)").unwrap_or_else(|_| r"C:\Program Files (x86)".into())).join(r"Windows Kits\10\bin");
    let mut versions: Vec<PathBuf> = std::fs::read_dir(&base).ok()?.filter_map(|e| e.ok()).map(|e| e.path().join(r"x64\rc.exe")).filter(|p| p.is_file()).collect();
    versions.sort();
    versions.pop()
}
