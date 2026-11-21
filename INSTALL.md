# AtomForge – Installation and Build Guide

AtomForge is written in Rust. One program, `AtomForge`, is both the desktop
application and, when given a mode such as `--build`, `--pipe` or `--science`,
the command line. The Python package (`python/atomforge`) uses two native
libraries built from the same sources: `atomforge_science_native` and
`atomforge_electronic`.

Prebuilt packages for Windows and Linux are on the
[releases page](https://github.com/albert-hzbn/AtomForge/releases/latest); this
guide is for building from source.

## Requirements

| Dependency | Version | Notes |
| --- | --- | --- |
| Rust toolchain | stable, 1.85 or later | Install with [rustup](https://rustup.rs/) |
| C/C++ linker | MSVC Build Tools (Windows) / GCC or Clang (Linux) | Rust links through the platform's toolchain |
| spglib | 2.x | Symmetry. Windows (MSVC): bundled in `rust/crates/af-spglib/vendor`. Linux: the system's `libsymspg` |
| GPU driver | Vulkan, Direct3D 12, Metal or OpenGL ES 3 | The desktop draws through [wgpu](https://wgpu.rs/) |
| Python | 3.8 or later | Only to run the packaging script and the Python package |

Every other dependency (egui/eframe, wgpu, glam, the file formats) is a Rust
crate that Cargo downloads on the first build. Structure files are read and
written by AtomForge's own readers and writers; Open Babel is no longer used.

---

## Windows

1. Install the [Visual Studio Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/)
   with the "Desktop development with C++" workload (MSVC and the Windows SDK).
2. Install Rust with `rustup-init.exe` (the default `x86_64-pc-windows-msvc` toolchain).
3. Build from a terminal in the repository:

```bash
cd rust
cargo build --release -p af-app -p af-ffi-science -p af-ffi-electronic
```

The program is `rust/target/release/AtomForge.exe`; `symspg.dll` is copied next to it.
The Windows SDK's `rc.exe` embeds the icon, manifest and version information
(set `RC` to its path if it is not found under `C:\Program Files (x86)\Windows Kits\10\bin`).

### Run

```bash
rust/target/release/AtomForge.exe
rust/target/release/AtomForge.exe structure.cif
rust/target/release/AtomForge.exe --build bulk --help
```

---

## Linux

### Install dependencies

```bash
# Debian / Ubuntu
sudo apt update
sudo apt install build-essential pkg-config curl libsymspg-dev \
                 libx11-dev libxkbcommon-dev libwayland-dev libvulkan1 mesa-vulkan-drivers
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### Build

```bash
cd rust
cargo build --release -p af-app -p af-ffi-science -p af-ffi-electronic
```

### Run

```bash
rust/target/release/AtomForge
rust/target/release/AtomForge structure.cif
```

---

## Tests

```bash
cd rust
cargo test --workspace --release          # unit and golden tests
cd ..
python rust/tools/run_ui_tests.py         # scripted desktop sessions (rust/tests/ui)
sh rust/tools/run_python_suite.sh         # the Python package against the native libraries
```

The UI scripts open real windows and drive the desktop through
`ATOMFORGE_UI_SCRIPT`; each writes screenshots to `rust/target/ui/`.

---

## Portable packages

```bash
python rust/tools/package.py
```

builds in release mode and writes, in `rust/target/package/`:

- Windows: `AtomForge-<version>-win64.zip` with `AtomForge/AtomForge.exe` and its libraries,
- Linux: `AtomForge-<version>-linux-x86_64.tar.gz` with `AtomForge/bin/AtomForge` and its
  libraries (spglib included; the program finds them beside it),

and in both the Python package (`AtomForge/python`), `README.md`, the manual
and the example structures. Extract the archive and run `AtomForge` directly.
`--no-build` packages the existing release build.

---

## Cleaning the build

```bash
cd rust
cargo clean
```

---

## Troubleshooting

### Choosing the renderer

The desktop draws through wgpu, which picks the best available backend
(Direct3D 12 or Vulkan on Windows, Vulkan or OpenGL on Linux, Metal on macOS).
To choose explicitly, in order of precedence:

- `--renderer vulkan|opengl|dx12|metal` on the command line,
- the environment variable `ATOMFORGE_RENDERER=vulkan|opengl|dx12|metal`,
- **Settings > Renderer** (saved as `renderer=` in `~/.atomforge_preferences.txt`, applied at the next start).

On remote sessions or virtual machines without a GPU driver, `--renderer opengl`
(with Mesa's software rasteriser on Linux) usually works.

### spglib is not found (Linux)

Install the development package (`libsymspg-dev` on Debian/Ubuntu,
`spglib-devel` on Fedora), or point the linker at a local build:

```bash
RUSTFLAGS="-L /path/to/spglib/lib" cargo build --release -p af-app
```

### The Python package cannot find the native libraries

Set `ATOMFORGE_SCIENCE_LIBRARY` and `ATOMFORGE_ELECTRONIC_LIBRARY` to the
libraries' full paths, or keep the `python` folder next to the program as in
the portable package. `ATOMFORGE_PATH` names the `AtomForge` program for the
functions that run it.
