# AtomForge – Installation and Build Guide

## Requirements

| Dependency | Version | Notes |
| --- | --- | --- |
| CMake | ≥ 3.16 | Build system |
| C++ compiler | C++17 | GCC / Clang / MinGW-w64 |
| OpenGL | 3.x | GPU rendering |
| GLFW3 | any current | Window and input |
| GLEW | any current | OpenGL extension loader |
| GLM | any current | Math library |
| Open Babel | 3.x | Structure file I/O |
| spglib / symspg | any current | **Optional** – symmetry features |

---

## Linux

### Install dependencies

```bash
sudo apt update
sudo apt install build-essential cmake pkg-config \
                 libglfw3-dev libglew-dev libglm-dev \
                 libopenbabel-dev

# optional: symmetry support
sudo apt install libsymspg-dev
```

### Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### Run

```bash
./build/AtomForge
./build/AtomForge structure.cif
```

### CMake options

| Option | Default | Description |
| --- | --- | --- |
| `CMAKE_BUILD_TYPE` | `Release` | `Debug` or `Release` |
| `ATOMFORGE_ENABLE_SPGLIB` | auto-detected | Enable spglib symmetry features |
| `BUILD_PORTABLE` | `OFF` | Bundle runtime libraries for redistribution |

---

## Windows (MSYS2 UCRT64)

> **Important:** Build and run from the MSYS2 UCRT64 shell.  
> PowerShell often misses required toolchain and DLL paths.

### Install MSYS2

Download and install from <https://www.msys2.org/>, then open the **UCRT64** shell.

### Update and install packages

```bash
pacman -Syu
# Close and reopen the UCRT64 shell, then run:
pacman -Su

pacman -S --needed \
    mingw-w64-ucrt-x86_64-toolchain \
    mingw-w64-ucrt-x86_64-cmake \
    mingw-w64-ucrt-x86_64-glfw \
    mingw-w64-ucrt-x86_64-glew \
    mingw-w64-ucrt-x86_64-glm \
    mingw-w64-ucrt-x86_64-openbabel \
    mingw-w64-ucrt-x86_64-pkgconf

# optional: symmetry support
pacman -S --needed mingw-w64-ucrt-x86_64-spglib
```

### Build

```bash
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### Run

```bash
./build/AtomForge.exe
./build/AtomForge.exe structure.cif
```

---

## Cleaning the build

Remove compiled objects while keeping the CMake configuration:

```bash
cmake --build build --target clean
```

Remove the entire build directory:

```bash
# Linux / MSYS2
rm -rf build
```

---

## Portable builds

Portable packaging bundles runtime libraries and Open Babel plugins into a redistributable archive.

### Windows portable ZIP

```bash
# In the MSYS2 UCRT64 shell:
cmake -S . -B build-mingw -G "MinGW Makefiles" \
      -DCMAKE_BUILD_TYPE=Release -DBUILD_PORTABLE=ON
cmake --build build-mingw -- -j1

# Create the package:
cpack --config build-mingw/CPackConfig.cmake -B build-mingw/package
```

Output: `build-mingw/package/AtomForge-1.0.0-win64.zip`

Contents:
- `AtomForge.exe`
- Required runtime DLLs
- Open Babel plugins under `openbabel/<version>/`
- `README.md`

Extract and run `AtomForge.exe` directly.

> Do not reuse a build folder configured with a different generator (e.g., NMake).  
> Use a dedicated folder such as `build-mingw`.

### Linux portable tarball

```bash
cmake -S . -B build-portable -DCMAKE_BUILD_TYPE=Release -DBUILD_PORTABLE=ON
cmake --build build-portable -j
cpack --config build-portable/CPackConfig.cmake -B build-portable/package
```

Output: `build-portable/package/AtomForge-*.tar.gz`

Contents:
- `AtomForge/bin/AtomForge`
- Required `.so` libraries copied next to the executable
- Open Babel plugins under `AtomForge/bin/openbabel/`
- `AtomForge/README.md`

Run:

```bash
tar xzf build-portable/package/AtomForge-*.tar.gz
cd AtomForge/bin
./AtomForge
```

### System installation (Linux)

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
sudo cmake --install build
```

Installs to `/usr/local/bin` by default.

---

## Troubleshooting

### Open Babel plugin path missing

If AtomForge reports missing Open Babel plugins, set `BABEL_LIBDIR` before launch:

```bash
# MSYS2 UCRT64
export BABEL_LIBDIR=/ucrt64/lib/openbabel/3.1.1/

# Linux (Debian/Ubuntu)
export BABEL_LIBDIR=/usr/lib/x86_64-linux-gnu/openbabel/3.1.1/

# Linux (from source)
export BABEL_LIBDIR=/usr/local/lib/openbabel/3.1.1/
```

Verify the path:

```bash
echo $BABEL_LIBDIR
ls "$BABEL_LIBDIR"
```

### Missing DLLs at startup (Windows)

Run from the UCRT64 shell, or ensure `C:/msys64/ucrt64/bin` is on your `PATH`.

### CMake cannot find dependencies

Verify that `pkgconf` and all required MSYS2 packages are installed in the same UCRT64 environment used for building and running.
