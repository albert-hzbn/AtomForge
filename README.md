# Atoms Editor

A small OpenGL-based molecular visualization demo using ImGui and Open Babel.

## Prerequisites

On Debian/Ubuntu (or derivatives):

```bash
sudo apt update
sudo apt install build-essential cmake libglfw3-dev libglew-dev libglm-dev \
                 libopenbabel-dev libopenbabel3
```

> **Note:** This project uses OpenGL 3.0+ (GLSL 130) and GLFW.

## Build

From the `src/` directory:

```bash
make
```

Then run:

```bash
./atoms_editor
```

## Usage

- Use **File → Open...** to pick a structure file.
- Supported file types: `.cif`, `.mol`, `.pdb`, `.xyz`, `.sdf`.
- Use mouse:
  - **Left drag** to rotate
  - **Scroll** to zoom

## Project layout

- `src/`: main application sources
- `src/ui/`: ImGui-based UI components
- `src/io/`: structure loading code (Open Babel)
- `imgui/`: bundled ImGui sources
- `backends/`: ImGui GLFW/OpenGL backends