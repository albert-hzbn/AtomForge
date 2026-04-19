# AtomForge – GitHub Copilot Instructions

## Project overview

AtomForge is a C++11 desktop application for interactive atomic structure building and visualization, targeting metallurgical simulation and atomistic modeling workflows (MD / DFT / first-principles). It runs on Linux and Windows (MSYS2 UCRT64 / MinGW-w64).

Core capabilities implemented so far:
- Bulk crystal generation with full space-group support (via spglib when available).
- CSL grain boundary bicrystal builder with sigma, rotation axis, and plane control.
- Nanocrystal builder supporting 8 shape primitives (sphere, ellipsoid, box, cylinder, octahedron, truncated octahedron, cuboctahedron, OBJ/STL mesh fill).
- Polycrystalline microstructure builder via Voronoi tessellation with Bunge Euler orientation control.
- Heterogeneous interface builder with 2D supercell strain matching.
- Common Neighbour Analysis (CNA) and Radial Distribution Function (RDF) analysis dialogs.
- Interactive atom editing: select, delete, grab/translate (Blender-style), box select.
- Measurement overlays: distance, angle, atom info, element labels, orientation axes.
- Lattice-plane and Voronoi cell overlays drawn directly into the viewport.
- Image export to PNG, JPEG, and SVG including high-resolution scaling.
- Full undo/redo stack.
- Three rendering tiers auto-selected by atom count.

---

## Build system

- **CMake ≥ 3.16**, C++17 (`CMAKE_CXX_STANDARD 17`, extensions OFF).
- Required libraries: OpenGL 3.x, GLFW3, GLEW, GLM, Open Babel 3 — all resolved via `pkg-config`.
- Optional: **spglib/symspg** (toggled by `ATOMFORGE_ENABLE_SPGLIB`); guarded with `#ifdef HAVE_SPGLIB`.
- Portable Windows packaging via `BUILD_PORTABLE=ON`; CPack produces `AtomForge-<ver>-win64.zip`.
- When adding a new `.cpp` source, add it to the `ATOMFORGE_SOURCES` list in `CMakeLists.txt` **and** include its directory in `target_include_directories`.
- Do not add new third-party libraries without a corresponding `pkg_check_modules` / `find_package` call and a note in `README.md`.

---

## Directory layout

```
src/
  main.cpp                       – entry point, calls runAtomsEditor()
  app/
    EditorApplication.cpp/.h     – main loop, all subsystem wiring, AdaptiveRenderState
    EditorState.h                – single struct holding all runtime mutable state
    EditorOps.cpp/.h             – pure editor operations (captureSnapshot, updateBuffers, applySnapshot, etc.)
    FileDropHandler.cpp/.h       – GLFW drop callback → StructureFileService
    ImageExport.cpp/.h           – off-screen render → PNG/JPEG/SVG
    InteractionHandlers.cpp/.h   – keyboard, mouse, picking, box-select, grab mode
    SceneView.cpp/.h             – FrameView computation (projection, view, light matrices)
    StructureFileService.cpp/.h  – load/save with optional supercell baked-in
    UndoRedo.h                   – UndoRedoManager + EditorSnapshot
    WindowSetup.cpp/.h           – GLFW window creation, camera callbacks
  algorithms/
    BulkCrystalBuilder.cpp/.h    – space-group crystal generation
    CSLComputation.cpp/.h        – CSL mathematics (M3 type, sigma search, rotation)
    InterfaceBuilder.cpp/.h      – 2D supercell strain matching, heterojunction construction
    MeshLoader.cpp/.h            – OBJ / STL triangle mesh parser
    NanoCrystalBuilder.cpp/.h    – shape-based cluster carving
    PolyCrystalBuilder.cpp/.h    – Voronoi polycrystal builder with IPF orientation
    VoronoiComputation.cpp/.h    – Voronoi diagram computation (PBC-aware)
  camera/
    Camera.cpp/.h                – arcball/orbit camera with pan, zoom, roll
  graphics/
    BillboardMesh.cpp/.h         – quad billboard geometry for BillboardImposters mode
    CustomStructureShaders.cpp/.h – mesh-fill VS/FS GLSL sources
    CylinderMesh.cpp/.h          – instanced bond cylinder geometry
    LowPolyMesh.cpp/.h           – 12-facet icosahedron for LowPolyInstancing mode
    Picking.cpp/.h               – ray-sphere CPU picker + ray unprojection
    Renderer.cpp/.h              – owns all GLSL programs, per-pass draw methods
    SceneBuffers.cpp/.h          – GPU VBOs, CPU pick caches, rendering mode
    Shader.cpp/.h                – compile/link helper for GLSL programs
    ShadowMap.cpp/.h             – depth FBO for shadow mapping
    SphereMesh.cpp/.h            – indexed sphere mesh for StandardInstancing mode
    StructureInstanceBuilder.cpp/.h – builds StructureInstanceData from Structure
  io/
    StructureLoader.cpp/.h       – Open Babel wrapper, AtomSite, Structure definitions
  math/
    StructureMath.h              – header-only: cell matrix helpers, min-image convention
  third_party/
    stb_image_write.h            – single-header PNG/JPEG writer
  ui/
    AtomContextMenu.cpp/.h       – right-click context menu for atoms
    BulkCrystalBuilderDialog.cpp/.h
    CommonNeighbourAnalysis.cpp/.h
    CSLGrainBoundaryDialog.cpp/.h
    CustomStructureDialog.cpp/.h – mesh-fill structure builder dialog
    EditMenuDialogs.cpp/.h       – element radii/colours/structure edit modals
    FileBrowser.cpp/.h           – unified menu bar (File/Edit/View/About) + all state flags
    ImGuiSetup.cpp/.h            – Dear ImGui init with HiDPI font loading
    InterfaceBuilderDialog.cpp/.h
    LatticePlaneOverlay.cpp/.h   – HKL plane overlay
    MeasurementOverlay.cpp/.h    – distance/angle/atom-info/element-label overlays, IPF legend
    NanoCrystalBuilderDialog.cpp/.h
    PeriodicTableDialog.cpp/.h   – interactive element picker
    PolyCrystalBuilderDialog.cpp/.h
    RadialDistributionAnalysis.cpp/.h
    StructureInfoDialog.cpp/.h   – structure metadata panel
    ThemeUtils.h                 – dark/light theme + HiDPI scale helpers
    TransformAtomsDialog.cpp/.h  – 3×3 supercell transform matrix UI
    VoronoiOverlay.cpp/.h        – Voronoi cell wire overlay
  util/
    ElementData.cpp/.h           – symbol, name, covalent radii, CPK colours (Z 1–118)
    PathUtils.cpp/.h             – cross-platform path helpers
imgui/                           – vendored Dear ImGui 1.x + GLFW/OpenGL3 backends
cmake/
  InstallPortableDependencies.cmake.in – Windows DLL bundling script
```

---

## Core data types

### `AtomSite` — `src/io/StructureLoader.h`
Per-atom record. Fields: `symbol` (string), `atomicNumber` (int), `x`/`y`/`z` (double, Cartesian Å), `r`/`g`/`b` (float, display colour).

### `Structure` — `src/io/StructureLoader.h`
```
atoms          std::vector<AtomSite>              canonical atom list (no PBC images)
hasUnitCell    bool
cellVectors    double[3][3]                        row-major lattice vectors (Å)
cellOffset     double[3]                           Cartesian origin of cell
pbcBoundaryTol float                               tolerance for boundary image insertion
grainColors    std::vector<std::array<float,3>>    IPF-Z colours (same size as atoms)
grainRegionIds std::vector<int>                    per-atom grain ID for GB coloring
ipfLoadStatus  std::string                         metadata note
```

### `EditorState` — `src/app/EditorState.h`
Single global mutable state struct. Fields of interest:
```
structure                Structure
fileBrowser              FileBrowser           holds all View/Edit flag state
editMenuDialogs          EditMenuDialogs       element radii/colour tables
sceneBuffers             SceneBuffers
selectedInstanceIndices  std::vector<int>      GPU instance indices (not atom indices)
contextMenu              AtomContextMenu
measurementState         MeasurementOverlayState
structureInfoDialog      StructureInfoDialogState
undoRedo                 UndoRedoManager
voronoiDiagram           VoronoiDiagram
voronoiDirty             bool                  set true whenever structure changes
grabState                GrabState             Blender-style grab mode
suppressHistoryCommit    bool
pendingDefaultViewReset  bool
pendingDroppedFiles      std::vector<std::string>
```

### `EditorSnapshot` — `src/app/UndoRedo.h`
```
structure        Structure
elementRadii     std::vector<float>
elementColors    std::vector<glm::vec3>
elementShininess std::vector<float>
```
`operator==` is defined for change-detection before committing.

### `SceneBuffers` — `src/graphics/SceneBuffers.h`
Owns all GPU VBOs for atoms and bonds; also stores CPU-side pick caches when atom count ≤ 100 k. Fields: `instanceVBO`, `colorVBO`, `scaleVBO`, `bondStart/EndVBO`, `lineVAO/VBO`, `atomCount`, `bondCount`, `renderMode`, `orbitCenter`, plus CPU vectors `atomPositions`, `atomColors`, `atomRadii`, `atomShininess`, `atomIndices`.

### `RenderingMode` — `src/graphics/SceneBuffers.h`
```
StandardInstancing   ≤ 100 k atoms   full sphere mesh (indexed draw)
LowPolyInstancing    100 k – 10 M    12-facet icosahedron
BillboardImposters   ≥ 10 M          quad billboards, sphere lighting in FS
```
Auto-selected by `Renderer::selectRenderingMode(atomCount)`. During runtime, `AdaptiveRenderState` (inside `EditorApplication.cpp`) re-evaluates the mode based on a 30-frame rolling average (downgrades if avg > 33 ms, upgrades if avg < 20 ms).

### `StructureInstanceData` — `src/graphics/StructureInstanceBuilder.h`
Per-instance CPU arrays fed to VBOs: `positions`, `colors`, `scales`, `bondRadii`, `shininess`, `atomicNumbers`, `boxLines`, `atomIndices`, `orbitCenter`. Built by `buildStructureInstanceData(structure, useTransformMatrix, transformMatrix, elementRadii, elementShininess)`. Also exposes `buildSupercell(structure, transformMatrix)` which bakes the supercell into a new `Structure`.

### `FrameView` — `src/app/SceneView.h`
Per-frame computed matrices: `framebufferWidth/Height`, `windowWidth/Height`, `projection`, `view`, `lightMVP`, `lightPosition`, `cameraPosition`. Built by `buildFrameView(camera, sceneBuffers, useOrthographicView, frame)`.

### `FrameActionRequests` — `src/app/InteractionHandlers.h`
One-shot boolean flags produced by `beginFrameActionRequests()` and populated by `applyKeyboardShortcuts()`. Consumed in the main loop. Flags: `doDeleteSelected`, `requestMeasureDistance`, `requestMeasureAngle`, `requestAtomInfo`, `requestStructureInfo`, `requestUndo`, `requestRedo`, `requestViewAxisX/Y/Z`, `requestViewLatticeA/B/C`, `requestRotateCrystalX/Y/Z`.

---

## Architecture conventions

### Application loop (`EditorApplication.cpp`)
`runAtomsEditor()` creates the GLFW window, initialises all subsystems (OpenGL, Dear ImGui, meshes, renderer, shadow map), then enters the per-frame loop:
1. `AdaptiveRenderState::update` — possibly switches rendering mode.
2. Shadow pass — depth-only render for all atoms and bonds.
3. Colour pass — atoms, bonds, bounding box lines.
4. ImGui frame — menu bar + all dialogs + overlays.
5. `applyKeyboardShortcuts`, `handlePendingAtomPick`, `handleBoxSelection`, `handleGrabMode`.
6. Consume `FrameActionRequests`, apply undo/redo, view requests, measurements.
7. Swap buffers.

### State management
- All editor mutations go through `EditorState`. Never hold long-lived pointers to `Structure` internals; re-index after any modification.
- The **only** correct way to update GPU data after a structure change is to call `updateBuffers(state)` or `updateBuffers(state, structure)` from `EditorOps.h`.
- Undo/redo: call `state.undoRedo.commit(captureSnapshot(state))` after every user-visible structural change. Guard with `if (!state.suppressHistoryCommit)`. Use `state.suppressHistoryCommit = true` when batching sub-steps.
- PBC images are secondary display data — never write them into `Structure::atoms`. Set `state.voronoiDirty = true` after structural changes so the Voronoi diagram is recomputed on the next frame.
- `EditorState::selectedInstanceIndices` stores GPU instance indices (1-to-1 with `SceneBuffers::atomIndices`), not raw `Structure::atoms` indices. Use `sceneBuffers.atomIndices[instanceIdx]` to go from instance to atom.
- `GrabState` is an ephemeral per-gesture struct. Do not persist it across frames after the grab is confirmed or cancelled.

### Graphics pipeline
- `Renderer` owns all GLSL programs: `atomProgram`, `atomLowPolyProgram`, `atomBillboardProgram`, `bondProgram`, `shadowProgram`, `shadowLowPolyProgram`, `shadowBillboardProgram`, `bondShadowProgram`, `lineProgram`.
- `SceneBuffers` owns GPU VBOs. Call `SceneBuffers::upload(...)` / `SceneBuffers::update(...)` after structural changes; do not write to VBOs directly outside that struct.
- Shadow mapping uses a dedicated depth FBO (`ShadowMap`). `beginShadowPass` / `endShadowPass` bracket the depth-only draw; the depth texture is then sampled in the colour pass.
- For atom picking: when atom count ≤ 100 k, CPU ray–sphere via `pickAtom()` in `Picking.h`. For larger structures the CPU cache is disabled; fall back to bounding-sphere heuristics or skip picking.
- `pickRayDir` / `pickRayOrigin` in `Picking.h` unproject window-space cursor to world-space ray.

### Camera (`Camera`)
- Arcball/orbit model: `yaw`, `pitch`, `roll` (degrees) + `distance` + `panOffset`.
- Static GLFW callbacks (`mouseButton`, `cursor`, `scroll`) write to `Camera::instance`.
- `pendingClick` / `pendingRightClick`: one-frame flags set when a short drag ends. Consumed by `handlePendingAtomPick` / `handleRightClick`.
- `dragAccum` / `rightDragAccum`: accumulated Manhattan distance; gates click vs drag disambiguation.

### UI / ImGui
- All dialogs live in `src/ui/`. Analysis dialogs (`CommonNeighbourAnalysis`, `RadialDistributionAnalysis`) use a member `drawMenuItem` + `drawDialog` split. Builder dialogs use a standalone `Render(EditorState&, ...)` free function.
- `FileBrowser` is the **unified menu bar host** — it draws File, Edit, View, and About menus and owns all boolean view-flag state (show bonds, element labels, orthographic view, box-select mode, lattice planes, Voronoi, atom colour mode, theme). Access flags via its accessor methods; never replicate them elsewhere.
- `EditMenuDialogs` owns the per-element radius and colour tables (`elementRadii`, `elementColors`, `elementShininess`). These are the ground truth for display properties, consumed by `buildStructureInstanceData`.
- Use `ImGui::BeginPopupModal` for dialogs that block the main window; use `ImGui::Begin` for non-modal tool panels.
- Do not call `ImGui::*` outside the `src/ui/` and `src/app/` layers.
- Theme utilities live in `src/ui/ThemeUtils.h`; use them for HiDPI-aware font/scale setup. Never hard-code pixel sizes.
- Toast notifications are drawn in `EditorApplication`; push `ToastNotification{message, expiresAt}` into the local queue to surface transient feedback.

### Algorithm modules (`src/algorithms/`)
- Each algorithm is **pure compute**: no ImGui calls, no OpenGL, no global state.
- Builders return a result struct with `success`, `message`, and output `Structure`. They do **not** mutate `EditorState` directly; callers apply the result and call `updateBuffers`.
- `CSLComputation.h` defines the internal `M3` matrix type (3×3 double, row-major), helpers `mul3`, `inv3`, `det3`, and the CSL search functions. Do not use `glm::mat3` inside algorithm math; use `M3` there and convert at call sites.
- `StructureMath.h` provides the bridge: `makeCellMatrix(structure)` → `glm::mat3`, `tryMakeCellMatrices`, `minimumImageDelta`, `tryCartesianToFractional`. Bond tolerance constants: `kBondToleranceFactor = 1.18f`, `kMinBondDistance = 0.10f`.
- `MeshLoader` parses OBJ and binary/ASCII STL into flat triangle lists; consumed by `NanoCrystalBuilder` (mesh-fill shape) and `CustomStructureDialog`.
- `VoronoiComputation`: PBC-aware Voronoi; result is `VoronoiDiagram` (vector of `VoronoiCell`, each containing convex polygon faces). Recompute when `EditorState::voronoiDirty` is true.

### I/O
- `StructureLoader` wraps Open Babel. Always go through `loadStructureFromFile(path, structure, errorMessage)` for file input. `isSupportedStructureFile(filename)` validates extensions before invoking Open Babel.
- `StructureFileService` provides `loadStructureFromPath` (with error text) and `saveStructureWithOptionalSupercell` (bakes the 3×3 transform matrix into atom coords before writing).
- `FileDropHandler` handles GLFW drag-and-drop events; it validates extensions and queues paths in `EditorState::pendingDroppedFiles`, which the main loop drains.
- Never call Open Babel APIs outside `src/io/`.

### Image export
- `ImageExport.h` exports `exportStructureImage(request, view, bgColor, showBonds, sceneBuffers, renderer, ...)`. It performs an off-screen render at `request.resolutionScale × window size` using `stb_image_write` (PNG/JPEG) or a custom SVG path builder.
- `ImageExportFormat`: `Png`, `Jpg`, `Svg`. `ImageExportRequest` carries `outputPath`, `format`, `includeBackground`, `resolutionScale`.

### Element data (`src/util/ElementData.h`)
- `elementSymbol(z)` / `elementName(z)` — Z 1–118, returns `"?"` / `"Unknown"` for out-of-range.
- `makeLiteratureCovalentRadii()` — Cordero 2008 covalent radii, vector size 119 (index 0 unused).
- `makeDefaultElementColors()` — CPK-style colours, same indexing.

---

## Coding standards

- **C++17** — use features like structured bindings, `std::optional`, `std::variant`, `if constexpr`, and fallthrough. Avoid C++20+ features (concepts, ranges, three-way comparison, designated initializers in most cases).
- Include guards: `#pragma once` on every header.
- GLM for all linear-algebra in rendering/camera code (`glm::vec3`, `glm::mat3`, `glm::mat4`, `glm::quat`). Do not introduce Eigen or another LA library.
- Algorithm math that is crystallographic (CSL, lattice) uses the local `M3` type or plain `double[3][3]` arrays; convert with `makeCellMatrix` / `glm::make_mat3` at the boundary.
- No exceptions in rendering or hot paths; use return-code booleans and result structs.
- OpenGL calls require GLEW-resolved pointers (`#include <GL/glew.h>` before any GL headers).
- No RTTI (`dynamic_cast`, `typeid`).
- Prefer free functions over class methods for stateless algorithms.
- Structs with only data fields and no invariants should not have constructors beyond member-initializers in the declaration.
- Raw arrays (`int mat[3][3]`) are used for transform/CSL matrices throughout; pass them as `const int (&m)[3][3]` to preserve dimensions.
- All per-element tables are indexed by atomic number (1-based, index 0 unused), size 119.

---

## Rendering mode details

| Mode | Trigger | Geometry | Draw call |
|---|---|---|---|
| `StandardInstancing` | ≤ 100 k | `SphereMesh` (full subdivided sphere) | `glDrawElementsInstanced` |
| `LowPolyInstancing` | 100 k – 10 M | `LowPolyMesh` (12-facet icosahedron) | `glDrawElementsInstanced` |
| `BillboardImposters` | ≥ 10 M | `BillboardMesh` (one quad per atom) | `glDrawArraysInstanced`, sphere lit in fragment shader |

`AdaptiveRenderState` (private to `EditorApplication.cpp`) adjusts mode at runtime. It has a 60-frame cool-down after each switch and a 30-frame rolling average window. Downgrade threshold: > 33 ms avg. Upgrade threshold: < 20 ms avg.

---

## Atom coloring modes

Defined as `AtomColorMode` in `src/ui/FileBrowser.h`:

| Value | Meaning |
|---|---|
| `ElementType` | CPK default colours from `EditMenuDialogs::elementColors` |
| `CrystalOrientation` | IPF-Z from `Structure::grainColors`; requires non-empty `grainColors` vector |
| `GrainBoundary` | Highlights grain-boundary atoms; uses `Structure::grainRegionIds` when present, else infers from CNA |

Coloring logic lives in `StructureInstanceBuilder.cpp`. Add a new mode by: (1) appending to `AtomColorMode`, (2) adding a branch in `buildStructureInstanceData`, (3) adding a menu item in `FileBrowser`.

---

## Builder patterns

### Bulk crystal builder (`BulkCrystalBuilder`)
- `CrystalSystem` enum: Triclinic, Monoclinic, Orthorhombic, Tetragonal, Trigonal, Hexagonal, Cubic.
- `LatticeParameters`: a/b/c (Å) + α/β/γ (°).
- `SpaceGroupRange` table (`kSpaceGroupRanges[7]`) maps system to first/last space group.
- `buildLatticeFromParameters` → `glm::mat3`. `validateParameters` / `applySystemConstraints` enforce crystal-system constraints.
- `hallBySpaceGroup()` returns reference to static vector of Hall numbers (used with spglib or internal tables).

### Nanocrystal builder (`NanoCrystalBuilder`)
- `NanoShape` enum (8 values). `NanoParams` carries all size/shape parameters plus `autoReplicate`, `autoCenterFromAtoms`, `setOutputCell`, `vacuumPadding`.
- `NanoBuildResult`: `success`, `message`, `inputAtoms`, `outputAtoms`, `shape`, `estimatedDiameter`, `repA/B/C`.

### Polycrystal builder (`PolyCrystalBuilder`)
- `PolyParams`: box dimensions (Å), `numGrains`, `seed`, `orientationMode` (`AllRandom` / `AllSpecified` / `PartialSpecified`), `specifiedOrientations` (Bunge Euler φ₁/Φ/φ₂ per grain).
- `buildPolycrystal` writes `Structure::grainColors` (IPF-Z) and `Structure::grainRegionIds`.
- References `glm::quat` for rotation; uses `glm/gtc/quaternion.hpp`.

### Interface builder (`InterfaceBuilder`)
- `Mat2` (2×2 double) for in-plane lattice representations.
- `generateUniqueSupercells` enumerates all coprime supercell matrices up to `nmax` with area ≤ `maxCells`.
- `strainComponents` / `meanAbsStrain` / `cubicElasticDensity` for strain matching scoring.
- `equivalentLatticeKey` deduplicates lattices with the same geometry.

### CSL builder (`CSLComputation`)
- `SigmaCandidate`: `sigma`, `m`/`n`, `thetaDeg`, `csl[3][3]`, `plane[3][3]`.
- `getRotateMatrix(axis, angleDeg)` → `M3`. `oLatticeToCsl`, `reduceCsl`, `orthogonalizeCsl`, `getCslMatrix` are the CSL math pipeline.

---

## Overlay and measurement system

All overlays are drawn via `ImDrawList*` in the ImGui render phase, projected with the same `projection` + `view` matrices used for the 3D scene.

| Overlay | Header | Function |
|---|---|---|
| Lattice planes | `LatticePlaneOverlay.h` | `drawLatticePlanesOverlay(drawList, ...)` |
| Voronoi cells | `VoronoiOverlay.h` | `drawVoronoiOverlay(drawList, ...)` |
| Distance / angle / atom info | `MeasurementOverlay.h` | `drawMeasurementOverlays(...)` + `drawMeasurementPopups(...)` |
| Element labels | `MeasurementOverlay.h` | `drawElementLabelsOverlay(...)` |
| Orientation axes | `MeasurementOverlay.h` | `drawOrientationAxesOverlay(...)` |
| IPF legend | `MeasurementOverlay.h` | `drawIPFTriangleLegend(...)` |

`MeasurementOverlayState` in `EditorState` holds all popup flags and message buffers. `processMeasurementRequests` fills them from the current selection.

`LatticePlane` struct: `h`, `k`, `l` (Miller), `offset`, `opacity`, `color[3]`, `visible`. Stored in `FileBrowser::latticePlanes`.

---

## Feature areas to extend

1. **Export formats** – add new writer in `src/io/`; hook into `StructureFileService::saveStructureWithOptionalSupercell` and the File → Save As menu in `FileBrowser`.
2. **New structure builders** – add `src/algorithms/FooBuilder.cpp/.h` + `src/ui/FooBuilderDialog.cpp/.h`; register in `EditorApplication.cpp`'s main loop; add `ATOMFORGE_SOURCES` entries in `CMakeLists.txt`.
3. **New analysis tool** – add a dialog struct in `src/ui/` with `drawMenuItem(bool enabled)` + `drawDialog(const Structure&)` members, following `CommonNeighbourAnalysisDialog`. Compute in `src/algorithms/`.
4. **New shader effect** – add `GLuint fooProgram = 0` to `Renderer`, compile in `Renderer::init()`, and add a `drawFoo(...)` method. Add matching shadow variant if the effect must cast/receive shadows.
5. **New coloring mode** – append to `AtomColorMode` → add branch in `StructureInstanceBuilder.cpp` → add menu item in `FileBrowser`.
6. **Keyboard shortcuts** – add a `bool requestFoo` field to `FrameActionRequests`, populate it in `applyKeyboardShortcuts`, handle it in the main loop. Mirror in the appropriate menu.
7. **New overlay** – implement a free function `drawFooOverlay(ImDrawList*, projection, view, width, height, ...)` in `src/ui/`; call it in the ImGui render phase in `EditorApplication.cpp`.
8. **spglib-dependent features** – guard with `#ifdef HAVE_SPGLIB`; test both code paths.
9. **New per-atom property** – add the field to `AtomSite`, extend `EditorSnapshot::operator==`, update `StructureLoader` to populate it, add display/edit UI in `EditMenuDialogs`.

---

## Patterns to follow when adding code

```cpp
// New algorithm result struct (src/algorithms/)
struct FooResult {
    bool        success = false;
    std::string message;
    Structure   output;
    int         generatedAtoms = 0;
};

FooResult buildFoo(const Structure& reference, const FooParams& params);

// Caller pattern in a dialog callback (src/ui/)
FooResult result = buildFoo(state.structure, params);
if (result.success)
{
    state.structure = result.output;
    updateBuffers(state);
    state.voronoiDirty = true;
    if (!state.suppressHistoryCommit)
        state.undoRedo.commit(captureSnapshot(state));
}

// New analysis dialog (src/ui/)
struct FooAnalysisDialog {
    void drawMenuItem(bool enabled);      // call inside an open ImGui menu
    void drawDialog(const Structure& s);  // call once per frame outside menus
private:
    bool m_openRequested = false;
};

// New non-modal overlay (src/ui/)
void drawFooOverlay(ImDrawList* drawList,
                    const glm::mat4& projection,
                    const glm::mat4& view,
                    int w, int h,
                    const Structure& structure,
                    bool enabled);

// New keyboard shortcut
// In FrameActionRequests: bool requestFoo = false;
// In applyKeyboardShortcuts:
if (ImGui::IsKeyPressed(ImGuiKey_F, false) && !io.WantCaptureKeyboard)
    requests.requestFoo = true;
// In main loop:
if (requests.requestFoo) { /* apply */ }

// New shader program
// In Renderer.h:  GLuint fooProgram = 0;
// In Renderer::init():
fooProgram = createProgram(fooVertSrc, fooFragSrc);
// In Renderer: add drawFoo(projection, view, ...) method.
```

---

## Version control

- Prefer small, focused commits that address one logical change at a time; avoid mixing refactors, feature work, formatting churn, and generated-file updates in the same commit.
- Before concluding a task, review the working tree (`git status`, `git diff --stat`, `git diff`) and ensure the change set matches the requested scope.
- Do not rewrite, discard, or revert user-authored changes unless explicitly asked. If the worktree contains unrelated edits, leave them intact and work around them.
- Treat `build/`, `build-portable/`, `package/`, and other generated outputs as non-source artifacts. Do not edit or commit them unless the task is explicitly about packaging, build-system output, or generated assets.
- When source changes alter behavior, keep documentation and citations in sync in the same branch or commit when practical.
- If a change requires multiple commits, split them by concern: infrastructure/setup first, core implementation next, and documentation or cleanup last.
- Prefer non-destructive git workflows. Avoid history rewrites (`reset --hard`, force pushes, rebases that drop commits) unless the user explicitly requests them.

---

## What NOT to do

- Do not add C++20 or later features (project targets C++17 strictly).
- Do not write directly to `Structure::atoms` from a UI callback; go through a builder or edit helper and then call `updateBuffers`.
- Do not store GPU state (VAOs, VBOs, shader IDs) anywhere except `SceneBuffers` or `Renderer`.
- Do not use `std::shared_ptr` or `std::unique_ptr` for scene objects; the ownership model is value-based structs.
- Do not call Open Babel APIs outside `src/io/`.
- Do not call `ImGui::*` outside `src/ui/` or `src/app/`.
- Do not replicate view-flag state (`showBonds`, `atomColorMode`, etc.) outside `FileBrowser`; read them via its accessors.
- Do not allocate CPU pick caches (`atomPositions`, `atomRadii`, etc.) for atom counts > 100 k; the code already skips this.
- Do not hard-code pixel sizes or font sizes; use `ThemeUtils.h` scale helpers.
- Do not add new third-party libraries without updating `CMakeLists.txt` and `README.md`.
- Do not write PBC image atoms into `Structure::atoms`; images are generated on the fly for display only.
- Do not store per-frame temporaries (matrices, draw lists) in `EditorState`; keep them local to the frame or in `FrameView`.
