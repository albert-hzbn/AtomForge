# ==============================================================================
# Configuration — edit these sections to add/remove files or change flags
# ==============================================================================

EXE   = AtomForge
CXX   = g++
IMGUI = ../
PKG_CONFIG ?= pkg-config
LINK_MODE ?= shared
STATIC_LIB_DIRS ?= /usr/lib /usr/lib/x86_64-linux-gnu /usr/local/lib /usr/local/lib/x86_64-linux-gnu
PORTABLE_DIR ?= dist/$(EXE)
PORTABLE_LIB_DIR := $(PORTABLE_DIR)/lib

# ------------------------------------------------------------------------------
# Project sources — one file per line for easy diff / addition
# ------------------------------------------------------------------------------

APP_SOURCES = \
    main.cpp                              \
    app/EditorApplication.cpp             \
    app/EditorOps.cpp                     \
    app/InteractionHandlers.cpp           \
    app/WindowSetup.cpp                   \
    ElementData.cpp

CAMERA_SOURCES = \
    camera/Camera.cpp

GRAPHICS_SOURCES = \
    graphics/CylinderMesh.cpp             \
    graphics/Renderer.cpp                 \
    graphics/Picking.cpp                  \
    graphics/SceneBuffers.cpp             \
    graphics/Shader.cpp                   \
    graphics/ShadowMap.cpp                \
    graphics/SphereMesh.cpp               \
    graphics/StructureInstanceBuilder.cpp

UI_SOURCES = \
    ui/AtomContextMenu.cpp                \
    ui/BulkCrystalBuilderDialog.cpp       \
    ui/CSLGrainBoundaryDialog.cpp         \
    ui/CommonNeighbourAnalysis.cpp        \
    ui/EditMenuDialogs.cpp                \
    ui/FileBrowser.cpp                    \
    ui/ImGuiSetup.cpp                     \
    ui/MeasurementOverlay.cpp             \
    ui/PeriodicTableDialog.cpp            \
    ui/RadialDistributionAnalysis.cpp     \
    ui/StructureInfoDialog.cpp            \
    ui/TransformAtomsDialog.cpp

IO_SOURCES = \
    io/StructureLoader.cpp

IMGUI_SOURCES = \
    $(IMGUI)/imgui/imgui.cpp              \
    $(IMGUI)/imgui/imgui_draw.cpp         \
    $(IMGUI)/imgui/imgui_tables.cpp       \
    $(IMGUI)/imgui/imgui_widgets.cpp      \
    $(IMGUI)/backends/imgui_impl_glfw.cpp \
    $(IMGUI)/backends/imgui_impl_opengl3.cpp

SOURCES = $(APP_SOURCES) $(CAMERA_SOURCES) $(GRAPHICS_SOURCES) \
          $(UI_SOURCES) $(IO_SOURCES) $(IMGUI_SOURCES)

# ------------------------------------------------------------------------------
# Compiler flags
# ------------------------------------------------------------------------------

CXXFLAGS = \
    -std=c++11 -O2 -g -Wall -Wformat \
    -I. -I$(IMGUI)/imgui -I$(IMGUI)/backends \
    -Icamera -Igraphics -Imath -Iui -Iio

LDFLAGS =

ifeq ($(LINK_MODE),static)
    PKG_CONFIG_LIBS_MODE = --static --libs
    LDFLAGS += -static -static-libgcc -static-libstdc++
else ifeq ($(LINK_MODE),portable)
    PKG_CONFIG_LIBS_MODE = --libs
    LDFLAGS += -Wl,-rpath,'$$ORIGIN/lib'
else
    PKG_CONFIG_LIBS_MODE = --libs
endif

# ------------------------------------------------------------------------------
# Platform-specific flags
# ------------------------------------------------------------------------------

# Detect Windows when running GNU Make under cmd.exe, MSYS2, MinGW, or Cygwin.
UNAME_S := $(shell uname -s 2>/dev/null)
IS_WINDOWS := 0

ifeq ($(OS),Windows_NT)
    IS_WINDOWS := 1
endif
ifneq (,$(findstring MINGW,$(UNAME_S)))
    IS_WINDOWS := 1
endif
ifneq (,$(findstring MSYS,$(UNAME_S)))
    IS_WINDOWS := 1
endif
ifneq (,$(findstring CYGWIN,$(UNAME_S)))
    IS_WINDOWS := 1
endif

ifeq ($(IS_WINDOWS),1)
    ECHO_MESSAGE = Windows

    # Keep the default target name executable on Windows.
    ifneq ($(suffix $(EXE)),.exe)
        EXE := $(EXE).exe
    endif

    # Prefer pkg-config in MSYS2/MinGW; fallback to common library names when .pc files are absent.
    HAVE_PKG_CONFIG := $(shell command -v $(PKG_CONFIG) >/dev/null 2>&1 && echo yes)
    ifeq ($(HAVE_PKG_CONFIG),yes)
        HAVE_WINDOWS_PKGS := $(shell $(PKG_CONFIG) --exists glfw3 glew openbabel-3 && echo yes)
    endif

    ifeq ($(HAVE_WINDOWS_PKGS),yes)
        CXXFLAGS += $(shell $(PKG_CONFIG) --cflags glfw3 glew openbabel-3)
        LIBS      = $(shell $(PKG_CONFIG) $(PKG_CONFIG_LIBS_MODE) glfw3 glew openbabel-3)
    else
        CXXFLAGS += -I/mingw64/include -I/mingw64/include/openbabel3
        LIBS      = -lglfw3 -lglew32 -lopenbabel
    endif

    LIBS += -lopengl32 -lgdi32 -luser32 -lkernel32 -lshell32 -lcomdlg32 -lole32 -luuid

    ifeq ($(HAVE_PKG_CONFIG),yes)
        # spglib package naming may differ across environments.
        SPGLIB_CFLAGS := $(shell $(PKG_CONFIG) --cflags spglib 2>/dev/null)
        SPGLIB_LIBS   := $(shell $(PKG_CONFIG) $(PKG_CONFIG_LIBS_MODE) spglib 2>/dev/null)
        ifeq ($(strip $(SPGLIB_LIBS)),)
            SPGLIB_CFLAGS := $(shell $(PKG_CONFIG) --cflags symspg 2>/dev/null)
            SPGLIB_LIBS   := $(shell $(PKG_CONFIG) $(PKG_CONFIG_LIBS_MODE) symspg 2>/dev/null)
        endif

        ifneq ($(strip $(SPGLIB_LIBS)),)
            CXXFLAGS += -DATOMS_ENABLE_SPGLIB $(SPGLIB_CFLAGS)
            LIBS += $(SPGLIB_LIBS)
        endif
    endif
else ifeq ($(UNAME_S), Linux)
  ECHO_MESSAGE = Linux
    CXXFLAGS    += $(shell $(PKG_CONFIG) --cflags glfw3 openbabel-3)
    LIBS         = $(shell $(PKG_CONFIG) $(PKG_CONFIG_LIBS_MODE) glfw3 openbabel-3) -lGLEW -lGL -ldl -lpthread

    # spglib package naming differs across distros:
    # some expose pkg-config as "spglib", others as "symspg".
        SPGLIB_CFLAGS := $(shell $(PKG_CONFIG) --cflags spglib 2>/dev/null)
        SPGLIB_LIBS   := $(shell $(PKG_CONFIG) $(PKG_CONFIG_LIBS_MODE) spglib 2>/dev/null)
    ifeq ($(strip $(SPGLIB_LIBS)),)
                SPGLIB_CFLAGS := $(shell $(PKG_CONFIG) --cflags symspg 2>/dev/null)
                SPGLIB_LIBS   := $(shell $(PKG_CONFIG) $(PKG_CONFIG_LIBS_MODE) symspg 2>/dev/null)
    endif

    # Fallback when pkg-config metadata is missing but libraries are present.
    ifeq ($(strip $(SPGLIB_LIBS)),)
        ifneq ($(wildcard /usr/lib/x86_64-linux-gnu/libsymspg.a),)
            SPGLIB_LIBS := -lsymspg
        else ifneq ($(wildcard /usr/lib/x86_64-linux-gnu/libspglib.a),)
            SPGLIB_LIBS := -lspglib
        else ifneq ($(wildcard /usr/lib/x86_64-linux-gnu/libsymspg.so),)
            SPGLIB_LIBS := -lsymspg
        else ifneq ($(wildcard /usr/lib/x86_64-linux-gnu/libspglib.so),)
            SPGLIB_LIBS := -lspglib
        endif
    endif

    ifneq ($(strip $(SPGLIB_LIBS)),)
        CXXFLAGS += -DATOMS_ENABLE_SPGLIB $(SPGLIB_CFLAGS)
        LIBS += $(SPGLIB_LIBS)
    endif
endif

# ==============================================================================
# Build rules — no need to edit below this line
# ==============================================================================

.RECIPEPREFIX := >

OBJS = $(SOURCES:.cpp=.o)

ifeq ($(LINK_MODE),static)
    PRE_LINK_CHECKS = check-static-libs
endif

%.o: %.cpp
> $(CXX) $(CXXFLAGS) -c $< -o $@

all: $(EXE)
> @echo Build complete for $(ECHO_MESSAGE)

ifeq ($(IS_WINDOWS),1)
portable:
> @echo "LINK_MODE=portable is Linux-only in this Makefile."
> @exit 1
else
portable:
> $(MAKE) LINK_MODE=portable portable-bundle
endif

portable-bundle: $(EXE)
> mkdir -p $(PORTABLE_LIB_DIR)
> cp -f $(EXE) $(PORTABLE_DIR)/$(EXE)
> ldd $(EXE) | awk '/=> \/.*\.so/ { print $$3 }' | while read -r lib; do \
>   case "$$(basename "$$lib")" in \
>     libc.so.*|libm.so.*|libpthread.so.*|libdl.so.*|librt.so.*|ld-linux*.so.*|libGL.so.*|libGLX.so.*|libGLdispatch.so.*) continue ;; \
>   esac; \
>   cp -Lf "$$lib" $(PORTABLE_LIB_DIR)/; \
> done
> printf '%s\n' '#!/usr/bin/env sh' 'set -eu' 'HERE=$$(CDPATH= cd -- "$$(dirname -- "$$0")" && pwd)' 'export LD_LIBRARY_PATH="$$HERE/lib$${LD_LIBRARY_PATH:+:$$LD_LIBRARY_PATH}"' 'exec "$$HERE/$(EXE)" "$$@"' > $(PORTABLE_DIR)/run.sh
> chmod +x $(PORTABLE_DIR)/run.sh
> echo "Portable bundle created in $(PORTABLE_DIR)"

check-static-libs:
> @missing=""; \
> find_archive() { \
>   pattern="$$1"; \
>   for dir in $(STATIC_LIB_DIRS); do \
>     if find "$$dir" -maxdepth 2 -name "$$pattern" -print -quit 2>/dev/null | grep -q .; then \
>       return 0; \
>     fi; \
>   done; \
>   return 1; \
> }; \
> if ! find_archive 'libglfw*.a'; then missing="$$missing glfw"; fi; \
> if ! find_archive 'libopenbabel*.a'; then missing="$$missing openbabel"; fi; \
> if ! find_archive 'libGL.a'; then missing="$$missing OpenGL"; fi; \
> if ! find_archive 'libsymspg*.a' && ! find_archive 'libspglib*.a'; then missing="$$missing spglib/symspg"; fi; \
> if [ -n "$$missing" ]; then \
>   echo "Static link requested, but these static libraries were not found:$$missing"; \
>   echo "Checked directories: $(STATIC_LIB_DIRS)"; \
>   echo "Install the matching *-dev packages that ship .a archives, or use LINK_MODE=shared."; \
>   exit 1; \
> fi

$(EXE): $(PRE_LINK_CHECKS) $(OBJS)
> $(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
> rm -f $(EXE) $(OBJS)
> rm -rf dist

