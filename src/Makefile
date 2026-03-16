# ==============================================================================
# Configuration — edit these sections to add/remove files or change flags
# ==============================================================================

EXE   = atoms_editor
CXX   = g++
IMGUI = ../

# ------------------------------------------------------------------------------
# Project sources — one file per line for easy diff / addition
# ------------------------------------------------------------------------------

APP_SOURCES = \
    main.cpp                              \
    ElementData.cpp

CAMERA_SOURCES = \
    camera/Camera.cpp

GRAPHICS_SOURCES = \
    graphics/Renderer.cpp                 \
    graphics/Picking.cpp                  \
    graphics/SceneBuffers.cpp             \
    graphics/Shader.cpp                   \
    graphics/ShadowMap.cpp                \
    graphics/SphereMesh.cpp               \
    graphics/StructureInstanceBuilder.cpp

MATH_SOURCES = \
    math/MathUtils.cpp

UI_SOURCES = \
    ui/AtomContextMenu.cpp                \
    ui/EditMenuDialogs.cpp                \
    ui/FileBrowser.cpp                    \
    ui/ImGuiSetup.cpp                     \
    ui/PeriodicTableDialog.cpp            \
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
          $(MATH_SOURCES) $(UI_SOURCES) $(IO_SOURCES) $(IMGUI_SOURCES)

# ------------------------------------------------------------------------------
# Compiler flags
# ------------------------------------------------------------------------------

CXXFLAGS = \
    -std=c++11 -O2 -g -Wall -Wformat \
    -I. -I$(IMGUI)/imgui -I$(IMGUI)/backends \
    -Icamera -Igraphics -Imath -Iui -Iio

# ------------------------------------------------------------------------------
# Platform-specific flags
# ------------------------------------------------------------------------------

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S), Linux)
  ECHO_MESSAGE = Linux
  CXXFLAGS    += $(shell pkg-config --cflags glfw3 openbabel-3)
  LIBS         = $(shell pkg-config --libs   glfw3 openbabel-3) -lGLEW -lGL -ldl -lpthread
endif

# ==============================================================================
# Build rules — no need to edit below this line
# ==============================================================================

OBJS = $(SOURCES:.cpp=.o)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

all: $(EXE)
	@echo Build complete for $(ECHO_MESSAGE)

$(EXE): $(OBJS)
	$(CXX) -o $@ $^ $(LIBS)

clean:
	rm -f $(EXE) $(OBJS)

