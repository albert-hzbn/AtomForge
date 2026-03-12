#
# Cross Platform Makefile
#

EXE = atoms_editor

IMGUI_DIR = ../

SRC_DIRS = . camera graphics math ui io

# --------------------------------------------------
# SOURCES
# --------------------------------------------------

SOURCES = main.cpp \
          camera/Camera.cpp \
          graphics/SphereMesh.cpp \
          graphics/Shader.cpp \
          graphics/ShadowMap.cpp \
          math/MathUtils.cpp \
          ui/Controls.cpp \
          io/StructureLoader.cpp

SOURCES += $(IMGUI_DIR)/imgui/imgui.cpp \
           $(IMGUI_DIR)/imgui/imgui_draw.cpp \
           $(IMGUI_DIR)/imgui/imgui_tables.cpp \
           $(IMGUI_DIR)/imgui/imgui_widgets.cpp

SOURCES += $(IMGUI_DIR)/backends/imgui_impl_glfw.cpp \
           $(IMGUI_DIR)/backends/imgui_impl_opengl3.cpp

OBJS = $(SOURCES:.cpp=.o)

UNAME_S := $(shell uname -s)

# --------------------------------------------------
# COMPILER FLAGS
# --------------------------------------------------
CXXFLAGS = -std=c++11 -O2
CXXFLAGS += -I. 
CXXFLAGS += -I$(IMGUI_DIR)/imgui
CXXFLAGS += -I$(IMGUI_DIR)/backends
CXXFLAGS += -Icamera -Igraphics -Imath -Iui -Iio
CXXFLAGS += -g -Wall -Wformat


# --------------------------------------------------
# PLATFORM FLAGS
# --------------------------------------------------

ifeq ($(UNAME_S), Linux)

ECHO_MESSAGE = "Linux"

CXXFLAGS += `pkg-config --cflags glfw3`
CXXFLAGS += `pkg-config --cflags openbabel-3`

LIBS += `pkg-config --libs openbabel-3`
LIBS += -lGLEW
LIBS += `pkg-config --libs glfw3`
LIBS += -lGL -ldl -lpthread

endif

# --------------------------------------------------
# BUILD RULES
# --------------------------------------------------

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


$(IMGUI_DIR)/imgui/%.o: $(IMGUI_DIR)/imgui/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(IMGUI_DIR)/backends/%.o: $(IMGUI_DIR)/backends/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
    
all: $(EXE)
	@echo Build complete for $(ECHO_MESSAGE)

$(EXE): $(OBJS)
	$(CXX) -o $@ $^ $(LIBS)



clean:
	rm -f $(EXE) $(OBJS) *.o