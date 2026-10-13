# eviVLP chemical-modification build.
#
# This file is tracked in the repo at chemical_modification/src/Makefile and
# synced into the flat working folder by sync_repo.py, so the build recipe
# lives under version control instead of in someone's shell history.
#
# The GUI must be rebuilt after any change to chem_GUI.cpp, FileProcessor.cpp
# or MoleculeViewer.cpp; the IC generators after any change to their own source.
# Build them on macOS - they link against the local wxWidgets and OpenGL.

CXX      := g++
CXXSTD   := -std=c++17
WX_CXX   := $(shell wx-config --cxxflags)
WX_LIBS  := $(shell wx-config --libs gl,core,base)

GUI_SRC  := chem_GUI.cpp FileProcessor.cpp MoleculeViewer.cpp

.PHONY: all gui ic clean
all: gui ic

# The GUI. Exactly the command that produced the working binary.
gui: GUI
GUI: $(GUI_SRC) FileProcessor.h MoleculeViewer.h
	$(CXX) $(CXXSTD) $(GUI_SRC) $(WX_CXX) $(WX_LIBS) -o GUI

# Internal-coordinate table generators. Standalone programs, each with its own
# main() and no wxWidgets dependency.
ic: IC_naa IC_ntrm_NPC IC_table_one
IC_naa: IC_table_naa.cpp
	$(CXX) $(CXXSTD) IC_table_naa.cpp -o IC_naa
IC_ntrm_NPC: IC_table_ntrm_NPC.cpp
	$(CXX) $(CXXSTD) IC_table_ntrm_NPC.cpp -o IC_ntrm_NPC
# The drug's IC table (RESI ONE). The binary used to be built by hand with
# -o IC_table_one from a source called IC_table_fixed.cpp, which is why nothing
# recorded how to rebuild it; the source is now named after its output.
IC_table_one: IC_table_one.cpp
	$(CXX) $(CXXSTD) IC_table_one.cpp -o IC_table_one

clean:
	rm -f GUI IC_naa IC_ntrm_NPC IC_table_one
