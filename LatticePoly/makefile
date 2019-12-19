# Compiler
CXX         := g++

# Linker
LD          := g++

# Target binary program
TARGET      := lat

# Version
VSN_MAJ     := 1
VSN_MIN     := 0

# VTK library installation path
VTK_PATH    := /Users/mtortora/apps/homebrew

# VTK version
VTK_VSN     := 8.2

# Directory paths
SRCDIR      := src
BUILDDIR    := build
TARGETDIR   := bin
INCDIR      := include
DATDIR      := data

# File extensions
SRCEXT      := cpp
HDREXT      := hpp
DEPEXT      := d
OBJEXT      := o

# Flags, libraries and includes
INC         := -I$(INCDIR) -isystem$(VTK_PATH)/include/vtk-$(VTK_VSN)
LIB         := -Wl,-rpath,$(VTK_PATH)/lib -L$(VTK_PATH)/lib -lvtkCommonCore-$(VTK_VSN) -lvtkCommonDataModel-$(VTK_VSN) -lvtkCommonExecutionModel-$(VTK_VSN) -lvtkFiltersSources-$(VTK_VSN) -lvtkIOXML-$(VTK_VSN)
FPATHS      := -D__DPATH__=$(CURDIR)/$(DATDIR)
FEXTRA      := -D__VERSION_MAJOR__=$(VSN_MAJ) -D__VERSION_MINOR__=$(VSN_MIN) -DNDEBUG
CXXFLAGS    := -std=c++0x -Wshadow -Wall -Wextra -O3 -ffast-math

# Object files lists
SOURCES     := $(shell find $(SRCDIR) -type f -iname "*.$(SRCEXT)")
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))

# Default make
all: dirs $(TARGET)
	@echo -e "\033[1;31mCompiled $(TARGET) v$(VSN_MAJ).$(VSN_MIN)\033[0m"
	@echo -e "\033[1;37m*****************\033[0m"

# Make directories
dirs: FORCE
	@echo -e "\033[1;37m*****************\033[0m"
	@mkdir -p $(BUILDDIR)
	@mkdir -p $(TARGETDIR)

# Run
run: clean all
	@LD_LIBRARY_PATH=$(VTK_PATH)/lib DYLD_LIBRARY_PATH=$(VTK_PATH)/lib ./$(TARGETDIR)/$(TARGET)

# Clean objects/files
clean:
	@echo -e "\033[1;31mCleaning object and data files...\033[0m"
	@$(RM) -rf $(OBJECTS)
	@$(RM) -rf $(DATDIR)/vtk/*
	
# Clean directories
cleaner: clean
	@echo -e "\033[1;31mCleaning directories...\033[0m"
	@$(RM) -rf $(BUILDDIR)
	@$(RM) -rf $(TARGETDIR)
	
# Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

# Link
$(TARGET): $(OBJECTS)
	@echo -e "\033[1;31mLinking object files...\033[0m"
	@$(LD) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)

# Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@echo -e "\033[1;34mBuilding: \033[0m" $<
	@mkdir -p $(dir $@)
	@$(CXX) $(CXXFLAGS) $(FPATHS) $(FEXTRA) $(INC) -c -o $@ $<
	@$(CXX) $(CXXFLAGS) $(FPATHS) $(FEXTRA) $(INC) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@$(RM) -f $(BUILDDIR)/$*.$(DEPEXT).tmp


# Force dependency - include for rules to be run at every call
FORCE:

# Non-file targets
.PHONY: all dirs run clean cleaner
