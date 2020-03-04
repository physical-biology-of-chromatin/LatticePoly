# Compiler
CXX        := g++

# Linker
LD         := g++

# Target binary program
TARGET     := lat

# Version
VSN_MAJ    := 1
VSN_MIN    := 3

# Directory paths
SRCDIR     := src
BUILDDIR   := build
TARGETDIR  := bin
INCDIR     := include
LIBDIR     := lib
DATDIR     := data

# File extensions
SRCEXT     := cpp
HDREXT     := hpp
DEPEXT     := d
OBJEXT     := o

# VTK library installation path
VTKDIR     := $(LIBDIR)/VTK/build

# Required VTK modules
MODULES    := CommonComputationalGeometry CommonCore CommonDataModel CommonExecutionModel CommonMath CommonMisc CommonSystem CommonTransforms FiltersCore FiltersGeneral FiltersSources IOCore IOXML IOXMLParser

VTKMODS    := $(foreach mod,$(MODULES),-DVTK_MODULE_ENABLE_VTK_$(mod)=YES)
VTKLIBS    := $(foreach mod,$(MODULES),-lvtk$(mod))

# Flags, libraries and includes
INC        := -I$(INCDIR) -isystem$(INCDIR)/vtk
LIB        := -Wl,-rpath,$(LIBDIR) -L$(LIBDIR) $(VTKLIBS) -lvtksys -lvtkzlib -lvtkdoubleconversion -lvtkexpat -lvtkloguru -lvtklz4 -lvtklzma
FPATHS     := -D__DPATH__=$(CURDIR)/$(DATDIR)
FEXTRA     := -D__VERSION_MAJOR__=$(VSN_MAJ) -D__VERSION_MINOR__=$(VSN_MIN) -DNDEBUG
CXXFLAGS   := -std=c++0x -Wshadow -Wall -Wextra -O3 -ffast-math
CMAKEFLAGS := -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$(CURDIR) -DVTK_VERSIONED_INSTALL=OFF -DVTK_BUILD_TESTING=OFF -DVTK_GROUP_ENABLE_Rendering=NO -DVTK_GROUP_ENABLE_StandAlone=NO $(VTKMODS)

# Object files lists
SOURCES    := $(shell find $(SRCDIR) -type f -iname "*.$(SRCEXT)")
OBJECTS    := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))

# VTK objects and directories for full clean
VTKCLEAN   := $(VTKDIR) $(LIBDIR)/cmake $(TARGETDIR)/vtk*

# Default make
all: dirs $(TARGET)
	@echo "\033[1;31mCompiled $(TARGET) v$(VSN_MAJ).$(VSN_MIN)\033[0m"
	@echo "\033[1;37m*****************\033[0m"

# Make directories
dirs: FORCE
	@echo "\033[1;37m*****************\033[0m"
	@mkdir -p $(BUILDDIR)
	@mkdir -p $(TARGETDIR)

# Run
run: all
	@$(RM) -r $(DATDIR)/output
	@./$(TARGETDIR)/$(TARGET) $(DATDIR)/input.cfg

# Clean objects/files
clean:
	@echo "\033[1;31mCleaning object files...\033[0m"
	@$(RM) -rf $(OBJECTS)
	
# Clean library files
cleanlib:
	@echo "\033[1;31mCleaning library files...\033[0m"
	@$(RM) -rf $(VTKCLEAN)

# Clean directories
cleaner: clean cleanlib
	
# Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

# Link
$(TARGET): $(OBJECTS)
ifeq ($(wildcard $(INCDIR)/vtk/.),)
	@echo "\033[1;31mCould not locate VTK library - please run 'make libvtk'\033[0m"
	@exit 1
else
	@echo "\033[1;31mLinking object files...\033[0m"
	@$(LD) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)
endif

# Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
ifeq ($(wildcard $(INCDIR)/vtk/.),)
	@echo "\033[1;31mCould not locate VTK library - please run 'make libvtk'\033[0m"
	@exit 1
else
	@echo "\033[1;34mBuilding: \033[0m" $<
	@mkdir -p $(dir $@)
	@$(CXX) $(CXXFLAGS) $(FPATHS) $(FEXTRA) $(INC) -c -o $@ $<
	@$(CXX) $(CXXFLAGS) $(FPATHS) $(FEXTRA) $(INC) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@$(RM) -f $(BUILDDIR)/$*.$(DEPEXT).tmp
endif

# Build (minimal) VTK library
libvtk: $(VTKDIR)/Makefile
	@$(MAKE) -C $(VTKDIR) -j4
	@$(MAKE) -C $(VTKDIR) install
	@echo "\033[1;34mSuccessfully compiled VTK library\033[0m"

$(VTKDIR)/Makefile:
	@echo "\033[1;34mBuilding VTK library...\033[0m"
	@mkdir -p $(VTKDIR)
	@cd $(VTKDIR) && cmake $(CMAKEFLAGS) ..

# Force dependency - include for rules to be run at every call
FORCE:

# Non-file targets
.PHONY: all dirs run clean cleanlib cleaner libvtk
