##
##  Makefile
##  LatticePoly
##
##  Created by mtortora on 02/01/2022.
##  Copyright © 2022 ENS Lyon. All rights reserved.
##


# Compiler
CXX        := g++

# Linker
LD         := g++

# CMake executable
CMAKE      := cmake

# Target binary program
TARGET     := lat

# Version
VSN_MAJ    := 2
VSN_MIN    := 0

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
# VTKDIR     := $(LIBDIR)/VTK/build
HDF5DIR    := $(LIBDIR)/hdf5/build

# Required VTK modules
MODULES    := CommonComputationalGeometry CommonCore CommonDataModel CommonExecutionModel CommonMath CommonMisc CommonSystem CommonTransforms FiltersCore FiltersGeneral FiltersSources IOCore IOXML IOXMLParser

# VTKMODS    := $(foreach mod,$(MODULES),-DVTK_MODULE_ENABLE_VTK_$(mod)=YES)
# VTKLIBS    := $(foreach mod,$(MODULES),-lvtk$(mod))

HDF5LIBS   := -lhdf5 -lhdf5_cpp


# Flags, libraries and includes
INC        := -I$(INCDIR) -isystem$(INCDIR)/hdf5 #-isystem$(INCDIR)/vtk
LIB        := -Wl,-rpath,$(CURDIR)/$(LIBDIR) -L$(LIBDIR) $(HDF5LIBS) #$(VTKLIBS) -lvtksys -lvtkzlib -lvtkdoubleconversion -lvtkexpat -lvtkloguru -lvtklz4 -lvtklzma
FPATHS     := -D__DPATH__=$(CURDIR)/$(DATDIR)
FEXTRA     := -D__VERSION_MAJOR__=$(VSN_MAJ) -D__VERSION_MINOR__=$(VSN_MIN) -DNDEBUG
CXXFLAGS   := -std=c++11 -Wshadow -Wall -Wextra -O3 -ffast-math -flto
CMAKEFLAGS := -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$(CURDIR)  #-DVTK_VERSIONED_INSTALL=OFF -DVTK_BUILD_TESTING=OFF -DVTK_GROUP_ENABLE_Rendering=NO -DVTK_GROUP_ENABLE_StandAlone=NO $(VTKMODS)
CMAKEFHDF5 := -G "Unix Makefiles" -DCMAKE_BUILD_TYPE:STRING=Release -DBUILD_SHARED_LIBS:BOOL=ON -DBUILD_TESTING:BOOL=OFF -DHDF5_BUILD_TOOLS:BOOL=ON -DHDF5_BUILD_CPP_LIB:BOOL=ON

# Object files lists
SOURCES    := $(shell find $(SRCDIR) -type f -iname "*.$(SRCEXT)")
OBJECTS    := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))

# VTK objects and directories for full clean
# VTKCLEAN   := $(VTKDIR) $(LIBDIR)/cmake $(TARGETDIR)/vtk*
# HDF5CLEAN  := $(HDF5DIR) $(TARGETDIR)/hdf5*
# Default make
all: dirs $(TARGET)
	@echo "\033[1;31mCompiled $(TARGET) v$(VSN_MAJ).$(VSN_MIN)\033[0m"
	@echo "\033[1;37m*****************\033[0m"

# Make directories
dirs:
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
	@$(RM) -rf $(HDF5CLEAN)
# @$(RM) -rf $(VTKCLEAN)

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
	@$(LD) -o $(TARGETDIR)/$(TARGET) $(CXXFLAGS) $^ $(LIB)
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

# Setup (minimal) VTK build
# $(VTKDIR)/Makefile:
# 	@echo "\033[1;34mBuilding VTK library...\033[0m"
# 	@mkdir -p $(VTKDIR)
# 	@cd $(VTKDIR) && $(CMAKE) $(CMAKEFLAGS) ..

# Install VTK headers/dynamic library files
# libvtk: $(VTKDIR)/Makefile
# 	@$(MAKE) -C $(VTKDIR) -j4
# 	@$(MAKE) -C $(VTKDIR) install
# 	@echo "\033[1;34mSuccessfully compiled VTK library\033[0m"

libhdf5:
	@echo "\033[1;34mBuilding HDF5 library...\033[0m"
#	mkdir -p $(HDF5DIR)
# cd $(HDF5DIR) && $(CMAKE) $(CMAKEFHDF5) ..
#	cd $(HDF5DIR) && $(CMAKE) --build . --config Release
#	cd $(HDF5DIR) && cpack -C Release CPackConfig.cmake
#	cd $(HDF5DIR) && ./HDF5-1.14.5-Linux.sh
	mv $(HDF5DIR)/HDF5-1.14.5-Linux/HDF_Group/HDF5/1.14.5/$(TARGETDIR)/* $(TARGETDIR)/
	mkdir -p $(INCDIR)/hdf5
	mv $(HDF5DIR)/HDF5-1.14.5-Linux/HDF_Group/HDF5/1.14.5/$(INCDIR)/* $(INCDIR)/hdf5/
	mkdir -p $(LIBDIR)/cmake/hdf5
	mv $(HDF5DIR)/HDF5-1.14.5-Linux/HDF_Group/HDF5/1.14.5/cmake/* $(LIBDIR)/cmake/hdf5/
	mkdir -p share/hdf5
	mv $(HDF5DIR)/HDF5-1.14.5-Linux/HDF_Group/HDF5/1.14.5/share/* share/hdf5/
	mkdir -p lib/pkgconfig
	mv $(HDF5DIR)/HDF5-1.14.5-Linux/HDF_Group/HDF5/1.14.5/$(LIBDIR)/pkgconfig/* $(LIBDIR)/pkgconfig/
	rm -rf $(HDF5DIR)/HDF5-1.14.5-Linux/HDF_Group/HDF5/1.14.5/$(LIBDIR)/pkgconfig
	mv $(HDF5DIR)/HDF5-1.14.5-Linux/HDF_Group/HDF5/1.14.5/$(LIBDIR)/* $(LIBDIR)/
	echo "\033[1;34mSuccessfully compiled HDF5 library\033[0m"


# Non-file targets
.PHONY: all dirs run clean cleanlib cleaner libhdf5 #libvtk
