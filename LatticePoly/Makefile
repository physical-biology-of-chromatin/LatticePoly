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

# HDF5 library installation path
HDF5DIR    := $(LIBDIR)/hdf5

HDF5LIBS   := -lhdf5 -lhdf5_cpp

# Flags, libraries and includes
INC        := -I$(INCDIR) -isystem$(INCDIR)/hdf5
LIB        := -Wl,-rpath,$(CURDIR)/$(LIBDIR) -L$(LIBDIR) $(HDF5LIBS)
FPATHS     := -D__DPATH__=$(CURDIR)/$(DATDIR)
FEXTRA     := -D__VERSION_MAJOR__=$(VSN_MAJ) -D__VERSION_MINOR__=$(VSN_MIN) -DNDEBUG
CXXFLAGS   := -std=c++11 -Wshadow -Wall -Wextra -O3 -ffast-math -flto
CMAKEFHDF5 := -G "Unix Makefiles" -DCMAKE_BUILD_TYPE:STRING=Release -DBUILD_SHARED_LIBS:BOOL=ON -DHDF5_BUILD_EXAMPLES:BOOL=OFF -DBUILD_TESTING:BOOL=OFF -DHDF5_BUILD_TOOLS:BOOL=ON -DHDF5_BUILD_CPP_LIB:BOOL=ON -DCMAKE_INSTALL_PREFIX=../../ -DHDF5_INSTALL_INCLUDE_DIR=include/hdf5 -DHDF5_INSTALL_CMAKE_DIR=$(HDF5DIR)/cmake

# Object files lists
SOURCES    := $(shell find $(SRCDIR) -type f -iname "*.$(SRCEXT)")
OBJECTS    := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))

# HDF5 objects and directories for full clean
HDF5CLEAN  := $(BUILDDIR)/hdf5 $(LIBDIR)/cmake $(TARGETDIR)/hdf5*

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

# Clean directories
cleaner: clean cleanlib

# Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

# Link
$(TARGET): $(OBJECTS)
ifeq ($(wildcard $(LIBDIR)/hdf5/.),)
	@echo "\033[1;31mCould not locate HDF5 library - please run 'make libhdf5'\033[0m"
	@exit 1
else
	@echo "\033[1;31mLinking object files...\033[0m"
	@$(LD) -o $(TARGETDIR)/$(TARGET) $(CXXFLAGS) $^ $(LIB)
endif

# Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
ifeq ($(wildcard $(LIBDIR)/hdf5/.),)
	@echo "\033[1;31mCould not locate HDF5 library - please run 'make libhdf5'\033[0m"
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


libhdf5:
	@echo "\033[1;34mBuilding HDF5 library...\033[0m"
	@cd $(HDF5DIR) && git checkout hdf5-1_14_3
	@cd $(HDF5DIR) && $(CMAKE) $(CMAKEFHDF5) -B ../../build/hdf5 .
	@cd $(HDF5DIR) && $(CMAKE) --build ../../build/hdf5
	@cd $(HDF5DIR) && $(CMAKE) --install ../../build/hdf5
	@echo "\033[1;34mSuccessfully compiled HDF5 library\033[0m"

# Non-file targets
.PHONY: all dirs run clean cleanlib cleaner libhdf5
