# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /lrde/home/movn/Documents/code/code_edwin

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /lrde/home/movn/Documents/code/code_edwin/program

# Include any dependencies generated for this target.
include apps/g2/CMakeFiles/satmaxtree-cli.dir/depend.make

# Include the progress variables for this target.
include apps/g2/CMakeFiles/satmaxtree-cli.dir/progress.make

# Include the compile flags for this target's objects.
include apps/g2/CMakeFiles/satmaxtree-cli.dir/flags.make

apps/g2/CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.o: apps/g2/CMakeFiles/satmaxtree-cli.dir/flags.make
apps/g2/CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.o: ../apps/g2/satmaxtree-cli.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /lrde/home/movn/Documents/code/code_edwin/program/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/g2/CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.o"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/g2 && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.o -c /lrde/home/movn/Documents/code/code_edwin/apps/g2/satmaxtree-cli.cpp

apps/g2/CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.i"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/g2 && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /lrde/home/movn/Documents/code/code_edwin/apps/g2/satmaxtree-cli.cpp > CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.i

apps/g2/CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.s"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/g2 && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /lrde/home/movn/Documents/code/code_edwin/apps/g2/satmaxtree-cli.cpp -o CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.s

apps/g2/CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.o.requires:
.PHONY : apps/g2/CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.o.requires

apps/g2/CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.o.provides: apps/g2/CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.o.requires
	$(MAKE) -f apps/g2/CMakeFiles/satmaxtree-cli.dir/build.make apps/g2/CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.o.provides.build
.PHONY : apps/g2/CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.o.provides

apps/g2/CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.o.provides.build: apps/g2/CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.o

# Object files for target satmaxtree-cli
satmaxtree__cli_OBJECTS = \
"CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.o"

# External object files for target satmaxtree-cli
satmaxtree__cli_EXTERNAL_OBJECTS = \
"/lrde/home/movn/Documents/code/code_edwin/program/apps/g2/CMakeFiles/g2-tools.dir/satmaxtree.cpp.o" \
"/lrde/home/movn/Documents/code/code_edwin/program/apps/g2/CMakeFiles/g2-tools.dir/compute_ctos.cpp.o" \
"/lrde/home/movn/Documents/code/code_edwin/program/apps/g2/CMakeFiles/g2-tools.dir/routines.cpp.o" \
"/lrde/home/movn/Documents/code/code_edwin/program/apps/g2/CMakeFiles/g2-tools.dir/compute_g2.cpp.o"

apps/g2/satmaxtree-cli: apps/g2/CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.o
apps/g2/satmaxtree-cli: apps/g2/CMakeFiles/g2-tools.dir/satmaxtree.cpp.o
apps/g2/satmaxtree-cli: apps/g2/CMakeFiles/g2-tools.dir/compute_ctos.cpp.o
apps/g2/satmaxtree-cli: apps/g2/CMakeFiles/g2-tools.dir/routines.cpp.o
apps/g2/satmaxtree-cli: apps/g2/CMakeFiles/g2-tools.dir/compute_g2.cpp.o
apps/g2/satmaxtree-cli: apps/g2/CMakeFiles/satmaxtree-cli.dir/build.make
apps/g2/satmaxtree-cli: /usr/lib/libtbb.so
apps/g2/satmaxtree-cli: /usr/lib/libtbbmalloc.so
apps/g2/satmaxtree-cli: /usr/lib/libfreeimage.so
apps/g2/satmaxtree-cli: apps/g2/CMakeFiles/satmaxtree-cli.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable satmaxtree-cli"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/g2 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/satmaxtree-cli.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/g2/CMakeFiles/satmaxtree-cli.dir/build: apps/g2/satmaxtree-cli
.PHONY : apps/g2/CMakeFiles/satmaxtree-cli.dir/build

apps/g2/CMakeFiles/satmaxtree-cli.dir/requires: apps/g2/CMakeFiles/satmaxtree-cli.dir/satmaxtree-cli.cpp.o.requires
.PHONY : apps/g2/CMakeFiles/satmaxtree-cli.dir/requires

apps/g2/CMakeFiles/satmaxtree-cli.dir/clean:
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/g2 && $(CMAKE_COMMAND) -P CMakeFiles/satmaxtree-cli.dir/cmake_clean.cmake
.PHONY : apps/g2/CMakeFiles/satmaxtree-cli.dir/clean

apps/g2/CMakeFiles/satmaxtree-cli.dir/depend:
	cd /lrde/home/movn/Documents/code/code_edwin/program && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lrde/home/movn/Documents/code/code_edwin /lrde/home/movn/Documents/code/code_edwin/apps/g2 /lrde/home/movn/Documents/code/code_edwin/program /lrde/home/movn/Documents/code/code_edwin/program/apps/g2 /lrde/home/movn/Documents/code/code_edwin/program/apps/g2/CMakeFiles/satmaxtree-cli.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/g2/CMakeFiles/satmaxtree-cli.dir/depend
