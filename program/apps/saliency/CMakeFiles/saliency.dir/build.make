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
include apps/saliency/CMakeFiles/saliency.dir/depend.make

# Include the progress variables for this target.
include apps/saliency/CMakeFiles/saliency.dir/progress.make

# Include the compile flags for this target's objects.
include apps/saliency/CMakeFiles/saliency.dir/flags.make

apps/saliency/CMakeFiles/saliency.dir/saliency.cpp.o: apps/saliency/CMakeFiles/saliency.dir/flags.make
apps/saliency/CMakeFiles/saliency.dir/saliency.cpp.o: ../apps/saliency/saliency.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /lrde/home/movn/Documents/code/code_edwin/program/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/saliency/CMakeFiles/saliency.dir/saliency.cpp.o"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/saliency.dir/saliency.cpp.o -c /lrde/home/movn/Documents/code/code_edwin/apps/saliency/saliency.cpp

apps/saliency/CMakeFiles/saliency.dir/saliency.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/saliency.dir/saliency.cpp.i"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /lrde/home/movn/Documents/code/code_edwin/apps/saliency/saliency.cpp > CMakeFiles/saliency.dir/saliency.cpp.i

apps/saliency/CMakeFiles/saliency.dir/saliency.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/saliency.dir/saliency.cpp.s"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /lrde/home/movn/Documents/code/code_edwin/apps/saliency/saliency.cpp -o CMakeFiles/saliency.dir/saliency.cpp.s

apps/saliency/CMakeFiles/saliency.dir/saliency.cpp.o.requires:
.PHONY : apps/saliency/CMakeFiles/saliency.dir/saliency.cpp.o.requires

apps/saliency/CMakeFiles/saliency.dir/saliency.cpp.o.provides: apps/saliency/CMakeFiles/saliency.dir/saliency.cpp.o.requires
	$(MAKE) -f apps/saliency/CMakeFiles/saliency.dir/build.make apps/saliency/CMakeFiles/saliency.dir/saliency.cpp.o.provides.build
.PHONY : apps/saliency/CMakeFiles/saliency.dir/saliency.cpp.o.provides

apps/saliency/CMakeFiles/saliency.dir/saliency.cpp.o.provides.build: apps/saliency/CMakeFiles/saliency.dir/saliency.cpp.o

apps/saliency/CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.o: apps/saliency/CMakeFiles/saliency.dir/flags.make
apps/saliency/CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.o: ../apps/attributes/MSERArgparser.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /lrde/home/movn/Documents/code/code_edwin/program/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/saliency/CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.o"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.o -c /lrde/home/movn/Documents/code/code_edwin/apps/attributes/MSERArgparser.cpp

apps/saliency/CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.i"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /lrde/home/movn/Documents/code/code_edwin/apps/attributes/MSERArgparser.cpp > CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.i

apps/saliency/CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.s"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /lrde/home/movn/Documents/code/code_edwin/apps/attributes/MSERArgparser.cpp -o CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.s

apps/saliency/CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.o.requires:
.PHONY : apps/saliency/CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.o.requires

apps/saliency/CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.o.provides: apps/saliency/CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.o.requires
	$(MAKE) -f apps/saliency/CMakeFiles/saliency.dir/build.make apps/saliency/CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.o.provides.build
.PHONY : apps/saliency/CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.o.provides

apps/saliency/CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.o.provides.build: apps/saliency/CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.o

apps/saliency/CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.o: apps/saliency/CMakeFiles/saliency.dir/flags.make
apps/saliency/CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.o: ../apps/attributes/meaningfullnessArgparser.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /lrde/home/movn/Documents/code/code_edwin/program/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/saliency/CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.o"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.o -c /lrde/home/movn/Documents/code/code_edwin/apps/attributes/meaningfullnessArgparser.cpp

apps/saliency/CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.i"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /lrde/home/movn/Documents/code/code_edwin/apps/attributes/meaningfullnessArgparser.cpp > CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.i

apps/saliency/CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.s"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /lrde/home/movn/Documents/code/code_edwin/apps/attributes/meaningfullnessArgparser.cpp -o CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.s

apps/saliency/CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.o.requires:
.PHONY : apps/saliency/CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.o.requires

apps/saliency/CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.o.provides: apps/saliency/CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.o.requires
	$(MAKE) -f apps/saliency/CMakeFiles/saliency.dir/build.make apps/saliency/CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.o.provides.build
.PHONY : apps/saliency/CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.o.provides

apps/saliency/CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.o.provides.build: apps/saliency/CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.o

# Object files for target saliency
saliency_OBJECTS = \
"CMakeFiles/saliency.dir/saliency.cpp.o" \
"CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.o" \
"CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.o"

# External object files for target saliency
saliency_EXTERNAL_OBJECTS =

apps/saliency/saliency: apps/saliency/CMakeFiles/saliency.dir/saliency.cpp.o
apps/saliency/saliency: apps/saliency/CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.o
apps/saliency/saliency: apps/saliency/CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.o
apps/saliency/saliency: apps/saliency/CMakeFiles/saliency.dir/build.make
apps/saliency/saliency: /usr/lib/libtbb.so
apps/saliency/saliency: /usr/lib/libtbbmalloc.so
apps/saliency/saliency: /usr/lib/libfreeimage.so
apps/saliency/saliency: /lrde/home/movn/local/include/boost/lib/libboost_program_options.so
apps/saliency/saliency: apps/saliency/CMakeFiles/saliency.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable saliency"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/saliency.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/saliency/CMakeFiles/saliency.dir/build: apps/saliency/saliency
.PHONY : apps/saliency/CMakeFiles/saliency.dir/build

apps/saliency/CMakeFiles/saliency.dir/requires: apps/saliency/CMakeFiles/saliency.dir/saliency.cpp.o.requires
apps/saliency/CMakeFiles/saliency.dir/requires: apps/saliency/CMakeFiles/saliency.dir/__/attributes/MSERArgparser.cpp.o.requires
apps/saliency/CMakeFiles/saliency.dir/requires: apps/saliency/CMakeFiles/saliency.dir/__/attributes/meaningfullnessArgparser.cpp.o.requires
.PHONY : apps/saliency/CMakeFiles/saliency.dir/requires

apps/saliency/CMakeFiles/saliency.dir/clean:
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && $(CMAKE_COMMAND) -P CMakeFiles/saliency.dir/cmake_clean.cmake
.PHONY : apps/saliency/CMakeFiles/saliency.dir/clean

apps/saliency/CMakeFiles/saliency.dir/depend:
	cd /lrde/home/movn/Documents/code/code_edwin/program && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lrde/home/movn/Documents/code/code_edwin /lrde/home/movn/Documents/code/code_edwin/apps/saliency /lrde/home/movn/Documents/code/code_edwin/program /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/CMakeFiles/saliency.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/saliency/CMakeFiles/saliency.dir/depend
