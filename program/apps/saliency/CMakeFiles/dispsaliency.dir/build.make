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
include apps/saliency/CMakeFiles/dispsaliency.dir/depend.make

# Include the progress variables for this target.
include apps/saliency/CMakeFiles/dispsaliency.dir/progress.make

# Include the compile flags for this target's objects.
include apps/saliency/CMakeFiles/dispsaliency.dir/flags.make

apps/saliency/__/__/mln/qt/moc_imageviewer.cxx: ../mln/qt/imageviewer.hpp
	$(CMAKE_COMMAND) -E cmake_progress_report /lrde/home/movn/Documents/code/code_edwin/program/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Generating __/__/mln/qt/moc_imageviewer.cxx"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/__/__/mln/qt && /usr/lib/x86_64-linux-gnu/qt4/bin/moc @/lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/__/__/mln/qt/moc_imageviewer.cxx_parameters

apps/saliency/moc_dispsaliency.cxx: ../apps/saliency/dispsaliency.hpp
	$(CMAKE_COMMAND) -E cmake_progress_report /lrde/home/movn/Documents/code/code_edwin/program/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Generating moc_dispsaliency.cxx"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/lib/x86_64-linux-gnu/qt4/bin/moc @/lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/moc_dispsaliency.cxx_parameters

apps/saliency/imageviewer.moc.cpp: ../mln/qt/imageviewer.hxx
	$(CMAKE_COMMAND) -E cmake_progress_report /lrde/home/movn/Documents/code/code_edwin/program/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Generating imageviewer.moc.cpp"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/lib/x86_64-linux-gnu/qt4/bin/moc @/lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/imageviewer.moc.cpp_parameters

apps/saliency/CMakeFiles/dispsaliency.dir/dispsaliency.cpp.o: apps/saliency/CMakeFiles/dispsaliency.dir/flags.make
apps/saliency/CMakeFiles/dispsaliency.dir/dispsaliency.cpp.o: ../apps/saliency/dispsaliency.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /lrde/home/movn/Documents/code/code_edwin/program/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/saliency/CMakeFiles/dispsaliency.dir/dispsaliency.cpp.o"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/dispsaliency.dir/dispsaliency.cpp.o -c /lrde/home/movn/Documents/code/code_edwin/apps/saliency/dispsaliency.cpp

apps/saliency/CMakeFiles/dispsaliency.dir/dispsaliency.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dispsaliency.dir/dispsaliency.cpp.i"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /lrde/home/movn/Documents/code/code_edwin/apps/saliency/dispsaliency.cpp > CMakeFiles/dispsaliency.dir/dispsaliency.cpp.i

apps/saliency/CMakeFiles/dispsaliency.dir/dispsaliency.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dispsaliency.dir/dispsaliency.cpp.s"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /lrde/home/movn/Documents/code/code_edwin/apps/saliency/dispsaliency.cpp -o CMakeFiles/dispsaliency.dir/dispsaliency.cpp.s

apps/saliency/CMakeFiles/dispsaliency.dir/dispsaliency.cpp.o.requires:
.PHONY : apps/saliency/CMakeFiles/dispsaliency.dir/dispsaliency.cpp.o.requires

apps/saliency/CMakeFiles/dispsaliency.dir/dispsaliency.cpp.o.provides: apps/saliency/CMakeFiles/dispsaliency.dir/dispsaliency.cpp.o.requires
	$(MAKE) -f apps/saliency/CMakeFiles/dispsaliency.dir/build.make apps/saliency/CMakeFiles/dispsaliency.dir/dispsaliency.cpp.o.provides.build
.PHONY : apps/saliency/CMakeFiles/dispsaliency.dir/dispsaliency.cpp.o.provides

apps/saliency/CMakeFiles/dispsaliency.dir/dispsaliency.cpp.o.provides.build: apps/saliency/CMakeFiles/dispsaliency.dir/dispsaliency.cpp.o

apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.o: apps/saliency/CMakeFiles/dispsaliency.dir/flags.make
apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.o: apps/saliency/__/__/mln/qt/moc_imageviewer.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /lrde/home/movn/Documents/code/code_edwin/program/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.o"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.o -c /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/__/__/mln/qt/moc_imageviewer.cxx

apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.i"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/__/__/mln/qt/moc_imageviewer.cxx > CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.i

apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.s"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/__/__/mln/qt/moc_imageviewer.cxx -o CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.s

apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.o.requires:
.PHONY : apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.o.requires

apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.o.provides: apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.o.requires
	$(MAKE) -f apps/saliency/CMakeFiles/dispsaliency.dir/build.make apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.o.provides.build
.PHONY : apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.o.provides

apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.o.provides.build: apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.o

apps/saliency/CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.o: apps/saliency/CMakeFiles/dispsaliency.dir/flags.make
apps/saliency/CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.o: apps/saliency/moc_dispsaliency.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /lrde/home/movn/Documents/code/code_edwin/program/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/saliency/CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.o"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.o -c /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/moc_dispsaliency.cxx

apps/saliency/CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.i"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/moc_dispsaliency.cxx > CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.i

apps/saliency/CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.s"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/moc_dispsaliency.cxx -o CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.s

apps/saliency/CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.o.requires:
.PHONY : apps/saliency/CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.o.requires

apps/saliency/CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.o.provides: apps/saliency/CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.o.requires
	$(MAKE) -f apps/saliency/CMakeFiles/dispsaliency.dir/build.make apps/saliency/CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.o.provides.build
.PHONY : apps/saliency/CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.o.provides

apps/saliency/CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.o.provides.build: apps/saliency/CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.o

apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.o: apps/saliency/CMakeFiles/dispsaliency.dir/flags.make
apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.o: ../mln/qt/qtimage.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /lrde/home/movn/Documents/code/code_edwin/program/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.o"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.o -c /lrde/home/movn/Documents/code/code_edwin/mln/qt/qtimage.cpp

apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.i"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /lrde/home/movn/Documents/code/code_edwin/mln/qt/qtimage.cpp > CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.i

apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.s"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /lrde/home/movn/Documents/code/code_edwin/mln/qt/qtimage.cpp -o CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.s

apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.o.requires:
.PHONY : apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.o.requires

apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.o.provides: apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.o.requires
	$(MAKE) -f apps/saliency/CMakeFiles/dispsaliency.dir/build.make apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.o.provides.build
.PHONY : apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.o.provides

apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.o.provides.build: apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.o

apps/saliency/CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.o: apps/saliency/CMakeFiles/dispsaliency.dir/flags.make
apps/saliency/CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.o: apps/saliency/imageviewer.moc.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /lrde/home/movn/Documents/code/code_edwin/program/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object apps/saliency/CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.o"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.o -c /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/imageviewer.moc.cpp

apps/saliency/CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.i"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/imageviewer.moc.cpp > CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.i

apps/saliency/CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.s"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && /usr/local/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/imageviewer.moc.cpp -o CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.s

apps/saliency/CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.o.requires:
.PHONY : apps/saliency/CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.o.requires

apps/saliency/CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.o.provides: apps/saliency/CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.o.requires
	$(MAKE) -f apps/saliency/CMakeFiles/dispsaliency.dir/build.make apps/saliency/CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.o.provides.build
.PHONY : apps/saliency/CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.o.provides

apps/saliency/CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.o.provides.build: apps/saliency/CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.o

# Object files for target dispsaliency
dispsaliency_OBJECTS = \
"CMakeFiles/dispsaliency.dir/dispsaliency.cpp.o" \
"CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.o" \
"CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.o" \
"CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.o" \
"CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.o"

# External object files for target dispsaliency
dispsaliency_EXTERNAL_OBJECTS =

apps/saliency/dispsaliency: apps/saliency/CMakeFiles/dispsaliency.dir/dispsaliency.cpp.o
apps/saliency/dispsaliency: apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.o
apps/saliency/dispsaliency: apps/saliency/CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.o
apps/saliency/dispsaliency: apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.o
apps/saliency/dispsaliency: apps/saliency/CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.o
apps/saliency/dispsaliency: apps/saliency/CMakeFiles/dispsaliency.dir/build.make
apps/saliency/dispsaliency: /usr/lib/libtbb.so
apps/saliency/dispsaliency: /usr/lib/libtbbmalloc.so
apps/saliency/dispsaliency: /usr/lib/libfreeimage.so
apps/saliency/dispsaliency: /usr/lib/x86_64-linux-gnu/libQtGui.so
apps/saliency/dispsaliency: /usr/lib/x86_64-linux-gnu/libQtCore.so
apps/saliency/dispsaliency: apps/saliency/CMakeFiles/dispsaliency.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable dispsaliency"
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dispsaliency.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
apps/saliency/CMakeFiles/dispsaliency.dir/build: apps/saliency/dispsaliency
.PHONY : apps/saliency/CMakeFiles/dispsaliency.dir/build

apps/saliency/CMakeFiles/dispsaliency.dir/requires: apps/saliency/CMakeFiles/dispsaliency.dir/dispsaliency.cpp.o.requires
apps/saliency/CMakeFiles/dispsaliency.dir/requires: apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/moc_imageviewer.cxx.o.requires
apps/saliency/CMakeFiles/dispsaliency.dir/requires: apps/saliency/CMakeFiles/dispsaliency.dir/moc_dispsaliency.cxx.o.requires
apps/saliency/CMakeFiles/dispsaliency.dir/requires: apps/saliency/CMakeFiles/dispsaliency.dir/__/__/mln/qt/qtimage.cpp.o.requires
apps/saliency/CMakeFiles/dispsaliency.dir/requires: apps/saliency/CMakeFiles/dispsaliency.dir/imageviewer.moc.cpp.o.requires
.PHONY : apps/saliency/CMakeFiles/dispsaliency.dir/requires

apps/saliency/CMakeFiles/dispsaliency.dir/clean:
	cd /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency && $(CMAKE_COMMAND) -P CMakeFiles/dispsaliency.dir/cmake_clean.cmake
.PHONY : apps/saliency/CMakeFiles/dispsaliency.dir/clean

apps/saliency/CMakeFiles/dispsaliency.dir/depend: apps/saliency/__/__/mln/qt/moc_imageviewer.cxx
apps/saliency/CMakeFiles/dispsaliency.dir/depend: apps/saliency/moc_dispsaliency.cxx
apps/saliency/CMakeFiles/dispsaliency.dir/depend: apps/saliency/imageviewer.moc.cpp
	cd /lrde/home/movn/Documents/code/code_edwin/program && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /lrde/home/movn/Documents/code/code_edwin /lrde/home/movn/Documents/code/code_edwin/apps/saliency /lrde/home/movn/Documents/code/code_edwin/program /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency /lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/CMakeFiles/dispsaliency.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : apps/saliency/CMakeFiles/dispsaliency.dir/depend

