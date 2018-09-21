# Install script for directory: /lrde/home/movn/Documents/code/code_edwin/apps

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/lrde/home/movn/Documents/code/code_edwin/program/apps/tests/cmake_install.cmake")
  include("/lrde/home/movn/Documents/code/code_edwin/program/apps/imview/cmake_install.cmake")
  include("/lrde/home/movn/Documents/code/code_edwin/program/apps/maxtree_comparison/cmake_install.cmake")
  include("/lrde/home/movn/Documents/code/code_edwin/program/apps/attributes/cmake_install.cmake")
  include("/lrde/home/movn/Documents/code/code_edwin/program/apps/simplification/cmake_install.cmake")
  include("/lrde/home/movn/Documents/code/code_edwin/program/apps/clattice/cmake_install.cmake")
  include("/lrde/home/movn/Documents/code/code_edwin/program/apps/g2/cmake_install.cmake")
  include("/lrde/home/movn/Documents/code/code_edwin/program/apps/mumford_shah_on_tree/cmake_install.cmake")
  include("/lrde/home/movn/Documents/code/code_edwin/program/apps/misc/cmake_install.cmake")
  include("/lrde/home/movn/Documents/code/code_edwin/program/apps/theo/cmake_install.cmake")
  include("/lrde/home/movn/Documents/code/code_edwin/program/apps/tosgui/cmake_install.cmake")
  include("/lrde/home/movn/Documents/code/code_edwin/program/apps/saliency/cmake_install.cmake")
  include("/lrde/home/movn/Documents/code/code_edwin/program/apps/supervised-gui/cmake_install.cmake")
  include("/lrde/home/movn/Documents/code/code_edwin/program/apps/hierachical_seg-gui/cmake_install.cmake")
  include("/lrde/home/movn/Documents/code/code_edwin/program/apps/hyperspec/cmake_install.cmake")

endif()

