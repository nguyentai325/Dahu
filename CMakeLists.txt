cmake_minimum_required(VERSION 2.8)
project(MILENA)

# Add FindTBB directory to CMake's module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

find_package(Boost 1.42 COMPONENTS program_options unit_test_framework REQUIRED)
find_package(FreeImage REQUIRED)
find_package(TBB)
find_package(Qt4)

if(TBB_FOUND)
        MESSAGE (STATUS "TBB found in ${TBB_INSTALL_DIR} + ${TBB_LIBRARY_DIRS}")
else()
        MESSAGE (FATAL "TBB not found.")
endif()



# add source directory to find mln path
include_directories("${PROJECT_SOURCE_DIR}")
include_directories(${TBB_INCLUDE_DIRS})
include_directories(SYSTEM ${Boost_INCLUDE_DIR})
link_directories(${TBB_LIBRARY_DIRS})
link_libraries(${TBB_LIBRARIES} rt)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wextra -Wall -std=c++0x")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -W -Wextra -Wall -g -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -W -Wextra -Wall -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -W -Wextra -Wall -O2 -g -ggdb")


# default buidling type
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()


# Testing
ENABLE_TESTING()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
macro(add_core_test Testname Exename)
add_test(${ARGV})
add_dependencies(check ${Exename})
endmacro(add_core_test)

#### Configuration of the Compiler ###
# set BOOST to use variadic macros
add_definitions(-D BOOST_PP_VARIADICS=1)



# Subdirectories
add_subdirectory(tests EXCLUDE_FROM_ALL)
add_subdirectory(bench EXCLUDE_FROM_ALL)
add_subdirectory(apps EXCLUDE_FROM_ALL)
add_subdirectory(cli EXCLUDE_FROM_ALL)

set(CPACK_SOURCE_GENERATOR "TBZ2")
set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_CURRENT_BINARY_DIR}/;/\\\\.git/;.*~;#.*#")
INCLUDE(CPack)
add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
