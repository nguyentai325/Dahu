find_package(Boost REQUIRED COMPONENTS program_options filesystem system)

include_directories(clattice ${Boost_INCLUDE_DIRS})

add_executable(clattice clattice.cpp)
target_link_libraries(clattice
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${FreeImage_LIBRARIES}
)