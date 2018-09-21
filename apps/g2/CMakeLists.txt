link_libraries(${FreeImage_LIBRARIES})


add_library(g2-tools OBJECT satmaxtree.cpp compute_ctos.cpp routines.cpp compute_g2.cpp)
add_executable(compute_ctos-cli  compute_ctos-cli.cpp $<TARGET_OBJECTS:g2-tools>)
target_link_libraries(compute_ctos-cli ${Boost_PROGRAM_OPTIONS_LIBRARY})

#add_library(g2-tools-16 OBJECT routines.cpp compute_g2.cpp)
add_executable(compute_ctos_16-cli compute_ctos-cli.cpp $<TARGET_OBJECTS:g2-tools>)
set_target_properties(compute_ctos_16-cli
  PROPERTIES COMPILE_DEFINITIONS MLN_INPUT_VALUE_TYPE=mln::rgb16)
target_link_libraries(compute_ctos_16-cli ${Boost_PROGRAM_OPTIONS_LIBRARY})



add_executable(satmaxtree-cli satmaxtree-cli.cpp $<TARGET_OBJECTS:g2-tools>)
