link_libraries(${FreeImage_LIBRARIES})

add_executable(mumford_shah_on_tree mumford_shah_on_tree.cpp)
add_executable(mumford_shah_on_tree_16 mumford_shah_on_tree.cpp)


add_executable(mumford_shah_on_tree_full mumford_shah_on_tree_full.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)


add_executable(dahudistance dahudistance.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)



add_executable(dahu_distance_updown_border dahu_distance_updown_border.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_distance_updown_4_border dahu_distance_updown_4_border.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(tree_color_test tree_color_test.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_distance_updown_4_borders_test dahu_distance_updown_4_borders_test.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_distance_updown_border_test dahu_distance_updown_border_test.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_slic_full dahu_slic_full.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)


add_executable(dahu_slic_full_4_bors dahu_slic_full_4_bors.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(test_path test_path.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(mumford_shah_slic mumford_shah_slic.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)


set_target_properties(mumford_shah_on_tree_16 PROPERTIES COMPILE_DEFINITIONS MLN_INPUT_VALUE_TYPE=mln::rgb16)


target_link_libraries(dahudistance ${TBB_LIBRARIES})
target_link_libraries(dahu_distance_updown_border ${TBB_LIBRARIES})
target_link_libraries(dahu_distance_updown_border_test ${TBB_LIBRARIES})
target_link_libraries(tree_color_test ${TBB_LIBRARIES})
target_link_libraries(dahu_distance_updown_4_border ${TBB_LIBRARIES})
target_link_libraries(dahu_distance_updown_4_borders_test ${TBB_LIBRARIES})
target_link_libraries(dahu_slic_full ${TBB_LIBRARIES})
target_link_libraries(dahu_slic_full_4_bors ${TBB_LIBRARIES})
target_link_libraries(test_path ${TBB_LIBRARIES})
target_link_libraries(mumford_shah_on_tree_full ${TBB_LIBRARIES})
target_link_libraries(mumford_shah_slic ${TBB_LIBRARIES})
