link_libraries(${FreeImage_LIBRARIES})

add_executable(segmentation_gradient segmentation_gradient.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(segment_node_tos segment_node_tos.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)



target_link_libraries(segmentation_gradient ${TBB_LIBRARIES})
target_link_libraries(segment_node_tos ${TBB_LIBRARIES})





