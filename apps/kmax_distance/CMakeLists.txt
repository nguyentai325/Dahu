link_libraries(${FreeImage_LIBRARIES})

add_executable(kmax_distance_border kmax_distance_border.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)





target_link_libraries(kmax_distance_border ${TBB_LIBRARIES})






