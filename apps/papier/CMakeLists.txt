link_libraries(${FreeImage_LIBRARIES})

add_executable(dahu_gray_tree dahu_gray_tree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(fast_mbd fast_mbd.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(fast_waterflow_mbd fast_waterflow_mbd.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(fast_dahu_gray fast_dahu_gray.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_color_tree dahu_color_tree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(fast_waterflow_mbd_color fast_waterflow_mbd_color.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(fast_mbd_color fast_mbd_color.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(fast_dahu_color fast_dahu_color.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(roi roi.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_mtos_vs_tos dahu_mtos_vs_tos.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_lab_tree dahu_lab_tree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_mtos_vs_tos_lab dahu_mtos_vs_tos_lab.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_bordure_vs_4_borders dahu_bordure_vs_4_borders.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_path_cost dahu_path_cost.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(distance_2_points distance_2_points.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(fast_water_flow_mbd_color_shortest_path fast_water_flow_mbd_color_shortest_path.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(shortest_path_tree_vs_image shortest_path_tree_vs_image.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(eval_speed eval_speed.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(inter_intra_distance inter_intra_distance.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_vs_fast_mbd dahu_vs_fast_mbd.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(mbd_mst mbd_mst.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(mbd_mst_color mbd_mst_color.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_color_tree_queue dahu_color_tree_queue.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_color_tree_priorityqueue dahu_color_tree_priorityqueue.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_color_tree_4borders_priority dahu_color_tree_4borders_priority.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(shortest_path_mbd_methods shortest_path_mbd_methods.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(synthetic_image_segmentation_errors synthetic_image_segmentation_errors.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(random_seeded_dahu random_seeded_dahu.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)


target_link_libraries(dahu_gray_tree ${TBB_LIBRARIES})
target_link_libraries(fast_mbd ${TBB_LIBRARIES})
target_link_libraries(fast_waterflow_mbd ${TBB_LIBRARIES})
target_link_libraries(fast_dahu_gray ${TBB_LIBRARIES})
target_link_libraries(dahu_color_tree ${TBB_LIBRARIES})
target_link_libraries(fast_waterflow_mbd_color ${TBB_LIBRARIES})
target_link_libraries(fast_mbd_color ${TBB_LIBRARIES})
target_link_libraries(fast_dahu_color ${TBB_LIBRARIES})
target_link_libraries(roi ${TBB_LIBRARIES})
target_link_libraries(dahu_mtos_vs_tos ${TBB_LIBRARIES})
target_link_libraries(dahu_lab_tree ${TBB_LIBRARIES})
target_link_libraries(dahu_mtos_vs_tos_lab ${TBB_LIBRARIES})
target_link_libraries(dahu_bordure_vs_4_borders ${TBB_LIBRARIES})
target_link_libraries(dahu_path_cost ${TBB_LIBRARIES})
target_link_libraries(distance_2_points ${TBB_LIBRARIES})
target_link_libraries(fast_water_flow_mbd_color_shortest_path ${TBB_LIBRARIES})	
target_link_libraries(shortest_path_tree_vs_image ${TBB_LIBRARIES})	
target_link_libraries(eval_speed ${TBB_LIBRARIES})	
target_link_libraries(inter_intra_distance ${TBB_LIBRARIES})	
target_link_libraries(dahu_vs_fast_mbd ${TBB_LIBRARIES})	
target_link_libraries(mbd_mst ${TBB_LIBRARIES})	
target_link_libraries(mbd_mst_color ${TBB_LIBRARIES})	
target_link_libraries(dahu_color_tree_queue ${TBB_LIBRARIES})	
target_link_libraries(dahu_color_tree_priorityqueue ${TBB_LIBRARIES})
target_link_libraries(dahu_color_tree_4borders_priority ${TBB_LIBRARIES})
target_link_libraries(shortest_path_mbd_methods ${TBB_LIBRARIES})
target_link_libraries(synthetic_image_segmentation_errors ${TBB_LIBRARIES})
target_link_libraries(random_seeded_dahu ${TBB_LIBRARIES})


