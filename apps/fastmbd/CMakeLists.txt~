link_libraries(${FreeImage_LIBRARIES})

add_executable(fastmbd fastmbd.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(waterflow_mbd waterflow_mbd.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)


add_executable(queue_waterflow_mbd queue_waterflow_mbd.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)



add_executable(fastmbd_folder fastmbd_folder.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(waterflow_mbd_folder waterflow_mbd_folder.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(fast_dahu fast_dahu.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)


add_executable(fast_dahu_border fast_dahu_border.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)


add_executable(fast_dahu_border_color fast_dahu_border_color.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)


add_executable(fast_dahu_border_color_folder fast_dahu_border_color_folder.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)


target_link_libraries(fastmbd ${TBB_LIBRARIES})
target_link_libraries(waterflow_mbd ${TBB_LIBRARIES})
target_link_libraries(fastmbd_folder ${TBB_LIBRARIES})
target_link_libraries(waterflow_mbd_folder ${TBB_LIBRARIES})
target_link_libraries(fast_dahu ${TBB_LIBRARIES})
target_link_libraries(fast_dahu_border ${TBB_LIBRARIES})
target_link_libraries(queue_waterflow_mbd ${TBB_LIBRARIES})
target_link_libraries(fast_dahu_border_color ${TBB_LIBRARIES})

target_link_libraries(fast_dahu_border_color_folder ${TBB_LIBRARIES})

