link_libraries(${FreeImage_LIBRARIES})

add_executable(seed_GOP seed_GOP.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_seed dahu_seed.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_seed_seg dahu_seed_seg.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_updown dahu_updown.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_updown_folder dahu_updown_folder.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)


add_executable(dahu_updown_folder2015 dahu_updown_folder2015.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)


add_executable(dahu_seed_interative dahu_seed_interative.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_seed_saliency dahu_seed_saliency.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(shortestpath shortestpath.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(test_shortest_path test_shortest_path.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(test_shortest_path1 test_shortest_path1.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_lab dahu_lab.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_seed_queue dahu_seed_queue.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

target_link_libraries(seed_GOP ${TBB_LIBRARIES})
target_link_libraries(dahu_seed ${TBB_LIBRARIES})
target_link_libraries(dahu_seed_seg ${TBB_LIBRARIES})
target_link_libraries(dahu_updown ${TBB_LIBRARIES})
target_link_libraries(dahu_updown_folder ${TBB_LIBRARIES})
target_link_libraries(dahu_updown_folder2015 ${TBB_LIBRARIES})
target_link_libraries(dahu_seed_interative ${TBB_LIBRARIES})
target_link_libraries(dahu_seed_saliency ${TBB_LIBRARIES})
target_link_libraries(shortestpath ${TBB_LIBRARIES})
target_link_libraries(test_shortest_path ${TBB_LIBRARIES})
target_link_libraries(test_shortest_path1 ${TBB_LIBRARIES})
target_link_libraries(dahu_lab ${TBB_LIBRARIES})
target_link_libraries(dahu_seed_queue ${TBB_LIBRARIES})


