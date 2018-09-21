link_libraries(${FreeImage_LIBRARIES})

add_executable(slic slic.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahudistance_slic dahudistance_slic.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(tos tos.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(dahu_slic dahu_slic.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)

add_executable(mtos mtos.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/satmaxtree.cpp
  ${CMAKE_SOURCE_DIR}/apps/g2/compute_g2.cpp)


target_link_libraries(slic ${TBB_LIBRARIES})
target_link_libraries(dahudistance_slic ${TBB_LIBRARIES})
target_link_libraries(tos ${TBB_LIBRARIES})
target_link_libraries(dahu_slic ${TBB_LIBRARIES})
target_link_libraries(mtos ${TBB_LIBRARIES})
