link_libraries(${FreeImage_LIBRARIES})
add_definitions(-DMLN_IMG_PATH="$(CMAKE_SOURCE_DIR)/img/")


add_executable(freeimage freeimage.cpp)
add_core_test(freeimage freeimage)

