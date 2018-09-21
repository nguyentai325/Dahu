add_executable(io		io.cpp)
add_core_test(io		io)
target_link_libraries(io ${FreeImage_LIBRARIES})
#set_target_properties(io PROPERTIES OUTPUT_NAME maxtree_old)
