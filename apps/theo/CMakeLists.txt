add_executable(ctos_grainfilter ctos_grainfilter.cpp $<TARGET_OBJECTS:g2-tools>)
target_link_libraries(ctos_grainfilter ${FreeImage_LIBRARIES})
