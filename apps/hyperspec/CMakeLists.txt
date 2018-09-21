add_executable(profiles_from_tos profiles_from_tos.cpp $<TARGET_OBJECTS:g2-tools>)
target_link_libraries(profiles_from_tos ${FreeImage_LIBRARIES})

add_executable(profiles_from_marginal_tos profiles_from_marginal_tos.cpp $<TARGET_OBJECTS:g2-tools>)
target_link_libraries(profiles_from_marginal_tos ${FreeImage_LIBRARIES})