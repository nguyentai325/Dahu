link_libraries(${FreeImage_LIBRARIES})
add_executable(simplify_gray simplify_gray.cpp)
add_executable(simplify_cluminance simplify_cluminance.cpp)
add_executable(simplify_color simplify_color.cpp ../tos/colorToSGrad.cpp)



