link_libraries(${FreeImage_LIBRARIES})


add_executable(stat stat.cpp)

add_executable(MSER-cli MSER-cli.cpp attributes.cpp)
add_executable(meaningfullness-cli meaningfullness-cli.cpp curvature.cpp attributes.cpp)
target_link_libraries(MSER-cli ${Boost_PROGRAM_OPTIONS_LIBRARY})
target_link_libraries(meaningfullness-cli ${Boost_PROGRAM_OPTIONS_LIBRARY})