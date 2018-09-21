link_libraries(${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_definitions(-DBOOST_TEST_DYN_LINK)

add_subdirectory(tos)
add_subdirectory(attributes)