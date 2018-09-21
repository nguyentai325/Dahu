project(tests)
# add link flag to test framework
link_libraries(${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_definitions(-DBOOST_TEST_DYN_LINK)
add_definitions(-DMLN_IMG_PATH="$(CMAKE_SOURCE_DIR)/img/")

SUBDIRS(core draw io colors morpho accu labeling graph graphcut transform)

