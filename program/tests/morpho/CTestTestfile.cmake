# CMake generated Testfile for 
# Source directory: /lrde/home/movn/Documents/code/code_edwin/tests/morpho
# Build directory: /lrde/home/movn/Documents/code/code_edwin/program/tests/morpho
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(maxtree "maxtree")
add_test(saturate "saturate")
add_test(dilate "dilate")
add_test(erode "erode")
add_test(gradient "gradient")
add_test(opening "opening")
add_test(extinction "extinction")
subdirs(datastruct)
subdirs(maxtree)
subdirs(alphatree)
subdirs(component_tree)
