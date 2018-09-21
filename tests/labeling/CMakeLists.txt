add_executable(blobs blobs.cpp)
add_executable(rag   rag.cpp)

add_core_test(blobs blobs)
add_core_test(rag   rag)
