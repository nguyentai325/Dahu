if (${QT4_FOUND})

  INCLUDE(${QT_USE_FILE})

  link_libraries(${FreeImage_LIBRARIES})


  QT4_WRAP_CPP(hierarchical_seg-gui_MOC
    ${CMAKE_HOME_DIRECTORY}/mln/qt/imageviewer.hpp
    brush.hpp)

  QT4_GENERATE_MOC(${CMAKE_HOME_DIRECTORY}/mln/qt/imageviewer.hxx
    ${CMAKE_CURRENT_BINARY_DIR}/imageviewer.moc.cpp)


  add_executable(hierarchical_seg-gui  hierachical_seg-gui.cpp brush.cpp
    ${hierarchical_seg-gui_MOC}
    ${CMAKE_HOME_DIRECTORY}/mln/qt/qtimage.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/imageviewer.moc.cpp)

  add_executable(hierarchical_seg-cli  hierachical_seg-cli.cpp)

  target_link_libraries(hierarchical_seg-gui ${QT_QTMAIN_LIBRARY} ${QT_LIBRARIES})

endif(${QT4_FOUND})

