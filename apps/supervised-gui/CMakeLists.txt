if(QT4_FOUND)

  INCLUDE(${QT_USE_FILE})

  link_libraries(${FreeImage_LIBRARIES})


  QT4_WRAP_CPP(supervised_MOC
    ${CMAKE_HOME_DIRECTORY}/mln/qt/imageviewer.hpp
    brush.hpp)

  QT4_GENERATE_MOC(${CMAKE_HOME_DIRECTORY}/mln/qt/imageviewer.hxx
    ${CMAKE_CURRENT_BINARY_DIR}/imageviewer.moc.cpp)

  add_executable(supervised-cli  supervised-cli.cpp)
  add_executable(supervised-cli-multi  supervised-cli-multi.cpp)
  add_executable(automatic  automatic.cpp)
  add_executable(nettoyage  nettoyage.cpp)
  add_executable(generate_res  generate_res.cpp)

  add_executable(supervised-gui-3c  supervised-gui-3c.cpp brush.cpp
    ${supervised_MOC}
    ${CMAKE_HOME_DIRECTORY}/mln/qt/qtimage.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/imageviewer.moc.cpp)


  target_link_libraries(supervised-gui-3c ${QT_QTMAIN_LIBRARY} ${QT_LIBRARIES})
  target_link_libraries(supervised-cli ${QT_QTMAIN_LIBRARY} ${QT_LIBRARIES})
  target_link_libraries(nettoyage ${Boost_PROGRAM_OPTIONS_LIBRARIES})

  set_property(TARGET supervised-cli-multi APPEND PROPERTY INCLUDE_DIRECTORIES "/work/carlinet/src/eigen/")

endif(QT4_FOUND)



