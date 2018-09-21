link_libraries(${FreeImage_LIBRARIES})

if(Qt4_FOUND)

  INCLUDE(${QT_USE_FILE})


  QT4_WRAP_CPP(dispsaliency_HEADERS_MOC
    ${CMAKE_HOME_DIRECTORY}/mln/qt/imageviewer.hpp
    dispsaliency.hpp
    )

  QT4_GENERATE_MOC(${CMAKE_HOME_DIRECTORY}/mln/qt/imageviewer.hxx
    ${CMAKE_CURRENT_BINARY_DIR}/imageviewer.moc.cpp)

  add_executable(dispsaliency dispsaliency.cpp
    ${dispsaliency_HEADERS_MOC}
    ${CMAKE_HOME_DIRECTORY}/mln/qt/qtimage.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/imageviewer.moc.cpp
    )

  target_link_libraries(dispsaliency ${QT_QTMAIN_LIBRARY} ${QT_LIBRARIES})

endif(Qt4_FOUND)

add_executable(saliency saliency.cpp
  ${CMAKE_HOME_DIRECTORY}/apps/attributes/MSERArgparser.cpp
  ${CMAKE_HOME_DIRECTORY}/apps/attributes/meaningfullnessArgparser.cpp)

target_link_libraries(saliency ${Boost_PROGRAM_OPTIONS_LIBRARY})

