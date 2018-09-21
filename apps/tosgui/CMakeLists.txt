if(Qt4_FOUND AND QWT_FOUND)

  INCLUDE(${QT_USE_FILE})
  include_directories(${QWT_INCLUDE_DIR})

  link_libraries(${FreeImage_LIBRARIES})


  QT4_WRAP_CPP(tosgui_HEADERS_MOC
    ${CMAKE_HOME_DIRECTORY}/mln/qt/imageviewer.hpp
    qattribute.hpp
    dispatcher.hpp
    plotwindow.hpp
    attribute.hpp
    )

  QT4_GENERATE_MOC(${CMAKE_HOME_DIRECTORY}/mln/qt/imageviewer.hxx
    ${CMAKE_CURRENT_BINARY_DIR}/imageviewer.moc.cpp)

  add_executable(tosgui tosgui.cpp qattribute.cpp dispatcher.cpp plotwindow.cpp attribute.cpp
    ${tosgui_HEADERS_MOC}
    ${CMAKE_HOME_DIRECTORY}/mln/qt/qtimage.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/imageviewer.moc.cpp
    )

  add_executable(toscgui toscgui.cpp qattribute.cpp dispatcher.cpp
    ${tosgui_HEADERS_MOC}
    ${CMAKE_HOME_DIRECTORY}/mln/qt/qtimage.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/imageviewer.moc.cpp
    )


  target_link_libraries(tosgui
    ${QWT_LIBRARY}
    ${QT_QTMAIN_LIBRARY}
    ${QT_LIBRARIES}
    )

  target_link_libraries(toscgui
    ${QWT_LIBRARY}
    ${QT_QTMAIN_LIBRARY}
    ${QT_LIBRARIES}
    )

endif(Qt4_FOUND AND QWT_FOUND)