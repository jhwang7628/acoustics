# =============================================================================
#   Tetmesh Viewer
# =============================================================================
if ( BUILD_TETVIEWER )
    # ===================
    # Dependencies
    # ===================
    
        # package for opengl and glut
        find_package(OpenGL REQUIRED)
        find_package(GLUT REQUIRED)

        if(USE_QT5)
            find_package(Qt5Widgets REQUIRED)
        else(USE_QT5)
            find_package(Qt4 REQUIRED)
        endif(USE_QT5)

        # check the existence of libQGLViewer
        find_package(QGLViewer 1.0 REQUIRED)

        set(TETVIEWER_APP_SRCS
            TetViewerFrame.h    TetViewerFrame.cpp
            TetViewerCanvas.h   TetViewerCanvas.cpp
            tetviewer.ui
            images.qrc
            )

    # ===================
    # Settings
    # ===================
        set(QT_USE_QTXML 1)
        set(QT_USE_QTOPENGL 1)

    # ===================
    # Target
    # ===================
        if(USE_QT5)
            include_directories(${QT_USE_FILE})
        else(USE_QT5)
            include(${QT_USE_FILE})
        endif(USE_QT5)
        
        set(CMAKE_INCLUDE_CURRENT_DIR ON)

        add_executable(tetviewer
            ${TETVIEWER_APP_SRCS}
            ${GENERIC_SRC}          ${IO_SRC}           ${DEFORMABLE_SRC}
            ${LINEARALGEBRA_SRC})

        target_include_directories(tetviewer
            # PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
            PUBLIC ${GLUT_INCLUDE_DIR}
            PUBLIC ${OPENGL_INCLUDE_DIR}
            PUBLIC ${QGLViewer_INCLUDE_DIR}
            PUBLIC ${Boost_INCLUDE_DIR}
            )

        target_link_libraries(tetviewer
            ${GLUT_LIBRARIES}       ${QT_LIBRARIES}         ${QGLViewer_LIBRARY}
            ${ICC_LIBS}             ${OPENGL_LIBRARIES}     rt
            )

        set_target_properties(tetviewer PROPERTIES AUTOMOC TRUE)
        set_target_properties(tetviewer PROPERTIES AUTOMOC_MOC_OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)
        set_target_properties(tetviewer PROPERTIES AUTOUIC TRUE)
        set_target_properties(tetviewer PROPERTIES AUTORCC TRUE)

        if(USE_QT5)
            qt5_use_modules(tetviewer Xml OpenGL)
        endif(USE_QT5)

        set_target_properties(tetviewer PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( BUILD_TETVIEWER ) 