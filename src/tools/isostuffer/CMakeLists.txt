# =============================================================================
#   Isostuffer GUI
# =============================================================================
if ( BUILD_ISOSTUFFER_GUI )
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

        set(ISOSTUFFER_GUI_SRCS
            GuiMain.cpp
            ../../ui/IsoStufferFrame.h      ../../ui/IsoStufferFrame.cpp
            ../../ui/IsoStufferFrame.ui

            ../../ui/IsoStufferCanvas.h     ../../ui/IsoStufferCanvas.cpp
            
            ../../ui/ParamsDialog.h
            ../../ui/Params.ui

            ../../io/TglMeshReader.hpp      ../../io/TglMeshReader.cpp 
            
            ../../ui/images.qrc
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

        add_executable(isostuffer-gui
            ${CONFIG_SOURCES}       ${ISOSTUFFER_GUI_SRCS}         ${ISOSTUFFER_SOURCES}
            ${UTILS_SRC}            ${LINEARALGEBRA_SRC}           ${GENERIC_SRC}
            ${IO_SRC})

        target_include_directories(isostuffer-gui
            # PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
            PUBLIC ${GLUT_INCLUDE_DIR}
            PUBLIC ${OPENGL_INCLUDE_DIR}
            PUBLIC ${QGLViewer_INCLUDE_DIR}
            )

        target_link_libraries(isostuffer-gui
            ${GLUT_LIBRARIES}       ${QT_LIBRARIES}         ${QGLViewer_LIBRARY}
            ${ICC_LIBS}             ${OPENGL_LIBRARIES}     rt
            ${MKL_LIBRARY})

        set_target_properties(isostuffer-gui PROPERTIES AUTOMOC TRUE)
        set_target_properties(isostuffer-gui PROPERTIES AUTOMOC_MOC_OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)
        set_target_properties(isostuffer-gui PROPERTIES AUTOUIC TRUE)
        set_target_properties(isostuffer-gui PROPERTIES AUTORCC TRUE)

        if(USE_QT5)
            qt5_use_modules(isostuffer-gui Xml OpenGL)
        endif(USE_QT5)

        set_target_properties(isostuffer-gui PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( BUILD_ISOSTUFFER_GUI ) 


# =============================================================================
#   Command-line isostuffer
# =============================================================================
if ( BUILD_ISOSTUFFER )
    # ===================
    # Target
    # ===================
        add_executable(isostuffer
            ../../io/TglMeshReader.hpp    ../../io/TglMeshReader.cpp
            main.cpp
            ${CONFIG_SOURCES}           ${ISOSTUFFER_SOURCES}       ${UTILS_SRC}
            ${LINEARALGEBRA_SRC}               ${GENERIC_SRC}              ${IO_SRC})

        target_link_libraries(isostuffer ${ICC_LIBS}        boost_program_options
                            rt)

        set_target_properties(isostuffer PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( BUILD_ISOSTUFFER )