add_definitions(-DUSE_RECORDER)

# =============================================================================
#   Builds a rigid-body solver GUI
# =============================================================================
if ( BUILD_RIGIDSIM_GUI )
    # ===================
    # Dependencies
    # ===================
    
        # package for opengl and glut
        find_package(OpenGL REQUIRED)
        find_package(GLUT REQUIRED)
        find_package(LibConfig 1.0 REQUIRED)

        if(USE_QT5)
            find_package(Qt5Widgets REQUIRED)
        else(USE_QT5)
            find_package(Qt4 REQUIRED)
        endif(USE_QT5)

        # check the existence of libQGLViewer
        find_package(QGLViewer 1.0 REQUIRED)

        set(RIGID_SIM_GUI_APP_SRCS
            RigidSim.h    		RigidSim.cpp
            RigidSim.ui 		images.qrc

            RigidCanvas.h 					RigidCanvas.cpp
            demo/DemoDropObjs.h 			demo/DemoDropObjs.cpp
            demo/DemoDropObjsWithFixed.h 	demo/DemoDropObjsWithFixed.cpp
            demo/DemoShakingPiggy.h 		demo/DemoShakingPiggy.cpp
            demo/DemoPlanarCollision.h      demo/DemoPlanarCollision.cpp
            )

    # ===================
    # Settings
    # ===================
        set(QT_USE_QTXML 1)
        set(QT_USE_QTOPENGL 1)

    # ===================
    # Target
    # ===================

        set(CMAKE_INCLUDE_CURRENT_DIR ON)

        add_executable(rigidsim-gui
        	${RIGID_SIM_GUI_APP_SRCS}
            ${GENERIC_SRC}                  ${DEFORMABLE_SRC}
			${IO_SRC}                       ${RIGID_SRC}
			${LOGGING_SRC}                  ${LINEARALGEBRA_SRC}
            ${RIGID_IO_SRC}
            ${MODAL_MODEL_SRC}
            ${SNDGEN_SRC}
            ${RIGID_PARSER_SRC}
            )

        target_include_directories(rigidsim-gui
            PUBLIC ${GLUT_INCLUDE_DIR}
            PUBLIC ${OPENGL_INCLUDE_DIR}
            PUBLIC ${QGLViewer_INCLUDE_DIR}
            PUBLIC ${Boost_INCLUDE_DIR}
            )

        target_link_libraries(rigidsim-gui
            ${GLUT_LIBRARIES}       ${QT_LIBRARIES}         
            ${QGLViewer_LIBRARY}	${OPENGL_LIBRARIES}
            ${LibConfig_LIBRARY}	${BLAS_LAPACK_LIBS}
            ${MKL_LIBRARY}
            )

        set_target_properties(rigidsim-gui PROPERTIES AUTOMOC TRUE)
        set_target_properties(rigidsim-gui PROPERTIES AUTOMOC_MOC_OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)
        set_target_properties(rigidsim-gui PROPERTIES AUTOUIC TRUE)
        set_target_properties(rigidsim-gui PROPERTIES AUTORCC TRUE)

        if(USE_QT5)
            qt5_use_modules(rigidsim-gui Xml OpenGL)
        endif(USE_QT5)

        set_target_properties(rigidsim-gui PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( BUILD_RIGIDSIM_GUI ) 
