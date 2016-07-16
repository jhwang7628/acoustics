# =============================================================================
#   Modal Viewer
# =============================================================================
# This file is adapted from tools/tetviewer/CMakeLists.txt
#
#
# ===================
# Dependencies
# ===================

    if (NOT USE_QT5) 
        message(FATAL_ERROR "Need to use Qt5 to build this target.")
    endif (NOT USE_QT5)
        
    # package for opengl and glut
    find_package(OpenGL REQUIRED)
    find_package(GLUT REQUIRED)
    find_package(Qt5Widgets REQUIRED)
    # check the existence of libQGLViewer
    find_package(QGLViewer 1.0 REQUIRED)

# ===================
# Settings
# ===================
    set(QT_USE_QTXML 1)
    set(QT_USE_QTOPENGL 1)

# ===================
# Definitions
# ===================
add_definitions(-DUSE_RECORDER) # needed for testing recorder

# ===================
# Target
# ===================
    set(CMAKE_INCLUDE_CURRENT_DIR ON)

    add_executable(modalviewer
        main.cpp
        ${UI_SRC}
        ${IO_SRC}
        ${GENERIC_SRC}
        ${GENERAL_TESTING_SRCS} 
        ${LINEARALGEBRA_SRC} 
        ${DEFORMABLE_SRC} 
        ${RIGID_SRC}
        ${RIGID_IO_SRC} 
        ${GEOMETRY_SRC}
        ${LOGGING_SRC}
        ${LINEARALGEBRA_BLAS_LAPACK_SRC}  
        ${SNDGEN_SRC}
        ${MODAL_MODEL_SRC}
        ${WAVESOLVER_SRC}
        ${MESH_UTILS_SRC}
        ${FIELD_SRC}                      
        ${DISTANCEFIELD_SRC}              
        ${PARSER_SRC} 
        ${GTS_GEOMETRY_SRC}               
        ${CGAL_GEOMETRY_SRC} 
        ${UTILS_SRC}                      
        ${TRANSFER_SRC}                   
        ${WAVELET_TRANSFER_SRC} 
        ${MATH_SRC} 
        )

    # target_include_directories(modalviewer
    #     PUBLIC ${GLUT_INCLUDE_DIR}
    #     PUBLIC ${OPENGL_INCLUDE_DIR}
    #     PUBLIC ${QGLViewer_INCLUDE_DIR}
    #     PUBLIC ${Boost_INCLUDE_DIR}
    #     )

    target_link_libraries(modalviewer
        ${GLUT_LIBRARIES}       
        ${QT_LIBRARIES}         
        ${QGLViewer_LIBRARY}
        ${OPENGL_LIBRARIES}     
        rt
        ${MKL_LIBRARY}
        ${BLAS_LAPACK_LIBS}             
        ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}               
        ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}                  
        ${LibConfig_LIBRARY}
        boost_thread
        boost_system                    
        boost_timer
        boost_filesystem
        boost_regex
        )

    set_target_properties(modalviewer PROPERTIES AUTOMOC TRUE)
    set_target_properties(modalviewer PROPERTIES AUTOMOC_MOC_OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)
    set_target_properties(modalviewer PROPERTIES AUTOUIC TRUE)
    set_target_properties(modalviewer PROPERTIES AUTORCC TRUE)
    qt5_use_modules(modalviewer Xml OpenGL)

    set_target_properties(modalviewer PROPERTIES LINK_FLAGS ${LINK_FLAGS})
