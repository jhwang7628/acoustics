# =============================================================================
#   FDTD Acoustic Simulator Viewer
# =============================================================================
# This file is adapted from tools/modalviewer/CMakeLists.txt
#
# ===================
# Dependencies
# ===================
if (BUILD_FDTD_ACOUSTIC_SIMULATOR_VIEWER)
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

    add_executable(fdtd-acoustic-simulator-viewer
        fdtd-acoustic-simulator-viewer.cpp
        ${DISTANCEFIELD_SRC}              
        ${WAVESOLVER_SRC} 
        ${MESH_UTILS_SRC}
        ${GENERAL_TESTING_SRCS} 
        ${FDTD_ACOUSTIC_SIMULATOR_TESTING_SRCS} 
        ${GENERIC_SRC}                    
        ${IO_SRC} 
        ${LOGGING_SRC}                    
        ${LINEARALGEBRA_SRC} 
        ${MATH_SRC} 
        ${LINEARALGEBRA_BLAS_LAPACK_SRC}  
        ${GEOMETRY_SRC} 
        ${FIELD_SRC}                      
        ${UTILS_SRC}                      
        ${UTILS_OPENGL_SRC}                      
        ${PARSER_SRC} 
        ${GTS_GEOMETRY_SRC}               
        ${CGAL_GEOMETRY_SRC} 
        ${TRANSFER_SRC}                   
        ${WAVELET_TRANSFER_SRC} 
        ${SPARSE_LINEAR_ALGEBRA_SRC} 
        ${MODAL_MODEL_SRC}
        ${MODAL_MATERIAL_SRC}
        ${DEFORMABLE_SRC} 
        ${RIGID_SRC}
        ${RIGID_IO_SRC} 
        ${IMPULSE_SERIES_SRC} 
        ${RIGID_PARSER_SRC}
        ${UI_SRC}
        ${BEM_SRC}
        ${BEM_IO_SRC}
        ${SND_IO_SRC}
        ${INTERP_SRC}
        )

    target_link_libraries(fdtd-acoustic-simulator-viewer
        ${GLUT_LIBRARIES}       
        ${QT_LIBRARIES}         
        ${QGLViewer_LIBRARY}
        ${OPENGL_LIBRARIES}     
        ${MKL_LIBRARY}
        ${BLAS_LAPACK_LIBS}             
        ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}               
        ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}                  
        ${LibConfig_LIBRARY}
        ${VLFEAT_LIBRARY}
        sndfile
        boost_thread
        boost_system                    
        boost_timer
        boost_filesystem
        boost_regex
        boost_program_options
        )

    set_target_properties(fdtd-acoustic-simulator-viewer PROPERTIES AUTOMOC TRUE)
    set_target_properties(fdtd-acoustic-simulator-viewer PROPERTIES AUTOMOC_MOC_OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)
    set_target_properties(fdtd-acoustic-simulator-viewer PROPERTIES AUTOUIC TRUE)
    set_target_properties(fdtd-acoustic-simulator-viewer PROPERTIES AUTORCC TRUE)
    qt5_use_modules(fdtd-acoustic-simulator-viewer Xml OpenGL)

    set_target_properties(fdtd-acoustic-simulator-viewer PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif (BUILD_FDTD_ACOUSTIC_SIMULATOR_VIEWER)
