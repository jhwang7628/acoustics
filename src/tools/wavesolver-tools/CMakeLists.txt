# =============================================================================
#   Builds tools for processing FDTD wavesolver data
# =============================================================================
if ( BUILD_WAVESOLVER_TOOLS )

    # ===================
    # Definitions
    # ===================
    add_definitions(-DUSE_RECORDER) # needed for testing recorder
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

    # ===================
    # Target
    # ===================
    add_executable(convert-dat-vtk
        ${GRID_SOURCES}
        ${GENERIC_SRC}                    ${IO_SRC}
        ${LOGGING_SRC}                    ${LINEARALGEBRA_SRC}
        ${DISTANCEFIELD_SRC}              ${MATH_SRC}
        ${LINEARALGEBRA_BLAS_LAPACK_SRC}  ${GEOMETRY_SRC}
        ${FIELD_SRC}                      
        ${UTILS_SRC}                      ${PARSER_SRC}
        ${RIGID_PARSER_SRC}               ${LIBCONFIG_PARSER_SRC}
        ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
        ${WAVESOLVER_SRC}
        ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC}
        ${SPARSE_LINEAR_ALGEBRA_SRC} 
        ${LEASTSQUARESURFACE_SOURCES}
        ${RIGID_SRC} 
        ${RIGID_IO_SRC}
        ${IMPULSE_SERIES_SRC} 
        ${MODAL_MODEL_SRC}
        ${DEFORMABLE_SRC}
        ${MESH_UTILS_SRC}
        convert_dat_vtk.cpp) 

    add_executable(convert-dat-vtk-nonuniform
        convert_dat_vtk_nonuniform.cpp) 

    target_link_libraries(convert-dat-vtk
        ${MKL_LIBRARY}
        ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}                  
        ${LibConfig_LIBRARY}
        ${VLFEAT_LIBRARY}
        boost_thread
        boost_timer
        boost_system                    
        boost_filesystem
        boost_regex) 

    target_link_libraries(convert-dat-vtk-nonuniform
        boost_thread
        boost_timer
        boost_system                    
        boost_filesystem
        boost_regex) 
endif ( BUILD_WAVESOLVER_TOOLS )
