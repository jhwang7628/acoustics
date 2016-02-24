# =============================================================================
#   Builds the unit testing tools
# =============================================================================
if (BUILD_UNIT_TESTING)
    # ===================
    # Dependencies
    # ===================

    set(GRID_TESTING_SRCS 
        test_GridWithObject.cpp)

    set(LEASTSQUARESURFACE_TESTING_SRCS 
        test_LeastSquareSurface.cpp)
    # set(TURBSOUND_INCLUDE_DIR
    #     /home/jui-hsien/code/turbsound_postfluid/src
    #     /home/jui-hsien/code/turbsound_postfluid/tools)


    # ===================
    # Settings
    # ===================
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


    # ===================
    # Target
    # ===================

    # include_directories(
    #     ${TURBSOUND_INCLUDE_DIR}
    #     )

    add_executable(unit_testing_leastSquareSurface
        ${LEASTSQUARESURFACE_TESTING_SRCS}
        ${LEASTSQUARESURFACE_SOURCES}
        )
    add_executable(unit_testing_grid
        ${GRID_TESTING_SRCS}
        ${GRID_SOURCES}
        ${GENERIC_SRC}                    ${IO_SRC}
        ${LOGGING_SRC}                    ${LINEARALGEBRA_SRC}
        ${DISTANCEFIELD_SRC}              ${MATH_SRC}
        ${LINEARALGEBRA_BLAS_LAPACK_SRC}  ${GEOMETRY_SRC}
        ${FIELD_SRC}                      ${WAVESOLVER_SRC}
        ${UTILS_SRC}                      ${PARSER_SRC}
        ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
        ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC}
        )

    target_link_libraries(unit_testing_grid
        ${MKL_LIBRARY}
        ${LEASTSQUARESURFACE_SOURCES}
        ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}                  
        boost_thread
        boost_system                    
        boost_regex) 

endif ( BUILD_UNIT_TESTING )
