# =============================================================================
#   Builds the tools for differentiating HRIR computed from wavesolver
# =============================================================================
if (BUILD_DIFFERENTIATE_INTERPOLATE_HRIR)
    # ===================
    # Dependencies
    # ===================
    set(DIFFERENTIATE_INTERPOLATE_HRIR_SRCS 
        DifferentiateInterpolateHRIR.cpp)
    set(TURBSOUND_INCLUDE_DIR
        /home/jui-hsien/code/turbsound_postfluid/src
        /home/jui-hsien/code/turbsound_postfluid/tools)


    # ===================
    # Settings
    # ===================
    set(CMAKE_CXX_FLAGS "-std=c++11")
    set(TARGET_NAME differentiate-interpolate-HRIR)


    # ===================
    # Target
    # ===================

    include_directories(
        ${TURBSOUND_INCLUDE_DIR}
        )

    add_executable(${TARGET_NAME}
        ${DIFFERENTIATE_INTERPOLATE_HRIR_SRCS}
        ${LEASTSQUARESURFACE_SOURCES}
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

    target_link_libraries(${TARGET_NAME}
        ${MKL_LIBRARY}
        ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}                  
        boost_thread
        boost_system                    
        boost_timer
        boost_regex) 

endif ( BUILD_DIFFERENTIATE_INTERPOLATE_HRIR )
