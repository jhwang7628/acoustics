# =============================================================================
#   Builds a tool for generating rigid-body sound, based on existing transfer
#   and impact data
# =============================================================================
if ( BUILD_RIGID_SNDGEN )
    add_executable(sound-generator
        RigidFBemSndGen.cpp
        ${SNDGEN_SRC}                     ${MATH_SRC}
        ${GENERIC_SRC}                    ${IO_SRC})

    target_link_libraries(sound-generator
        ${LibConfig_LIBRARY}              ${BLAS_LAPACK_LIBS}
        sndfile                           ${MKL_LIBRARY})
endif ( BUILD_RIGID_SNDGEN )

# =============================================================================
#   Builds a tool for estimating Hertz time scales for each impact in a
#   simulation and writing this data out to a set of files which can later
#   be used by sound-generator
# =============================================================================
if ( BUILD_ESTIMATE_TIMESCALE )
    add_executable(estimate-timescales
        estimate_timescales.cpp
        ${IO_SRC}                         ${TRANSFER_SRC}
        ${UTILS_SRC}                      ${GEOMETRY_SRC}
        ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
        ${PARSER_SRC}                     ${LINEARALGEBRA_SRC}
        ${LINEARALGEBRA_BLAS_LAPACK_SRC}         ${MATH_SRC}
        ${WAVELET_TRANSFER_SRC})

    target_link_libraries(estimate-timescales
        ${CGAL_LIBRARIES}                 ${GTS_LIBRARIES}
        ${BLAS_LAPACK_LIBS}               ${GSL_LIBRARY}
        boost_thread                      boost_system
        ${MKL_LIBRARY})
endif ( BUILD_ESTIMATE_TIMESCALE )