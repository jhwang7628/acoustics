# =============================================================================
#   Builds a command-line PAN wavesolver
# =============================================================================
if ( BUILD_RUN_PULSE_SOLVE )
    add_executable(run-pulse-solve
        run_pulse_solve.cpp
        ${GENERIC_SRC}                    ${IO_SRC}
        ${LOGGING_SRC}                    ${LINEARALGEBRA_SRC}
        ${DISTANCEFIELD_SRC}              ${MATH_SRC}
        ${LINEARALGEBRA_BLAS_LAPACK_SRC}         ${GEOMETRY_SRC}
        ${FIELD_SRC}                      ${WAVESOLVER_SRC}
        ${UTILS_SRC}                      ${PARSER_SRC}
        ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
        ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC})

    target_link_libraries(run-pulse-solve
        ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}                  boost_thread
        boost_system                    ${MKL_LIBRARY})
endif ( BUILD_RUN_PULSE_SOLVE )

# =============================================================================
#   Builds a program which reads the output of run-pulse-solve or
#   precompute-pan-gui and converts it in to a multi-term PAN approximation
# =============================================================================
if ( BUILD_MULTITERM_PAN )
    add_executable(build-multiterm-pan
        build_pulse_approximation_system.cpp
        ${GENERIC_SRC}                    ${IO_SRC}
        ${LOGGING_SRC}                    ${LINEARALGEBRA_SRC}
        ${DISTANCEFIELD_SRC}              ${MATH_SRC}
        ${LINEARALGEBRA_BLAS_LAPACK_SRC}         ${GEOMETRY_SRC}
        ${FIELD_SRC}                      ${WAVESOLVER_SRC}
        ${UTILS_SRC}                      ${PARSER_SRC}
        ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
        ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC})

    target_link_libraries(build-multiterm-pan
        ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}                  boost_thread
        boost_system                    ${MKL_LIBRARY})
endif ( BUILD_MULTITERM_PAN )

# =============================================================================
#   Builds a program which reads the output of run-pulse-solve or
#   precompute-pan-gui and converts it in to a multi-term PAN approximation
#   with wavelet-based compression
# =============================================================================
if ( BUILD_COMPRESSED_PAN )
    add_executable(build-compressed-pan
        build_compressed_pulse_approximation.cpp
        ${GENERIC_SRC}                    ${IO_SRC}
        ${LOGGING_SRC}                    ${LINEARALGEBRA_SRC}
        ${DISTANCEFIELD_SRC}              ${MATH_SRC}
        ${LINEARALGEBRA_BLAS_LAPACK_SRC}         ${GEOMETRY_SRC}
        ${FIELD_SRC}                      ${WAVESOLVER_SRC}
        ${UTILS_SRC}                      ${PARSER_SRC}
        ${GTS_GEOMETRY_SRC}               ${CGAL_GEOMETRY_SRC}
        ${TRANSFER_SRC}                   ${WAVELET_TRANSFER_SRC})

    target_link_libraries(build-compressed-pan
        ${BLAS_LAPACK_LIBS}             ${GTS_LIBRARIES}
        ${GLIB_LIBRARIES}               ${CGAL_LIBRARIES}
        ${GSL_LIBRARY}                  boost_thread
        boost_system                    ${MKL_LIBRARY})
endif ( BUILD_COMPRESSED_PAN )