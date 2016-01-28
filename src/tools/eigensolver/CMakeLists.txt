# =============================================================================
#   Modal eigensolver (eigensolver)
# =============================================================================
if ( BUILD_EIGENSOLVER )
    # ===================
    # Dependencies
    # ===================

        find_package(ARPACK REQUIRED)
        find_package(SuperLU REQUIRED)

    # ===================
    # Target: eigensolver
    # ===================
        add_executable(eigensolver
            eigensolver.cpp
            ${LOGGING_SRC}        ${UTILS_SRC}
            ${IO_SRC}             ${SUPERLU_INTERFACE_SRC}
            ${LINEARALGEBRA_BLAS_LAPACK_SRC}
            )

        target_include_directories(eigensolver
            PUBLIC ${ARPACKPP_INCLUDE_DIR}
            PUBLIC ${SUPERLU_INCLUDES})

        target_link_libraries(eigensolver
            ${SUPERLU_LIBRARIES}  ${ARPACK_LIBRARY}
            ${BLAS_LAPACK_LIBS}   ${MKL_LIBRARY}
            boost_program_options rt)

        set_target_properties(eigensolver PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( BUILD_EIGENSOLVER )

# =============================================================================
#   Modal eigensolver (arpack-eigensolver)
# =============================================================================
if ( BUILD_EIGENSOLVER )
    # ===================
    # Dependencies
    # ===================

        find_package(ARPACK REQUIRED)
        find_package(SuperLU REQUIRED)

    # ===================
    # Target: arpack-eigensolver
    # ===================
        add_executable(arpack-eigensolver
            arpack-eigensolver.cpp
            ${EIGENSOLVER_SRC}
            ${LOGGING_SRC}        ${UTILS_SRC}
            ${IO_SRC}             ${SUPERLU_INTERFACE_SRC}
            ${LINEARALGEBRA_BLAS_LAPACK_SRC}
            )

        target_include_directories(arpack-eigensolver
            PUBLIC ${ARPACKPP_INCLUDE_DIR}
            PUBLIC ${SUPERLU_INCLUDES})

        target_link_libraries(arpack-eigensolver
            ${SUPERLU_LIBRARIES}  ${ARPACK_LIBRARY}
            ${BLAS_LAPACK_LIBS}   ${MKL_LIBRARY}
            boost_program_options rt)

        set_target_properties(arpack-eigensolver PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( BUILD_EIGENSOLVER )