# =============================================================================
#   Builds tools for writing rigid-body data to disk (inertia and mass)
# =============================================================================
if ( BUILD_INIT_RIGID_TOOLS )
    add_executable(write-inertia
        write_inertia.cpp
        ${GENERIC_SRC}                    ${IO_SRC})

    add_executable(write-mass-center
        write_mass_center.cpp
        ${GENERIC_SRC}                    ${IO_SRC})

    add_executable(init-mass-centers
        init_mass_centers.cpp
        ${GENERIC_SRC}                    ${IO_SRC}
        ${GEOMETRY_SRC}                   ${LINEARALGEBRA_BLAS_LAPACK_SRC})

    target_link_libraries(init-mass-centers
        ${BLAS_LAPACK_LIBS}               ${MKL_LIBRARY})
endif ( BUILD_INIT_RIGID_TOOLS )