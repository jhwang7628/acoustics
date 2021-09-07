# =============================================================================
#   Elasticity solver (builds mass and stiffness matrices)
# =============================================================================
if ( BUILD_ELASTICITY_SOLVER )
    # ===================
    # Target
    # ===================
        add_executable(elasticity_solver
            elasticity_solver.cpp
            ${LOGGING_SRC}
            ${UTILS_SRC}
            ${DEFORMABLE_SRC}
            ${LINEARALGEBRA_SRC}
            ${GENERIC_SRC}
            ${IO_SRC}
            )

        target_link_libraries(elasticity_solver
            boost_program_options
            rt
            ${MKL_LIBRARY})

        set_target_properties(elasticity_solver PROPERTIES LINK_FLAGS ${LINK_FLAGS})
endif ( BUILD_ELASTICITY_SOLVER )
